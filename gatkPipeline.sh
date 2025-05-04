#!/bin/bash
set -euo pipefail

####################
# CONFIGURATION
####################
THREADS=$(nproc)                   # Uses up to 6 cores for parallel tasks, if available
if (( THREADS > 10 )); then
  THREADS=10
fi
REF_FASTA="Common/universe.fasta"  # Reference genome
TARGET_BED="Common/targetsPad100.bed" # Target regions
SAMPLE_IDS=("child" "father" "mother") # Sample IDs: child 0, father 1 and mother 2
DOCKER_GATK="sudo docker run -v $(pwd)/Common:/ref -v $(pwd)/out:/data broadinstitute/gatk"	# GATK via Docker, direct installation is possible but not recommended

####################
# PIPELINE START
####################

# 1. Alignment (BWA-MEM + multithreading)
echo "=== STEP 1/6: Alignment ==="
# Check if the reference genome is indexed, otherwise it creates the indexes for bwa
if [[ ! -f "Common/universe.fasta.bwt" ]]; then
  echo "Indexing the reference genome..."
  bwa index Common/universe.fasta
else
  echo "Reference genome is already indexed."
fi
# Check if there's the dictionary file, otherwise it creates the indexes for bwa
if [[ ! -f "Common/universe.dict" ]]; then
  echo "Creating the dictionary file of the reference genome..."
  $DOCKER_GATK gatk CreateSequenceDictionary \
    -R /ref/universe.fasta \
    -O /ref/universe.dict
else
  echo "Dictionary file is already present."
fi
for SAMPLE in "${SAMPLE_IDS[@]}"; do		# For each patient
	if [[ ! -f "out/${SAMPLE}.bam" ]]; then	# If no corresponding BAM file is in the out folder proceed
		echo "Aligning ${SAMPLE}..."			# Alignment is done via BWA MEM, which is more precise than Bowtie2 for clinical cases
	# Here a read group (@RG) header is added to the alignment. ID:${SAMPLE} is the unique identifier, SM:${SAMPLE} is the sample name.
	# The output is piped to samtools to be sorted and saved as a .bam
	# Then the .bam is indexed by samtools that creates the corresponding .bam.bai index file
		bwa mem -t $THREADS -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
		$REF_FASTA "in/${SAMPLE}.fq.gz" | \
		samtools sort -@ $THREADS -o "out/${SAMPLE}.bam" -
		samtools index -@ $THREADS "out/${SAMPLE}.bam"
	else
		echo "Skipping ${SAMPLE}, BAM file already exists."
	fi
done

# 2. Generate GVCFs (Parallelized)
echo "=== STEP 2/6: GVCF Generation ==="
for SAMPLE in "${SAMPLE_IDS[@]}"; do
  {
    echo "Generating GVCF for ${SAMPLE}..."
	# After generating GVCFs for all samples, you can combine them into a multi-sample GVCF and perform joint genotyping. This way, you can detect rare variants in a population by leveraging information across multiple samples at once.
    $DOCKER_GATK gatk HaplotypeCaller \
      -R /ref/$(basename $REF_FASTA) \
      -I /data/${SAMPLE}.bam \
      -O /data/${SAMPLE}.g.vcf.gz \
      -ERC GVCF \
      -L /ref/$(basename $TARGET_BED) \
      --native-pair-hmm-threads $((THREADS/2))	# Paired Hidden Markov Model it's resource intensive, half of the available threads are dedicated
  } &	# The & allows for concurrent operation
done

# Wait for all background jobs to finish
wait


# 3. Joint Genotyping
echo "=== STEP 3/6: Joint Genotyping ==="
#	GenomicsDBImport combines the 3 individual GVCFs into a GenomicsDB workspace, saved in the folder trio_db
$DOCKER_GATK gatk GenomicsDBImport \
  -V /data/child.g.vcf.gz \
  -V /data/father.g.vcf.gz \
  -V /data/mother.g.vcf.gz \
  --genomicsdb-workspace-path /data/trio_db \
  -L /ref/$(basename $TARGET_BED)
#	GenotypeGVCFs is used to call variants across all samples at once from the GenomicsDB created earlier
$DOCKER_GATK gatk GenotypeGVCFs \
  -R /ref/$(basename $REF_FASTA) \
  -V gendb:///data/trio_db \
  -O /data/trio_joint.vcf.gz \

# 4. Variant Filtering (Trio-aware)
echo "=== STEP 4/6: Variant Filtering ==="
#	Filter out low quality reads and only keeps variants with at least 10 read depths
bcftools view out/trio_joint.vcf.gz | \
  bcftools filter -i 'QUAL>=10 && INFO/DP>=10' | \
  bcftools view -T $TARGET_BED -Oz -o out/trio_filtered.vcf.gz
bcftools index out/trio_filtered.vcf.gz

# 5. Inheritance Patterns
# Recall that Child is 0, Father 1 and Mother 2 (see upstream)
echo "=== STEP 5/6: Inheritance Analysis ==="
# Autosomal Recessive (Compound Het/Hom)
#	add  || (GT[0]="het" && GT[1]="het" && GT[2]="het") to add healthy carriers children, which may be relevant for some diseases
bcftools view out/trio_filtered.vcf.gz | \
  bcftools filter -i '(GT[0]="hom" && GT[1]="het" && GT[2]="het")' \
  -Oz -o out/AR_candidates.vcf.gz

# Autosomal Dominant (De Novo Variant assumption)
bcftools view out/trio_filtered.vcf.gz | \
bcftools filter -i ‘(GT[0]=”het”||GT[0]=”homalt”)&&GT[1]=”homref”&& GT[2]=”homref”’ \
  -Oz -o out/de_novo_candidates.vcf.gz

# 6. Cleanup & Report
echo "=== STEP 6/6: Cleanup ==="
rm -rf out/trio_db  # Remove temp database
echo "=== Pipeline Complete ==="
echo "AR candidates: out/AR_candidates.vcf.gz"
echo "De novo candidates: out/de_novo_candidates.vcf.gz"
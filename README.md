# Trio Analysis Pipeline
Genomic analysis pipeline for identifying rare disease variants in familial trios. The goal is to determine wheter the child is affected by a rare genetic desease in chr16.
  ## Students
- [Oberti Gabriele](https://github.com/QuellObii)  
- [Rosa Michelangelo](https://github.com/michelangelorosa)

## Quick Start

  

### 1. Prerequisites

  

Software:

bash

sudo apt-get install -y bwa samtools bcftools docker

docker: broadinstitute/gatk:4.1.3.0  

# Input Files:

  

Common/

  

├── universe.fasta # Reference genome

  

├── targetsPad100.bed # Target regions

  

in/

  

├── child.fq.gz # Child's reads

  

├── father.fq.gz # Father's reads

  

├── mother.fq.gz # Mother's reads

### 2. Run Pipeline

chmod +x gatkPipeline.sh

Place trio files in the "in" folder  

./gatkPipeline.sh

## Output Files

Output files will be saved in the out folder

Autosomal recessive candidates: out/AR_candidates.vcf.gz

De novo mutation candidates :out/de_novo_candidates.vcf.gz

## Configuration

Edit these variables in the script if needed:

  

REF_FASTA="Common/universe.fasta"

  

TARGET_BED="Common/targetsPad100.bed"

## Runtime
⏱️ ~30 minutes on typical hardware (16GB RAM, 4 cores)
## ⚠️ Important Legal Notice
This pipeline is for RESEARCH USE ONLY. It is not intended for and must not be used for diagnostic purposes. Clinical interpretation requires review by qualified geneticists and additional validation.

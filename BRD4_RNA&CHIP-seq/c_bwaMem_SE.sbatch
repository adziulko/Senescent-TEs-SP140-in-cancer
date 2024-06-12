#!/bin/bash

## Script for running bwa
## Date: 17 Sept 2020

## Example usage:
## sbatch --array 0-21 c_bwaMem_SE.sbatch

## bwaIndexDir options: hg38, mm10, bosTau8, canFam3
## bwaIndex options: hg38.main.fa, mm10.fa, bosTau8.fa, canFam3.fa

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set variables
inDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/b_trimmed_fastq
outDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/c_bwa_bam
bwaIndexDir=/scratch/Shares/chuong/genomes/hg38
bwaIndex=hg38.main.fa

# Load bwa and samtools
module load bwa/0.7.15 samtools/1.8

# Define query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run bwa
pwd; hostname; date

echo "bwa version: 0.7.15-r1140"
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")
echo $(date +"[%b %d %H:%M:%S] Converting bams to sams, filtering out low quality reads, and sorting...")

bwa mem -t 8 \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools view -Sb -q 10 - \
| samtools sort -@ 8 - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <int' - number of threads to use
# '-Sb' - convert sam (bwa output) to bam
# '-q <int>' - skip alignments with MAPQ (i.e. quality) score smaller than <int>
# '-@ <int>' - number of threads to use

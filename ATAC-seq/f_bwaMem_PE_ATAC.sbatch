#!/bin/bash

## Script for running bwa
## Date: 26 Sept 2020
## for d_p300

## Example usage:
## sbatch --array 0-81 c_bwaMem_SE.sbatch

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
inDir=/Shares/CL_Shared/data/adam/k_RS_WI38/atac_seq/b_bbduk
outDir=/Shares/CL_Shared/data/adam/k_RS_WI38/atac_seq/c_bwa_bam
bwaIndexDir=/scratch/Shares/chuong/genomes/hg38
bwaIndex=hg38.main.fa

# Load bwa and samtools
module load bwa/0.7.15 samtools/1.8

# Define query files
queries=($(ls ${inDir}/*bbduk.fq.gz | xargs -n 1 basename | sed 's/_1_bbduk.fq.gz//g' | \
sed 's/_2_bbduk.fq.gz//g' | uniq))

# Run bwa
pwd; hostname; date

mkdir -p ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

echo "bwa version: 0.7.15-r1140"
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")
echo $(date +"[%b %d %H:%M:%S] Converting bams to sams, filtering out low quality reads, and sorting...")

bwa mem -t 8 \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_bbduk.fq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_bbduk.fq.gz \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

# Generate flagstat files from unfiltered sam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from unfiltered sam...")

samtools flagstat \
-@ 8 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_flagstat.txt

# Convert sam files to bam files, filtering out low quality (<10) and unmapped (through '-F 4' option)
echo $(date +"[%b %d %H:%M:%S] Converting sam to unfiltered bam & sorting")

samtools view \
-@ 8 \
-Sb \
-q 10 \
-F 4 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
| samtools sort \
-@ ${numThreads} \
-T ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} \
- \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted bam...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing sam...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <int' - specifies number of threads to use
# '-Sb' - convert sam (bwa output) to bam
# '-q <int>' - skip alignments with MAPQ (i.e. quality) score smaller than <int>
# '-@ <int>' - specifies number of threads to use

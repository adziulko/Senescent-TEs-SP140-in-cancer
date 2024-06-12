#!/bin/bash

## Script for running hisat2 (RNA-seq read mapping to reference genome)
## Date: 23 Aug 2021

## Example usage:
## cd /Shares/CL_Shared/data/rotations/adam/aging_TE_rotation_2020/sbatch_scripts
## inDir=/Shares/CL_Shared/data/rotations/adam/aging_TE_rotation_2020/RNA_Seq/trim_fastq \
## outDir=/Shares/CL_Shared/data/rotations/adam/aging_TE_rotation_2020/RNA_Seq/hisat2 \
## hisatIndexDir=/Shares/CL_Shared/db/genomes/mm10/index/hisat2 \
## hisatIndex=mm10 \
## rnaStrandness=FR \
## sbatch --array 0-2 hisat2_PE.sbatch

## hisatIndexDir options: hg38, mm10, bosTau8
## hisatIndex: hg38.main, mm10, bosTau8
## rnaStrandness options: FR (forward) or RF (reverse stranded)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J hisat2_PE
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set constant variables
numThreads=8
inDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/b_trimmomatic/trinity
outDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/g_hisat2
hisatIndexDir=/Shares/CL_Shared/db/genomes/hg38/index/hisat2
hisatIndex=hg38
rnaStrandness=rf

# Set query files
# note: remove 1/2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*paired.fastq.gz | xargs -n 1 basename | sed 's/_1paired.fastq.gz//g' \
	| sed 's/_2paired.fastq.gz//g' | uniq))

# Load hisat2 and samtools
module load hisat2 samtools

# Run hisat2
pwd; hostname; date

echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running hisat2...")

hisat2 \
-p ${numThreads} \
--rna-strandness ${rnaStrandness} \
--no-softclip \
-x ${hisatIndexDir}/${hisatIndex} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1paired.fastq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2paired.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

# Create .sorted.uniq.bai index files for upload to UCSC browser
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-p <int>' - number of threads
# '--rna-strandness <rf/fr>' - specify strand-specific information; library kit typically has this information
# '--no-softclip' - disallows soft clipping
# '-x <file.basename>' - hisat2 index basename; include everything up to .main.1.ht2
# '-1 <file_1_paired.fastq.gz' - mate file 1
# '-2 <file_2_paired.fastq.gz' - mate file 2
# '-q <int>' - filter out mapped reads by quality threshold
# '-Sb' - convert sam to bam
# '-@' - number of threads

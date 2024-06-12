#!/bin/bash

## Script for merging and sorting
## Date: 21 July 2021

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=24GB

# Job name and output
#SBATCH -J merge_sort.Bam
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Define key variables
workDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/sorted_bam \


# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bam file...")


samtools merge \
-@ 8 \
${workDir}/sen_123_merge.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/sorted_bam/sen_br1_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/sorted_bam/sen_br2_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/sorted_bam/sen_br3_sorted.bam \


echo "Samtools version: see stdout" 'samtools --version'


echo "sorting merged bam file"
samtools sort \
-@ 8 \
-o ${workDir}/sen_123_merge.sorted.bam \
${workDir}/sen_123_merge.bam


echo "index resorted bams"
samtools index \
-@ 8 \
${workDir}/sen_123_merge.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

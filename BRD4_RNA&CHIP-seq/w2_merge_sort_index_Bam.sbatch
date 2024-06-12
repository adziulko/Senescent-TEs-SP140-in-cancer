#!/bin/bash

## Script for merging, sorting, and indexing bam files
## Date: 16 July 2021

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=24GB

# Job name and output
#SBATCH -J merge_sort_index_trinity.Bam
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Define key variables
workDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity \

# define query bam files
# queries=($(ls ${workDir}/*br*bam | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bam file...")


samtools merge \
-@ 8 \
${workDir}/qui_sen_merge.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/qui_br1_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/qui_br2_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/qui_br3_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/sen_br1_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/sen_br2_sorted.bam \
/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/sen_br3_sorted.bam \


echo "Samtools version: see stdout" 'samtools --version'


echo "sorting merged bam file"
samtools sort \
-@ 8 \
-o ${workDir}/qui_sen_merge.sorted.bam \
${workDir}/qui_sen_merge.bam


echo "index resorted bams"
samtools index \
-@ 8 \
${workDir}/qui_sen_merge.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

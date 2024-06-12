#!/bin/bash

## Script indexing bam files
## I only indexed pro/qui/sen H3 bam files
## Date: 10 Dec 2020
## Example usage: sbatch --array 0-5 s_indexBam.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J indexBam
#SBATCH -o /Users/%u/slurm_out/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Define key variables
inDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/c_bwa_bam \

# define query bam files
queries=($(ls ${inDir}/*H3*_r*bam | xargs -n 1 basename))

# load modules
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bam file...")


samtools index ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

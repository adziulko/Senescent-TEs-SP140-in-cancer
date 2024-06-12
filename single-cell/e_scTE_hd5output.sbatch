#!/bin/bash

## Example usage:
## sbatch --array 0-10 scTE_hd5output.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J scTE
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err

# load modules
module load samtools

# set query files
queries=($(ls *_clean.bam | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting scTE...")

scTE \
-i ${queries[$SLURM_ARRAY_TASK_ID]} \
-o ${queries[$SLURM_ARRAY_TASK_ID]%_clean.bam}_scTE \
-x /Shares/CL_Shared/db/genomes/hg38/index/scTE/hg38.exclusive.idx \
-p 8 -CB CB -UMI UB --hdf5 True

# can also add --hdf5 True to produce a .h5 file instead of .csv

echo $(date +"[%b %d %H:%M:%S] Done!")

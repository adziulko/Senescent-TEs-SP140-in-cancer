#!/bin/bash

## Example usage:
## sbatch --array 0-x gzip.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J gzip
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err

## set variables
inDir=/scratch/Users/addz6536/q_single_cell_ovary_lily/rna_seq/a_lily_raw_fastq

# set query files
queries=($(ls $inDir/*.fastq | xargs -n 1 basename))

# run the thing
pwd; hostname; date

# keep original files, for now
cd $inDir
gzip -c ${queries[$SLURM_ARRAY_TASK_ID]} > ${queries[$SLURM_ARRAY_TASK_ID]}.gz
cd ..
# only delete original uncompressed files after all files have finished gzipping

echo $(date +"[%b %d %H:%M:%S] Done!")

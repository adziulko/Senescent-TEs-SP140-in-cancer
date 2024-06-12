#!/bin/bash

## Script to run cellranger count on single cell fastq files
##
## Example usage:
## inDir=C109 sbatch --array 0-9 cellranger_count.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J cellranger
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err


# increase ulimit to what cellranger recommends
ulimit -n 16000

# echo the new ulimit (for troubleshooting later)
bash -c 'ulimit -a'

# define query files
inDir=/scratch/Users/addz6536/q_single_cell_ovary_lily/rna_seq/a_lily_raw_fastq
outDir=/scratch/Users/addz6536/q_single_cell_ovary_lily/rna_seq/b_cell_ranger
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename | cut -d "_" -f 1,2,3 | sort | uniq))

# try running cellranger count on one sample
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting cellranger counts...")

cd $outDir

cellranger count \
--id=${queries[$SLURM_ARRAY_TASK_ID]} \
--fastqs=$inDir \
--sample=${queries[$SLURM_ARRAY_TASK_ID]} \
--transcriptome=/Shares/CL_Shared/db/genomes/hg38/index/CellRanger/refdata-gex-GRCh38-2020-A

cd ..

echo $(date +"[%b %d %H:%M:%S] Done!")

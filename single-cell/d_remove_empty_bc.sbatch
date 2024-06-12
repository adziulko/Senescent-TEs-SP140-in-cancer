#!/bin/bash

## Example usage:
## sbatch --array 0-10 remove_empty_CB_reads.sbatch

#General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J remove_empty_CB_reads
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err

# load modules
module load samtools

# define query files
inDir=/scratch/Users/addz6536/q_single_cell_ovary_lily/rna_seq/c_bams

# set query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Removing CB reads…")

cd $inDir

samtools view ${queries[$SLURM_ARRAY_TASK_ID]} -h \
| awk '/^@/ || /CB:/' \
| samtools view -h -b \
> ${queries[$SLURM_ARRAY_TASK_ID]}.clean

echo $(date +"[%b %d %H:%M:%S] Indexing bam...")

samtools index ${queries[$SLURM_ARRAY_TASK_ID]}.clean

cd ..

echo $(date +"[%b %d %H:%M:%S] Done!")

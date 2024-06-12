#!/bin/bash

## Script for running fastqc using an array
## Date: 02 Nov 2020

## Example usage:
##sbatch --array 0-71 o_STAR_map.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 2
#SBATCH -n 24
#SBATCH --time=1:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J star
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Load FastQC
module load STAR

# Define constant variables
numThreads=16
inDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/b_trimmomatic
outDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star
genDir=/Shares/CL_Shared/db/genomes/hg38/index/STAR/gencode.v34_refSeq

#Define folder names
queries=($(ls ${inDir}/*_paired.fastq.gz | xargs -n 1 basename | sed 's/_1_paired.fastq.gz//g' \
	| sed 's/_2_paired.fastq.gz//g' | uniq))

# Run FastQC
pwd; hostname; date

echo "Running STAR..."
echo "FastQC version: "$(STAR --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

STAR --runThreadN ${numThreads} \
--genomeDir ${genDir} --readFilesIn ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_paired.fastq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_paired.fastq.gz \
--readFilesCommand zcat \
--alignIntronMin 20 --alignIntronMax 1000000 \
--outFileNamePrefix ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 \
--outSAMtype BAM Unsorted SortedByCoordinate


echo $(date +"[%b %d %H:%M:%S] Done")

#STAR mapping options optimized for TEtranscript (multimapping)

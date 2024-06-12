#!/bin/bash

## Script for converting bedgraphs to bigwigs
## Date: 6 Dec 2023

## Example usage:
## sbatch f_bdg_to_bigwig.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bedGraphToBigWig
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set variables
inDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/h_meme/JUND_bed_output
outDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/h_meme/JUND_bed_output
chromSizes=/Shares/CL_Shared/db/genomes/hg38/fa/hg38.chrom.sizes

# Define query files
queries=($(ls ${inDir}/MA0492.1.bdg | xargs -n 1 basename))

# Convert bedgraphs to bigwigs
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort \
-k1,1 -k2,2n \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")
echo "bedGraphToBigWig v 4"

/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg \
${chromSizes} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.bw

echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Done")

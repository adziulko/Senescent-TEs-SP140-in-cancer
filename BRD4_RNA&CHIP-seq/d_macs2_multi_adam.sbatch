#!/bin/bash

## Example usage:
## cd /Shares/CL_Shared/data/adam/c_BRD4/sbatch_scripts
## sbatch --array 0-2 d_H3_macs2_multi.sbatch
## for c_BRD4

## General settings
#SBATCH -p long # Partition or queue. In this case, short!
#SBATCH --nodes=1 # Only use a single node
#SBATCH --cpus-per-task=32 # cpus (max = 64)
#SBATCH --mem=64gb # Memory limit (max = 512 gb per node, but only request as needed)
#SBATCH --time=4:00:00 # Time limit hrs:min:sec

## Job name and output
#SBATCH --job-name=macs2
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

## Email Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

## Load modules
module load python/2.7.14/MACS/2.1.1

## Code to execute

inDir="/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/c_bwa_bam" # sorted alignment of actual ChIP-seq library
outDir="/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/d_macs2_hs"
MACSGENOME="hs" # either the size of mappable genome eg 3e9 or a species abbreviation, hs = human, mm = mouse.

queries=($(ls ${inDir}/*Brd*_r2* | xargs -n 1 basename | sed 's/.sorted.uniq.bam//g' | uniq)) #use *K4me3* for other job
controlFiles=($(ls ${inDir}/*none* | xargs -n 1 basename | sed 's/.sorted.uniq.bam//g' | uniq))


# Call peaks with macs2
pwd; hostname; date

echo "Calling peaks with macs2"
echo "macs2 version: "$(macs2 --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

# Run MACS2
macs2 callpeak \
-t ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam \
-c ${inDir}/${controlFiles[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam \
-n ${queries[$SLURM_ARRAY_TASK_ID]} \
-g $MACSGENOME \
--outdir $outDir \
--SPMR -B

echo $(date +"[%b %d %H:%M:%S] Done")

# --SPMR and -B are to make a normalized bedgraph file
# -g hs sets genome size to homo sapiens preset (2.7e9)
# Effective genome size. It can be 1.0e+9 or 1000000000,
#                       or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
#                       (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
#                       fruitfly (1.2e8), Default:hs
# use "-f BAMPE" if using paired end data

# Using effective genome size of 3049315783 for hg38
# aligns all fastq files in a folder into bam files (unique + multimap)

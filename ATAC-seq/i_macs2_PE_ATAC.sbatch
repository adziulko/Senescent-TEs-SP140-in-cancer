#!/bin/bash

## Script for running macs2
## Date: 5 May 2022

## Example usage:
## sbatch --array=0-4 macs2_PE.sbatch


## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=12:00:00
#SBATCH --mem=16G

# Job name and output
#SBATCH -J macs2
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Define variables
inDir=/Shares/CL_Shared/data/adam/k_RS_WI38/atac_seq/c_bowtie2
outDir=/Shares/CL_Shared/data/adam/k_RS_WI38/atac_seq/d_macs2_hs
MACSGENOME="hs"

# Module load macs2
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Call peaks with macs2
pwd; hostname; date

echo "Calling peaks with macs2"
echo "macs2 version: "$(macs2 --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

macs2 callpeak \
--treatment ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam} \
--outdir ${outDir} \
-g $MACSGENOME \
--SPMR -B --keep-dup all --format BAMPE --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '--treatment <file>' - input file
# '--name <name>' - file name; will be used to name output
# '-B' - whether or not to save extended fragment pileup and local lambda tracks at every bp into a bedGraph file (default = false)
# '--SPMR' - save signal per million reads for fragment pileup profiles; requires -B to be set (default = false)
# '--gsize <int>' - effective genome size; can be 1.0e+9 or 1000000000; can also use shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
# '--format <AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE>' - format of tag file (default = AUTO); use '-f BAMPE' if using paired end data
# '--outdir <outDir>' - output will be written here
# '--keep-dup' - specifies whether macs2 will retain or discard duplicates

#!/bin/bash

## Script to download a bunch of sras
## Date: 10 Sep 2020

## Example usage:
## sra_file= SraRunTable.txt
## outDir=/Shares/CL_Shared/data/adam/senescent_L1_IFN/rna_seq/c_raw_fastq
## sbatch --array 0-21 a_sraDownload_multi_adam.sbatch      (chip:x=22)

## This was input into command line (in sbatch directory):
## sbatch --array 0-8 faster.sbatch

## General settings
#SBATCH -p short                 #Partition or queue (short = <24 hours)
#SBATCH --nodes=1                #Only use a single node
#SBATCH --cpus-per-task=8        #cpus (max = 64)
#SBATCH --mem=16gb               #Memory limit (max = 512 gb per node, but only request as needed)
#SBATCH --time=5:00:00           #Time limit hrs:min:sec

## Job name and output
#SBATCH --job-name=sraDownload                        #Job name
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out  #Standard output and error log
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err   #%A is job number while %a provides number for each sra

## Email Settings
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=addz6536@colorado.edu

## Load modules
module load sra

## set variables
sra_file=SraRunTable.txt
outDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/a_raw_fastq
inDir=/Shares/CL_Shared/data/adam/c_BRD4

# set query files
sra_list=($(cat $inDir/$sra_file | awk '(NR>1){print $1}'))

## Code to execute
pwd; hostname; date
echo $(date +"[%b %d %H:%M:%S] Start download...")

## download SRR files
fasterq-dump -O ${outDir} --split-files --gzip ${sra_list[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

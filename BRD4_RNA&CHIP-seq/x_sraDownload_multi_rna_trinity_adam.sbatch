#!/bin/bash

## Script to download a bunch of sras
## Date: 12 July 2021

## Example usage:

## General settings
#SBATCH -p short                 #Partition or queue (short = <24 hours)
#SBATCH --nodes=1                #Only use a single node
#SBATCH --cpus-per-task=8        #cpus (max = 64)
#SBATCH --mem=16gb               #Memory limit (max = 512 gb per node, but only request as needed)
#SBATCH --time=2:00:00           #Time limit hrs:min:sec

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
sra_file=SraRunTable2_RNAseq.txt
outDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/a_raw_fastq/trinity
inDir=/Shares/CL_Shared/data/adam/c_BRD4

# set query files
sra_list=($(cat $inDir/$sra_file | awk '(NR>1){print $1}'))

## Code to execute
pwd; hostname; date
echo $(date +"[%b %d %H:%M:%S] Start download...")

## download SRR files
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ${sra_list[$SLURM_ARRAY_TASK_ID]} --gzip  -O ${outDir}


echo $(date +"[%b %d %H:%M:%S] Done!")

#!/bin/bash

## Script to download a bunch of ftp links
## Date: 17 April 2023

## Example usage:
## sbatch a_CCLE_Download_multi_adam.sbatch

## General settings
#SBATCH -p long                 #Partition or queue (short = <24 hours)
#SBATCH --nodes=1                #Only use a single node
#SBATCH --cpus-per-task=8        #cpus (max = 64)
#SBATCH --mem=16gb               #Memory limit (max = 512 gb per node, but only request as needed)
#SBATCH --time=12-00:00           #Time limit hrs:min:sec

## Job name and output
#SBATCH --job-name=CCLE_Download                        #Job name
#SBATCH --output=/Users/%u/slurm_out/slurm-%j.out  #Standard output and error log
#SBATCH --error=/Users/%u/slurm_err/slurm-%j.err   #%j inserts job number

## Email Settings
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=addz6536@colorado.edu


pwd; hostname; date
echo $(date +"[%b %d %H:%M:%S] Start download...")


# download files
xargs -n 1 curl -O < fwd_rev_fastq_download.txt



echo $(date +"[%b %d %H:%M:%S] Done!")

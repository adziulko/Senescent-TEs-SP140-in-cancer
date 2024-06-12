#!/bin/bash

## Script for making md5sum file
## Example usage:

## General settings
#SBATCH -p short
#SBATCH -N 2
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J md5sum
#SBATCH --output=/Users/%u/slurm_out/slurm-%j.out
#SBATCH --error=/Users/%u/slurm_err/slurm-%j.err


## fasta format of reference genome
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on 'hostname'

echo $(date +"[%b %d %H:%M:%S] Running code...")

md5sum *gz > md5sum.txt

echo $(date +"[%b %d %H:%M:%S] Done")

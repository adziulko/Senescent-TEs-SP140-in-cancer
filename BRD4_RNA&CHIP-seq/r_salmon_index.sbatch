#!/bin/bash

## Script for running salmon index
## Date: 21 April 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## fastaFile=/Users/CL_Shared/db/genomes/mm10/fa/gencode.vM18.transcripts.fa \
## outDir=/Users/CL_Shared/db/genomes/mm10/index/salmon/gencode.vM18.transcripts_salmon_index \
## sbatch salmon_index.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J salmon_index
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8
kmerLen=31

# Load salmon
module load salmon

# Run salmon index
pwd; hostname; date

echo "Salmon version: "$(salmon --version)
echo "Processing file: "${fastaFile}
echo "Output location: "${outDir}
echo "Starting salmon index..."

salmon index -p ${numThreads} -k ${kmerLen} --gencode -t ${fastaFile} -i ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-p <int>' - number of threads
# '-k <int>' - build index using auxilliary k-mer hash over k-mers of <int>
# '--gencode' - tells salmon to expect gencode fasta format as input
# '-t <fastaFile>' - transcript fasta file
# '-i <outDir>' - index output directory

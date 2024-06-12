#!/bin/bash

## Script for running trimmomatic
## Date: 16 Sept 2020
## for c_BRD4

## Example usage (chip):
## sbatch --array 0-21 b_trimmomatic_multi_SE.sbatch

## adapterFile options: NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq3-PE-2.fa, TruSeq3-PE.fa

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J trimmomatic
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set constant variables
trimExe=/opt/trimmomatic/0.36/trimmomatic-0.36.jar
numThreads=8
inDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/a_raw_fastq
outDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/b_trimmed_fastq
logOut=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/b_trimmed_fastq/log
adapterFile=/opt/trimmomatic/0.36/adapters/TruSeq3-SE.fa

# Define query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run Trimmomatic
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Trimmomatic version: "$(java -jar ${trimExe} --version)
echo "Starting trimmomatic..."

java -jar ${trimExe} SE -threads ${numThreads} \
-trimlog ${logOut}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}_trimlog.txt \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}_trimmed.fastq.gz \
ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/${adapterFile}:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# 'SE' - single-end input
# '-threads <int>' - number of threads to use
# '-trimlog <file.txt>' - creates a log of all read trimmings
# 'ILLUMINACLIP <file.fa>' - cut adapter and other illumina-specific sequences from the read
# 'SLIDINGWINDOW <int:int>' - performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold.
# 'LEADING <int>' - cut bases off the start of a read, if below a threshold quality
# 'TRAILING <int>' - cut bases off the end of a read, if below a threshold quality
# 'MINLEN <int>' - retain only reads of some length (in bp)

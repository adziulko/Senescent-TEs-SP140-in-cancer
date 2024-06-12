#!/bin/bash

## Script for running bbduk trimming
## Date: 25 May 2023

## Example usage:
## sbatch --array 0-x bbduk_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bbduk
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu



# Set constant variables
numThreads=1
inDir=/scratch/Users/addz6536/r_CCLE/rna_seq/a_raw_fastq/batch3
outDir=/scratch/Users/addz6536/r_CCLE/rna_seq/j_bbduk/batch2
adapterFile=/opt/bbmap/38.05/resources/adapters.fa

#change fq file names to remove unnecessary string:
#for file in *gz; do mv "$file" "${file/_CKD*L3/}" ; done

# Load module
module load bbmap

# Define query files
# Removes R1/R2 to generate a unique identifier for each pair of files
# Modify to fit input file name if necessary
queries=($(ls ${inDir}/*fq.gz | xargs -n 1 basename | sed --expression='s/_1.fq.gz//g' | \
sed --expression='s/_2.fq.gz//g' | uniq))


# Run bbduk
pwd; hostname; date

echo "bbduk version: "$(bbduk.sh --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running bbduk...")

bbduk.sh \
-Xmx32g \
t=${numThreads} \
in1=${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1.fq.gz \
in2=${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2.fq.gz \
out1=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_bbduk.fq.gz \
out2=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_bbduk.fq.gz \
ref=${adapterFile} \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
maq=10 \
qtrim=r \
trimq=10 \
tpe \
tbo \
literal=AAAAAAAAAAAAAAAAAAAAAAA \
stats=${outDir}/log/${queries[$SLURM_ARRAY_TASK_ID]}.trimstats.txt \
refstats=${outDir}/log/${queries[$SLURM_ARRAY_TASK_ID]}.refstats.txt \
ehist=${outDir}/log/${queries[$SLURM_ARRAY_TASK_ID]}.ehist.txt

echo $(date +"[%b %d %H:%M:%S] Done")

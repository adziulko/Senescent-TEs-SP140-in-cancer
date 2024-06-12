#!/bin/bash

## Script for running salmon quant with paired end short reads
## Date: 26 Apr 2023

## Example usage:
# inDir=first_1000_files_fastq outDir=salmonoutput_IL18R1_FLAIR \
# gtfFile=/scratch/Users/ativ2716/data/FLAIR_GTF/flair_filter_transcripts_IL18R1modified.gtf \
# indexDir=/scratch/Users/ativ2716/data/FLAIR_GTF/FLAIR_index_IL18R1 \
# sbatch --array 0-999 salmonQuantPE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J salmonQuantPE
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Load modules/programs
# Note: Don't load the Fiji version of Salmon
# Updated version salmon 1.8 is located on CL_Shared at: /Shares/CL_Shared/programs/salmon-1.8.0_linux_x86_64/bin
# Put this in your bash profile PATH before running this script

# Set constant variables
indexDir=/Shares/CL_Shared/data/adam/aa_misc/salmon_GTF_SP140/SP140_4_index
gtfFile=/Shares/CL_Shared/data/adam/aa_misc/salmon_GTF_SP140/gencode.v39.4_SP140.transcripts.gtf
inDir=/scratch/Users/addz6536/r_CCLE/rna_seq/a_raw_fastq/batch5
outDir=/scratch/Users/addz6536/r_CCLE/rna_seq/d_salmon_quant_4iso

# Define query files
queries=($(ls ${inDir}/*.gz | xargs -n 1 basename | sed 's/_1.fq.gz//g' | sed 's/_2.fq.gz//g' | uniq))

# Run salmon
pwd; hostname; date

echo "Salmon version: "$(salmon --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting salmon pseudo-alignment...")

salmon quant --libType A -p 8 \
-i ${indexDir} --geneMap ${gtfFile} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1.fq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2.fq.gz \
--validateMappings --rangeFactorizationBins 4 --gcBias -w 10000 \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_transcripts_quant

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-libType <arg>' - library type; specify '-libType A' for automatic detection or see manual for other options
# '-p <int>' - number of threads
# '-i <indexDir>' - location of salmon index directory
# '--geneMap <gtfFile>' - gtf file containing mapping of transcripts to genes; allows for gene-level abundance estimates
# '-1 <file>' - file containing #1 mates
# '-2 <file>' - file containing #2 mates
# '--validateMappings' - validate mappings using alignment-based verification (highly recommended)
# '--rangeFactorizationBins <int>' - factorizes the likelihood used in quantification; common value to select is 4
# '--gcBias' - perform GC bias correction; typically recommended
# '-o <dir>' - output directory
# '-w <num>' - reads mapping to more than this many places won't be considered

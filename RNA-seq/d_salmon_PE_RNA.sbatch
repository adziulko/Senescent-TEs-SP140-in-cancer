#!/bin/bash

## Script for running salmon quant with PE reads
## Date: 10 June 2022

## Example usage:
## sbatch --array 0-17 salmon_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J salmon_quant_PE
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=addz6536@colorado.edu

# Set constant variables
numThreads=8
indexDir=/Shares/CL_Shared/db/genomes/hg38/annotations/gencode.v34.comp.transcripts.fa.salmon_quasi.idx
gtfFile=/Shares/CL_Shared/db/genomes/hg38/annotations/gencode.v34.annotation.gtf
inDir=/Shares/CL_Shared/data/adam/l_LTR8B_KD_IFNB_HT1080/rna_seq/b_bbduk
outDir=/Shares/CL_Shared/data/adam/l_LTR8B_KD_IFNB_HT1080/rna_seq/d_salmon 

# Set query files
# note: remove 1/2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*bbduk.fq.gz | xargs -n 1 basename | sed 's/_1_bbduk.fq.gz//g' \
	| sed 's/_2_bbduk.fq.gz//g' | uniq))

# Load module
module load salmon

# Run salmon
pwd; hostname; date

echo "salmon version: "$(salmon --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Starting salmon pseudo-alignment..."

salmon quant -i ${indexDir} --libType A -p ${numThreads} --geneMap ${gtfFile} \
-1 <(gunzip -c ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_bbduk.fq.gz) \
-2 <(gunzip -c ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_bbduk.fq.gz) \
--validateMappings --rangeFactorizationBins 4 --gcBias \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_transcripts_quant

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-i <indexDir>' - location of salmon index directory
# '-libType <arg>' - library type; specify '-libType A' for automatic detection or see manual for other options
# '-p <int>' - number of threads
# '--geneMap <gtfFile>' - gtf file containing mapping of transcripts to genes; allows for gene-level abundance estimates
# '-1 <file>' - file containing #1 mates
# '-2 <file>' - file containing #2 mates
# '--validateMappings' - validate mappings using alignment-based verification (highly recommended)
# '--rangeFactorizationBins <int>' - factorizes the likelihood used in quantification; common value to select is 4
# '--gcBias' - perform GC bias correction; typically recommended
# '-o <dir>' - ooutput directory

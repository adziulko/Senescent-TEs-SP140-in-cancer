#!/bin/bash

## Script for running deeptools bamCoverage
## Date: 29 Sep 2023

## Example usage:
## sbatch --array 0-8 u_deeptools_bam_to_bigwig.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J deeptools_bamCov
#SBATCH -o /Users/%u/slurm_out/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# load modules
module load singularity
module load samtools

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
inDir=/scratch/Users/addz6536/r_CCLE/rna_seq/k_hisat2
outDir=/scratch/Users/addz6536/r_CCLE/rna_seq/m_bigwig
ccle_file=/scratch/Users/addz6536/r_CCLE/rna_seq/m_bigwig/BW_samples_TEcounts.txt


# set query files
sra_list=($(cat $ccle_file | awk '(NR>1){print $1}'))
#queries=($(ls ${inDir}/*sort*bam | xargs -n 1 basename ))


# Define deeptools parameters
binSize=1
numCPU=16

# Use "bamCoverage" to generate continuous coverage signal from bam files
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Convert bam into bigwig with bamCoverage...")

# Start by generating bigwig from forward strand data
echo $(date +"[%b %d %H:%M:%S] Generate forward strand bigwigs...")

singularity exec --bind /scratch/Users ${deeptools} \
bamCoverage \
-b ${inDir}/${sra_list[$SLURM_ARRAY_TASK_ID]}.k100.bam \
-o ${outDir}/${sra_list[$SLURM_ARRAY_TASK_ID]}_fwd.bw \
--binSize=${binSize} \
-p ${numCPU} \
--normalizeUsing CPM \
--filterRNAstrand forward \
--ignoreForNormalization chrX chrM

# Start by generating bigwig from reverse strand data
echo $(date +"[%b %d %H:%M:%S] Generate reverse strand bigwigs...")

singularity exec --bind /scratch/Users ${deeptools} \
bamCoverage \
-b ${inDir}/${sra_list[$SLURM_ARRAY_TASK_ID]}.k100.bam \
-o ${outDir}/${sra_list[$SLURM_ARRAY_TASK_ID]}_rev.bw \
--binSize=${binSize} \
-p ${numCPU} \
--normalizeUsing CPM \
--filterRNAstrand reverse \
--ignoreForNormalization chrX chrM

echo $(date +"[%b %d %H:%M:%S] Done!")

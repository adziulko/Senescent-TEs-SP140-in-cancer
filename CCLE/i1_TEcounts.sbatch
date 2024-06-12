#!/bin/bash

## Script for running TEcounts
## Date: 28 Sept 2023
## Example usage:
## sbatch p_TEtranscript.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J TEtranscript
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=addz6536@colorado.edu


# variables
ccle_file=/scratch/Users/addz6536/r_CCLE/rna_seq/l_TEtranscripts/samples_TEcounts.txt
inDir=/scratch/Users/addz6536/r_CCLE/rna_seq/k_hisat2
outDir=/scratch/Users/addz6536/r_CCLE/rna_seq/l_TEtranscripts

# set query files
queries=($(cat $ccle_file | awk '(NR>1){print $1}'))


#Run TEcounts
module load samtools
module load python/2.7.14/tetranscripts/2.1.4


TEcount \
-b $inDir/${queries[$SLURM_ARRAY_TASK_ID]}.k100.bam \
--project $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.k100.bam} \
--sortByPos \
--GTF /Shares/CL_Shared/db/genomes/hg38/annotations/hg38.Gencode.V32.main.gtf \
--TE /Shares/CL_Shared/db/genomes/hg38/annotations/GRCh38_GENCODE_rmsk_TE.gtf

echo $(date +"[%b %d %H:%M:%S] Done")
#--stranded reverse

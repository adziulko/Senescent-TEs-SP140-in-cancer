#!/bin/bash

## Script for running TEtranscript
## Date: 4 Nov 2020
## Example usage:
## sbatch p_TEtranscript.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 2
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

# Load TEtranscript
module load python/2.7.14/tetranscripts/2.1.4

# variables
inDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star \
outDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/d_TEtranscript \
geneFile=/Shares/CL_Shared/db/genomes/hg38/annotations/hg38.Gencode.V32.main.gtf \
TEgtf=/Shares/CL_Shared/db/genomes/hg38/annotations/GRCh38_GENCODE_rmsk_TE.gtf

#Run TEtranscript

TEtranscripts --format BAM --sortByPos --mode multi --stranded reverse -n DESeq_default \
-c ${inDir}/qui_br1_sorted.bam \
${inDir}/qui_br2_sorted.bam \
${inDir}/qui_br3_sorted.bam \
-t ${inDir}/sen_br1_sorted.bam \
${inDir}/sen_br2_sorted.bam \
${inDir}/sen_br3_sorted.bam \
--GTF ${geneFile} --TE ${TEgtf} --project ${outDir}/c_qui_t_sen_TEtran_out

##Notes
##--sortByPos - input bam files are sorted ; -i number of iterations to optimize multimap read assignment ; --mode multi - distribute TE count over multiple alignments
##--stranded  - "no" for PE or no stranded in general; "forward" e.g. QIAseq stranded; "reverse" e.g., Illumina TruSeq stranded - makes a huge difference if incorrectly specified
##-n - normalization: DESeq1 by default; TC: total annotated read count; quant: quantile normalization - DON'T bother. Anyway it won't affect the output .cntTable that you will import and analyse in DESeq2.
## NB = job may fail if normalization fails (e.g., no R module loaded ..). No panic, you have the .cntTable anyway.

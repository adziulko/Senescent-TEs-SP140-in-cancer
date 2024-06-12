#!/bin/bash

## Script for running TEtranscript
## Date: 10 June 2022
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
inDir=/Shares/CL_Shared/data/adam/l_LTR8B_KD_IFNB_HT1080/rna_seq/c_star
outDir=/Shares/CL_Shared/data/adam/l_LTR8B_KD_IFNB_HT1080/rna_seq/e_TEtranscripts
geneFile=/Shares/CL_Shared/db/genomes/hg38/annotations/hg38.Gencode.V32.main.gtf
TEgtf=/Shares/CL_Shared/db/genomes/hg38/annotations/GRCh38_GENCODE_rmsk_TE.gtf

#Run TEtranscript

TEtranscripts --format BAM --sortByPos --mode multi --stranded forward -n DESeq_default \
-c ${inDir}/13_Gw2L_S55Aligned.sortedByCoord.out.bam \
${inDir}/2G1L_S44Aligned.sortedByCoord.out.bam \
-t ${inDir}/14_3w3L_S56Aligned.sortedByCoord.out.bam \
${inDir}/18_3_1_L_S60Aligned.sortedByCoord.out.bam \
${inDir}/3_3w1L_S45Aligned.sortedByCoord.out.bam \
${inDir}/4_3w2L_S46Aligned.sortedByCoord.out.bam \
--GTF ${geneFile} --TE ${TEgtf} --project ${outDir}/c_GFPia_t_LTR8KD_LPS_TEtran_out

##Notes
##--sortByPos - input bam files are sorted ; -i number of iterations to optimize multimap read assignment ; --mode multi - distribute TE count over multiple alignments
##--stranded  - "no" for PE or no stranded in general; "forward" e.g. QIAseq stranded; "reverse" e.g., Illumina TruSeq stranded - makes a huge difference if incorrectly specified
##you can figure out strandedness by seeing 'logs' result of salmon. ISR=forward; ISF=reverse; IU=no
##-n - normalization: DESeq1 by default; TC: total annotated read count; quant: quantile normalization - DON'T bother. Anyway it won't affect the output .cntTable that you will import and analyse in DESeq2.
## NB = job may fail if normalization fails (e.g., no R module loaded ..). No panic, you have the .cntTable anyway.

#!/bin/bash

## Script for building a gmap genome index from trinity fasta output. This output will be used to upload to UCSC genome browser.
## Date: 21 July 2021


# General settings
#SBATCH -p short
#SBATCH -N 2
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J gmap
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu


#Define constant varibles
numThreads=14
GenomeDir=/Shares/CL_Shared/programs/gmap/share
DB=hg38.main 
refFile=/scratch/Users/addz6536/c_BRD4/trinity/denovo_trinity/sen_samples/trinity/Trinity.fasta   #trintiy fasta output
outDir=/scratch/Users/addz6536/c_BRD4/trinity/denovo_trinity/sen_samples/trinity/gmap
outName=BRD4_senOnly_denovo_trinity_UCSC_gmap.gff


# Run gmap-build
pwd; hostname; date
echo $(date +"[%b %d %H:%M:%S] Running gmap-build...")

/Shares/CL_Shared/programs/gmap/bin/gmap  -D ${GenomeDir} -d ${DB} \
-A -S -t ${numThreads} \
-b -f 2 \
${refFile} > ${outDir}/${outName}



echo $(date +"[%b %d %H:%M:%S] Done")

#gmap flags:
#-A run both alignment and mapping
#-S show summary alignment
#-t multithread (2 standard + numThread specified)
#-b map bothstrands of genome
#-O keep order of cDNA as in input file
#-f output format: gff3_gene (or 2)= GFF3 gene format, gff3_match_cdna (or 3)= GFF3 cDNA_match format, gff3_match_est (or 4) = GFF3 EST_match format,

#!/bin/bash

## Script for running bwa
## Date: 27 Oct 2021

## bwaIndexDir options: hg38, mm10, bosTau8, bosTau9, canFam3
## bwaIndex options: hg38.main.fa, mm10.fa, bosTau8.fa, bosTau9.main.fa, canFam3.fa

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=8:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set constant variables
genomeChrFile=/Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.chrom.sizes
numThreads=4
nonChrM=$(cat ${genomeChrFile} | awk '{print $1}' | grep -v chrM | tr '\n' ' ')
inDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/b_trimmed_fastq
outDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/c_bwa_bam/bwa2
bwaIndexDir=/scratch/Shares/chuong/genomes/hg38
bwaIndex=hg38.main.fa

# Load modules
module load bwa
module load samtools
module load bedtools

# Define query files
queries=($(ls ${inDir}/*fq.gz | xargs -n 1 basename | sed --expression='s/_1_bbduk.fq.gz//g' | sed --expression='s/_2_bbduk.fq.gz//g' | uniq))

# Run bwa
pwd; hostname; date

echo "Aligning reads to the genome & outputting sam files"
echo "bwa version: "$(bwa)
echo "Samtools version: "$(samtools --version)
echo "Bedtools version: "$(bedtools --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bwa mem \
-t ${numThreads} \
${bwaIndexDir}/${bwaIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_1_bbduk.fq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_2_bbduk.fq.gz \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

# Generate flagstat files from unfiltered sam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from unfiltered sam...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
> ${outDir}/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_flagstat.txt

# Convert sam files to bam files, filtering out low quality (<10) and unmapped (through '-F 4' option)
echo $(date +"[%b %d %H:%M:%S] Converting sam to unfiltered bam & sorting")

samtools view \
-@ ${numThreads} \
-Sb \
-q 10 \
-F 4 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam \
| samtools sort \
-@ ${numThreads} \
- \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index unfiltered sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing unfiltered, sorted bams...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Filter out chrM reads & sort
echo $(date +"[%b %d %H:%M:%S] Removing chrM reads and sorting...")

samtools view \
-b ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
${nonChrM} \
| samtools sort \
-@ ${numThreads} \
- \
> ${outDir}/filtered/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam

# Index filtered sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing filtered, sorted bams...")

samtools index \
${outDir}/filtered/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
> ${outDir}/filtered/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam.bai

# Generate flagstat files from unfiltered, sorted bam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from filtered bam...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/filtered/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
> ${outDir}/filtered/flagstat/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_flagstat.txt

# Determine percentage of total reads that mapped to mitochondrial genes
echo $(date +"[%b %d %H:%M:%S] Calculating percentage of total reads mapping to mitochondrial genes")

chrMreads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam chrM`
totalReads=`samtools view -c ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam`
fractionMreads=`echo "100 * ${chrMreads} / ${totalReads}" | bc -l`
touch ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${queries[$SLURM_ARRAY_TASK_ID]} >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${totalReads} 'total mapped reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${chrMreads} 'mitochondrial reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${fractionMreads} 'percentage of mitochondrial reads from total mapped reads' >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt

# Remove intermediate files
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sam

echo $(date +"[%b %d %H:%M:%S] Done")

### Explanation of bwa options
### '-t <int' - specifies number of threads to use

### Explanation of samtools options
### '-@ <int>' - number of threads
### '-Sb' - convert sam to bam
### '-q <int>' - quality threshold (filter out reads with quality below int)
### '-F 4' - exclude unmapped reads
### '-c' - print only the count of matching records

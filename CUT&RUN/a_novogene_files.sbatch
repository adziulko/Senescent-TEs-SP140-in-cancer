#!/bin/bash

# Date: 12 January 2022

# Define path to directory downloaded from Novogene
workingDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/a_raw_fastq/usftp21.novogene.com # directory that was downloaded from Novogene

# Define path to directory to host all fastq files (after copying)
fastqDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/a_raw_fastq2

# Check downloaded files
pwd; hostname; date

echo "Checking downloaded files..."

cd ${workingDir}/

md5sum -c ${workingDir}/MD5.txt > ${workingDir}/md5_check.txt

# Copy downloaded files to fastqDir
echo $(date +"[%b %d %H:%M:%S] Copying downloaded files to fastqDir")

cd ${workingDir}/raw_data

for i in $(ls -d -- */)
do
	cd ${workingDir}/raw_data/${i}/
	echo "Copying files for sample: "$(ls *.fq.gz | sed 's/_1.fq.gz//g' | sed 's/_2.fq.gz//g' | uniq)
	cp *_1.fq.gz ${fastqDir}/
	cp *_2.fq.gz ${fastqDir}/
	cd ${workingDir}/raw_data/
done

# Check downloaded files after copying
echo $(date +"[%b %d %H:%M:%S] Checking downloaded files after copying...")

grep -v ".zip" ${workingDir}/MD5.txt | awk '{print $1}' > ${fastqDir}/MD5.txt.tmp1
grep -v ".zip" ${workingDir}/MD5.txt | awk '{print $2}' | xargs -n 1 basename | awk -v dir="$fastqDir" '{print dir "/" $1}' > ${fastqDir}/MD5.txt.tmp2
paste ${fastqDir}/MD5.txt.tmp1 ${fastqDir}/MD5.txt.tmp2 > ${fastqDir}/MD5.txt
rm MD5.txt.tmp*

md5sum -c ${fastqDir}/MD5.txt > ${fastqDir}/md5_check.txt

echo $(date +"[%b %d %H:%M:%S] Done")

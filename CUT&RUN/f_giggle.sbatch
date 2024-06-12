#!/bin/bash

## Script to search for repeat enrichment using giggle
## Date: 11 Oct 2023
## for h_U2_HT_CR

## Example usage:
## sbatch --array 0-x e_giggle.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --propagate=NONE

# Job name and output
#SBATCH -J giggleSearch
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Declare variables
inDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/d_macs2_hs
outDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/f_giggle
repeatIndexDir=/Shares/CL_Shared/db/giggle/hg38/repeats/indexed
repeatPath=/Shares/CL_Shared/db/giggle/hg38/repeats/sorted/
genomeSize=3209286105
suffix=narrowPeak

# define query files
queries=($(ls ${inDir}/*${suffix} | xargs -n 1 basename))

# load modules
module load samtools/1.13

# run the thing
pwd; hostname; date

echo "samtools version: "$(samtools --version)
echo "giggle version: "$(giggle)
echo "Processing peak file: "${queries[$SLURM_ARRAY_TASK_ID]}

# add '.bed' to the end of each file
# giggle requires all input files to have '.bed' as a suffix
# I add a step below to remove the '.bed' suffix
echo $(date +"[%b %d %H:%M:%S] Adding suffix...")

rename ${suffix} ${suffix}.bed ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

# bgzip each bed file
echo $(date +"[%b %d %H:%M:%S] Zipping files...")

bgzip ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

# run giggle search
echo $(date +"[%b %d %H:%M:%S] Starting giggle search...")

giggle search \
-i ${repeatIndexDir} \
-q ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed.gz \
-s -g ${genomeSize} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp

# remove full index name from output
echo $(date +"[%b %d %H:%M:%S] Removing full repeat name...")

awk '(NR>1){print $0}' ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp \
| sed 's/.bed.gz//g' \
| sed "s|$repeatPath||g" \
| sed "s|sorted/||g" \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats

# sort output by descending giggle score
# note: only use unsorted output as input into makeGiggleMatrix.sbatch!
echo $(date +"[%b %d %H:%M:%S] Sorting by descending combo score...")

sort -k8,8rn \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.sorted

# remove temporary file
echo $(date +"[%b %d %H:%M:%S] Removing intermediate file...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.giggleStats.tmp

# decompress bgzip'd bed files
echo $(date +"[%b %d %H:%M:%S] Decompressing bgzip'd bed files...")

bgzip -d ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed.gz

# remove '.bed' suffix added on earlier
echo $(date +"[%b %d %H:%M:%S] Adding suffix...")

rename ${suffix}.bed ${suffix} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed

echo $(date +"[%b %d %H:%M:%S] Done")

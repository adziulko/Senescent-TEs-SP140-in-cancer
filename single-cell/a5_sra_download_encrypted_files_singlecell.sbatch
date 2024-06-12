#!/bin/bash

## Example usage:
##sbatch --array 0-x sra_download_encrypted_files.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J sraDownload
#SBATCH --output=/Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH --error=/Users/%u/slurm_err/slurm_%A_%a.err


## set variables
#sra_file=Sra2.txt
outDir=/scratch/Users/addz6536/q_single_cell_ovary_lily/rna_seq/a4
#inDir=/scratch/Users/addz6536/q_single_cell_ovary_lily

# set query files
#queries=($(cat $inDir/$sra_file | grep -v "#" | awk '{print $1}'))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

#note: need to specify --include-technical and --split-files
# in order to download all three files for single cell RNAseq
# also, don't need to module load sra, since the newest version (sra3.0.0) is in my ~/.bash_profile PATH

fasterq-dump -e 8 -O $outDir --include-technical --split-files -v -v \
-t /scratch/Users/addz6536/aaa_fasterq-dump-temp \
--ngc /scratch/Users/addz6536/prj_31282_D37419.ngc \
SRR13914005

echo $(date +"[%b %d %H:%M:%S] Done!")

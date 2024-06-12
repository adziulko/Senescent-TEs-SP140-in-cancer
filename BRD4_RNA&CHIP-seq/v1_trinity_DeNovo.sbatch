#!/bin/bash
#SBATCH -p long
#SBATCH --time=48:00:00
#SBATCH --mem=500gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

# Job name and output
#SBATCH -J BRD4_denovo_trinity
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on 'hostname'
echo Job started 'date'

# Define constant variables
trinityFiles=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/a_raw_fastq/trinity/trinity_samples.txt
outDir=/scratch/Users/addz6536/c_BRD4/trinity


module purge
module load python/2.7.14/numpy/1.11.0
module load salmon/0.9.1 # newly installed salmon/0.12.0 in fiji isn't compatible with trinity in miniconda as of 3/3/2019
module load bowtie

pwd; hostname; date


Trinity --trimmomatic --seqType fq --max_memory 500G --CPU 48 --samples_file ${trinityFiles} --SS_lib_type RF --output ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

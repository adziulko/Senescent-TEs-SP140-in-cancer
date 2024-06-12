#!/bin/bash
#SBATCH -p long
#SBATCH --time=28:00:00
#SBATCH --mem=200gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

# Job name and output
#SBATCH -J BRD4_GenomeGuided_trinity
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on 'hostname'
echo Job started 'date'

# Define constant variables
inFile=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/qui_sen_merge.sorted.bam
outDir=/scratch/Users/addz6536/c_BRD4/trinity/genome_guided_trinity/attempt_2_v2.2/trinity


module purge
module load python/2.7.14/numpy/1.11.0
module load salmon/0.9.1 # newly installed salmon/0.12.0 in fiji isn't compatible with trinity in miniconda as of 3/3/2019
module load bowtie

pwd; hostname; date


Trinity --genome_guided_bam ${inFile} --genome_guided_max_intron 30000 \
--seqType fq --max_memory 200G --CPU 24 --SS_lib_type RF --output ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

#!/bin/bash
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --mem=75gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

# Job name and output
#SBATCH -J BRD4_stringtie
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu


# Define constant variables
inFile=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/c_star/trinity/qui_sen_merge.sorted.bam
outDir=/Shares/CL_Shared/data/adam/c_BRD4/rna_seq/h_stringtie
ref=/Shares/CL_Shared/db/genomes/hg38/annotations/gencode.v36.annotation.gtf

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on 'hostname'
echo Job started 'date'



module load stringtie/1.3.3b


pwd; hostname; date


stringtie ${inFile} -p 12 -G ${ref} --rf -o ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

#!/bin/bash

## Script for running deeptools
## Date: 6 Dec 2023

## Example usage:
## sbatch g_hip_map.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J deeptools
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=addz6536@colorado.edu

# load modules
module load singularity

# define key variables
fileName=LTR8Bfullrepbase_JUND
bedFile=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/d_macs2_hs/LTR8B_fullrepbase_HT_H3K27ac.narrowPeaks.bed
wigFile=HT_H3K27ac_JUND_bw.txt
refDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/e_bigwig/no_control
outDir=/Shares/CL_Shared/data/adam/h_U2_HT_CR/cutrun/g_deeptools_heatmap

deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
wigs=$(cat ${refDir}/${wigFile})
samplesLabel=$(cat ${refDir}/${wigFile} | xargs -n 1 basename | sed 's/_treat_pileup.bw//g' | sed 's/.......$//')

# run the thing
pwd; hostname; date

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Compute matrix...")

# Use "computeMatrix" to generate data underlying heatmap
windowLeft=4000
windowRight=4000
binSize=10
numCPU=16

singularity exec --bind /Shares/CL_Shared ${deeptools} \
computeMatrix reference-point \
--referencePoint TSS \
--scoreFileName ${wigs} \
--regionsFileName ${bedFile} \
--beforeRegionStartLength ${windowLeft} \
--afterRegionStartLength ${windowRight} \
--binSize ${binSize} \
--missingDataAsZero \
-o ${outDir}/${fileName}.mat.gz \
-p ${numCPU}

echo $(date +"[%b %d %H:%M:%S] Plot heatmap...")

# Use "plotHeatmap" to create a png or pdf
zMin=0
yMin=0

singularity exec --bind /Shares/CL_Shared ${deeptools} \
plotHeatmap \
-m ${outDir}/${fileName}.mat.gz \
--outFileName ${outDir}/${fileName}.png \
--outFileSortedRegions ${outDir}/${fileName}.dt.bed \
--outFileNameMatrix ${outDir}/${fileName}.matrix.tab \
--sortRegions descend \
--colorMap Blues \
--zMin ${zMin} \
--yMin ${yMin} \
--regionsLabel ${fileName}.bed \
--samplesLabel ${samplesLabel} \
--plotTitle ${fileName}_heatmap \

echo $(date +"[%b %d %H:%M:%S] Done")

# Oncogene Induced Senescent TEs & a novel SP140 isoform in cancer

### Scripts and files used in:

Adam Dziulko, Holly Allen, Edward B. Chuong (2024) "An endogenous retrovirus regulates tumor-specific expression of the immune transcriptional regulator SP140".
DOI: [10.1093/hmg/ddae084](https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddae084/7673981)

### GEO:
High-throughput sequencing data (RNA-seq and CUT&RUN) have been deposited in the Gene Expression Omnibus (GEO) with the accession code [GSE256012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256012)

### Programs used:
- GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
- bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
- deepTools v3.0.1 (https://deeptools.readthedocs.io/en/develop/index.html)
- samtools v1.16.1 (http://www.htslib.org/)
- BBDuk/BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
- Trimmomatic v0.36 (http://www.usadellab.org/cms/?page=trimmomatic)
- hisat2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
- subread/featureCounts v1.6.2 (http://subread.sourceforge.net/)
- TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)
- MACS2 v2.1.1 (https://pypi.org/project/MACS2/)
- BWA v0.7.15 (https://github.com/lh3/bwa)
- DESeq2 v1.38.3 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- Recount2 ‘snapcount’ package in R (v4.2.2) (https://academic.oup.com/bioinformatics/article/34/1/114/4101942)
- fastq-dump v2.10.5 from SRA Toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
- MEME Suite v5.4.1 (https://meme-suite.org/meme/)
- MUSCLE v3.8.1551 (https://www.drive5.com/muscle/)
- singularity v3.1.1 (https://github.com/hpcng/singularity)

### Public databases:
- The Cancer Genome Atlas (TCGA), via the Genomic Data Commons (https://gdc.cancer.gov/)
- Cancer Cell Line Encyclopedia (CCLE) (https://www.ebi.ac.uk/ena/browser/view/PRJNA523380)
- RJunBase (http://www.rjunbase.org/)
- Cistrome (http://cistrome.org/)
- ENCODE (https://www.encodeproject.org/)
- JASPAR (http://jaspar.genereg.net/)
- Human Endogenous Retrovirus Database (https://herv.img.cas.cz/)
- UCSC (https://genome.ucsc.edu/)
- GENCODE Release 34 (https://www.gencodegenes.org/human/)
- Dfam v2.0 (https://dfam.org/home)
- GREAT  (http://great.stanford.edu/public/html/)

### Contents
Scripts used for analysis of Next-generation sequencing (NGS) genomic data.

This repository contains scripts for the following dataset types:

- CUT&RUN
- ChIP-seq
- RNA-seq
- ATAC-seq
- Single-cell analysis
- Cancer Cell Line Encyclopedia (CCLE) RNA-seq
- Tasdemir et al (BRD4 connects enhancer remodeling to senescence immune surveillance) reanalysis (RNA&CHIP-seq)

These scripts are designed to work on CU Boulder's computing system, Fiji, which uses Slurm as a resource manager and job scheduler.

## Overrepresented transposable element (TE) families that overlap with H3K27ac mark in Oncogene Induced Senescent (OIS) cells
Started with [Chip-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74238) (H3K27ac - active enhancers and transcription) fastq files from OIS IMR90 cells (from Tasdemir et al [(BRD4 connects enhancer remodeling to senescence immune surveillance)](https://aacrjournals.org/cancerdiscovery/article/6/6/612/5661/BRD4-Connects-Enhancer-Remodeling-to-Senescence)):

**ChIP-seq Workflow:**
1) Download fastq files from GEO
    - [a_sraDownload_multi_adam.sbatch](BRD4_RNA&CHIP-seq/a_sraDownload_multi_adam.sbatch)
2)  Trim adapters on raw fastq files using trimmomatic
    - [b_trimmomatic_multi_SE.sbatch](BRD4_RNA&CHIP-seq/b_trimmomatic_multi_SE.sbatch)
3) Take trimmed fastq files and run through BWA-MEM to map reads to hg38 genome 
    - [c_bwaMem_SE.sbatch](BRD4_RNA&CHIP-seq/c_bwaMem_SE.sbatch)
4) Take BAM files from BWA-MEM output and convert to peak files using MACS2
   - [d_macs2_multi_adam.sbatch](BRD4_RNA&CHIP-seq/d_macs2_multi_adam.sbatch)
5) Take peak files from MACS2 output and run through colacalization analysis (GIGGLE) to output list of overrepresented TEs that overlap with H3K27ac in OIS
   - [e_giggle.sbatch](BRD4_RNA&CHIP-seq/e_giggle.sbatch)
6) Input TE giggle results into R, then make volcano plot
   - [chip_deseq_cleanedup.R](BRD4_RNA&CHIP-seq/chip_deseq_cleanedup.R)
  
## OIS specific TE loci within enhancer acting distance (~150kb) of differentially regulated genes
Started with OIS TE list from above & [RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74324) fastq files from OIS IMR90 cells (from Tasdemir et al [(BRD4 connects enhancer remodeling to senescence immune surveillance)](https://aacrjournals.org/cancerdiscovery/article/6/6/612/5661/BRD4-Connects-Enhancer-Remodeling-to-Senescence)):

**RNA-seq Workflow:**
1) Download fastq files from GEO
    - [m_sraDownload_multi_rna_adam.sbatch](BRD4_RNA&CHIP-seq/m_sraDownload_multi_rna_adam.sbatch)
2)  Trim adapters on raw fastq files using trimmomatic
    - [n_trimmomatic_multi_PE.sbatch](BRD4_RNA&CHIP-seq/n_trimmomatic_multi_PE.sbatch)
3) Take trimmed fastq files and run through STAR to map reads to hg38 genome 
    - [o_STAR_map.sbatch](BRD4_RNA&CHIP-seq/o_STAR_map.sbatch)
4) Take BAM files from STAR output and convert to bigwig file (compresses the file/map and makes a more readable format to look at genome map on UCSC genome browser)
   - [u_deeptools_bam_to_bigwig.sbatch](BRD4_RNA&CHIP-seq/u_deeptools_bam_to_bigwig.sbatch)
5) Take trimmed fastq files and run through Salmon to quantify the expression of transcripts
    - [q_salmon_PE.sbatch](BRD4_RNA&CHIP-seq/q_salmon_PE.sbatch)
6) Input Salmon quantification files into R, then run DeSeq2 for differential expression (tximport version to convert transcripts to gene level) 
    - [DESeq2_c_BRD4_salmon_rna_cleanedup.R](BRD4_RNA&CHIP-seq/DESeq2_c_BRD4_salmon_rna_cleanedup.R)

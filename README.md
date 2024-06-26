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

## Overrepresented transposable element (TE) families in Oncogene Induced Senescent (OIS) cells
Started with [Chip-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74238) & [RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74324) fastq files from OIS IMR90 cells (from Tasdemir et al [(BRD4 connects enhancer remodeling to senescence immune surveillance)](https://aacrjournals.org/cancerdiscovery/article/6/6/612/5661/BRD4-Connects-Enhancer-Remodeling-to-Senescence)):



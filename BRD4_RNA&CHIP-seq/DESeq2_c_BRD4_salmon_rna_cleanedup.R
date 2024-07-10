#good link explaining DESeq2: 
#http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#p-values-and-adjusted-p-values

#clear global environment
rm(list = ls())

#Packages
install.packages(c("tidyverse", "scales", "WriteXLS", "BiocManager", "hexbin"))
install.packages("plotly")
install.packages("htmlwidgets")
install.packages("flexdashboard")
BiocManager::install("DESeq2")
library(DESeq2)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(flexdashboard) #for markdown (saving multiple html in one place)
library(plyr) #for binding tables with different column names



setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/c_BRD4/rna_seq_rscript/")


###read in table
pro_br1_quant.sf <- read.table("pro_br1_quant.sf", sep="\t", header = TRUE)
pro_br2_quant.sf <- read.table("pro_br2_quant.sf", sep="\t", header = TRUE)
pro_br3_quant.sf <- read.table("pro_br3_quant.sf", sep="\t", header = TRUE)
qui_br1_quant.sf <- read.table("qui_br1_quant.sf", sep="\t", header = TRUE)
qui_br2_quant.sf <- read.table("qui_br2_quant.sf", sep="\t", header = TRUE)
qui_br3_quant.sf <- read.table("qui_br3_quant.sf", sep="\t", header = TRUE)
sen_br1_quant.sf <- read.table("sen_br1_quant.sf", sep="\t", header = TRUE)
sen_br2_quant.sf <- read.table("sen_br2_quant.sf", sep="\t", header = TRUE)
sen_br3_quant.sf <- read.table("sen_br3_quant.sf", sep="\t", header = TRUE)

###clean transcript Name column
pro_br1_quant.sf <- pro_br1_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
pro_br2_quant.sf <- pro_br2_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
pro_br3_quant.sf <- pro_br3_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
qui_br1_quant.sf <- qui_br1_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
qui_br2_quant.sf <- qui_br2_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
qui_br3_quant.sf <- qui_br3_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
sen_br1_quant.sf <- sen_br1_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
sen_br2_quant.sf <- sen_br2_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
sen_br3_quant.sf <- sen_br3_quant.sf %>% separate(Name, sep = "\\|", into = "transcript_e", remove = TRUE)
#quant_gene_sf <- read.table("pro_br1_quant.genes.sf", sep="\t", header = TRUE)

###remove unncessary col
pro_br1_quant.sf = pro_br1_quant.sf[c(-2, -3, -4)]
colnames(pro_br1_quant.sf) = c('transcript_e', 'pro_br1')
pro_br2_quant.sf = pro_br2_quant.sf[c(-2, -3, -4)]
colnames(pro_br2_quant.sf) = c('transcript_e', 'pro_br2')
pro_br3_quant.sf = pro_br3_quant.sf[c(-2, -3, -4)]
colnames(pro_br3_quant.sf) = c('transcript_e', 'pro_br3')
qui_br1_quant.sf = qui_br1_quant.sf[c(-2, -3, -4)]
colnames(qui_br1_quant.sf) = c('transcript_e', 'qui_br1')
qui_br2_quant.sf = qui_br2_quant.sf[c(-2, -3, -4)]
colnames(qui_br2_quant.sf) = c('transcript_e', 'qui_br2')
qui_br3_quant.sf = qui_br3_quant.sf[c(-2, -3, -4)]
colnames(qui_br3_quant.sf) = c('transcript_e', 'qui_br3')
sen_br1_quant.sf = sen_br1_quant.sf[c(-2, -3, -4)]
colnames(sen_br1_quant.sf) = c('transcript_e', 'sen_br1')
sen_br2_quant.sf = sen_br2_quant.sf[c(-2, -3, -4)]
colnames(sen_br2_quant.sf) = c('transcript_e', 'sen_br2')
sen_br3_quant.sf = sen_br3_quant.sf[c(-2, -3, -4)]
colnames(sen_br3_quant.sf) = c('transcript_e', 'sen_br3')

###merge
pro_vs_sen = merge(sen_br1_quant.sf, sen_br2_quant.sf, by = "transcript_e")
pro_vs_sen = merge(pro_vs_sen, sen_br3_quant.sf, by = "transcript_e")
pro_vs_sen = merge(pro_vs_sen, pro_br1_quant.sf, by = "transcript_e")
pro_vs_sen = merge(pro_vs_sen, pro_br2_quant.sf, by = "transcript_e")
pro_vs_sen = merge(pro_vs_sen, pro_br3_quant.sf, by = "transcript_e")

qui_vs_sen = merge(sen_br1_quant.sf, sen_br2_quant.sf, by = "transcript_e")
qui_vs_sen = merge(qui_vs_sen, sen_br3_quant.sf, by = "transcript_e")
qui_vs_sen = merge(qui_vs_sen, qui_br1_quant.sf, by = "transcript_e")
qui_vs_sen = merge(qui_vs_sen, qui_br2_quant.sf, by = "transcript_e")
qui_vs_sen = merge(qui_vs_sen, qui_br3_quant.sf, by = "transcript_e")

pro_vs_sen = pro_vs_sen %>% remove_rownames %>% column_to_rownames(var="transcript_e")
qui_vs_sen = qui_vs_sen %>% remove_rownames %>% column_to_rownames(var="transcript_e")


###import ensemble table (will need later for merging with deseq output)
ensemble3 = read.table("ensemble_transcript_gene_update.txt", sep=",", header = TRUE)
colnames(ensemble3) = c('Gene', 'transcript_name', 'gene_name')



###Convert the countdata table into matrix format.
pro_vs_sen <- as.matrix(pro_vs_sen)
storage.mode(pro_vs_sen) = "integer"

qui_vs_sen <- as.matrix(qui_vs_sen)
storage.mode(qui_vs_sen) = "integer"


### Create a "conditions" vector, with each entry corresponding to your samples based on the sample order in the ''coldata'' table.
condition2 <- c("senescence", "senescence", "senescence", 
               "quiescence", "quiescence", "quiescence")

condition <- c("senescence", "senescence", "senescence", 
               "proliferate", "proliferate", "proliferate")


coldata <- data.frame(row.names=colnames(pro_vs_sen), condition)

coldata2 <- data.frame(row.names=colnames(qui_vs_sen), condition2)


### dds
dds <- DESeqDataSetFromMatrix( countData = pro_vs_sen,
                                   colData = coldata,
                                   design = ~ condition )

dds2 <- DESeqDataSetFromMatrix( countData = qui_vs_sen,
                               colData = coldata2,
                               design = ~ condition2 )


### Set the reference level
dds$condition <- relevel( dds$condition, ref = "proliferate" )
dds_ds <- DESeq(dds)
resultsNames(dds_ds)
res_sen_vs_pro <- results(dds_ds, contrast=c("condition", "senescence", "proliferate"))

dds2$condition2 <- relevel( dds2$condition2, ref = "quiescence" )
dds_ds2 <- DESeq(dds2)
resultsNames(dds_ds2)
res_sen_vs_qui <- results(dds_ds2, contrast=c("condition2", "senescence", "quiescence"))


### Then, filter out all "N/A" entries using na.omit() and sort your res table by ascending p-value. 
### List the number of entries that have an adjusted p-value less than 0.05.
res_sen_vs_qui <- na.omit(res_sen_vs_qui)
res_sen_vs_qui <- res_sen_vs_qui[order(res_sen_vs_qui$padj), ]
table(res_sen_vs_qui$padj<0.05)

res_sen_vs_pro <- na.omit(res_sen_vs_pro)
res_sen_vs_pro <- res_sen_vs_pro[order(res_sen_vs_pro$padj), ]
table(res_sen_vs_pro$padj<0.05)

### Optional: Make an MA plot plotting FA (-10 to 10) on the y-axis and mean of normalized counts on the x-axis
#plotMA(res_sen_vs_pro, ylim=c(-10,10))

### Merge your results table with the normalized count data table. 
### To clarify, this adds additional columns to your results table to reflect the normalized 
### counts for your samples - it does not change the values in your results table in any way.
resdata_sen_vs_qui <- merge(as.data.frame(res_sen_vs_qui), 
                                   as.data.frame(counts(dds_ds2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(resdata_sen_vs_qui)[1] <- "Gene"

resdata_sen_vs_pro <- merge(as.data.frame(res_sen_vs_pro), 
                            as.data.frame(counts(dds_ds, normalized = TRUE)), by = "row.names", sort = FALSE)
names(resdata_sen_vs_pro)[1] <- "Gene"

### calculate number of hits left over after filtering padj (false discovery rate, so padj < .05 is a false discovery rate of 5%)
sum(resdata_sen_vs_qui$padj < 0.005, na.rm=TRUE )

sum(resdata_sen_vs_pro$padj < 0.005, na.rm=TRUE )

### Filter out entires that have an adjusted p-value greater than 0.05.
resdata_sen_vs_qui_padj <- subset(resdata_sen_vs_qui, padj <= 0.05)

resdata_sen_vs_pro_padj <- subset(resdata_sen_vs_pro, padj <= 0.05)

### Separate FC>0 and FC<0 entries into their own objects and sort by FC.
resdata_sen_vs_qui_padj_ISGs <- subset(resdata_sen_vs_qui_padj, log2FoldChange >= 0)
resdata_sen_vs_qui_padj_ISGs <- resdata_sen_vs_qui_padj_ISGs[order(-resdata_sen_vs_qui_padj_ISGs$log2FoldChange), ]
resdata_sen_vs_qui_padj_IDGs <- subset(resdata_sen_vs_qui_padj, log2FoldChange <= 0)
resdata_sen_vs_qui_padj_IDGs <- resdata_sen_vs_qui_padj_IDGs[order(resdata_sen_vs_qui_padj_IDGs$log2FoldChange), ]

resdata_sen_vs_pro_padj_ISGs <- subset(resdata_sen_vs_pro_padj, log2FoldChange >= 0)
resdata_sen_vs_pro_padj_ISGs <- resdata_sen_vs_pro_padj_ISGs[order(-resdata_sen_vs_pro_padj_ISGs$log2FoldChange), ]
resdata_sen_vs_pro_padj_IDGs <- subset(resdata_sen_vs_pro_padj, log2FoldChange <= 0)
resdata_sen_vs_pro_padj_IDGs <- resdata_sen_vs_pro_padj_IDGs[order(resdata_sen_vs_pro_padj_IDGs$log2FoldChange), ]

### diff orders
order_resdata_sen_vs_qui_padj_IDGs <- resdata_sen_vs_qui_padj_IDGs[order(resdata_sen_vs_qui_padj_IDGs$padj), ]
order_resdata_sen_vs_qui_padj_ISGs <- resdata_sen_vs_qui_padj_ISGs[order(resdata_sen_vs_qui_padj_ISGs$padj), ]
foldchange_resdata_sen_vs_qui_padj_ISGs <- resdata_sen_vs_qui_padj_ISGs[order(-resdata_sen_vs_qui_padj_ISGs$log2FoldChange), ]
pval_resdata_sen_vs_qui_padj_ISGs <- resdata_sen_vs_qui_padj_ISGs[order(resdata_sen_vs_qui_padj_ISGs$pvalue), ]

order_resdata_sen_vs_pro_padj_IDGs <- resdata_sen_vs_pro_padj_IDGs[order(resdata_sen_vs_pro_padj_IDGs$padj), ]
order_resdata_sen_vs_pro_padj_ISGs <- resdata_sen_vs_pro_padj_ISGs[order(resdata_sen_vs_pro_padj_ISGs$padj), ]
foldchange_resdata_sen_vs_pro_padj_ISGs <- resdata_sen_vs_pro_padj_ISGs[order(-resdata_sen_vs_pro_padj_ISGs$log2FoldChange), ]
pval_resdata_sen_vs_pro_padj_ISGs <- resdata_sen_vs_pro_padj_ISGs[order(resdata_sen_vs_pro_padj_ISGs$pvalue), ]


### merge ensemble to get gene name and isoform
ensem_resdata_sen_vs_qui_padj_ISGs_named = merge(ensemble3, resdata_sen_vs_qui_padj_ISGs, by = "Gene")
ensem_resdata_sen_vs_qui_padj_IDGs_named = merge(ensemble3, resdata_sen_vs_qui_padj_IDGs, by = "Gene")

ensem_resdata_sen_vs_pro_padj_ISGs_named = merge(ensemble3, resdata_sen_vs_pro_padj_ISGs, by = "Gene")
ensem_resdata_sen_vs_pro_padj_IDGs_named = merge(ensemble3, resdata_sen_vs_pro_padj_IDGs, by = "Gene")

### ensem diff orders
ensem_order_resdata_sen_vs_qui_padj_IDGs <- ensem_resdata_sen_vs_qui_padj_IDGs_named[order(ensem_resdata_sen_vs_qui_padj_IDGs_named$padj), ]
ensem_order_resdata_sen_vs_qui_padj_ISGs <- ensem_resdata_sen_vs_qui_padj_ISGs_named[order(ensem_resdata_sen_vs_qui_padj_ISGs_named$padj), ]
ensem_foldchange_resdata_sen_vs_qui_padj_ISGs <- ensem_resdata_sen_vs_qui_padj_ISGs_named[order(-ensem_resdata_sen_vs_qui_padj_ISGs_named$log2FoldChange), ]


ensem_order_resdata_sen_vs_pro_padj_IDGs <- ensem_resdata_sen_vs_pro_padj_IDGs_named[order(ensem_resdata_sen_vs_pro_padj_IDGs_named$padj), ]
ensem_order_resdata_sen_vs_pro_padj_ISGs <- ensem_resdata_sen_vs_pro_padj_ISGs_named[order(ensem_resdata_sen_vs_pro_padj_ISGs_named$padj), ]
ensem_foldchange_resdata_sen_vs_pro_padj_ISGs <- ensem_resdata_sen_vs_pro_padj_ISGs_named[order(-ensem_resdata_sen_vs_pro_padj_ISGs_named$log2FoldChange), ]

### Merge all tables together so that I can have one file containing all info
ensem_order_resdata_sen_vs_qui_padj_IDGs$data.source <- paste0("down_senVSqui", ensem_order_resdata_sen_vs_qui_padj_IDGs$data.source)
ensem_order_resdata_sen_vs_qui_padj_ISGs$data.source <- paste0("up_senVSqui", ensem_order_resdata_sen_vs_qui_padj_ISGs$data.source)

ensem_order_resdata_sen_vs_pro_padj_IDGs$data.source <- paste0("down_senVSpro", ensem_order_resdata_sen_vs_pro_padj_IDGs$data.source)
ensem_order_resdata_sen_vs_pro_padj_ISGs$data.source <- paste0("up_senVSpro", ensem_order_resdata_sen_vs_pro_padj_ISGs$data.source)

all_significant_up_down_genes_senVSpro.qui <- rbind.fill(ensem_order_resdata_sen_vs_qui_padj_IDGs,
                                                    ensem_order_resdata_sen_vs_qui_padj_ISGs,
                                                    ensem_order_resdata_sen_vs_pro_padj_IDGs,
                                                    ensem_order_resdata_sen_vs_pro_padj_ISGs)









##The following is to create list of senescent specific upregulated genes within enhancer distance of senescent specific TE family, LTR8B
##########
### read in curated list / add column names / remove duplicate genes
setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Papers/curated_files/")
sen_age_inflam_TF_epig_apop_genes <- read.table("aging_senescence_genes2.txt", sep="\t", header = FALSE) #list of genes taken from multiple sources associated with inflammation/aging/epigenome
setwd("/Users/adamdziulko/Documents/Chuong Lab/Projects/c_BRD4/rna_seq_rscript/")
### set column names
colnames(sen_age_inflam_TF_epig_apop_genes) = c('gene_name', 'gene_type')
### remove duplicate genes in list
dedup_sen_age_inflam_TF_epig_apop_genes = distinct(sen_age_inflam_TF_epig_apop_genes, gene_name, .keep_all = TRUE)

### make table of upregulated/downregulated sen genes vs qui and pro
all_sig_up_genes_senVSpro.qui <- rbind.fill(ensem_order_resdata_sen_vs_qui_padj_ISGs,
                                            ensem_order_resdata_sen_vs_pro_padj_ISGs)
#write.table(all_sig_up_genes_senVSpro.qui, file="all_sig_up_genes_senVSpro.qui.txt", row.names = FALSE, sep ="\t")
all_sig_up_genes_senVSpro.qui = read.table("all_sig_up_genes_senVSpro.qui.txt", sep="\t", header = TRUE)
tximp_all_sig_up_genes_senVSpro.qui = read.table("all_sig_up_GENES2_senVSpro.qui.txt", sep="\t", header = TRUE)



### write in table of transcription start sites of TEs/clean
tss_qui = read.table("all_TEs_senPeaks_tss_sen_vs_qui_up_wb.bed", sep="\t", header = FALSE)
# tss_pro = read.table("all_TEs_senPeaks_tss_sen_vs_pro_up_wb.bed", sep="\t", header = FALSE) #adds nothing new

### remove unnecessary columns
tss_qui = tss_qui[c(-1, -2, -3, -4, -5, -6, -8, -9, -10, -12, -13)]
# tss_pro = tss_pro[c(-1, -2, -3, -4, -5, -6, -8, -9, -10, -12, -13)] #adds nothing new

tss_qui = tss_qui[tss_qui$V11== "LTR8B" | tss_qui$V11== "LTR8"| tss_qui$V11== "LTR10A"| tss_qui$V11== "LTR18A"| tss_qui$V11== "MER44B",]
##tss_pro = tss_pro[tss_pro$V11== "LTR8B" | tss_pro$V11== "LTR8"| tss_pro$V11== "LTR10A"| tss_pro$V11== "LTR18A"| tss_pro$V11== "MER44B",] #add nothing new

tss_qui = distinct(tss_qui, V7, .keep_all = TRUE)
tss_qui = tss_qui[-2]
colnames(tss_qui) = 'gene_name'


### read in table of genes within xkb of senescence specific H3k27ac peaks
#LTR8B_250kb_qui_senescent_gene <- read.table("250KBwindowTSS_LTR8B_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
#colnames(LTR8B_250kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
#LTR8B_250kb_qui_senescent_gene = LTR8B_250kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR8B_150kb_qui_senescent_gene <- read.table("150KBwindowTSS_LTR8B_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR8B_150kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR8B_150kb_qui_senescent_gene = LTR8B_150kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

#LTR8_250kb_qui_senescent_gene <- read.table("250KBwindowTSS_LTR8_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE, quote = "")
#colnames(LTR8_250kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
#LTR8_250kb_qui_senescent_gene = LTR8_250kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR8_150kb_qui_senescent_gene <- read.table("150KBwindowTSS_LTR8_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR8_150kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR8_150kb_qui_senescent_gene = LTR8_150kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR18A_150kb_qui_senescent_gene <- read.table("150KBwindowTSS_LTR18A_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR18A_150kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR18A_150kb_qui_senescent_gene = LTR18A_150kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR10A_150kb_qui_senescent_gene <- read.table("150KBwindowTSS_LTR10A_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR10A_150kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR10A_150kb_qui_senescent_gene = LTR10A_150kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

MER44B_150kb_qui_senescent_gene <- read.table("150KBwindowTSS_MER44B_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(MER44B_150kb_qui_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
MER44B_150kb_qui_senescent_gene = MER44B_150kb_qui_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]


#LTR8B_250kb_pro_senescent_gene <- read.table("250KBwindowTSS_LTR8B_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
#colnames(LTR8B_250kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
#LTR8B_250kb_pro_senescent_gene = LTR8B_250kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR8B_150kb_pro_senescent_gene <- read.table("150KBwindowTSS_LTR8B_senPeaks_minus_proPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR8B_150kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR8B_150kb_pro_senescent_gene = LTR8B_150kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

#LTR8_250kb_pro_senescent_gene <- read.table("250KBwindowTSS_LTR8_senPeaks_minus_quiPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE, quote = "")
#colnames(LTR8_250kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
#LTR8_250kb_pro_senescent_gene = LTR8_250kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR8_150kb_pro_senescent_gene <- read.table("150KBwindowTSS_LTR8_senPeaks_minus_proPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR8_150kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR8_150kb_pro_senescent_gene = LTR8_150kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR18A_150kb_pro_senescent_gene <- read.table("150KBwindowTSS_LTR18A_senPeaks_minus_proPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR18A_150kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR18A_150kb_pro_senescent_gene = LTR18A_150kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

LTR10A_150kb_pro_senescent_gene <- read.table("150KBwindowTSS_LTR10A_senPeaks_minus_proPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(LTR10A_150kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
LTR10A_150kb_pro_senescent_gene = LTR10A_150kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

MER44B_150kb_pro_senescent_gene <- read.table("150KBwindowTSS_MER44B_senPeaks_minus_proPeaks_H3K27ac.narrowPeaks.bed", sep="\t", header = FALSE, fill = TRUE)
colnames(MER44B_150kb_pro_senescent_gene) = c('chr_TE', 'TE_start', 'TE_end', '4', '5', '6', '7', '8', '9', 'chr_gene', 'gene_start', 'gene_end', 'ens_ID', '14', '15', 'gene_name', 'gene_descrip', 'gene_class')
MER44B_150kb_pro_senescent_gene = MER44B_150kb_pro_senescent_gene[c(-4, -5, -6, -7, -8, -9, -14, -15)]

#add in TSS genes/TEs to LTR8B
LTR8B_150kb_qui_senescent_gene = dplyr::bind_rows(LTR8B_150kb_qui_senescent_gene,tss_qui)

#dedupTRAN_all_sig_up_genes_senVSpro.qui = distinct(all_sig_up_genes_senVSpro.qui, transcript_name, .keep_all = TRUE)
#dedupGENE_all_sig_up_genes_senVSpro.qui = distinct(all_sig_up_genes_senVSpro.qui, gene_name, .keep_all = TRUE)


### merge genes within LTR8B 250kb window with upregulated genes in oncogene-induced senescent cells
#LTR8B_250kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8B_250kb_qui_senescent_gene, by = "gene_name")              
#dedup_LTR8B_250kb_qui_upregulated_genes = ddply(LTR8B_250kb_qui_upregulated_genes,.(gene_name),nrow)
#colnames(dedup_LTR8B_250kb_qui_upregulated_genes)[2] = 'num_dup'

LTR8B_150kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8B_150kb_qui_senescent_gene, by = "gene_name")
dedup_LTR8B_150kb_qui_upregulated_genes = distinct(LTR8B_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

LTR8B_150kb_pro_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8B_150kb_pro_senescent_gene, by = "gene_name")
dedup_LTR8B_150kb_pro_upregulated_genes = distinct(LTR8B_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

### tximport version
tximp_LTR8B_150kb_qui_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR8B_150kb_qui_senescent_gene, by = "gene_name")
tximp_dedup_LTR8B_150kb_qui_upregulated_genes = distinct(tximp_LTR8B_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR8B_150kb_pro_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR8B_150kb_pro_senescent_gene, by = "gene_name")
tximp_dedup_LTR8B_150kb_pro_upregulated_genes = distinct(tximp_LTR8B_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

#LTR8_250kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8_250kb_qui_senescent_gene, by = "gene_name")              
#dedup_LTR8_250kb_qui_upregulated_genes = ddply(LTR8_250kb_qui_upregulated_genes,.(gene_name),nrow)
#colnames(dedup_LTR8_250kb_qui_upregulated_genes)[2] = 'num_dup'

LTR8_150kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8_150kb_qui_senescent_gene, by = "gene_name")              
dedup_LTR8_150kb_qui_upregulated_genes = distinct(LTR8_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

LTR8_150kb_pro_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR8_150kb_pro_senescent_gene, by = "gene_name")              
dedup_LTR8_150kb_pro_upregulated_genes = distinct(LTR8_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

LTR18A_150kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR18A_150kb_qui_senescent_gene, by = "gene_name")              
dedup_LTR18A_150kb_qui_upregulated_genes = distinct(LTR18A_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

LTR18A_150kb_pro_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR18A_150kb_pro_senescent_gene, by = "gene_name")              
dedup_LTR18A_150kb_pro_upregulated_genes = distinct(LTR18A_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

LTR10A_150kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR10A_150kb_qui_senescent_gene, by = "gene_name")              
dedup_LTR10A_150kb_qui_upregulated_genes = distinct(LTR10A_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

LTR10A_150kb_pro_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, LTR10A_150kb_pro_senescent_gene, by = "gene_name")              
dedup_LTR10A_150kb_pro_upregulated_genes = distinct(LTR10A_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

MER44B_150kb_qui_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, MER44B_150kb_qui_senescent_gene, by = "gene_name")              
dedup_MER44B_150kb_qui_upregulated_genes = distinct(MER44B_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

MER44B_150kb_pro_upregulated_genes = merge(all_sig_up_genes_senVSpro.qui, MER44B_150kb_pro_senescent_gene, by = "gene_name")              
dedup_MER44B_150kb_pro_upregulated_genes = distinct(MER44B_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

### tximport version
tximp_LTR8_150kb_qui_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR8_150kb_qui_senescent_gene, by = "gene_name")              
tximp_dedup_LTR8_150kb_qui_upregulated_genes = distinct(tximp_LTR8_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR8_150kb_pro_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR8_150kb_pro_senescent_gene, by = "gene_name")              
tximp_dedup_LTR8_150kb_pro_upregulated_genes = distinct(tximp_LTR8_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR18A_150kb_qui_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR18A_150kb_qui_senescent_gene, by = "gene_name")              
tximp_dedup_LTR18A_150kb_qui_upregulated_genes = distinct(tximp_LTR18A_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR18A_150kb_pro_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR18A_150kb_pro_senescent_gene, by = "gene_name")              
tximp_dedup_LTR18A_150kb_pro_upregulated_genes = distinct(tximp_LTR18A_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR10A_150kb_qui_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR10A_150kb_qui_senescent_gene, by = "gene_name")              
tximp_dedup_LTR10A_150kb_qui_upregulated_genes = distinct(tximp_LTR10A_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_LTR10A_150kb_pro_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, LTR10A_150kb_pro_senescent_gene, by = "gene_name")              
tximp_dedup_LTR10A_150kb_pro_upregulated_genes = distinct(tximp_LTR10A_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_MER44B_150kb_qui_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, MER44B_150kb_qui_senescent_gene, by = "gene_name")              
tximp_dedup_MER44B_150kb_qui_upregulated_genes = distinct(tximp_MER44B_150kb_qui_upregulated_genes, gene_name, .keep_all = TRUE)

tximp_MER44B_150kb_pro_upregulated_genes = merge(tximp_all_sig_up_genes_senVSpro.qui, MER44B_150kb_pro_senescent_gene, by = "gene_name")              
tximp_dedup_MER44B_150kb_pro_upregulated_genes = distinct(tximp_MER44B_150kb_pro_upregulated_genes, gene_name, .keep_all = TRUE)

### testing out different ways to quantify dups for LTR8B_150KB ##wont be using due to split TEs and isoform calling
#TRAN_LTR8B_150kb_qui_upregulated_genes = merge(dedupTRAN_all_sig_up_genes_senVSpro.qui, LTR8B_150kb_qui_senescent_gene, by = "gene_name")              
#dedupTRAN_LTR8B_150kb_qui_upregulated_genes = ddply(TRAN_LTR8B_150kb_qui_upregulated_genes,.(gene_name),nrow)
#colnames(dedupTRAN_LTR8B_150kb_qui_upregulated_genes)[2] = 'num_dup'

#GENE_LTR8B_150kb_qui_upregulated_genes = merge(dedupGENE_all_sig_up_genes_senVSpro.qui, LTR8B_150kb_qui_senescent_gene, by = "gene_name")              
#dedupGENE_LTR8B_150kb_qui_upregulated_genes = ddply(GENE_LTR8B_150kb_qui_upregulated_genes,.(gene_name),nrow)
#colnames(dedupGENE_LTR8B_150kb_qui_upregulated_genes)[2] = 'num_dup'


### merge curated list to upregulated genes withing x kb to TE / dedup curated list
#curated_LTR8B_250kb_qui_upregulated_genes = merge(dedup_LTR8B_250kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR8B_150kb_qui_upregulated_genes = merge(dedup_LTR8B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR8B_150kb_qui_upregulated_genes = curated_LTR8B_150kb_qui_upregulated_genes[ , c(1, 5, 9, 28)] 
curated_LTR8B_150kb_qui_upregulated_genes$TE = "LTR8B"

curated_LTR8B_150kb_pro_upregulated_genes = merge(dedup_LTR8B_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR8B_150kb_pro_upregulated_genes = curated_LTR8B_150kb_pro_upregulated_genes[ , c(1, 5, 9, 28)] 
curated_LTR8B_150kb_pro_upregulated_genes$TE = "LTR8B"

r_curated_LTR8B_150kb_pro.qui_upregulated_genes = rbind(curated_LTR8B_150kb_qui_upregulated_genes, curated_LTR8B_150kb_pro_upregulated_genes) 
r_curated_LTR8B_150kb_pro.qui_upregulated_genes = distinct(r_curated_LTR8B_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


### tximport version
tximp_curated_LTR8B_150kb_qui_upregulated_genes = merge(tximp_dedup_LTR8B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR8B_150kb_qui_upregulated_genes = tximp_curated_LTR8B_150kb_qui_upregulated_genes[ , c(1, 4, 27)] 
tximp_curated_LTR8B_150kb_qui_upregulated_genes$TE = "LTR8B"

tximp_curated_LTR8B_150kb_pro_upregulated_genes = merge(tximp_dedup_LTR8B_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR8B_150kb_pro_upregulated_genes = tximp_curated_LTR8B_150kb_pro_upregulated_genes[ , c(1, 4, 27)] 
tximp_curated_LTR8B_150kb_pro_upregulated_genes$TE = "LTR8B"

tximp_r_curated_LTR8B_150kb_pro.qui_upregulated_genes = rbind(tximp_curated_LTR8B_150kb_qui_upregulated_genes, tximp_curated_LTR8B_150kb_pro_upregulated_genes) 
tximp_r_curated_LTR8B_150kb_pro.qui_upregulated_genes = distinct(tximp_r_curated_LTR8B_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)

### Below makes files counting number of TE loci per genes / how many transcripts each gene is affecting ##not using it due to split TEs and isoform calling
#TRAN_curated_LTR8B_150kb_qui_upregulated_genes = merge(dedupTRAN_LTR8B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
#GENE_curated_LTR8B_150kb_qui_upregulated_genes = merge(dedupGENE_LTR8B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")


#curated_LTR8_250kb_qui_upregulated_genes = merge(dedup_LTR8_250kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")

curated_LTR8_150kb_qui_upregulated_genes = merge(dedup_LTR8_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR8_150kb_qui_upregulated_genes = curated_LTR8_150kb_qui_upregulated_genes[ , c(1, 5, 9, 28)]
curated_LTR8_150kb_qui_upregulated_genes$TE = "LTR8"

curated_LTR8_150kb_pro_upregulated_genes = merge(dedup_LTR8_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR8_150kb_pro_upregulated_genes = curated_LTR8_150kb_pro_upregulated_genes[ , c(1, 5, 9, 28)]
curated_LTR8_150kb_pro_upregulated_genes$TE = "LTR8"

r_curated_LTR8_150kb_pro.qui_upregulated_genes = rbind(curated_LTR8_150kb_qui_upregulated_genes, curated_LTR8_150kb_pro_upregulated_genes) 
r_curated_LTR8_150kb_pro.qui_upregulated_genes = distinct(r_curated_LTR8_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


curated_LTR18A_150kb_qui_upregulated_genes = merge(dedup_LTR18A_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR18A_150kb_qui_upregulated_genes = curated_LTR18A_150kb_qui_upregulated_genes[ , c(1, 5, 28)]
curated_LTR18A_150kb_qui_upregulated_genes$TE = "LTR18A"

curated_LTR18A_150kb_pro_upregulated_genes = merge(dedup_LTR18A_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR18A_150kb_pro_upregulated_genes = curated_LTR18A_150kb_pro_upregulated_genes[ , c(1, 5, 28)]
curated_LTR18A_150kb_pro_upregulated_genes$TE = "LTR18A"

r_curated_LTR18A_150kb_pro.qui_upregulated_genes = rbind(curated_LTR18A_150kb_qui_upregulated_genes, curated_LTR18A_150kb_pro_upregulated_genes) 
r_curated_LTR18A_150kb_pro.qui_upregulated_genes = distinct(r_curated_LTR18A_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


curated_LTR10A_150kb_qui_upregulated_genes = merge(dedup_LTR10A_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR10A_150kb_qui_upregulated_genes = curated_LTR10A_150kb_qui_upregulated_genes[ , c(1, 5, 28)]
curated_LTR10A_150kb_qui_upregulated_genes$TE = "LTR10A"

curated_LTR10A_150kb_pro_upregulated_genes = merge(dedup_LTR10A_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_LTR10A_150kb_pro_upregulated_genes = curated_LTR10A_150kb_pro_upregulated_genes[ , c(1, 5, 28)]
curated_LTR10A_150kb_pro_upregulated_genes$TE = "LTR10A"

r_curated_LTR10A_150kb_pro.qui_upregulated_genes = rbind(curated_LTR10A_150kb_qui_upregulated_genes, curated_LTR10A_150kb_pro_upregulated_genes) 
r_curated_LTR10A_150kb_pro.qui_upregulated_genes = distinct(r_curated_LTR10A_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


curated_MER44B_150kb_qui_upregulated_genes = merge(dedup_MER44B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_MER44B_150kb_qui_upregulated_genes = curated_MER44B_150kb_qui_upregulated_genes[ , c(1, 5, 28)]
curated_MER44B_150kb_qui_upregulated_genes$TE = "MER44B"

curated_MER44B_150kb_pro_upregulated_genes = merge(dedup_MER44B_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
curated_MER44B_150kb_pro_upregulated_genes = curated_MER44B_150kb_pro_upregulated_genes[ , c(1, 5, 28)]
curated_MER44B_150kb_pro_upregulated_genes$TE = "MER44B"

r_curated_MER44B_150kb_pro.qui_upregulated_genes = rbind(curated_MER44B_150kb_qui_upregulated_genes, curated_MER44B_150kb_pro_upregulated_genes) 
r_curated_MER44B_150kb_pro.qui_upregulated_genes = distinct(r_curated_MER44B_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


### tximport vetrsion
tximp_curated_LTR8_150kb_qui_upregulated_genes = merge(tximp_dedup_LTR8_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR8_150kb_qui_upregulated_genes = tximp_curated_LTR8_150kb_qui_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_LTR8_150kb_qui_upregulated_genes$TE = "LTR8"

tximp_curated_LTR8_150kb_pro_upregulated_genes = merge(tximp_dedup_LTR8_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR8_150kb_pro_upregulated_genes = tximp_curated_LTR8_150kb_pro_upregulated_genes[ , c(1,  4, 27)]
tximp_curated_LTR8_150kb_pro_upregulated_genes$TE = "LTR8"

tximp_r_curated_LTR8_150kb_pro.qui_upregulated_genes = rbind(tximp_curated_LTR8_150kb_qui_upregulated_genes, tximp_curated_LTR8_150kb_pro_upregulated_genes) 
tximp_r_curated_LTR8_150kb_pro.qui_upregulated_genes = distinct(tximp_r_curated_LTR8_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


tximp_curated_LTR18A_150kb_qui_upregulated_genes = merge(tximp_dedup_LTR18A_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR18A_150kb_qui_upregulated_genes = tximp_curated_LTR18A_150kb_qui_upregulated_genes[ , c(1,  4, 27)]
tximp_curated_LTR18A_150kb_qui_upregulated_genes$TE = "LTR18A"

tximp_curated_LTR18A_150kb_pro_upregulated_genes = merge(tximp_dedup_LTR18A_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR18A_150kb_pro_upregulated_genes = tximp_curated_LTR18A_150kb_pro_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_LTR18A_150kb_pro_upregulated_genes$TE = "LTR18A"

tximp_r_curated_LTR18A_150kb_pro.qui_upregulated_genes = rbind(tximp_curated_LTR18A_150kb_qui_upregulated_genes, tximp_curated_LTR18A_150kb_pro_upregulated_genes) 
tximp_r_curated_LTR18A_150kb_pro.qui_upregulated_genes = distinct(tximp_r_curated_LTR18A_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


tximp_curated_LTR10A_150kb_qui_upregulated_genes = merge(tximp_dedup_LTR10A_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR10A_150kb_qui_upregulated_genes = tximp_curated_LTR10A_150kb_qui_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_LTR10A_150kb_qui_upregulated_genes$TE = "LTR10A"

tximp_curated_LTR10A_150kb_pro_upregulated_genes = merge(tximp_dedup_LTR10A_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_LTR10A_150kb_pro_upregulated_genes = tximp_curated_LTR10A_150kb_pro_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_LTR10A_150kb_pro_upregulated_genes$TE = "LTR10A"

tximp_r_curated_LTR10A_150kb_pro.qui_upregulated_genes = rbind(tximp_curated_LTR10A_150kb_qui_upregulated_genes, tximp_curated_LTR10A_150kb_pro_upregulated_genes) 
tximp_r_curated_LTR10A_150kb_pro.qui_upregulated_genes = distinct(tximp_r_curated_LTR10A_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


tximp_curated_MER44B_150kb_qui_upregulated_genes = merge(tximp_dedup_MER44B_150kb_qui_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_MER44B_150kb_qui_upregulated_genes = tximp_curated_MER44B_150kb_qui_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_MER44B_150kb_qui_upregulated_genes$TE = "MER44B"

tximp_curated_MER44B_150kb_pro_upregulated_genes = merge(tximp_dedup_MER44B_150kb_pro_upregulated_genes, dedup_sen_age_inflam_TF_epig_apop_genes, by = "gene_name")
tximp_curated_MER44B_150kb_pro_upregulated_genes = tximp_curated_MER44B_150kb_pro_upregulated_genes[ , c(1, 4, 27)]
tximp_curated_MER44B_150kb_pro_upregulated_genes$TE = "MER44B"

tximp_r_curated_MER44B_150kb_pro.qui_upregulated_genes = rbind(tximp_curated_MER44B_150kb_qui_upregulated_genes, tximp_curated_MER44B_150kb_pro_upregulated_genes) 
tximp_r_curated_MER44B_150kb_pro.qui_upregulated_genes = distinct(tximp_r_curated_MER44B_150kb_pro.qui_upregulated_genes, gene_name, .keep_all = TRUE)


### concatene all r_curated files
all_r = rbind(r_curated_LTR8B_150kb_pro.qui_upregulated_genes, r_curated_LTR8_150kb_pro.qui_upregulated_genes, r_curated_LTR18A_150kb_pro.qui_upregulated_genes,
              r_curated_LTR10A_150kb_pro.qui_upregulated_genes, r_curated_MER44B_150kb_pro.qui_upregulated_genes)

dist_all_r = distinct(all_r, gene_name, .keep_all = TRUE)

### concatene LTR8/B r_curated files
all_LTR8.B_r = rbind(r_curated_LTR8B_150kb_pro.qui_upregulated_genes, r_curated_LTR8_150kb_pro.qui_upregulated_genes)

dist_all_LTR8.B_r = distinct(all_LTR8.B_r, gene_name, .keep_all = TRUE)
### remove unnecssary columns
dist_all_LTR8.B_r = dist_all_LTR8.B_r[c(-3)]
colnames(dist_all_LTR8.B_r) = c('Gene name', 'Log2 fold_change', 'TE')
dist_all_LTR8.B_r$`Log2 fold_change` = round(dist_all_LTR8.B_r$`Log2 fold_change`, digits = 4)

name.width <- max(sapply(names(dist_all_LTR8.B_r), nchar))
dist_all_LTR8.B_r$`Log2 fold_change` = format(dist_all_LTR8.B_r$`Log2 fold_change`, width = name.width, justify = "centre")

#colnames(pro_br1_quant.sf) = c('transcript_e', 'pro_br1')


### tximport version
tximp_all_r = rbind(tximp_r_curated_LTR8B_150kb_pro.qui_upregulated_genes, tximp_r_curated_LTR8_150kb_pro.qui_upregulated_genes, tximp_r_curated_LTR18A_150kb_pro.qui_upregulated_genes,
                    tximp_r_curated_LTR10A_150kb_pro.qui_upregulated_genes, tximp_r_curated_MER44B_150kb_pro.qui_upregulated_genes)

tximp_dist_all_r = distinct(tximp_all_r, gene_name, .keep_all = TRUE)

### write table
#write.table(r_curated_LTR8B_150kb_pro.qui_upregulated_genes, file="r_curated_LTR8B_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")
#write.table(r_curated_LTR8_150kb_pro.qui_upregulated_genes, file="r_curated_LTR8_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")
#write.table(r_curated_LTR18A_150kb_pro.qui_upregulated_genes, file="r_curated_LTR18A_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")
#write.table(r_curated_LTR10A_150kb_pro.qui_upregulated_genes, file="r_curated_LTR10A_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")
#write.table(r_curated_MER44B_150kb_pro.qui_upregulated_genes, file="r_curated_MER44B_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")
#write.table(all_r, file="r_curated_LTR8B.LTR8.LTR18A.LTR10A.MER44B_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")

#write.table(tximp_dist_all_r, file="tximp_dist_all_r_LTR8B.LTR8.LTR18A.LTR10A.MER44B_150kb_pro.qui_upregulated_genes", row.names = FALSE, sep ="\t")

tximp_dist_all_r = read.table("tximp_dist_all_r_LTR8B.LTR8.LTR18A.LTR10A.MER44B_150kb_pro.qui_upregulated_genes", sep="\t", header = TRUE)



setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/c_BRD4/rna_seq_rscript/")
all_r = read.table("r_curated_LTR8B.LTR8.LTR18A.LTR10A.MER44B_150kb_pro.qui_upregulated_genes", sep="\t", header = TRUE)
all_r = all_r[c(-2,-3)]
all_r = all_r[c(1,8,14,23,48,58),]
colnames(all_r) = c('Gene', 'ERV')
all_r = all_r[c(3,2,1,4,5,6),]




### write table
#write.table(ensem_order_resdata_sen_vs_qui_padj_IDGs, file="sen_vs_qui_padj_down.txt", row.names = FALSE, sep ="\t")
#write.table(ensem_order_resdata_sen_vs_qui_padj_ISGs, file="sen_vs_qui_padj_up.txt", row.names = FALSE, sep ="\t")
#write.table(ensem_foldchange_resdata_sen_vs_qui_padj_ISGs, file="sen_vs_qui_foldchange_up.txt", row.names = FALSE, sep ="\t")

#write.table(ensem_order_resdata_sen_vs_pro_padj_IDGs, file="sen_vs_pro_padj_down.txt", row.names = FALSE, sep ="\t")
#write.table(ensem_order_resdata_sen_vs_pro_padj_ISGs, file="sen_vs_pro_padj_up.txt", row.names = FALSE, sep ="\t")
#write.table(ensem_foldchange_resdata_sen_vs_pro_padj_ISGs, file="sen_vs_pro_foldchange_up.txt", row.names = FALSE, sep ="\t")

#write.table(all_significant_up_down_genes_senVSpro.qui, 
 #           file="all_significant_up_down_genes_senVSpro.qui.txt", row.names = FALSE, sep ="\t")



### how many genes lost due to ensemble merge?
nrow(ensem_resdata_sen_vs_qui_padj_ISGs_named)
nrow(resdata_sen_vs_qui_padj_ISGs)
nrow(ensem_resdata_sen_vs_qui_padj_IDGs_named)
nrow(resdata_sen_vs_qui_padj_IDGs)
nrow(ensem_resdata_sen_vs_pro_padj_ISGs_named)
nrow(resdata_sen_vs_pro_padj_ISGs)
nrow(ensem_resdata_sen_vs_pro_padj_IDGs_named)
nrow(resdata_sen_vs_pro_padj_IDGs)

### remove unnecessary columns
#ensem_resdata_sen_vs_qui_padj_ISGs_named = ensem_resdata_sen_vs_qui_padj_ISGs_named[c(-1, -2)]
#ensem_resdata_sen_vs_qui_padj_IDGs_named = ensem_resdata_sen_vs_qui_padj_IDGs_named[c(-1, -2)]

#ensem_resdata_sen_vs_pro_padj_ISGs_named = ensem_resdata_sen_vs_pro_padj_ISGs_named[c(-1, -2)]
#ensem_resdata_sen_vs_pro_padj_IDGs_named = ensem_resdata_sen_vs_pro_padj_IDGs_named[c(-1, -2)]

### Create table of uniqe genes 
### import gene transcripts bed file
gene_transcripts_tss_bed <- read.table("gencode.v34.transcripts.tss.bed", sep="\t", header = FALSE)
colnames(gene_transcripts_tss_bed) = c('chr', 'bp1', 'bp2', 'Gene', '?', 'strand')
### merge ISG tables by ENS (gene)
tss_sen_vs_qui_up = merge(gene_transcripts_tss_bed, ensem_resdata_sen_vs_qui_padj_ISGs_named, by = "Gene")
tss_sen_vs_pro_up = merge(gene_transcripts_tss_bed, ensem_resdata_sen_vs_pro_padj_ISGs_named, by = "Gene")
### reorder columns to fit bed files
tss_sen_vs_qui_up <- tss_sen_vs_qui_up[, c(2, 3, 4, 5, 6, 1, 8)]
tss_sen_vs_pro_up <- tss_sen_vs_pro_up[, c(2, 3, 4, 5, 6, 1, 8)]
### order by chromosome
tss_sen_vs_qui_up = tss_sen_vs_qui_up[order(tss_sen_vs_qui_up[,1], tss_sen_vs_qui_up[,2] ),]
tss_sen_vs_pro_up = tss_sen_vs_pro_up[order(tss_sen_vs_pro_up[,1], tss_sen_vs_pro_up[,2] ),]
### write files
#write.table(tss_sen_vs_qui_up, file="tss_sen_vs_qui_up.bed", row.names = FALSE, sep ="\t", col.names = FALSE)
#write.table(tss_sen_vs_pro_up, file="tss_sen_vs_pro_up.bed", row.names = FALSE, sep ="\t", col.names = FALSE)

### uniqe genes upregulated 
length(unique(ensem_resdata_sen_vs_qui_padj_ISGs_named$gene_name)) #org = 5588
length(unique(ensem_resdata_sen_vs_pro_padj_ISGs_named$gene_name)) #org = 4497
### number of overlapping of pro and qui
pro_qui_up_overlap = merge(ensem_resdata_sen_vs_qui_padj_ISGs_named, ensem_resdata_sen_vs_pro_padj_ISGs_named, by = "Gene")
length(unique(pro_qui_up_overlap$gene_name.x)) #org = 3166


### order the new tables by fold change
#gene_resdata_sen_vs_qui_padj_ISGs_named <- gene_resdata_sen_vs_qui_padj_ISGs_named[order(-gene_resdata_sen_vs_qui_padj_ISGs_named$log2FoldChange), ]
#gene_resdata_sen_vs_qui_padj_IDGs_named <- gene_resdata_sen_vs_qui_padj_IDGs_named[order(gene_resdata_sen_vs_qui_padj_IDGs_named$log2FoldChange), ]

#order_gene_resdata_sen_vs_qui_padj_IDGs_named <- gene_resdata_sen_vs_qui_padj_IDGs_named[order(gene_resdata_sen_vs_qui_padj_IDGs_named$padj), ]
#order_gene_resdata_sen_vs_qui_padj_ISGs_named <- gene_resdata_sen_vs_qui_padj_ISGs_named[order(gene_resdata_sen_vs_qui_padj_ISGs_named$padj), ]

### write the named tables
#write.table(gene_resdata_sen_vs_qui_padj_ISGs_named, 
#            file="gene_resdata_sen_vs_qui_padj_increase_named.txt", row.names = FALSE, sep ="\t")
#write.table(gene_resdata_sen_vs_qui_padj_IDGs_named, 
#            file="gene_resdata_sen_vs_qui_padj_decrease_named.txt", row.names = FALSE, sep ="\t")




### whole code for making barplots comparing senescence vs quiesence counts
######
#create mean of biological replicates
order_TE_resdata_sen_vs_qui_padj_ISGs$quiesence_mean = rowMeans(order_TE_resdata_sen_vs_qui_padj_ISGs[, c(11,12,13)])
order_TE_resdata_sen_vs_qui_padj_ISGs$senescence_mean = rowMeans(order_TE_resdata_sen_vs_qui_padj_ISGs[, c(8,9,10)])
### rearrange the data to fit ggplot (used tidyr: https://uc-r.github.io/tidyr)
mean_barplot_long_data = order_TE_resdata_sen_vs_qui_padj_ISGs %>% gather(cell_stage, mean_count, quiesence_mean, senescence_mean)
### make TE and cell_type factors to fit ggplot
mean_barplot_long_data$TE = as.factor(mean_barplot_long_data$TE)
mean_barplot_long_data$mean_type = as.factor(mean_barplot_long_data$cell_stage)
### calling ggplot
mean_graph = ggplot(mean_barplot_long_data, aes(x=TE, y=mean_count, fill=cell_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
mean_graph
######


### p-value/subset (first 25) code for making barplots comparing senescence vs quiesence counts
######
### create mean of biological replicates
pval_TE_resdata_sen_vs_qui_padj_ISGs$quiesence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(11,12,13)])
pval_TE_resdata_sen_vs_qui_padj_ISGs$senescence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(8,9,10)])
### subset data
subset_pval_barplot_TE_data = subset(pval_TE_resdata_sen_vs_qui_padj_ISGs, senescence_mean > 50)
#subset_pval_barplot_TE_data = head(subset_pval_barplot_TE_data,25)
### rearrange the data to fit ggplot (used tidyr: https://uc-r.github.io/tidyr)
pval_mean_barplot_long_data = subset_pval_barplot_TE_data %>% gather(cell_stage, mean_count, quiesence_mean, senescence_mean)
### make TE and cell_type factors to fit ggplot
pval_mean_barplot_long_data$TE = as.factor(pval_mean_barplot_long_data$TE)
pval_mean_barplot_long_data$mean_type = as.factor(pval_mean_barplot_long_data$cell_stage)
### reorder the long data by pvalue
pval_mean_barplot_long_data <- pval_mean_barplot_long_data[order(pval_mean_barplot_long_data$pvalue), ]
### calling ggplot
pval_mean_graph = ggplot(pval_mean_barplot_long_data, aes(x=reorder(TE, pvalue), y=mean_count, fill=cell_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("Transposable Elements & Satellites") +
  ylab("Mean Count") +
  ### Add the text 
  geom_text(aes(label = round(mean_count, 0), y = mean_count + 20), position = position_dodge(0.9), size = 2) + 
  ### add arrow
  annotate("segment", x = 2, xend = 18, y = 1500, yend = 1500, colour = "black", size=1, arrow=arrow()) + 
  ### add "increasing p-value"
  geom_text(x=10, y=1540, label="Increasing p-value")
pval_mean_graph
### save graph
pdf(file = "pval_all_TE_senescenceVSquiesence_bargraph", width = 16, height = 12)
pval_mean_graph
dev.off()
######



### p-value/subset (25-75) code for making barplots comparing senescence vs quiesence counts
######
### create mean of biological replicates
pval_TE_resdata_sen_vs_qui_padj_ISGs$quiesence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(11,12,13)])
pval_TE_resdata_sen_vs_qui_padj_ISGs$senescence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(8,9,10)])
### subset data
subset_25to75_pval_barplot_TE_data = subset(pval_TE_resdata_sen_vs_qui_padj_ISGs, senescence_mean > 50)
subset_25to75_pval_barplot_TE_data = subset_25to75_pval_barplot_TE_data[26:75,]
### rearrange the data to fit ggplot (used tidyr: https://uc-r.github.io/tidyr)
pval_25to75_mean_barplot_long_data = subset_25to75_pval_barplot_TE_data %>% gather(cell_stage, mean_count, quiesence_mean, senescence_mean)
### make TE and cell_type factors to fit ggplot
pval_25to75_mean_barplot_long_data$TE = as.factor(pval_25to75_mean_barplot_long_data$TE)
pval_25to75_mean_barplot_long_data$mean_type = as.factor(pval_25to75_mean_barplot_long_data$cell_stage)
### reorder the long data by pvalue
pval_25to75_mean_barplot_long_data <- pval_25to75_mean_barplot_long_data[order(pval_25to75_mean_barplot_long_data$pvalue), ]
### calling ggplot
pval_25to75_mean_graph = ggplot(pval_25to75_mean_barplot_long_data, aes(x=reorder(TE, pvalue), y=mean_count, fill=cell_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("Transposable Elements & Satellites") +
  ylab("Mean Count") +
  ### Add the text 
  geom_text(aes(label = round(mean_count, 0), y = mean_count + 100), position = position_dodge(0.9), size = 2) + 
  ### add arrow
  annotate("segment", x = 2, xend = 18, y = 15750, yend = 15750, colour = "black", size=1, arrow=arrow()) + 
  ### add "increasing p-value"
  geom_text(x=10, y=16000, label="Increasing p-value")
pval_25to75_mean_graph
### save graph
pdf(file = "pval_25to75_TE_senescenceVSquiesence_bargraph", width = 16, height = 12)
pval_25to75_mean_graph
dev.off()
######


### p-value ordered/subset (top 100) then foldchange ordered (1-25) code for making barplots comparing senescence vs quiesence counts
######
### create mean of biological replicates
pval_TE_resdata_sen_vs_qui_padj_ISGs$quiesence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(11,12,13)])
pval_TE_resdata_sen_vs_qui_padj_ISGs$senescence_mean = rowMeans(pval_TE_resdata_sen_vs_qui_padj_ISGs[, c(8,9,10)])
### subset data
pval.subset_top100_barplot_TE_data = subset(pval_TE_resdata_sen_vs_qui_padj_ISGs, senescence_mean > 50)
pval.subset_top100_barplot_TE_data = pval.subset_top100_barplot_TE_data[1:100,]
foldchange_pval.subset_top100_barplot_TE_data <- pval.subset_top100_barplot_TE_data[order(-pval.subset_top100_barplot_TE_data$log2FoldChange), ]
### rearrange the data to fit ggplot (used tidyr: https://uc-r.github.io/tidyr)
foldchange_pval_top25_mean_barplot_long_data = foldchange_pval.subset_top100_barplot_TE_data %>% gather(cell_stage, mean_count, quiesence_mean, senescence_mean)
### make TE and cell_type factors to fit ggplot
foldchange_pval_top25_mean_barplot_long_data$TE = as.factor(foldchange_pval_top25_mean_barplot_long_data$TE)
foldchange_pval_top25_mean_barplot_long_data$mean_type = as.factor(foldchange_pval_top25_mean_barplot_long_data$cell_stage)
### reorder the long data by foldchange
foldchange_pval_top25_mean_barplot_long_data <- foldchange_pval_top25_mean_barplot_long_data[order(-foldchange_pval_top25_mean_barplot_long_data$log2FoldChange), ]
### calling ggplot
pval.subset_top100_barplot = ggplot(foldchange_pval_top25_mean_barplot_long_data, aes(x=reorder(TE, -log2FoldChange), y=mean_count, fill=cell_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("Transposable Elements & Satellites") +
  ylab("Mean Count") +
  ggtitle("[1]Filtered (late senescence mean > 50); [2]Subsetted (lowest 100 padj values (4.96e-156 to 1.46e-24)); [3]Ordered (Log2 Fold-Change (6 to 1.3))") +
  coord_cartesian(ylim=c(0, 17500)) +
  ### Add the text 
  geom_text(aes(label = round(mean_count, 0), y = mean_count + 100), position = position_dodge(0.9), size = 2) + 
  ### add arrow
  annotate("segment", x = 2, xend = 18, y = 15750, yend = 15750, colour = "black", size=1, arrow=arrow()) + 
  ### add "increasing p-value"
  geom_text(x=10, y=16000, label="Decreasing log2 Fold-Change") +
  geom_text(x=97, y=17500, label="64341", size = 3) +
  geom_text(x=99, y=16000, label="25341", size = 3)
pval.subset_top100_barplot
### save graph
pdf(file = "foldchangeORDER_top100pval_TE_senescenceVSquiesence_bargraph", width = 30, height = 12)
pval.subset_top100_barplot
dev.off()
######



### test code (subset data) for making barplots comparing senescence vs quiesence counts
######
### create mean of biological replicates
order_TE_resdata_sen_vs_qui_padj_ISGs$quiesence_mean = rowMeans(order_TE_resdata_sen_vs_qui_padj_ISGs[, c(11,12,13)])
order_TE_resdata_sen_vs_qui_padj_ISGs$senescence_mean = rowMeans(order_TE_resdata_sen_vs_qui_padj_ISGs[, c(8,9,10)])
### subset data
mean_barplot_subset_data = head(order_TE_resdata_sen_vs_qui_padj_ISGs,5)
### rearrange the data to fit ggplot (used tidyr: https://uc-r.github.io/tidyr)
mean_barplot_subset_long_data = mean_barplot_subset_data %>% gather(cell_stage, mean_count, quiesence_mean, senescence_mean)
### make TE and cell_type factors to fit ggplot
mean_barplot_subset_long_data$TE = as.factor(mean_barplot_subset_long_data$TE)
mean_barplot_subset_long_data$mean_type = as.factor(mean_barplot_subset_long_data$cell_stage)
### calling ggplot
mean_subset_graph = ggplot(mean_barplot_subset_long_data, aes(x=TE, y=mean_count, fill=cell_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
mean_subset_graph
#########




# below is code for making boxplots comparing senescence vs quiesence counts
#########
boxplot_test_data = head(order_TE_resdata_sen_vs_qui_padj_ISGs,5)
boxplot_long_data = boxplot_test_data %>% gather(replicate, count, senescence.br3:senescence.br1, quiesence.br3:quiesence.br1)
separate_boxplot_long_data <- boxplot_long_data %>% separate(replicate, c("cell_type", "replicate_number"), sep = "\\.")
separate_boxplot_long_data$count = as.integer(separate_boxplot_long_data$count)
###TE and cell_type columns need to be 'factors'
separate_boxplot_long_data$TE = as.factor(separate_boxplot_long_data$TE)
separate_boxplot_long_data$cell_type = as.factor(separate_boxplot_long_data$cell_type)
###ggplot graph code below
box_plot_graph = ggplot(separate_boxplot_long_data, aes(x=TE, y=count, fill=cell_type)) + 
  geom_boxplot()
box_plot_graph
#########



### code for adding TE year to table
##########
### import table
TE_age_table <- read.table("hg38.repeats.info", sep="\t", header = TRUE)
### remove simple repeats
TE_age_table = subset(TE_age_table, !grepl("Simple_repeat", TE_age_table$repFamily) )
### separate data
separate_TE_up_count <- TE_resdata_sen_vs_qui_padj_ISGs %>% separate(TE, c("repName", "Family", "Class"), sep = "\\:")
### merge
merge_separate_TE_up_count <- merge(separate_TE_up_count, TE_age_table, by = "repName", sort = FALSE )
### filter out larger pvalues (necessary because y scale too large for the very small values)
merge_separate_TE_up_count = subset(merge_separate_TE_up_count, pvalue <= 0.000005)
### log pvalue column
merge_separate_TE_up_count$log2_pvalue <- log(merge_separate_TE_up_count$pvalue, 2)
### scatter plot of repdivergence vs pvalue 
pvalue_repDiv_scatter = ggplot(merge_separate_TE_up_count, aes(x=repDiv, y=log2_pvalue, label=repName, label2=Family, color=repClass)) +
  ggtitle("Transposable Elements/Satellites Divergence by pavlue scatterplot") +
  geom_smooth(method = "lm", se = FALSE, aes(group=repClass)) +
  geom_point()
### make plot interactive
pvalue_repDiv_scatter = ggplotly(pvalue_repDiv_scatter)
### save interactive plot as html
htmlwidgets::saveWidget(as_widget(pvalue_repDiv_scatter), "pvalue_repDiv_scatter.html")
### view plot
pvalue_repDiv_scatter

### scatter plot fo repdivergence vs foldchange
foldchange_repDiv_scatter = ggplot(merge_separate_TE_up_count, aes(x=repDiv, y=log2FoldChange, label=repName, label2=Family, color=repClass)) +
  ggtitle("Transposable Elements/Satellites Divergence by log2-foldchange scatterplot") +
  geom_smooth(method = "lm", se = FALSE, aes(group=repClass)) +
  geom_point()
### make plot interactive
foldchange_repDiv_scatter = ggplotly(foldchange_repDiv_scatter)
### save interactive plot as html
htmlwidgets::saveWidget(as_widget(foldchange_repDiv_scatter), "foldchange_repDiv_scatter.html")
### view plot
foldchange_repDiv_scatter





#_______________________________________________________________________________
### Volcano plot #3 (contrast1.2_D_SPvsB6; pairwise; batched)

dist_all_LTR8.B_r = subset(dist_all_LTR8.B_r,  padj > 5e-70)

### Define cutoffs for labeling genes
### Tweak these until you're happy with how many genes are labeled
value4_ltr8s = subset(dist_all_LTR8.B_r, padj< 5e-2 & log2FoldChange > 1 & padj > 5e-70)
nrow(value4_ltr8s)
value5_ltr8s = subset(dist_all_LTR8.B_r, padj< 5e-2 & log2FoldChange < -1)
nrow(value5_ltr8s)


### Prepare your volcano plot
DEG_150kb_LTR8B_volcano = ggplot(dist_all_LTR8.B_r, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("DeSeq2; 150kb; Ltr8/B") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10 adjusted pvalue", x = "log2 fold change") +
  ### Set all dots color to grey
  geom_point(data=dist_all_LTR8.B_r, colour = "grey") + 
  ### If pvalue < 0.05, change dot color to green
  #geom_point(data=new[which(new$padj<0.05),], colour = "grey") + 
  ### If log2FC > 1, change dot color to orange
  #geom_point(data=new[which(abs(new$log2FoldChange)>1),], colour = "grey") +
  ### If both, change dot color to blue
  geom_point(data = value4_ltr8s, colour = "#A8DDB5", size = 4) +
  geom_point(data = value5_ltr8s, colour = "blue", size = 4) +
  ### Add text label for interesting outliers
  ### Add text label for interesting outliers
  geom_text_repel(data = value4_ltr8s, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 5, force = 8) +
  geom_text_repel(data = value5_ltr8s, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 5, force = 8) +
  theme(axis.text.x = element_text(face="bold", size=22),
        axis.text.y = element_text(face="bold", size=22)) +
  theme(axis.title=element_text(size=28), plot.title = element_text(size=28, face = "bold")) 
DEG_150kb_LTR8B_volcano
#dev.copy(pdf, height = 12, width = 12, pointsize=4, file.path(pngDir, "volcano_bIFNg_2h.pdf"))
#dev.off()
figure = ggarrange(DEG_150kb_LTR8B_volcano, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_DEG_150kb_LTR8B_volcano3", ".pdf"), width = 14, height = 10)
figure
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()
#_______________________________________________________________________________







#####
### Create dataset for TEs and for Genes only
res_dds <- results(dds_ds)
res_dds_Gene <- res_dds[1:67638,]
res_dds_TE <- res_dds[67639:68659,]





### Make volcano plots.
library(ggrepel)
library(dplyr)

### Adjusted p-value
data <- res_sen_vs_qui[ ,c( 2,6 ) ] %>% as.data.frame( )
data <- na.omit(data)
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes(x= log2FoldChange, y= -log10( padj ) ), color = "grey" ) +
  theme_bw( ) +
  ### Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggtitle( "Volcano plot: DE Genes/TEs - quiesence vs late senescence" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 Adj.pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.05 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=2, shape=5,stroke=1  )






num1 = ggplot(data = data, aes(x= log2FoldChange, y= -log10(padj))) +
  ### Change the 'theme' or appearance of the plot such that the background is white
  theme_bw() +
  ### Give your plot a title
  ggtitle("Volcano plot: DE Genes/TEs - quiesence vs late senescence") +
  ### Name your axes
  labs( y = "-log10 Adj.pval ", x = "log2 fold change" ) +
  ### Adjust the size and placement (hjust & vjust) of your labels
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15, hjust = 0.5, vjust = 0.5)) +
  ### Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ### Plot 'insigificant' (subset above) repeats as dark grey dots (color = "#4D4D4D") that are smaller (size = 4) and more transparent (alpha = 0.5) than significant repeats
  geom_point(data = insig, color = "#4D4D4D", size = 4, alpha = 0.5, stroke = 0) +
  ### Plot positively enriched repeats as light blue dots that are larger and more opaque than insignificant repeats
  geom_point(data = pos_enrich, color = "#4EB3D3", size = 8, stroke = 0, alpha = 0.8) +
  ### Plot negatively enriched repeats as light green dots that are larger and more opaque than insigificant repeats
  geom_point(data = neg_enrich, color = "#A8DDB5", size = 8, stroke = 0, alpha = 0.8) +
  ### Tell ggplot to make the x-axis scale from a log2(obs / exp) of -6 to 6 with breaks at intervals of 2
  ### Change this to fit your data!
  scale_x_continuous(limits = c(-6,6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  ### Label your positively enriched repeats such that the labels are (ideally) not overlapping
  ### Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  geom_text_repel(data = pos_enrich, 
                  aes(x = log2(obs / exp), y = -log10(fishers_two_tail)), 
                  label = pos_enrich$name, 
                  size = 6, 
                  force = 10,
                  nudge_y = 0.5,
                  nudge_x = 0.5,
                  point.padding = 0.25) +
  ### Label your negatively enriched repeats such that the labels are (ideally) not overlapping
  ### Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  geom_text_repel(data = neg_enrich,
                  aes(x = log2(obs / exp), y = -log10(fishers_two_tail)),
                  label = neg_enrich$name,
                  size = 6,
                  force = 15,
                  nudge_y = 0.5,
                  nudge_x = 0.5,
                  point.padding = 0.25)



#######################################
### Pvalue             
data <- res_UNvsK9M[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE Genes/TEs - K9M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=1, shape=5,stroke=1  )

data <- res_UNvsK36M[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE Genes/TEs - K36M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=1, shape=5,stroke=1  )



res_K9M_TE <- res_UNvsK9M[67639:68659,]
data <- res_K9M_TE[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE TEs - K9M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pvalue ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=2, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) > 1 ), ], color = "#f4a548", size=2, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),], color = "#862a5c", size=2, shape=5,stroke=1  ) +
  geom_text_repel(data = data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),],  
                  aes(x = log2FoldChange, y = -log10( pvalue )), 
                  label = rownames(sig.data),
                  size = 2,
                  nudge_y = -0.1,
                  #nudge_x = 0.5,
                  segment.color = NA)












#####################







### Remove regions with zero counts (this is a sanity check - there shouldn't be any regions with zero counts across all samples)
dds <- dds[ rowSums( counts( dds ) ) > 10, ]
nrow( dds)

#write.table(counts(dds), sep="\t", file="rawcounts_filtered.txt")

### Normalization to export data and make custom HeatMaps
dds_N <- estimateSizeFactors(dds)
dds_N <- estimateDispersions(dds_N)
Ndds <- counts(dds_N, normalized=TRUE)
#write.csv(as.data.frame(Ndds), file="Normalized.values.csv")

### Run DESeq2
dds_ds <- DESeq(dds)

### DE genes: genes that have differences in the effect of treatment between sexes
resultsNames( dds_ds)
### Plot dispersion estimates
plotDispEsts( dds_ds)

res_sen_vs_qui <- results(dds, contrast=c("condition", "IFNG_2h", "untreated"))







### Create dataset for TEs and for Genes only
res_dds <- results(dds_ds)
res_dds_Gene <- res_dds[1:67638,]
res_dds_TE <- res_dds[67639:68659,]

res_UNvsK9M <- results(dds_ds, contrast=c("treat","K9M","un"))
write.table(as.data.frame(res_UNvsK9M), sep="\t", file="Contrast_K9MvsUN.txt")

res_UNvsK36M <- results(dds_ds, contrast=c("treat","K36M","un"))
write.table(as.data.frame(res_UNvsK36M), sep="\t", file="Contrast_K36MvsUN.txt")

res_dds_TE_K9<- results(dds_ds[67639:68659,], contrast=c("treat","K9M","un"))
res_dds_TE_K36<- results(dds_ds[67639:68659,], contrast=c("treat","K36M","un"))
### Adj pvalue
res_dds_TE_K9M <- na.omit( res_UNvsK9M[67639:68659,] )
res_dds_TE_K9M <- res_dds_TE_K9M[ order( res_dds_TE_K9M$padj ), ]
table( res_dds_TE$padj<0.05 )
#13

res_dds_TE_K36M <- na.omit( res_UNvsK36M[67639:68659,] )
res_dds_TE_K36M <- res_dds_TE_K36M[ order( res_dds_TE_K36M$padj ), ]
table( res_dds_TE_K36M$padj<0.05 )
#81

### pvalue
res_dds_TE_K9M <- na.omit( res_UNvsK9M[67639:68659,] )
res_dds_TE_K9M <- res_dds_TE_K9M[ order( res_dds_TE_K9M$pvalue ), ]
table( res_dds_TE$pvalue<0.05 )
#47

res_dds_TE_K36M <- na.omit( res_UNvsK36M[67639:68659,] )
res_dds_TE_K36M <- res_dds_TE_K36M[ order( res_dds_TE_K36M$pvalue ), ]
table( res_dds_TE_K36M$pvalue<0.05 )
#166

### Make MA plots.
resultsNames(dds_ds)
plotMA( res_UNvsK9M , ylim = c( -10, 10 ), alpha = 0.05 )
plotMA( res_UNvsK36M , ylim = c( -10, 10 ), alpha = 0.05 )


plotMA( res_dds_TE_K9, ylim = c( -10, 10 ), alpha = 0.05 )
plotMA( res_dds_TE_K36, ylim = c( -10, 10 ), alpha = 0.05 )

### Merge the results tables with the normalized count data.
### Extract normalized counts
dds_ds_tot <- as.data.frame( counts(dds_ds), normalized = TRUE )
colnames( dds_ds_tot )

### Merge results tables with normalized counts tables
resdata_K9vsUN <- merge( as.data.frame( res_UNvsK9M  ), as.data.frame(dds_ds_tot), by = "row.names", sort = FALSE )
names( resdata_K9vsUN)[ 1 ] <- "Gene/TE"
#write.table( resdata_K9vsUN, "K9vsUN_Normalized.txt", row.names = FALSE, sep = "\t" )

resdata_K36vsUN <- merge( as.data.frame( res_UNvsK36M ), as.data.frame(dds_ds_tot  ), by = "row.names", sort = FALSE )
names( resdata_K36vsUN)[ 1 ] <- "Gene/TE"
#write.table( resdata_K36vsUN, "K36vsUN_Normalized.txt", row.names = FALSE, sep = "\t" )


### Filter out entries with a padj value greater than 0.05 and log2FoldChange less than 1 (down) or greater than 1 (up)
resdata_K9vsUN_padj <- subset( resdata_K9vsUN, padj <= 0.05 )
resdata_K9vsUN_padj_log2FC <-subset(resdata_K9vsUN_padj, abs(log2FoldChange)>= 1)
resdata_K36vsUN_padj <- subset( resdata_K36vsUN, padj <= 0.05 )
resdata_K36vsUN_padj_log2FC <-subset(resdata_K36vsUN_padj, abs(log2FoldChange)>= 1)

### Filter out entries with a pvalue greater than 0.05 and log2FoldChange less than 1 (down) or greater than 1 (up)
resdata_K9vsUN_pvalue <- subset( resdata_K9vsUN, pvalue <= 0.05 )
resdata_K9vsUN_pvalue_log2FC <-subset(resdata_K9vsUN_pvalue, abs(log2FoldChange)>= 1)
resdata_K36vsUN_pvalue <- subset( resdata_K36vsUN, pvalue <= 0.05 )
resdata_K36vsUN_pvalue_log2FC <-subset(resdata_K36vsUN_pvalue, abs(log2FoldChange)>= 1)

resdata_LRT_TE_log2FC <-subset(resdata_TE, abs(log2FoldChange)>= 2)
resdata_LRT_TE_pvalue <-subset(resdata_TE, pvalue <= 0.05)

### Write DESeq2 output txt files
#write.table( resdata_K9vsUN_padj, "K9vsUN.DE.padj.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K9vsUN_padj_log2FC, "K9vsUN.DE.padj.Log2FC.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K36vsUN_padj, "K36vsUN.DE.padj.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K36vsUN_padj_log2FC, "K36vsUN.DE.padj.Log2FC.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K9vsUN_pvalue, "K9vsUN.DE.pval.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K9vsUN_pvalue_log2FC, "K9vsUN.DE.pval.Log2FC.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K36vsUN_pvalue, "K36vsUN.DE.pval.txt", quote = FALSE, row.names = FALSE, sep = "\t" )
#write.table( resdata_K36vsUN_pvalue_log2FC, "K36vsUN.DE.pval.Log2FC.txt", quote = FALSE, row.names = FALSE, sep = "\t" )

### Make volcano plots.
library(ggrepel)
library(dplyr)

### Adjusted p-value
data <- res_UNvsK9M[ ,c( 2,6 ) ] %>% as.data.frame( )
data <- na.omit(data)
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE Genes/TEs - K9M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 Adj.pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=2, shape=5,stroke=1  )

### Pvalue             
data <- res_UNvsK9M[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE Genes/TEs - K9M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=1, shape=5,stroke=1  )

data <- res_UNvsK36M[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) >= 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE Genes/TEs - K36M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pval ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=1, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) >=1), ], color = "#f4a548", size=1, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = sig.data, color = "#862a5c", size=1, shape=5,stroke=1  )



res_K9M_TE <- res_UNvsK9M[67639:68659,]
data <- res_K9M_TE[ ,c( 2,5 ) ] %>% as.data.frame( )
pvalue <- data[ ,c( 2 ) ]
log2FoldChange <- data[ ,c( 1 ) ]
sig.data<- data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),]
ggplot( data, aes( log2FoldChange, -log10( pvalue ) ), color = "grey" ) +
  theme_bw( ) +
  ggtitle( "Volcano plot: DE TEs - K9M vs wt hematopoietic SC" ) +
  theme( plot.title = element_text( hjust = 0.5 ) ) +
  labs( y = "-log10 pvalue ", x = "log2 fold change" ) +
  geom_point( color = "#4D4D4D", size = 1, alpha = 0.5, shape=4) +
  geom_point( data = data[ which( pvalue < 0.01 ), ], color = "#00909e", size=2, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) > 1 ), ], color = "#f4a548", size=2, alpha = 0.5, shape=1, stroke=1 ) +
  geom_point( data = data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),], color = "#862a5c", size=2, shape=5,stroke=1  ) +
  geom_text_repel(data = data[ which( abs( log2FoldChange ) > 1 & pvalue < 0.01 ),],  
                  aes(x = log2FoldChange, y = -log10( pvalue )), 
                  label = rownames(sig.data),
                  size = 2,
                  nudge_y = -0.1,
                  #nudge_x = 0.5,
                  segment.color = NA)

### Load required packages
library( "vsn" )
library( "dplyr" )
library( "ggplot2" )
library( "hexbin" )
library( "pheatmap" )
library( "RColorBrewer" )

### Transform count data via VST
vsd_val <- vst( dds, blind = FALSE )
vsd_data <- (assay(vsd_val))
write.table(as.data.frame(vsd_data), sep = "\t",file="vsd_All.txt")

vsd_G <- vsd_val[1:67638,]
vsd_TE <- vsd_val[67639:68659,]
### Generate scatterplot using VST
df_val <- bind_rows( as_data_frame( assay( vsd_val )[, 1:2] ) %>% mutate( transformation = "vst" ) )
colnames( df_val )[1:2] <- c( "x", "y" )
ggplot( df_val, aes( x = x, y = y ) ) + 
  geom_hex( bins = 80 ) +
  coord_fixed( ) + 
  facet_grid( . ~ transformation )

df_G <- bind_rows( as_data_frame( assay( vsd_G )[, 1:2] ) %>% mutate( transformation = "vst" ) )
colnames( df_G )[1:2] <- c( "x", "y" )
ggplot( df_G, aes( x = x, y = y ) ) + 
  geom_hex( bins = 80 ) +
  coord_fixed( ) + 
  facet_grid( . ~ transformation )

df_TE <- bind_rows( as_data_frame( assay( vsd_TE )[, 1:2] ) %>% mutate( transformation = "vst" ) )
colnames( df_TE )[1:2] <- c( "x", "y" )
ggplot( df_TE, aes( x = x, y = y ) ) + 
  geom_hex( bins = 80 ) +
  coord_fixed( ) + 
  facet_grid( . ~ transformation )

### Calculate distance
sampleDists_val <- dist( t( assay( vsd_val ) ) )
#sampleDists_val
sampleDistMatrix_val <- as.matrix( sampleDists_val )
#rownames( sampleDistMatrix_OCR ) <- paste( "BL3.1", vsd_OCR$condition_OCR, vsd_OCR$replicate_OCR, sep = "_" )
#colnames( sampleDistMatrix_OCR ) <- paste( "BL3.1", vsd_OCR$condition_OCR, vsd_OCR$replicate_OCR, sep = "_" )
colors <- colorRampPalette(  rev( brewer.pal( 9, "Blues" ) ) )( 255 )
pheatmap( sampleDistMatrix_val,
          clustering_distance_rows = sampleDists_val,
          clustering_distance_cols = sampleDists_val,
          col = colors )

### Genes
sampleDists_val_G <- dist( t( assay( vsd_G ) ) )
sampleDistMatrix_val_G <- as.matrix( sampleDists_val_G )
colors <- colorRampPalette(  rev( brewer.pal( 9, "Blues" ) ) )( 255 )
pheatmap( sampleDistMatrix_val_G,
          clustering_distance_rows = sampleDists_val_G,
          clustering_distance_cols = sampleDists_val_G,
          col = colors )

### TEs
sampleDists_val_TE <- dist( t( assay( vsd_TE ) ) )
sampleDistMatrix_val_TE <- as.matrix( sampleDists_val_TE )
colors <- colorRampPalette(  rev( brewer.pal( 9, "Blues" ) ) )( 255 )
pheatmap( sampleDistMatrix_val_TE,
          clustering_distance_rows = sampleDists_val_TE,
          clustering_distance_cols = sampleDists_val_TE,
          col = colors )

### Visualize VST distances in a PCA plot.
### Load required packages
library( "affy" )

### Make PCA plots
pcaData_valG <- plotPCA( vsd_G, intgroup = c( "treat" ), returnData = TRUE )
percentVar_val <- round( 100 * attr( pcaData_valG, "percentVar" ) )
ggplot( pcaData_valG, aes( x = PC1, y = PC2, color="treat") ) +
  geom_point( size = 3 ) +
  xlab( paste0( "PC1: ", percentVar_val[1], "% variance" ) ) +
  ylab( paste0( "PC2: ", percentVar_val[2], "% variance" ) ) +
  coord_fixed( ) + 
  labs( title = "Hematopoietic SC", subtitle = "PCA - VST transformed" )

plotPCA( vsd_TE, intgroup = c( "treat" ))

### Generate heatmaps showing hierarchal clustering of top 100 or 20 genes with the highest variance.

library( "genefilter" )
library( "data.table" )


### Get top 100 genes showing greatest variability & make heatmap
topVarRegions_A <- head( order( rowVars( assay( vsd_val ) ), decreasing = TRUE ), 100 )
mat_A  <- assay( vsd_val )[ topVarRegions_A, ]
mat_A  <- mat_A - rowMeans( mat_A )
colData( vsd_val )
anno_A <- as.data.frame( colData( vsd_val )[, "treat"] )
rownames( anno_A) <- colnames( mat_A )
colnames( anno_A )[1] <- "condition group"
pheatmap( mat_A, annotation_col = anno_A, fontsize_row = 7, fontsize_col = 10 )


topVarRegions_TE <- head( order( rowVars( assay( vsd_TE ) ), decreasing = TRUE ), 20 )
mat_TE  <- assay( vsd_TE )[ topVarRegions_TE, ]
mat_TE  <- mat_TE - rowMeans( mat_TE )
anno_TE <- as.data.frame( colData( vsd_TE )[, "treat"] )
rownames( anno_TE) <- colnames( mat_TE )
colnames( anno_TE )[1] <- "condition group"
pheatmap( mat_TE, annotation_col = anno_TE, fontsize_row = 7, fontsize_col = 10 )


### Get top 20 genes showing greatest variability & make heatmap
topVarRegions_A <- head( order( rowVars( assay( vsd_val ) ), decreasing = TRUE ), 20 )
mat_A  <- assay( vsd_val )[ topVarRegions_A, ]
mat_A  <- mat_A - rowMeans( mat_A )
anno_A <- as.data.frame( colData( vsd_val )[, "treat"] )
rownames( anno_A) <- colnames( mat_A )
colnames( anno_A )[1] <- "condition group"
pheatmap( mat_A, annotation_col = anno_A, fontsize_row = 9, fontsize_col = 10 )
###################





#### Define cutoffs for coloring dots
value1 = subset(new, padj < 0.05 & abs(log2FoldChange) < 1)
value2 = subset(new, padj > 0.05 & abs(log2FoldChange) > 1)
value3 = subset(new, padj < 0.05 & abs(log2FoldChange) > 1)

### Assess how many genes you have with padj < 0.05 and abs(log2FoldChange) > 1
nrow(value3)
head(value3)

### Define cutoffs for labeling genes
### Tweak these until you're happy with how many genes are labeled
value4 = subset(new, padj< 5e-5 & log2FoldChange > 1.3)
nrow(value4)
value5 = subset(new, padj< 5e-5 & log2FoldChange < -1.3)
nrow(value5)

### Prepare your volcano plot
num1 = ggplot(new, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("Differential expressed genes in HT1080 \n SP140 KD vs control") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10 adjusted pvalue", x = "log2 fold change") +
  ## Set all dots color to grey
  geom_point(data=new, colour = "grey") + 
  ## If pvalue < 0.05, change dot color to green
  #geom_point(data=new[which(new$padj<0.05),], colour = "grey") + 
  ## If log2FC > 1, change dot color to orange
  #geom_point(data=new[which(abs(new$log2FoldChange)>1),], colour = "grey") +
  ## If both, change dot color to blue
  geom_point(data = value4, colour = "red", size = 6) +
  geom_point(data = value5, colour = "blue", size = 6) +
  # Add text label for interesting outliers
  geom_text_repel(data = value4, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 6, force = 5) +
  geom_text_repel(data = value5, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 6, force = 5) +
  theme(axis.text.x = element_text(face="bold", size=22),
        axis.text.y = element_text(face="bold", size=22)) +
  theme(axis.title=element_text(size=28), plot.title = element_text(size=28, face = "bold")) 
num1
#dev.copy(pdf, height = 12, width = 12, pointsize=4, file.path(pngDir, "volcano_bIFNg_2h.pdf"))
#dev.off()
figure = ggarrange(num1, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_rnaseq_HT1080_SP140KDvscontrol_1.3fc_poster", ".pdf"), width = 8, height = 10)
figure
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()

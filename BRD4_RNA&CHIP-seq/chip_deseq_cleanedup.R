#################################################
#clear global environment
rm(list = ls())

# LOAD REQUIRED PACKAGES
BiocManager::install("ggpubr")

library("DESeq2")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("ggpubr")

#################################################
### PREP INPUT TABLE

### set your working dir
setwd("/Users/adamdziulko/Documents/Chuong Lab/Projects/c_BRD4/CHiP_r_script/")

### read in the tab-delimited count table
countdata <- read.table("bamCountsWithinPeaks.tab", sep="\t", header = FALSE)

### check the table looks correct
head(countdata)

### create a new column which is the chr and coordinates combined
countdata$chr <- paste(countdata$V1,countdata$V2, sep=":")
countdata$coord <- paste(countdata$chr,countdata$V3, sep="-")

### set countdata$coord to be the rownames
rownames(countdata) <- countdata$coord
head(countdata)

### remove the first four columns (chr, start, stop, regionlabel)
### and the last two columns (coord, chr) since we don't need them anymore
### retaining only the counts
countdata <- countdata[, c(5, 6, 7, 8, 9, 10)] 
head(countdata)
### pro count data (1)
countdata1 <- countdata[, c(-3,-4)] 
head(countdata1)
### qui count data (2)
countdata2 <- countdata[, c(-1,-2)] 
head(countdata2)

### add bam names as column names
colnames(countdata) <- c("pro_r1","pro_r2","qui_r1","qui_r2", "sen_r1","sen_r2")
head(countdata)

### convert table to matrix format
countdata <- as.matrix(countdata)
head(countdata)

### assign control vs treated samples
condition <- factor(c(rep("control", 2), rep("treated", 2)))

### create a "coldata" table containing the sample names with their appropriate condition (e.g. control versus cancer sample)
coldata1 <- data.frame(row.names=c("pro_r1","pro_r2","sen_r1","sen_r2"), condition)
coldata1
coldata2 <- data.frame(row.names=c("qui_r1","qui_r2","sen_r1","sen_r2"), condition)
coldata2

#################################################
### RUN DESEQ2

### construct a DESeqDataSet
dds1 <- DESeqDataSetFromMatrix(countData = countdata1, colData = coldata1, design = ~ condition)
dds1
dds2 <- DESeqDataSetFromMatrix(countData = countdata2, colData = coldata2, design = ~ condition)
dds2

### relevel to set the controls as the reference levels
dds1$condition <- relevel(dds1$condition, ref = "control")
dds1

dds2$condition <- relevel(dds2$condition, ref = "control")
dds2

### run deseq2
dds1 <- DESeq(dds1)
resultsNames(dds1)

dds2 <- DESeq(dds2)
resultsNames(dds2)

### call results
res_treated_vs_control1 <- results(dds1, contrast=c("condition", "treated", "control"))
head(res_treated_vs_control1)

res_treated_vs_control2 <- results(dds2, contrast=c("condition", "treated", "control"))
head(res_treated_vs_control2)

### omit rows with counts of "N/A" 
res_treated_vs_control1 <- na.omit(res_treated_vs_control1)
head(res_treated_vs_control1)

res_treated_vs_control2 <- na.omit(res_treated_vs_control2)
head(res_treated_vs_control2)

### sort results by ascending adjusted pvalue
res_treated_vs_control1 <- res_treated_vs_control1[order(res_treated_vs_control1$padj), ]
head(res_treated_vs_control1)

res_treated_vs_control2 <- res_treated_vs_control2[order(res_treated_vs_control2$padj), ]
head(res_treated_vs_control2)

### report the number of rows with an adjusted pvalue less than 0.05
table(res_treated_vs_control1$padj<0.05)

table(res_treated_vs_control2$padj<0.05)

#################################################
# SAVE RESULTS TABLES

### save a table of the raw counts and normalized counts
#raw_counts = counts(dds)
#normalized_counts <- counts(dds, normalized=TRUE)
#write.table(raw_counts, file="raw_counts.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")
#write.table(normalized_counts, file="normalized_counts.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

# extract regions that are significantly different in treated samples compared to control
# separate into two files based on whether they are significantly up (treated > control) or down (control > treated)
sigUp1 = subset(res_treated_vs_control1, padj<0.05 & log2FoldChange>1)
head(sigUp1)
nrow(sigUp1)
sigUp1 <- as.data.frame(sigUp1)
#write.table(sigUp1, file="pro_sen_sigUp.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

sigUp2 = subset(res_treated_vs_control2, padj<0.05 & log2FoldChange>1)
head(sigUp2)
nrow(sigUp2)
sigUp2 <- as.data.frame(sigUp2)
#write.table(sigUp2, file="qui_sen_sigUp.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

sigDown1 = subset(res_treated_vs_control1, padj<0.05 & log2FoldChange<(-1))
head(sigDown1)
nrow(sigDown1)
sigDown1 <- as.data.frame(sigDown1)
#write.table(sigDown1, file="pro_sen_sigDown.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

sigDown2 = subset(res_treated_vs_control2, padj<0.05 & log2FoldChange<(-1))
head(sigDown2)
nrow(sigDown2)
sigDown2 <- as.data.frame(sigDown2)
#write.table(sigDown2, file="qui_sen_sigDown.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

#keep overlaps between pro and qui (both = 3)
sigUp3 = merge(sigUp1, sigUp2, by =0)
length(unique(sigUp3$Row.names)) 
rownames(sigUp3) <- sigUp3$Row.names
sigUp3$baseMean <- rowMeans(sigUp3[,c('baseMean.x', 'baseMean.y')], na.rm=TRUE)
sigUp3$log2FoldChange <- rowMeans(sigUp3[,c('log2FoldChange.x', 'log2FoldChange.y')], na.rm=TRUE)
sigUp3$lfcSE <- rowMeans(sigUp3[,c('lfcSE.x', 'lfcSE.y')], na.rm=TRUE)
sigUp3$stat <- rowMeans(sigUp3[,c('stat.x', 'stat.y')], na.rm=TRUE)
sigUp3$pvalue <- rowMeans(sigUp3[,c('pvalue.x', 'pvalue.y')], na.rm=TRUE)
sigUp3$padj <- rowMeans(sigUp3[,c('padj.x', 'padj.y')], na.rm=TRUE)
sigUp3 <- sigUp3[, c(14, 15, 16, 17, 18, 19)] 
sigUp3 <- sigUp3[order(sigUp3$padj), ]
#write.table(sigUp3, file="proqui_sen_sigUp.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")


sigDown3 = merge(sigDown1, sigDown2, by =0)
length(unique(sigDown3$Row.names)) #org = 3166
rownames(sigDown3) <- sigDown3$Row.names
sigDown3$baseMean <- rowMeans(sigDown3[,c('baseMean.x', 'baseMean.y')], na.rm=TRUE)
sigDown3$log2FoldChange <- rowMeans(sigDown3[,c('log2FoldChange.x', 'log2FoldChange.y')], na.rm=TRUE)
sigDown3$lfcSE <- rowMeans(sigDown3[,c('lfcSE.x', 'lfcSE.y')], na.rm=TRUE)
sigDown3$stat <- rowMeans(sigDown3[,c('stat.x', 'stat.y')], na.rm=TRUE)
sigDown3$pvalue <- rowMeans(sigDown3[,c('pvalue.x', 'pvalue.y')], na.rm=TRUE)
sigDown3$padj <- rowMeans(sigDown3[,c('padj.x', 'padj.y')], na.rm=TRUE)
sigDown3 <- sigDown3[, c(14, 15, 16, 17, 18, 19)] 
sigDown3 <- sigDown3[order(sigDown3$padj), ]
#write.table(sigDown3, file="proqui_sen_sigDown.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

#################################################
### PLOT STUFF

### make an MA plot
plotMA(res_treated_vs_control1, ylim=c(-10,10))

plotMA(res_treated_vs_control2, ylim=c(-10,10))

### save a copy of the MA plot 
#dev.copy(png,'treated_vs_control_MAplot.png')
#dev.off()

### make a volcano plot using ggplot2
### first, make the results table a data frame
res_treated_vs_control1 <- as.data.frame(res_treated_vs_control1)
res_treated_vs_control1

res_treated_vs_control2 <- as.data.frame(res_treated_vs_control2)
res_treated_vs_control2

ggplot(res_treated_vs_control1, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  labs(y="-log10 adjusted pvalue", x = "log2 fold change") +
  ### set all dots to be grey
  geom_point(data=res_treated_vs_control1, colour = "grey") + 
  ### if pvalue<0.05, change dot color to green
  geom_point(data=res_treated_vs_control1[which(res_treated_vs_control1 $padj <0.05),], colour = "springgreen2") + 
  ### if log2FC >1, change dot color to orange
  geom_point(data=res_treated_vs_control1[which(abs(res_treated_vs_control1 $log2FoldChange)>1),], colour = "darkgoldenrod1") +
  ### if both, change dot color to blue
  geom_point(data=res_treated_vs_control1[which(abs(res_treated_vs_control1 $log2FoldChange)>1 & res_treated_vs_control1$padj<0.05),], colour = "royalblue1") 
  ### add text labels to the most significant regions
  #geom_text_repel(data =res_treated_vs_control1[which(res_treated_vs_control1 $padj <0.000005),], mapping = aes(log2FoldChange, -log10(padj), label = rownames(res_treated_vs_control1[which(res_treated_vs_control1 $padj <0.000005),])),size = 4,force = 1)

### save a copy of the volcano plot
dev.copy(png, res=200, height = 1000, width = 1000, pointsize=4, 'sen.treated_vs_pro.control_volcano.png')
dev.off()


 ggplot(res_treated_vs_control2, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  labs(y="-log10 adjusted pvalue", x = "log2 fold change") +
  ### set all dots to be grey
  geom_point(data=res_treated_vs_control2, colour = "grey") + 
  ### if pvalue<0.05, change dot color to green
  geom_point(data=res_treated_vs_control2[which(res_treated_vs_control2 $padj <0.05),], colour = "springgreen2") + 
  ### if log2FC >1, change dot color to orange
  geom_point(data=res_treated_vs_control2[which(abs(res_treated_vs_control2 $log2FoldChange)>1),], colour = "darkgoldenrod1") +
  ### if both, change dot color to blue
  geom_point(data=res_treated_vs_control2[which(abs(res_treated_vs_control2 $log2FoldChange)>1 & res_treated_vs_control2$padj<0.05),], colour = "royalblue1") +
  ### add text labels to the most significant regions
  geom_text_repel(data =res_treated_vs_control2[which(res_treated_vs_control2 $padj <0.000005),], mapping = aes(log2FoldChange, -log10(padj), label = rownames(res_treated_vs_control2[which(res_treated_vs_control2 $padj <0.000005),])),size = 4,force = 1)

### save a copy of the volcano plot
dev.copy(png, res=200, height = 1000, width = 1000, pointsize=4, 'sen.treated_vs_qui.control_volcano.png')
dev.off()







#_______________________________________________________________________________
### Volcano for giggle files
chip_deseq_UP_GIGGLErepeats <- read.table("proqui_sen_sigUp_vs_repeats_sorted.giggle",
                   sep = "",
                   header = FALSE)

chip_deseq_DOWN_GIGGLErepeats <- read.table("proqui_sen_sigDown_vs_repeats_sorted.giggle",
                                          sep = "",
                                          header = FALSE)

### Set the column names
colnames(chip_deseq_UP_GIGGLErepeats) <- c('repeat', 'repeat.total.loci', 'overlaps',
                    'odds.ratio', 'fisher.two.tail', 'fisher.left',
                    'fisher.right', 'giggle.score')

colnames(chip_deseq_DOWN_GIGGLErepeats) <- c('repeat', 'repeat.total.loci', 'overlaps',
                                           'odds.ratio', 'fisher.two.tail', 'fisher.left',
                                           'fisher.right', 'giggle.score')

### add column calculating percent of senescent loci vs total loci
chip_deseq_UP_GIGGLErepeats <- transform(chip_deseq_UP_GIGGLErepeats, 
                                         percent.sen.loci = (overlaps / repeat.total.loci) * 100)

chip_deseq_DOWN_GIGGLErepeats <- transform(chip_deseq_DOWN_GIGGLErepeats, 
                                         percent.sen.loci = (overlaps / repeat.total.loci) * 100)

### Calculate the expected overlap by dividing the observed overlap by the odds ratio & set this as the fourth column
chip_deseq_UP_GIGGLErepeats <- transform(chip_deseq_UP_GIGGLErepeats, exp = overlaps / odds.ratio) %>% .[,c(1,2,3,10,4,5,6,7,8,9)]

chip_deseq_DOWN_GIGGLErepeats <- transform(chip_deseq_DOWN_GIGGLErepeats, exp = overlaps / odds.ratio) %>% .[,c(1,2,3,10,4,5,6,7,8,9)]

### Calculate log2FC
chip_deseq_UP_GIGGLErepeats <- transform(chip_deseq_UP_GIGGLErepeats, log2FC = log2(overlaps / exp)) %>% .[,c(1,2,3,10,11,4,5,6,7,8,9)]

chip_deseq_DOWN_GIGGLErepeats <- transform(chip_deseq_DOWN_GIGGLErepeats, log2FC = log2(overlaps / exp)) %>% .[,c(1,2,3,10,11,4,5,6,7,8,9)]

### Reorder columns
chip_deseq_UP_GIGGLErepeats = chip_deseq_UP_GIGGLErepeats[,c(1,5,11,4,2,3,6,7,8,9,10)]

chip_deseq_DOWN_GIGGLErepeats = chip_deseq_DOWN_GIGGLErepeats[,c(1,5,11,4,2,3,6,7,8,9,10)]

### Subset your data for insignificant, positively enriched, and negatively enriched repeats
### The p-val and fold change cutoffs here are very arbitrary - play around with them to see what fits your data well
insig_giggleUP <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail > 1e-20 | overlaps / exp < 4)
pos_enrich_giggleUP <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp > 4)
neg_enrich_giggleUP <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp < .5)
### below is different cutoff (for figure 5 (all data))
insig_giggleUP2 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail > 1e-20 | overlaps / exp < 1)
pos_enrich_giggleUP2 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp > 1)
neg_enrich_giggleUP2 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp < 1)
### below us cutoff for odds.ratio in place of obs/exp
insig_giggleUP3 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail > 1e-20 | odds.ratio < 4)
pos_enrich_giggleUP3 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & odds.ratio > 4)
neg_enrich_giggleUP3 <- subset(chip_deseq_UP_GIGGLErepeats, fisher.two.tail < 1e-20 & odds.ratio < .5)

insig_giggleDOWN3 <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail > 1e-20 | odds.ratio < 4)
pos_enrich_giggleDOWN3 <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail < 1e-20 & odds.ratio > 4)
neg_enrich_giggleDOWN3 <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail < 1e-20 & odds.ratio < .5)

insig_giggleDOWN <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail > 1e-20 | overlaps / exp < 4)
pos_enrich_giggleDOWN <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp > 4)
neg_enrich_giggleDOWN <- subset(chip_deseq_DOWN_GIGGLErepeats, fisher.two.tail < 1e-20 & overlaps / exp < .5)



### pos/neg only in sen
### pos only in sen
pos_enrich_giggleUP$state = "sen"
pos_enrich_giggleDOWN$state = "pro.qui"
#
pos_enrich_giggleBOTH <- rbind(pos_enrich_giggleUP,pos_enrich_giggleDOWN)
#
pos_enrich_giggleBOTHminPROQUI = pos_enrich_giggleBOTH %>% group_by(repeat.) %>% filter(n() == 1)
#
pos_enrich_giggleBOTHminPROQUI = pos_enrich_giggleBOTHminPROQUI[grepl("sen", pos_enrich_giggleBOTHminPROQUI$state),]


### neg only in sen
neg_enrich_giggleUP$state = "sen"
neg_enrich_giggleDOWN$state = "pro.qui"
#
neg_enrich_giggleBOTH <- rbind(neg_enrich_giggleUP, neg_enrich_giggleDOWN)
#
neg_enrich_giggleBOTHminPROQUI = neg_enrich_giggleBOTH %>% group_by(repeat.) %>% filter(n() == 1)
#
neg_enrich_giggleBOTHminPROQUI = neg_enrich_giggleBOTHminPROQUI[grepl("sen", neg_enrich_giggleBOTHminPROQUI$state),]


### pos only in sen (different subset)
pos_enrich_giggleUP3$state = "sen"
pos_enrich_giggleDOWN3$state = "pro.qui"
#
pos_enrich_giggleBOTH3 <- rbind(pos_enrich_giggleUP3,pos_enrich_giggleDOWN3)
#
pos_enrich_giggleBOTHminPROQUI3 = pos_enrich_giggleBOTH3 %>% group_by(repeat.) %>% filter(n() == 1)
#
pos_enrich_giggleBOTHminPROQUI3 = pos_enrich_giggleBOTHminPROQUI3[grepl("sen", pos_enrich_giggleBOTHminPROQUI3$state),]


### neg only in sen (different subset)
neg_enrich_giggleUP3$state = "sen"
neg_enrich_giggleDOWN3$state = "pro.qui"
#
neg_enrich_giggleBOTH3 <- rbind(neg_enrich_giggleUP3, neg_enrich_giggleDOWN3)
#
neg_enrich_giggleBOTHminPROQUI3 = neg_enrich_giggleBOTH3 %>% group_by(repeat.) %>% filter(n() == 1)
#
neg_enrich_giggleBOTHminPROQUI3 = neg_enrich_giggleBOTHminPROQUI3[grepl("sen", neg_enrich_giggleBOTHminPROQUI3$state),]



### Plots
########

### Call ggplot & make your volcano plot
### Tell ggplot to plot log2(obs / exp) on the x-axis and -log10(p value) on the y-axis
num2 = ggplot(data = chip_deseq_UP_GIGGLErepeats, aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail))) +
  ### Change the 'theme' or appearance of the plot such that the background is white
  theme_bw() +
  ### Give your plot a title
  ggtitle(paste0("TE families in senescent ", "lung fibroblast")) +
  ### Name your axes
  labs(x = "log2 fold change (observed/expected)", y = "-log10 p value") +
  ### Adjust the size and placement (hjust & vjust) of your labels
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 35, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 35, hjust = 0.5, vjust = 0.5)) +
  ### Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ### Plot 'insigificant' (subset above) repeats as dark grey dots (color = "#4D4D4D") that are smaller (size = 4) and more transparent (alpha = 0.5) than significant repeats
  geom_point(data = insig_giggleUP, color = "#4D4D4D", size = 4, alpha = 0.5, stroke = 0) +
  ### Plot positively enriched repeats as light blue dots that are larger and more opaque than insignificant repeats
  geom_point(data = pos_enrich_giggleBOTHminPROQUI, color = "#A8DDB5", size = 8, stroke = 0, alpha = 0.8) +
  ### Plot negatively enriched repeats as light green dots that are larger and more opaque than insigificant repeats
  geom_point(data = neg_enrich_giggleBOTHminPROQUI, color = "red", size = 8, stroke = 0, alpha = 0.8) +
  ### Tell ggplot to make the x-axis scale from a log2(obs / exp) of -6 to 6 with breaks at intervals of 2
  ### Change this to fit your data!
  scale_x_continuous(limits = c(-6,6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(0,65),
                     breaks = c(0, 20, 40, 60))
  ### Label your positively enriched repeats such that the labels are (ideally) not overlapping
  ### Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  #geom_text_repel(data = pos_enrich_giggleBOTHminPROQUI, 
   #               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)), 
    #              label = pos_enrich_giggleBOTHminPROQUI$repeat., 
     #             size = 6, 
      #            force = 2,
       #           nudge_y = 0.5,
        #          nudge_x = 0.5,
         #         point.padding = 0.25) +
  ### Label your negatively enriched repeats such that the labels are (ideally) not overlapping
  ### Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  #geom_text_repel(data = neg_enrich_giggleBOTHminPROQUI,
   #               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)),
    #              label = neg_enrich_giggleBOTHminPROQUI$repeat.,
     #             size = 6,
      #            force = 2,
       #           nudge_y = 0.5,
        #          nudge_x = 0.5,
 3        #         point.padding = 0)
## Save your plot
figure2 = ggarrange(num2, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_giggle_sen_volcano_NONannotated", ".pdf"), width = 14, height = 15)
figure2
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()




# Call ggplot & make your volcano plot
## Tell ggplot to plot log2(obs / exp) on the x-axis and -log10(p value) on the y-axis
num3 = ggplot(data = chip_deseq_UP_GIGGLErepeats, aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail))) +
  ## Change the 'theme' or appearance of the plot such that the background is white
  theme_bw() +
  ## Give your plot a title
  ggtitle(paste0("TE families in senescent ", "lung fibroblast")) +
  ## Name your axes
  labs(x = "log2 fold change (observed/expected)", y = "-log10 p value") +
  ## Adjust the size and placement (hjust & vjust) of your labels
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 35, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 35, hjust = 0.5, vjust = 0.5)) +
  ## Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ## Plot 'insigificant' (subset above) repeats as dark grey dots (color = "#4D4D4D") that are smaller (size = 4) and more transparent (alpha = 0.5) than significant repeats
  geom_point(data = insig_giggleUP, color = "#4D4D4D", size = 4, alpha = 0.5, stroke = 0) +
  ## Plot positively enriched repeats as light blue dots that are larger and more opaque than insignificant repeats
  geom_point(data = pos_enrich_giggleBOTHminPROQUI, color = "#A8DDB5", size = 8, stroke = 0, alpha = 0.8) +
  ## Plot negatively enriched repeats as light green dots that are larger and more opaque than insigificant repeats
  geom_point(data = neg_enrich_giggleBOTHminPROQUI, color = "red", size = 8, stroke = 0, alpha = 0.8) +
  ## Tell ggplot to make the x-axis scale from a log2(obs / exp) of -6 to 6 with breaks at intervals of 2
  ## Change this to fit your data!
  scale_x_continuous(limits = c(-6,6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(0,65),
                     breaks = c(0, 20, 40, 60)) +
## Label your positively enriched repeats such that the labels are (ideally) not overlapping
## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
geom_text_repel(data = pos_enrich_giggleBOTHminPROQUI, 
               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)), 
              label = pos_enrich_giggleBOTHminPROQUI$repeat., 
             size = 8, 
            force = 2,
           nudge_y = 0.5,
          nudge_x = 0.5,
         point.padding = 0.25) +
## Label your negatively enriched repeats such that the labels are (ideally) not overlapping
## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
geom_text_repel(data = neg_enrich_giggleBOTHminPROQUI,
               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)),
              label = neg_enrich_giggleBOTHminPROQUI$repeat.,
             size = 8,
            force = 2,
           nudge_y = -1.5,
          nudge_x = 0.5,
         point.padding = 0)
## Save your plot
figure3 = ggarrange(num3, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_giggle_sen_volcano_comps_insigcorrection", ".pdf"), width = 14, height = 15)
figure3
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()


num6 = ggplot(data = chip_deseq_UP_GIGGLErepeats, aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail))) +
  ## Change the 'theme' or appearance of the plot such that the background is white
  theme_bw() +
  ## Give your plot a title
  ggtitle(paste0("TE families in senescent ", "lung fibroblast")) +
  ## Name your axes
  labs(x = "log2 fold change (observed/expected)", y = "-log10 p value") +
  ## Adjust the size and placement (hjust & vjust) of your labels
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 35, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 35, hjust = 0.5, vjust = 0.5)) +
  ## Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ## Plot 'insigificant' (subset above) repeats as dark grey dots (color = "#4D4D4D") that are smaller (size = 4) and more transparent (alpha = 0.5) than significant repeats
  geom_point(data = insig_giggleUP, color = "#4D4D4D", size = 4, alpha = 0.5, stroke = 0) +
  ## Plot positively enriched repeats as light blue dots that are larger and more opaque than insignificant repeats
  geom_point(data = pos_enrich_giggleBOTHminPROQUI, color = "#A8DDB5", size = 8, stroke = 0, alpha = 0.8) +
  ## Plot negatively enriched repeats as light green dots that are larger and more opaque than insigificant repeats
  geom_point(data = neg_enrich_giggleBOTHminPROQUI, color = "red", size = 8, stroke = 0, alpha = 0.8) +
  ## Tell ggplot to make the x-axis scale from a log2(obs / exp) of -6 to 6 with breaks at intervals of 2
  ## Change this to fit your data!
  scale_x_continuous(limits = c(-6,6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(0,65),
                     breaks = c(0, 20, 40, 60)) 
  ## Label your positively enriched repeats such that the labels are (ideally) not overlapping
  ## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  #geom_text_repel(data = pos_enrich_giggleBOTHminPROQUI, 
   #               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)), 
    #              label = pos_enrich_giggleBOTHminPROQUI$repeat., 
     #             size = 8, 
      #            force = 2,
       #           nudge_y = 0.5,
        #          nudge_x = 0.5,
         #         point.padding = 0.25) +
  ## Label your negatively enriched repeats such that the labels are (ideally) not overlapping
  ## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  #geom_text_repel(data = neg_enrich_giggleBOTHminPROQUI,
   #               aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)),
    #              label = neg_enrich_giggleBOTHminPROQUI$repeat.,
     #             size = 8,
      #            force = 2,
       #           nudge_y = -1.5,
        #          nudge_x = 0.5,
         #         point.padding = 0)
## Save your plot
figure6 = ggarrange(num6, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_giggle_sen_volcano_comps_insigcorrection_NONannotated", ".pdf"), width = 14, height = 15)
figure6
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()


###########
###testing out only plotting sig up
num4 = ggplot(data = chip_deseq_UP_GIGGLErepeats, aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail))) +
  ## Change the 'theme' or appearance of the plot such that the background is white
  theme_bw() +
  ## Give your plot a title
  ggtitle(paste0("TE families in senescent ", "lung fibroblast")) +
  ## Name your axes
  labs(x = "log2 fold change (observed/expected)", y = "-log10 p value") +
  ## Adjust the size and placement (hjust & vjust) of your labels
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 35, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 40, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 35, hjust = 0.5, vjust = 0.5)) +
  ## Draw a vertical dotted line at x = 0
  geom_vline(xintercept = 0, linetype = "dotted") +
  ## Plot 'insigificant' (subset above) repeats as dark grey dots (color = "#4D4D4D") that are smaller (size = 4) and more transparent (alpha = 0.5) than significant repeats
  geom_point(data = insig_giggleUP, color = "#4D4D4D", size = 4, alpha = 0.5, stroke = 0) +
  ## Plot positively enriched repeats as light blue dots that are larger and more opaque than insignificant repeats
  geom_point(data = pos_enrich_giggleUP, color = "#A8DDB5", size = 8, stroke = 0, alpha = 0.8) +
  ## Plot negatively enriched repeats as light green dots that are larger and more opaque than insigificant repeats
  geom_point(data = neg_enrich_giggleUP, color = "red", size = 8, stroke = 0, alpha = 0.8) +
  ## Tell ggplot to make the x-axis scale from a log2(obs / exp) of -6 to 6 with breaks at intervals of 2
  ## Change this to fit your data!
  scale_x_continuous(limits = c(-6,6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(0,65),
                     breaks = c(0, 20, 40, 60)) +
  ## Label your positively enriched repeats such that the labels are (ideally) not overlapping
  ## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  geom_text_repel(data = pos_enrich_giggleUP, 
                  aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)), 
                  label = pos_enrich_giggleUP$repeat., 
                  size = 8, 
                  force = 2,
                  nudge_y = 0.5,
                  nudge_x = 0.5,
                  point.padding = 0.25) +
  ## Label your negatively enriched repeats such that the labels are (ideally) not overlapping
  ## Note: the effectiveness of ggrepel in distributing your labels such that they are not overlapping is highly dependent on the size of your plot
  geom_text_repel(data = neg_enrich_giggleUP,
                  aes(x = log2(overlaps / exp), y = -log10(fisher.two.tail)),
                  label = neg_enrich_giggleUP$repeat.,
                  size = 8,
                  force = 2,
                  nudge_y = -1.5,
                  nudge_x = 0.5,
                  point.padding = 0)
## Save your plot
figure4 = ggarrange(num4, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_giggle_sen_removedBoth", ".pdf"), width = 14, height = 15)
figure4
#annotate_figure(figure, top = text_grob("test volcano", color = "red", face = "bold", size = 28))
dev.off()






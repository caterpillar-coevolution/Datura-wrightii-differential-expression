#DEG analysis of LEma induced Datura wrightii plants
#prepare R packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("DEGreport")
BiocManager::install("ComplexHeatmap")
BiocManager::install("pheatmap", force = TRUE)
install.packages("tidyverse")

library("pheatmap")
library("tidyverse")
library("DEGreport")
library( "DESeq2" )
library( "EnhancedVolcano" )
library(dplyr)
library(tibble)
#########################################################################
#begin pairwise analysis
#load data
metaData <- read.csv("Dwri_metadata.csv")
geneData <- read.csv("GeneCounts_from_STAR.csv")

#run model and store results
dds <- DESeqDataSetFromMatrix(countData=geneData, 
                              colData=metaData, 
                              design=~experiment+treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)

#order output by pvalue
res <- res[order(res$padj),]
head(res)

#set significance cutoff
padj.cutoff <- 0.05

#make a subset of the data with only significant genes
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

write.csv(as.data.frame(res),
          file="Datura_corrected_deseq_output.csv")
write.csv(as.data.frame(sig_res),
          file="Datura_corrected_deseq_output_SigGenes.csv") #supplemental table 3

#make PCA plot
vsdata <- vst(dds, blind=FALSE)
PCAplot <- plotPCA(vsdata, intgroup="treatment")
PCAplot + theme_classic()

#make volcano plot
EnhancedVolcano(res,
                lab = NA,
                labSize = 2.0,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Control vs. lema induced',
                pointSize = 1.5,
                pCutoff = 0.0075, #this is raw p-value, not adjusted
                FCcutoff = 2,
                col=c('black', 'black', 'red3' , 'red3'),
                colAlpha = 1.)
##################################################################################
# begin timecourse analysis
metaDataTC <- read.csv("Dwri_metadata_TConly.csv")
geneDataTC <- read.csv("GeneCounts_from_STAR_TConly.csv")

ddsTC <- DESeqDataSetFromMatrix(countData=geneDataTC, 
                              colData=metaDataTC, 
                              design=~treatment + timepoint + treatment:timepoint, tidy = TRUE)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~treatment + timepoint)
resTC <- results(ddsTC)

#sort results by p-value
resTC <- resTC[order(resTC$padj),]
head(resTC)

# Subset the LRT results to return genes with padj < 0.05
sig_resTC <- resTC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

#save results to excel files
write.csv(as.data.frame(resTC),
          file="Datura_corrected_deseq_output_TConly.csv")
write.csv(as.data.frame(sig_resTC),
          file="Datura_corrected_deseq_output_TC_SigGenes.csv") #supplemental table 4

# make a heatmap
betas <- coef(ddsTC)
colnames(betas)

topGenes <- head(order(resTC$padj),50)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

#make PCA plot
vsdataTC <- vst(ddsTC, blind=FALSE)
PCAplotTC <- plotPCA(vsdataTC, intgroup="timepoint")
PCAplotTC + theme_classic()

#DEG analysis of LEma induced Datura wrightii plants
#prepare R packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("DEGreport")
BiocManager::install("ComplexHeatmap")
BiocManager::install("pheatmap")
install.packages("tidyverse")
install.packages("cluster")

library ("cluster")
library("pheatmap")
library("tidyverse")
library("DEGreport")
library( "DESeq2" )
library( "EnhancedVolcano" )
library(dplyr)
library(tibble)
#############################################################
# begin timecourse analysis
metaDataTC <- read.csv("Dwri_metadata_TConly.csv", fileEncoding="UTF-8-BOM",row.names = 1)
geneDataTC <- read.csv("Datura_timecourse_geneCounts.csv", fileEncoding="UTF-8-BOM")

metaDataTC$timepoint <- as.factor(metaDataTC$timepoint)


#set significance cutoff
padj.cutoff <- 0.05

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
          file="Datura_timecourse_results.csv")
write.csv(as.data.frame(sig_resTC),
          file="Datura_timecourse_sigGenes.csv") #supplemental table 4

#make PCA plot
vsdataTC <- vst(ddsTC, blind=FALSE)
PCAplotTC <- plotPCA(vsdataTC, intgroup=c("treatment", "timepoint"))
PCAplotTC + theme_classic()
#############################################################################
#CLUSTER ANALYSIS USING DEGREPORT
# make matrix of transformed values
rld_mat <- assay(vsdataTC)
clustering_sig_genes <- sig_resTC %>%
  arrange(padj) %>%
  head(n=1000)

cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
clusters <- degPatterns(cluster_rlog, metadata = metaDataTC, time = "timepoint", col = "treatment")
degPlotCluster(clusters$normalized,"timepoint", "treatment")

#export results as a spreadsheet with gene IDs and their cluster
cluster_groups <- clusters$df
write.csv(as.data.frame(cluster_groups),
          file="Datura_timecourse_cluster_results.csv")

#plot specific groups
group1 <- subset(clusters$normalized, clusters$normalized$genes %in% group1$genes)
degPlotCluster(group1$normalized,"timepoint", "treatment")
##################################################################################
#join cluster results with main results table
sigClusters <- read.csv("Datura_timecourse_cluster_results.csv")
sigGenes <- read.csv("Datura_timecourse_sigGenes.csv")

sigJoined <- left_join(sigGenes, sigClusters, by = "gene")
write.csv(as.data.frame(sigJoined),
          file="Datura_timecourse_joined_geneList.csv")


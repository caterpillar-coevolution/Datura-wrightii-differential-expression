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
install.packages("cluster")


require("clusterProfiler")
library ("cluster")
library("pheatmap")
library("tidyverse")
library("DEGreport")
library( "DESeq2" )
library( "EnhancedVolcano" )
library(dplyr)
library(tibble)
##############################################################
#begin pairwise analysis of MAy 2012 only dataset
#load data
metaData <- read.csv("Dwri_metadata.csv")
geneData <- read.csv("pairwise_geneCounts.csv")

metaData$replicate <- as.factor(metaData$replicate)

#run model and store results
dds <- DESeqDataSetFromMatrix(countData=geneData, 
                              colData=metaData, 
                              design=~replicate+treatment, tidy = TRUE)
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
          file="Datura_pairwise_results.csv") #for GSEA later on
write.csv(as.data.frame(sig_res),
          file="Datura_pairwise_siggenes.csv") #supplemental table 3

#make PCA plot
vsdata <- vst(dds, blind=FALSE)
PCAplot <- plotPCA(vsdata, intgroup="treatment")
PCAplot + theme_classic()

#make volcano plot
EnhancedVolcano(res,
                lab = NA,
                labSize = 2.0,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Control vs. lema induced',
                pointSize = 1.5,
                pCutoff = 0.05, #this is raw p-value, not adjusted
                FCcutoff = 2,
                col=c('black', 'black', 'red3' , 'red3'),
                colAlpha = 1.)
##########################################################################
#GSEA with results from previous analysis

#load the data. separate gene lists for up and down regulated
geneList_upreg <- read.csv("Datura_pairwise_results_upreg.csv", fileEncoding="UTF-8-BOM")
geneList_downreg <- read.csv("Datura_pairwise_results_downreg.csv", fileEncoding="UTF-8-BOM")


#load the TERM2GENE file (connects gene IDs to their corresponding Go terms)
bp_terms <- read.csv("gene2bp.csv", fileEncoding="UTF-8-BOM")
mf_terms <- read.csv("gene2mf.csv", fileEncoding="UTF-8-BOM")
cc_terms <- read.csv("gene2cc.csv", fileEncoding="UTF-8-BOM")

#load the GO2NAME files
bp_names <- read.csv("bp2name.csv", fileEncoding="UTF-8-BOM")
mf_names <- read.csv("mf2name.csv", fileEncoding="UTF-8-BOM")
cc_names <- read.csv("cc2name.csv", fileEncoding="UTF-8-BOM")


#subset upregulated genes to only include genes with functional annotations
bp_upreg_joined <- inner_join(geneList_upreg, bp_terms, by = "gene_ID") 
mf_upreg_joined <- inner_join(geneList_upreg, mf_terms, by = "gene_ID") 
cc_upreg_joined <- inner_join(geneList_upreg, cc_terms, by = "gene_ID") 

#subset downregulated genelists like above
bp_downreg_joined <- inner_join(geneList_downreg, bp_terms, by = "gene_ID") 
mf_downreg_joined <- inner_join(geneList_downreg, mf_terms, by = "gene_ID") 
cc_downreg_joined <- inner_join(geneList_downreg, cc_terms, by = "gene_ID") 

#write joined datasets to csvs for later reference
write.csv(as.data.frame(bp_upreg_joined),
          file="bp_upreg_joined.csv")
write.csv(as.data.frame(mf_upreg_joined),
          file="mf_upreg_joined.csv")
write.csv(as.data.frame(cc_upreg_joined),
          file="cc_upreg_joined.csv")

#do the same for downregulated datasets
write.csv(as.data.frame(bp_downreg_joined),
          file="bp_downreg_joined.csv")
write.csv(as.data.frame(mf_downreg_joined),
          file="mf_downreg_joined.csv")
write.csv(as.data.frame(cc_downreg_joined),
          file="cc_downreg_joined.csv")

#QC BP dataset
bp_genes = bp_upreg_joined[,2]
names(bp_genes) = as.character(bp_upreg_joined[,1])
bp_genes = sort(bp_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
y <- GSEA(bp_genes, TERM2GENE = bp_terms, TERM2NAME = bp_names)
write.csv(y,
          file="bp_upreg_results.csv")
dotplot(y, showCategory=30)

#QC MF dataset
mf_genes = mf_upreg_joined[,2]
names(mf_genes) = as.character(mf_upreg_joined[,1])
mf_genes = sort(mf_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
x <- GSEA(mf_genes, TERM2GENE = mf_terms, TERM2NAME = mf_names)
write.csv(x,
          file="mf_upreg_results.csv")
dotplot(x, showCategory=30)

#QC MF dataset
cc_genes = cc_upreg_joined[,2]
names(cc_genes) = as.character(cc_upreg_joined[,1])
cc_genes = sort(cc_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
z <- GSEA(cc_genes, TERM2GENE = cc_terms, TERM2NAME = cc_names)
write.csv(z,
          file="cc_upreg_results.csv")
dotplot(z, showCategory=30)

#QC BP baseline dataset
bp_based = bp_downreg_joined[,2]
names(bp_based) = as.character(bp_downreg_joined[,1])
bp_based = sort(bp_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
a <- GSEA(bp_based, TERM2GENE = bp_terms, TERM2NAME = bp_names)
write.csv(a,
          file="bp_downreg_results.csv")
dotplot(a, showCategory=30)

#QC MF baseline dataset
mf_based = mf_downreg_joined[,2]
names(mf_based) = as.character(mf_downreg_joined[,1])
mf_based = sort(mf_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
b <- GSEA(mf_based, TERM2GENE = mf_terms, TERM2NAME = mf_names)
write.csv(b,
          file="mf_downreg_results.csv")
dotplot(b, showCategory=30)

#QC CC baseline dataset
cc_based = cc_downreg_joined[,2]
names(cc_based) = as.character(cc_downreg_joined[,1])
cc_based = sort(cc_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
c <- GSEA(cc_based, TERM2GENE = cc_terms, TERM2NAME = cc_names)
write.csv(c,
          file="cc_downreg_results.csv")
dotplot(c, showCategory=30)



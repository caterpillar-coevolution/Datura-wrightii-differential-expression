#enrichment analyses to find the role of WGD in herbivore interactions
BiocManager::install("clusterProfiler")
library("tidyverse")
require("clusterProfiler")

#load the datasets
term2gene <- read.csv("gene2mf.csv", fileEncoding="UTF-8-BOM")
term2name <- read.csv("mf2name.csv", fileEncoding="UTF-8-BOM")
WGD_genes <- read.csv("geneList_WGD.csv", fileEncoding="UTF-8-BOM")
DIS_genes <- read.csv("geneList_dispersed.csv", fileEncoding="UTF-8-BOM")
PROX_genes <- read.csv("geneList_proximal.csv", fileEncoding="UTF-8-BOM")
TAN_genes <- read.csv("geneList_tandems.csv", fileEncoding="UTF-8-BOM")
TC_genes <- read.csv("geneList_timecourseDEGs.csv", fileEncoding="UTF-8-BOM")
PW_genes <- read.csv("geneList_pairwiseDEGs.csv", fileEncoding="UTF-8-BOM")

#joining functions to make appropriate datasets
DE_genes <- full_join(PW_genes, TC_genes)
write.csv(as.data.frame(DE_genes),
          file="geneList_joinedDEGs.csv")

DEG_w_GO <- inner_join(DE_genes, term2gene)
write.csv(as.data.frame(DEG_w_GO),
          file="geneList_DEGs_w_goterms.csv")

WGD_w_GO <- inner_join(WGD_genes, term2gene)
write.csv(as.data.frame(WGD_w_GO),
          file="geneList_WGD_w_goterms.csv")

DIS_W_GO <- inner_join(DIS_genes, term2gene)
write.csv(as.data.frame(DIS_W_GO),
          file="geneList_dispersed_w_goterms.csv")

PROX_W_GO <- inner_join(PROX_genes, term2gene)
write.csv(as.data.frame(PROX_W_GO),
          file="geneList_dispersed_w_goterms.csv")

TAN_W_GO <- inner_join(TAN_genes, term2gene)
write.csv(as.data.frame(TAN_W_GO),
          file="geneList_dispersed_w_goterms.csv")

#try enrichment analysis with "all" genes
WGD_list <- as.vector(WGD_w_GO$gene)
x <- enricher(WGD_list, TERM2GENE = term2gene, TERM2NAME = term2name)
write.csv(x,
          file="WGD_enrichment_results.csv")
dotplot(x, showCategory=30)

DIS_list <- as.vector(DIS_W_GO$gene)
y <- enricher(DIS_list, TERM2GENE = term2gene, TERM2NAME = term2name)
write.csv(y,
          file="DIS_enrichment_results.csv")
dotplot(y, showCategory=30)

PROX_list <- as.vector(PROX_W_GO$gene)
z <- enricher(PROX_list, TERM2GENE = term2gene, TERM2NAME = term2name)
write.csv(z,
          file="PROX_enrichment_results.csv")
dotplot(z, showCategory=30)

TAN_list <- as.vector(TAN_W_GO$gene)
w <- enricher(TAN_list, TERM2GENE = term2gene, TERM2NAME = term2name)
write.csv(w,
          file="TAN_enrichment_results.csv")
dotplot(w, showCategory=30)

###################################
# load datasets of only DEGs of various dup types
WGD_DEGs <- read.csv("geneList_WGD_DEGs.csv", fileEncoding="UTF-8-BOM")
DIS_DEGs <- read.csv("geneList_DIS_DEGs.csv", fileEncoding="UTF-8-BOM")
PROX_DEGs <- read.csv("geneList_PROX_DEGs.csv", fileEncoding="UTF-8-BOM")
TAN_DEGs <- read.csv("geneList_TAN_DEGs.csv", fileEncoding="UTF-8-BOM")

#and now with just DEGS
DEG_list <- as.vector(DEG_w_GO$gene)
WGD_DEG_list <- as.vector(WGD_DEGs$gene)
a <- enricher(WGD_DEG_list, TERM2GENE = term2gene, TERM2NAME = term2name, universe = DEG_list)
write.csv(a,
          file="WGD_DEG_enrichment_results.csv")
dotplot(a, showCategory=30)

DIS_DEG_list <- as.vector(DIS_DEGs$gene)
b <- enricher(DIS_DEG_list, TERM2GENE = term2gene, TERM2NAME = term2name, universe = DEG_list)
write.csv(b,
          file="DIS_DEG_enrichment_results.csv")
dotplot(b, showCategory=30)

PROX_DEG_list <- as.vector(PROX_DEGs$gene)
c <- enricher(PROX_DEG_list, TERM2GENE = term2gene, TERM2NAME = term2name, universe = DEG_list)
write.csv(c,
          file="PROX_DEG_enrichment_results.csv")
dotplot(c, showCategory=30)

TAN_DEG_list <- as.vector(TAN_DEGs$gene)
d <- enricher(TAN_DEG_list, TERM2GENE = term2gene, TERM2NAME = term2name, universe = DEG_list)
write.csv(d,
          file="TAN_DEG_enrichment_results.csv")
dotplot(d, showCategory=30)

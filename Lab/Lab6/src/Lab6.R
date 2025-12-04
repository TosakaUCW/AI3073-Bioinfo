options(warn = -1)
rm(list = ls())
BiocManager::install("genefilter")
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(RColorBrewer))

rawCounts <- read.delim("E-GEOD-50760-raw-counts.tsv", check.names = FALSE)

head(rawCounts)[, 1:6]
dim(rawCounts)

sampleData <- read.delim("E-GEOD-50760-experiment-design.tsv", check.names = FALSE)
head(sampleData)

countData <- rawCounts
rownames(countData) <- countData[["Gene ID"]]
countData <- countData[, !(colnames(countData) %in% c("Gene ID", "Gene Name"))]

all(colnames(countData) %in% sampleData$Run)

sampleData <- sampleData[match(colnames(countData), sampleData$Run), ]
stopifnot(all(sampleData$Run == colnames(countData)))
condition_full <- sampleData[["Sample Characteristic[biopsy site]"]]
condition_full <- factor(condition_full)

table(condition_full)

keepSamples <- condition_full %in% c("normal", "primary tumor")
countData_sub <- countData[, keepSamples]
condition <- droplevels(condition_full[keepSamples])

table(condition)

colData <- data.frame(
  row.names = colnames(countData_sub),
  condition = condition
)

#
dds <- DESeqDataSetFromMatrix(
  countData = countData_sub,
  colData   = colData,
  design    = ~ condition
)

# dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "primary tumor", "normal"))

resOrdered <- res[order(res$padj), ]

summary(res)
resAll_df <- as.data.frame(resOrdered)
write.csv(resAll_df, file = "Lab6_DESeq2_all_genes_by_padj.csv")

resSig <- subset(resOrdered, padj < 0.1)
resSig_df <- as.data.frame(resSig)
write.csv(resSig_df, file = "Lab6_DESeq2_sig_genes_padj_lt_0.1.csv")
gene_id_list <- rownames(resSig)
length(gene_id_list)
head(gene_id_list)
GO_result <- enrichGO(
  gene          = gene_id_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

GO_df <- as.data.frame(GO_result)
write.csv(GO_df, file = "Lab6_GO_enrichment_sig_DEGs.csv", row.names = FALSE)

head(GO_df[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")])
pdf("Lab6_GO_dotplot.pdf", width = 8, height = 6)
dotplot(
  GO_result,
  x            = "GeneRatio",
  showCategory = 20,
  color        = "p.adjust",
  title        = "GO enrichment of DE genes"
)
dev.off()

pdf("Lab6_GO_barplot.pdf", width = 8, height = 6)
barplot(
  GO_result,
  x            = "GeneRatio",
  showCategory = 20,
  color        = "p.adjust",
  title        = "GO enrichment of DE genes"
)
dev.off()
vsd <- vst(dds, blind = FALSE)

library(genefilter)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)

mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
annotation_col <- data.frame(
  condition = colData(vsd)$condition
)
rownames(annotation_col) <- colnames(mat)
dist_genes <- dist(mat)
hc_genes <- hclust(dist_genes, method = "complete")

pdf("Lab6_gene_dendrogram.pdf", width = 8, height = 6)
plot(hc_genes, labels = FALSE, main = "Hierarchical clustering of top variable genes")
dev.off()
mycols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)

pdf("Lab6_gene_heatmap.pdf", width = 8, height = 10)
pheatmap(
  mat,
  color             = mycols,
  show_rownames     = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col    = annotation_col
)
dev.off()
cor_mat <- cor(assay(vsd))

pdf("Lab6_sample_correlation_heatmap.pdf", width = 8, height = 6)
pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  cluster_cols   = FALSE
)
dev.off()

head(GO_df$Description, 10)

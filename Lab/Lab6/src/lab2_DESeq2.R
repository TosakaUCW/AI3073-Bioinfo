# 1. 加载数据
# check.names = FALSE 非常重要，它能保留列名中的空格和括号，防止变成点
rawCounts <- read.delim("E-GEOD-50760-raw-counts.tsv", check.names = FALSE)
sampleData <- read.delim("E-GEOD-50760-experiment-design.tsv", stringsAsFactors = FALSE, check.names = FALSE)

# --- 调试检查开始 ---
print("检查列名是否正确:")
print(colnames(sampleData)[14:15]) # 打印包含 Factor Value 的列名看看
# --- 调试检查结束 ---

# 2. 数据预处理
# 提取数值部分，设置行名
counts_matrix <- rawCounts[, 3:ncol(rawCounts)]
rownames(counts_matrix) <- rawCounts[, 1]

# 3. 处理样本信息
# 使用反引号 `` 包裹带有特殊字符的列名
target_col <- "Factor Value[biopsy site]"

# 筛选只包含 primary tumor 和 normal 的行
# 这里的 which 确保即使有 NA 也不会报错
target_samples <- sampleData[which(sampleData[[target_col]] %in% c("primary tumor", "normal")), ]

# --- 关键检查：如果这里行数为0，说明列名或者筛选值不对 ---
if (nrow(target_samples) == 0) {
  stop("错误：筛选后的样本数为 0！请检查 sampleData 中的列名或值是否包含 'primary tumor' 和 'normal'。")
} else {
  print(paste("成功筛选出样本数:", nrow(target_samples)))
}

# 对应筛选 count 矩阵
counts_matrix_subset <- counts_matrix[, target_samples$Run]

# 确保列名顺序一致
if (!all(colnames(counts_matrix_subset) == target_samples$Run)) {
  stop("错误：Count矩阵的列名与样本ID不匹配！")
}

# 4. 设置分组因子
target_samples$condition <- factor(target_samples[[target_col]])

# --- 调试检查：打印当前的因子水平 ---
print("当前的分组水平:")
print(table(target_samples$condition))

# 重新设置参考水平 (现在应该能成功了，因为我们确保了数据存在)
target_samples$condition <- relevel(target_samples$condition, ref = "normal")

print("第一步完成，数据准备就绪。")

library(DESeq2)

print("--- 开始运行 DESeq2 差异表达分析 ---")

# 1. 构建 DESeqDataSet 对象
# design = ~ condition 表示我们要根据 'condition' 列来寻找差异
dds <- DESeqDataSetFromMatrix(countData = counts_matrix_subset,
                              colData = target_samples,
                              design = ~ condition)

# 2. 运行 DESeq 主流程 (标准化、离散度估计、统计检验)
dds <- DESeq(dds)

# 3. 提取结果
# contrast 参数指定比较方向：分母是 normal，分子是 primary tumor
res <- results(dds, contrast = c("condition", "primary tumor", "normal"))

# 4. 按 p-value 从小到大排序
resOrdered <- res[order(res$pvalue), ]

# 5. 打印一些统计摘要
print("差异表达结果摘要:")
summary(res)

# 6. 保存结果到 CSV (Task 2 要求)
# 将结果转换为数据框，并将行名(Gene ID) 变成单独的一列
res_df <- as.data.frame(resOrdered)
res_df$GeneID <- rownames(res_df)
# 调整列顺序，把 GeneID 放在第一列
res_df <- res_df[, c("GeneID", colnames(res_df)[1:ncol(res_df)-1])]

write.csv(res_df, file = "DE_genes_by_pvalue.csv", row.names = FALSE)
print("已保存差异基因列表为: DE_genes_by_pvalue.csv")


library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因组注释包
library(ggplot2)

print("--- 开始运行富集分析 (Enrichment Analysis) ---")

# 1. 筛选显著差异基因
# 标准通常是：矫正后P值 (padj) < 0.05 且 差异倍数绝对值 (|log2FoldChange|) > 1
sig_res <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
sig_gene_ids <- rownames(sig_res)

print(paste("筛选出的显著差异基因数量:", length(sig_gene_ids)))

# 如果基因太少（比如少于10个），后面的分析可能会报错，做个判断
if (length(sig_gene_ids) < 5) {
  warning("显著差异基因过少，富集分析可能无法产生结果。")
}

# 2. ID 转换 (Ensembl -> Entrez ID)
# 富集分析通常需要 Entrez ID
gene_entrez <- bitr(sig_gene_ids, 
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# 3. 运行 GO 富集分析 (Biological Process)
ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",          # BP = Biological Process (生物学过程)
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)          # 将结果中的 ID 转换回基因名方便阅读

# 4. 保存富集结果表
if (!is.null(ego)) {
  write.csv(as.data.frame(ego), "Enrichment_GO_Results.csv")
  print("已保存富集分析结果表: Enrichment_GO_Results.csv")
  
  # 5. 生成可视化图片并保存 (PDF)
  # Task 3 要求描述结论并附图，这里我们生成气泡图和条形图
  pdf("Enrichment_Figures.pdf", width = 10, height = 8)
  
  # 气泡图
  p1 <- dotplot(ego, showCategory=20) + 
    ggtitle("GO Enrichment Analysis (Dotplot)")
  print(p1)
  
  # 条形图
  p2 <- barplot(ego, showCategory=20) + 
    ggtitle("GO Enrichment Analysis (Barplot)")
  print(p2)
  
  dev.off()
  print("已保存富集分析图片: Enrichment_Figures.pdf")
  
} else {
  print("未找到显著富集的通路。")
}

### 第四步：层次聚类与热图绘制 (Task 4)
library(pheatmap)
library(ggplot2)

print("--- 开始绘制热图 (Heatmap) ---")

# 1. 数据标准化 (VST: Variance Stabilizing Transformation)
# 这一步是为了消除测序深度的差异，让样本之间可比，适合做聚类
# 这里的 dds 是你在第二步生成的对象
vsd <- vst(dds, blind = FALSE)

# 2. 选择用于绘图的基因
# 我们只画差异最显著的前 50 个基因，否则热图会密密麻麻看不清
# resOrdered 是你在第二步生成的排序后的结果
topGenes <- head(rownames(resOrdered), 50)

# 提取这50个基因的标准化表达矩阵
mat <- assay(vsd)[topGenes, ]

# 对数据进行“中心化” (Row-centering)
# 这样颜色代表的是：相对于该基因平均水平，是高(红)还是低(蓝)
mat <- mat - rowMeans(mat)

# 3. 准备注释条 (Annotation)
# 这会让热图上方显示每个样本属于 Normal 还是 Tumor
annotation_col <- data.frame(
  Condition = target_samples$condition,
  row.names = colnames(vsd) # 确保行名是样本ID
)

# 4. 开始绘图并保存
pdf("Heatmap_Clustering.pdf", width = 10, height = 12)

# 绘制热图
pheatmap(mat, 
         annotation_col = annotation_col,
         show_rownames = TRUE,  # 显示基因ID (如果是Ensembl ID可能不太直观，但符合要求)
         show_colnames = FALSE, # 样本名太多，不显示以免拥挤
         main = "Heatmap of Top 50 DE Genes",
         clustering_distance_rows = "correlation", # 基因按相关性聚类
         clustering_distance_cols = "euclidean",   # 样本按欧氏距离聚类
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) # 蓝-白-红 配色

# 额外绘制：样本距离聚类树 (单独展示样本是否聚成两堆)
sampleDists <- dist(t(assay(vsd)))
plot(hclust(sampleDists), 
     main = "Sample Clustering Dendrogram", 
     xlab = "", sub = "",
     labels = target_samples$condition) # 把树的标签设为 Condition，方便看有没有混杂

dev.off()
print("已保存热图和聚类树: Heatmap_Clustering.pdf")

print("=== 所有分析任务全部完成！ ===")
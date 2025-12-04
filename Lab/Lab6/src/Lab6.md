# Lab 6

![GO enrichment dotplot](Lab6_GO_dotplot.pdf)

![GO enrichment barplot](Lab6_GO_barplot.pdf)

![Gene dendrogram](Lab6_gene_dendrogram.pdf)

![Gene expression heatmap](Lab6_gene_heatmap.pdf)

![Sample correlation heatmap](Lab6_sample_correlation_heatmap.pdf)

### 3. Gene set enrichment analysis

```R
head(GO_df$Description, 10)
```
```
 [1] "ribosome biogenesis"                 "rRNA metabolic process"
 [3] "rRNA processing"                     "DNA replication"
 [5] "DNA-templated DNA replication"       "RNA localization"
 [7] "mitotic nuclear division"            "cell cycle DNA replication"
 [9] "connective tissue development"       "regulation of chromosome separation"
```

The dominant functional themes that distinguish **primary tumor** from **normal** tissue are mainly related to **ribosome biogenesis, RNA metabolism, and cell cycle / DNA replication**. 

These enriched functions indicate that tumor samples show increased activity in protein synthesis machinery (ribosome and rRNA-related processes) and strong up-regulation of cell cycle and DNA replication pathways. 

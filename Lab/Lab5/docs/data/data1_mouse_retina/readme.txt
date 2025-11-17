Dataset: Mouse retina (Macosko et al.) 
The dataset from Macosko et al. 38 contains scRNA-seq data of 44808 cells (after pruning singletons) obtained 
through Drop-seq from mouse retina (14 days post-natal). The expression matrix was obtained from GEO (GSE63472), 
while the cluster information was obtained from the host laboratory webpage (http://mccarrolllab.com/dropseq/). We 
used the normalized expression matrix [given as log((UMI counts per gene in a cell/Total UMI counts in cell)*10000)+1)]. 

In order to reduce the computational cost of the analysis, the dataset was down-sampled into a smaller set in which 
all the given cell types are represented. In mouse retina, the majority of the cells are rods 71, which according to the 
authors, in this data set correspond to more than 29000 cells. Since rods are the smallest cell type in mouse retina 72 and 
express fewer genes, they also contain higher levels of noise. In order to take a representative sample not overtaken by 
the rods content, Macosko et al. selected cells which express more than 900 genes. We used this same down-sampling 
approach, which resulted in a selection of 11020 cells, to build the gene regulatory network. Running GENIE3 (and 
RcisTarget) on the 12953 genes with more than 55.1 normalized counts (0.5 in 1% of the population) and detected in 
more than ~55 cells (0.5% of the population). This network was then evaluated on all the cells in the dataset, which led to 
an activity matrix including 123 regulons.
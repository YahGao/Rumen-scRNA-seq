setwd("/yhgao/rumen/")

# Seurat2monocle
library(monocle)
library(Seurat)

BW_AW_rds <- readRDS('BW_AW_tutorial.rds')
data <- as(as.matrix(BW_AW_rds@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = BW_AW_rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# Quality Control
BW_AW_cds <- monocle_cds
BW_AW_cds <- estimateSizeFactors(BW_AW_cds)
BW_AW_cds <- estimateDispersions(BW_AW_cds)
BW_AW_cds <- detectGenes(BW_AW_cds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(BW_AW_cds), num_cells_expressed >= 10))

# Cluster
disp_table <- dispersionTable(BW_AW_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
BW_AW_cds <- setOrderingFilter(BW_AW_cds, unsup_clustering_genes$gene_id)
BW_AW_cds <- reduceDimension(BW_AW_cds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
BW_AW_cds <- clusterCells(BW_AW_cds, num_clusters = 2)

# Trajectory
diff_test_res <- differentialGeneTest(BW_AW_cds[expressed_genes,], fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) 
BW_AW_cds <- setOrderingFilter(BW_AW_cds, ordering_genes)
BW_AW_cds <- reduceDimension(BW_AW_cds, max_components = 2, method = 'DDRTree')
BW_AW_cds <- orderCells(BW_AW_cds)

tiff(file="BW_AW_fig4a.tiff", width = 14.5, height = 11, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot_cell_trajectory(BW_AW_cds, color_by = "seurat_clusters") +
  labs(x = "Componet 1", y="Componet 2") +
  theme(legend.text = element_text(size=15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(0.5,'cm'),
        legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

tiff(file="0-BW_AW_fig4b.tiff", width = 14.5, height = 11, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot_cell_trajectory(BW_AW_cds, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 2) +
  theme(legend.position = "none")
dev.off()

tiff(file="0-BW_AW_fig4c.tiff", width = 14, height = 11, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot_cell_trajectory(BW_AW_cds, color_by = "Pseudotime") +
  labs(x = "Componet 1", y="Componet 2") +
  theme(legend.text = element_text(size=15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(0.5,'cm'),
        legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

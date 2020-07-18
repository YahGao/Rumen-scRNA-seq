setwd("/yhgao/rumen/")

library(Seurat)
library(scRNAseq)
library(scater)
library(cowplot)

#==================#
#   Seurat Step    #
#==================#
# Setup the Seurat objects
BW.data <- Read10X(data.dir = "/yhgao/rumen/BW/")
BW <- CreateSeuratObject(counts = BW.data, project = "BW", min.cells = 3, min.features = 200)
BW[["percent.mt"]] <- PercentageFeatureSet(BW, pattern = "^MT-")
BW <- subset(BW, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 30)
BW <- NormalizeData(BW, normalization.method = "LogNormalize", scale.factor = 10000)
BW <- FindVariableFeatures(BW, selection.method = "vst", nfeatures = 2000)

AW.data <- Read10X(data.dir = "/yhgao/rumen/AW")
AW <- CreateSeuratObject(counts = AW.data, project = "AW", min.cells = 3, min.features = 200)
AW[["percent.mt"]] <- PercentageFeatureSet(AW, pattern = "^MT-")
AW <- subset(AW, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 30)
AW <- NormalizeData(AW, normalization.method = "LogNormalize", scale.factor = 10000)
AW <- FindVariableFeatures(AW, selection.method = "vst", nfeatures = 2000)

BW_AW <- merge(BW, y = AW, add.cell.ids = c("BW","AW"), project = "rumen")
BW_AW.list <- SplitObject(BW_AW, split.by = "orig.ident")
BW_AW.list <- lapply(X = BW_AW.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Perform integration
rumen.anchors <- FindIntegrationAnchors(object.list = BW_AW.list, dims = 1:20)
rumen.combined <- IntegrateData(anchorset = rumen.anchors, dims = 1:20)

# Perform an integrated analysis
DefaultAssay(rumen.combined) <- "integrated"
rumen.combined <- ScaleData(rumen.combined, verbose = FALSE)
rumen.combined <- RunPCA(rumen.combined, npcs = 30, verbose = FALSE)
rumen.combined <- RunUMAP(rumen.combined, reduction = "pca", dims = 1:10)
rumen.combined <- FindNeighbors(rumen.combined, reduction = "pca", dims = 1:10)
rumen.combined <- FindClusters(rumen.combined, resolution = 0.24)

#======================#
#   Cell_cycle Step    #
#======================#
cell_cycle_g1_s <- read.table("cell_cycle_g1_s.txt", header = F)
cell_cycle_g2_m <- read.table("cell_cycle_g2_m.txt", header = F)
cell_cycle_genes <- rbind(cell_cycle_g1_s, cell_cycle_g2_m)
cell_cycle_genes <- unique(cell_cycle_genes)
colnames(cell_cycle_genes) <- "genes"

# Step 1. By cluster
# Cell cycle heatmap by cluster
rumen.combined_umi <- as.matrix(rumen.combined@assays$RNA@data)
rumen.combined_umi <- as.data.frame(rumen.combined_umi)
rumen.combined_umi$genes <- rownames(rumen.combined_umi)

rumen.combined_umi_cycle <- merge(cell_cycle_genes, rumen.combined_umi)
rownames(rumen.combined_umi_cycle) <- rumen.combined_umi_cycle$genes
rumen.combined_umi_cycle <- rumen.combined_umi_cycle[,-1]
rumen.combined_umi_cycle <- t(rumen.combined_umi_cycle)
rumen.combined_umi_cycle <- as.data.frame(rumen.combined_umi_cycle)
rumen.combined_umi_cycle$cells <- rownames(rumen.combined_umi_cycle)

rumen.combined_cluster <- as.data.frame(rumen.combined$seurat_clusters)
rumen.combined_cluster$cells <- rownames(rumen.combined_cluster)
colnames(rumen.combined_cluster) <- c("cluster", "cells")

rumen.combined_umi_cycle_cluster <- merge(rumen.combined_cluster, rumen.combined_umi_cycle)
rumen.combined_umi_cycle_cluster <- rumen.combined_umi_cycle_cluster[order(rumen.combined_umi_cycle_cluster[,2]),]
rownames(rumen.combined_umi_cycle_cluster) <- rumen.combined_umi_cycle_cluster$cells
rumen.combined_umi_cycle_cluster <- rumen.combined_umi_cycle_cluster[,-1]
rumen.combined_umi_cycle_cluster1 <- as.data.frame(rumen.combined_umi_cycle_cluster$cluster)
rumen.combined_umi_cycle_cluster2 <- rumen.combined_umi_cycle_cluster[,-1]
rumen.combined_umi_cycle_cluster3 <- apply(rumen.combined_umi_cycle_cluster2, 1, mean)
rumen.combined_umi_cycle_cluster3 <- as.data.frame(rumen.combined_umi_cycle_cluster3)
rumen.combined_umi_cycle_cluster4 <- cbind(rumen.combined_umi_cycle_cluster1, rumen.combined_umi_cycle_cluster3)
colnames(rumen.combined_umi_cycle_cluster4) <- c("cluster", "mean_cc_genes")

sdata <- split(rumen.combined_umi_cycle_cluster4,rumen.combined_umi_cycle_cluster4$cluster)
result <- lapply(sdata,function(x) x[order(x[,2], decreasing = T),])
result0 <- result$`0`
result1 <- result$`1`
result2 <- result$`2`
result3 <- result$`3`
result4 <- result$`4`
result5 <- result$`5`
result_all <- rbind(result0,result1,result2,result3,result4,result5)
result_all$cells <- rownames(result_all)

rumen.combined_umi_cycle_cluster$cells <- rownames(rumen.combined_umi_cycle_cluster)
result_all_umi <- merge(result_all, rumen.combined_umi_cycle_cluster, sort = F)
result_all_umi_cluster <- data.frame(result_all_umi$cluster)
rownames(result_all_umi_cluster) <- result_all_umi$cells
colnames(result_all_umi_cluster) <- "cluster"

rownames(result_all_umi) <- result_all_umi$cells
result_all_umi <- result_all_umi[,-c(1:3)]
result_all_umi <- t(result_all_umi)
result_all_umi <- log2(result_all_umi + 1)
result_all_umi <- result_all_umi - 1
result_all_umi <- as.data.frame(result_all_umi)
result_all_umi$genes <- rownames(result_all_umi)
result_all_umi_cc <- merge(cell_cycle_genes, result_all_umi, sort = F)
rownames(result_all_umi_cc) <- result_all_umi_cc$genes
result_all_umi_cc <- result_all_umi_cc[,-1]

library(pheatmap)
tiff(file="BW_AW_fig2a.tiff", width = 24, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
pheatmap(result_all_umi_cc, show_colnames =F, show_rownames = T, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "red"))(100),
         annotation_col=result_all_umi_cluster,annotation_legend = F,annotation_names_col=F,
         fontsize_row = 5, legend_breaks = c(-1:1), legend_labels = c("-1","0","1"))
dev.off()


# Cell cycle index by cluster
cc_rumen.combined_umi <- merge(cell_cycle_genes, rumen.combined_umi)
rownames(cc_rumen.combined_umi) <- cc_rumen.combined_umi[,1]
cc_rumen.combined_umi <- cc_rumen.combined_umi[,-1]
cc_rumen.combined_umi <- t(cc_rumen.combined_umi)
cc_rumen.combined_umi <- as.data.frame(cc_rumen.combined_umi)
cc_rumen.combined_umi$cells <- rownames(cc_rumen.combined_umi)

rumen.combined_cluster <- as.data.frame(rumen.combined$seurat_clusters)
colnames(rumen.combined_cluster) <- "cluster"
rumen.combined_cluster$cells <- rownames(rumen.combined_cluster)
rumen.combined_cluster_cells <- merge(rumen.combined_cluster, cc_rumen.combined_umi)
rumen.combined_cluster_cells <- rumen.combined_cluster_cells[order(rumen.combined_cluster_cells[,2]),]
rownames(rumen.combined_cluster_cells) <- rumen.combined_cluster_cells$cells
rumen.combined_cluster_cells <- rumen.combined_cluster_cells[,-1]

rumen.combined_cluster_cells1 <- rumen.combined_cluster_cells
rumen.combined_cluster_cells1$mean <- rowMeans(rumen.combined_cluster_cells1[,2:93])
rumen.combined_cluster_cells1 <- rumen.combined_cluster_cells1[,-c(2:93)]
rumen.combined_cluster_cells2 <- mean(rumen.combined_cluster_cells1$mean)
rumen.combined_cluster_cells2_1 <- rep(rumen.combined_cluster_cells2, times = 6436)
rumen.combined_cluster_cells2_1 <- as.data.frame(rumen.combined_cluster_cells2_1)
rumen.combined_cluster_cells3 <- cbind(rumen.combined_cluster_cells1, rumen.combined_cluster_cells2_1)
colnames(rumen.combined_cluster_cells3) <- c("cluster", "mean", "mean1")
rumen.combined_cluster_cells3$judge <- rumen.combined_cluster_cells3$mean > rumen.combined_cluster_cells3$mean1
rumen.combined_cluster_cells4 <- aggregate(rumen.combined_cluster_cells3[,4], list(rumen.combined_cluster_cells3[,1]), table)
rumen.combined_cluster_cells4_1 <- c(27, 1618, 218, 158, 222, 167)  
rumen.combined_cluster_cells5 <- aggregate(rumen.combined_cluster_cells3[,1], list(rumen.combined_cluster_cells3[,1]), length)
rumen.combined_cluster_cells5 <- cbind(rumen.combined_cluster_cells5, rumen.combined_cluster_cells4_1)
colnames(rumen.combined_cluster_cells5) <- c("cluster", "total", "ture")
rumen.combined_cluster_cells5$index <- rumen.combined_cluster_cells5$ture/rumen.combined_cluster_cells5$total
write.table(rumen.combined_cluster_cells5, "BW_AW_fig2a.txt", quote = F, row.names = F, sep = "\t")


# Step 2. By ident
# Cell cycle heatmap by ident

rumen.combined_umi_cycle_cells <- rumen.combined_umi_cycle
rumen.combined_umi_cycle_cells <- rumen.combined_umi_cycle_cells[,-93]
rumen.combined_umi_cycle_cells1 <- apply(rumen.combined_umi_cycle_cells, 1, mean)
rumen.combined_umi_cycle_cells1 <- as.data.frame(rumen.combined_umi_cycle_cells1)
rumen.combined_umi_cycle_cells1$cells <- rownames(rumen.combined_umi_cycle_cells1)
colnames(rumen.combined_umi_cycle_cells1) <- c("cells_mean", "cells")

rumen.combined_ident <- as.data.frame(rumen.combined$orig.ident)
rumen.combined_ident$cells <- rownames(rumen.combined_ident)
colnames(rumen.combined_ident) <- c("ident", "cells")

rumen.combined_umi_cycle_cells2 <- merge(rumen.combined_ident, rumen.combined_umi_cycle_cells1)
rumen.combined_umi_cycle_cells3 <- merge(rumen.combined_umi_cycle_cells2, rumen.combined_umi_cycle)
rownames(rumen.combined_umi_cycle_cells3) <- rumen.combined_umi_cycle_cells3$cells
rumen.combined_umi_cycle_cells3 <- rumen.combined_umi_cycle_cells3[,-1]

sdata1 <- split(rumen.combined_umi_cycle_cells3,rumen.combined_umi_cycle_cells3$ident)
result1 <- lapply(sdata1,function(x) x[order(x[,2], decreasing = T),])
result1_BW <- result1$BW
result1_AW <- result1$AW
result1_combined <- rbind(result1_BW, result1_AW)
result1_combined_ident <- as.data.frame(result1_combined$ident)
rownames(result1_combined_ident) <- rownames(result1_combined)
colnames(result1_combined_ident) <- "ident"

result1_combined <- result1_combined[,-c(1:2)]
result1_combined <- t(result1_combined)
result1_combined <- as.data.frame(result1_combined)
result1_combined$genes <- rownames(result1_combined)

result1_combined <- merge(cell_cycle_genes, result1_combined, sort = F)
rownames(result1_combined) <- result1_combined$genes
result1_combined <- result1_combined[,-1]
result1_combined <- log2(result1_combined + 1)
result1_combined <- result1_combined -1

library(pheatmap)
tiff(file="BW_AW_figs2.tiff", width = 24, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
pheatmap(result1_combined, show_colnames =F, show_rownames = T, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "red"))(100),
         annotation_col=result1_combined_ident,annotation_legend = F,annotation_names_col=F,
         fontsize_row = 5, legend_breaks = c(-1:1), legend_labels = c("-1","0","1"))
dev.off()

# Cell cycle index by ident
rumen.combined_ident <- as.data.frame(rumen.combined$orig.ident)
colnames(rumen.combined_ident) <- "ident"
rumen.combined_ident$cells <- rownames(rumen.combined_ident)
rumen.combined_ident_cells <- merge(rumen.combined_ident, cc_rumen.combined_umi)
rumen.combined_ident_cells <- rumen.combined_ident_cells[order(rumen.combined_ident_cells[,2]),]
rownames(rumen.combined_ident_cells) <- rumen.combined_ident_cells$cells
rumen.combined_ident_cells <- rumen.combined_ident_cells[,-1]

rumen.combined_ident_cells1 <- rumen.combined_ident_cells
rumen.combined_ident_cells1$mean <- rowMeans(rumen.combined_ident_cells1[,2:94])
rumen.combined_ident_cells1 <- rumen.combined_ident_cells1[,-c(2:94)]
rumen.combined_ident_cells2 <- mean(rumen.combined_ident_cells1$mean)
rumen.combined_ident_cells2_1 <- rep(rumen.combined_ident_cells2, times = 6436)
rumen.combined_ident_cells2_1 <- as.data.frame(rumen.combined_ident_cells2_1)
rumen.combined_ident_cells3 <- cbind(rumen.combined_ident_cells1, rumen.combined_ident_cells2_1)
colnames(rumen.combined_ident_cells3) <- c("ident", "mean", "mean1")
rumen.combined_ident_cells3$judge <- rumen.combined_ident_cells3$mean > rumen.combined_ident_cells3$mean1
rumen.combined_ident_cells4 <- aggregate(rumen.combined_ident_cells3[,4], list(rumen.combined_ident_cells3[,1]), table)
rumen.combined_ident_cells4_1 <- c(1025, 3001)  
rumen.combined_ident_cells5 <- aggregate(rumen.combined_ident_cells3[,1], list(rumen.combined_ident_cells3[,1]), length)
rumen.combined_ident_cells5 <- cbind(rumen.combined_ident_cells5, rumen.combined_ident_cells4_1)
colnames(rumen.combined_ident_cells5) <- c("ident", "total", "ture")
rumen.combined_ident_cells5$index <- rumen.combined_ident_cells5$ture/rumen.combined_ident_cells5$total
write.table(rumen.combined_ident_cells5, "BW_AW_figs2.txt", quote = F, row.names = F, sep = "\t")

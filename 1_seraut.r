setwd("/yhgao/rumen/")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)

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

tiff(file="BW_AW_fig1a.tiff", width = 11, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(rumen.combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.5) +
  scale_color_hue(labels=c("AW","BW"))+
  labs(x = "", y="") +
  theme(legend.text = element_text(size=11),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width=unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=1.5)))  
dev.off()

tiff(file="BW_AW_fig1b.tiff", width = 9.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(rumen.combined, reduction = "umap", label = T, pt.size = 0.5) +
  labs(x = "", y="") +
  theme(legend.text = element_text(size=11),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=1.5)))
dev.off()

tiff(file="BW_AW_fig1c.tiff", width = 14.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(rumen.combined, reduction = "umap", split.by = "orig.ident", pt.size = 0.5) +
  labs(x = "", y="") +
  theme(legend.text = element_text(size=11),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=1.2)))
dev.off()


# Identify conserved cell type markers
rumen.combined <- ScaleData(rumen.combined, verbose = FALSE, assay = "RNA") 
DefaultAssay(rumen.combined) <- "RNA"
cluster0.markers <- FindConservedMarkers(rumen.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
cluster1.markers <- FindConservedMarkers(rumen.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
cluster2.markers <- FindConservedMarkers(rumen.combined, ident.1 = 2, grouping.var = "orig.ident", verbose = FALSE)
cluster3.markers <- FindConservedMarkers(rumen.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
cluster4.markers <- FindConservedMarkers(rumen.combined, ident.1 = 4, grouping.var = "orig.ident", verbose = FALSE)
cluster5.markers <- FindConservedMarkers(rumen.combined, ident.1 = 5, grouping.var = "orig.ident", verbose = FALSE)

saveRDS(rumen.combined, file = "BW_AW_tutorial.rds")


# Plot Fig. 3
tf <- c("MKI67", "HMMR", "EZH2", "BRCA2")
ep <- c("TGFBI", "TGFB2", "TGFBR2") # Epithelial
ke <- c("KRT6A", "KRT7", "KRT8", "KRT17", "KRT18", "KRT19") # Keratins
ud <- c("FOSB", "IRF2BP2", "IRF7", "FGF2", "CDH13") # up-down
all <- c(tf, ep, ke, ud)

tiff(file="BW_AW_fig3a.tiff", width = 24, height = 16, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DoHeatmap(rumen.combined, features = all, 
          slot = "scale.data", assay = "RNA") + 
  theme(legend.text = element_text(size=12),
        legend.key.size=unit(0.5,'cm'),
        legend.position = "right",
        axis.text.y=element_text(size=12))
dev.off()

tiff(file="BW_AW_fig3b.tiff", width = 30, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(rumen.combined, features = all, pt.size = 0, log = T, ncol = 6)
dev.off()

tiff(file="BW_AW_fig3c.tiff", width = 39, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(rumen.combined, features = all, ncol = 6)
dev.off()



# Plot top DEG
cluster0.markers$cluster <- "0"
cluster1.markers$cluster <- "1"
cluster2.markers$cluster <- "2"
cluster3.markers$cluster <- "3"
cluster4.markers$cluster <- "4"
cluster5.markers$cluster <- "5"
cluster0.markers$genes <- rownames(cluster0.markers)
cluster1.markers$genes <- rownames(cluster1.markers)
cluster2.markers$genes <- rownames(cluster2.markers)
cluster3.markers$genes <- rownames(cluster3.markers)
cluster4.markers$genes <- rownames(cluster4.markers)
cluster5.markers$genes <- rownames(cluster5.markers)
cluster0.markers$order <- c(1:length(rownames(cluster0.markers)))
cluster1.markers$order <- c(1:length(rownames(cluster1.markers)))
cluster2.markers$order <- c(1:length(rownames(cluster2.markers)))
cluster3.markers$order <- c(1:length(rownames(cluster3.markers)))
cluster4.markers$order <- c(1:length(rownames(cluster4.markers)))
cluster5.markers$order <- c(1:length(rownames(cluster5.markers)))

cluster.markers <- rbind(cluster0.markers, cluster1.markers, cluster2.markers, cluster3.markers,
                         cluster4.markers, cluster5.markers)

top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = order)
tiff(file="BW_AW_figs3a1.tiff", width = 24, height = 16, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DoHeatmap(rumen.combined, features = top10$genes, 
          slot = "scale.data", assay = "RNA") + 
  theme(legend.text = element_text(size=8),
        legend.key.size=unit(0.25,'cm'),
        legend.position = "right",
        axis.text.y=element_text(size=6))
dev.off()


top1 <- cluster.markers %>% group_by(cluster) %>% top_n(n = -1, wt = order)
tiff(file="BW_AW_figs3a2.tiff", width = 15, height = 10, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(rumen.combined, features = top1$genes, pt.size = 0, slot = "counts", log = T, ncol = 3)
dev.off()

tiff(file="BW_AW_figs3a3.tiff", width = 18, height = 10, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(rumen.combined, features = top1$genes, ncol = 3)
dev.off()


# Plot top TFs
tf <- c("ATF4", "BCLAF1", "BRCA1", "CEBPZ", "DDIT3", "E2F8", "EP300", "ETS1", "EZH2", 
	    "GATA2", "HMMR", "KDM5A", "MAFG", "MKI67", "POLR2A", "SMARCA4" ,"SREBF2", "YY1")

tiff(file="BW_AW_figs3b1.tiff", width = 24, height = 16, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DoHeatmap(rumen.combined, features = tf, 
          slot = "scale.data", assay = "RNA") + 
  theme(legend.text = element_text(size=8),
        legend.key.size=unit(0.25,'cm'),
        legend.position = "right",
        axis.text.y=element_text(size=9))
dev.off()

tiff(file="BW_AW_figs3b2.tiff", width = 30, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(rumen.combined, features = tf, pt.size = 0, log = T, ncol = 6)
dev.off()

tiff(file="BW_AW_figs3b3.tiff", width = 39, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(rumen.combined, features = tf, ncol = 6)
dev.off()



# Percentage of cells
rumen.combined_cluster <- rumen.combined$seurat_clusters
write.table(rumen.combined_cluster, "rumen.combined_cluster.txt", col.names = F, quote = F, sep = "\t")

rumen.combined1 <- read.table("BW_AW_combine.txt", header = T, row.names = 1) # Prepare the input file in EXCEL or R
colnames(rumen.combined1) <- c("BW", "AW")

tiff(file="BW_AW_fig1d.tiff", width = 4, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(as.matrix(100*rumen.combined1), legend = rownames(rumen.combined1), border = NA,
        col = c('#F8766D','#E38900','#C49A00','#99A800','#53B400'),
        cex.axis = 1, cex.names = 0.7, ylim = c(0, 100),
        las = 1, width = 0.8, space = 0.25, beside = FALSE, 
        args.legend = list(x = 'right', bty = 'n', inset = -1.2, cex = 1, y.intersp = -0.8, 
                           x.intersp = 0.7, text.width = 1.8, ncol=1))
mtext('Percentage of cells (%)', cex = 1, side = 2, line = 2.3)
dev.off()

setwd("/yhgao/rumen/")

library(Seurat)
library(SingleR)
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
DefaultAssay(rumen.combined) <- "RNA"


#===================#
#   SingleR Step    #
#===================#
# Annotation

blue.se <- BlueprintEncodeData()
BW_AW_blue <- SingleR(test = as.SingleCellExperiment(rumen.combined), ref = blue.se, labels = blue.se$label.fine)

rumen.combined$SingleR.pruned.calls <- BW_AW_blue$pruned.labels
rumen.combined$SingleR.calls <- BW_AW_blue$labels
to.remove <- pruneScores(BW_AW_blue)
new.pruned <- BW_AW_blue$labels
new.pruned[pruneScores(BW_AW_blue, nmads=5)] <- NA
all.markers <- metadata(BW_AW_blue)$de.genes
rumen.combined$labels_blue <- BW_AW_blue$labels

tiff(file="BW_AW_s1a.tiff", width = 20, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(rumen.combined, reduction = "umap", label = F, pt.size = 2, group.by = "labels_blue") +
  theme(legend.text = element_text(size=15),
        legend.position = "right",
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 0.6),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))
dev.off()

tiff(file="BW_AW_s1b.tiff", width = 25, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(rumen.combined, reduction = "umap", label = T, pt.size = 0.5, split.by = "labels_blue", 
                  ncol = 4, label.size = 4) +  NoLegend() +
  theme(legend.text = element_text(size=15),
        legend.position = "right",
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 0.3),
        axis.text.x = element_text(size = 10, angle = 0, hjust = 0.6),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))
dev.off()


# Percentage of each cell type
unique(rumen.combined$labels_blue)
 [1] "Mesangial cells"      "Epithelial cells"     "Erythrocytes"         "Keratinocytes"       
 [5] "Fibroblasts"          "Pericytes"            "MEP"                  "Myocytes"            
 [9] "mv Endothelial cells" "Chondrocytes"         "Astrocytes"           "CLP"    

Mesangial <- subset(rumen.combined, subset=labels_blue=="Epithelial cells")
Epithelial_ident_c <- data.frame(Epithelial$orig.ident, Epithelial$seurat_clusters)
write.table(Epithelial_ident_c, "BW_AW_mesangial.txt", quote = F, sep = "\t")

Epithelial <- subset(rumen.combined, subset=labels_blue=="Epithelial cells")
Epithelial_ident_c <- data.frame(Epithelial$orig.ident, Epithelial$seurat_clusters)
write.table(Epithelial_ident_c, "BW_AW_epithelial.txt", quote = F, sep = "\t")

Erythrocytes <- subset(rumen.combined, subset=labels_blue=="Erythrocytes")
Erythrocytes_ident_c <- data.frame(Erythrocytes$orig.ident, Erythrocytes$seurat_clusters)
write.table(Erythrocytes_ident_c, "BW_AW_erythrocytes.txt", quote = F, sep = "\t")

Keratinocytes <- subset(rumen.combined, subset=labels_blue=="Keratinocytes")
Keratinocytes_ident_c <- data.frame(Keratinocytes$orig.ident, Keratinocytes$seurat_clusters)
write.table(Keratinocytes_ident_c, "BW_AW_Keratinocytes.txt", quote = F, sep = "\t")

Fibroblasts <- subset(rumen.combined, subset=labels_blue=="Fibroblasts")
Fibroblasts_ident_c <- data.frame(Fibroblasts$orig.ident, Fibroblasts$seurat_clusters)
write.table(Fibroblasts_ident_c, "BW_AW_Fibroblasts.txt", quote = F, sep = "\t")

Pericytes <- subset(rumen.combined, subset=labels_blue=="Pericytes")
Pericytes_ident_c <- data.frame(Pericytes$orig.ident, Pericytes$seurat_clusters)
write.table(Pericytes_ident_c, "BW_AW_Pericytes.txt", quote = F, sep = "\t")

MEP <- subset(rumen.combined, subset=labels_blue=="MEP")
MEP_ident_c <- data.frame(MEP$orig.ident, MEP$seurat_clusters)
write.table(MEP_ident_c, "BW_AW_MEP.txt", quote = F, sep = "\t")

Myocytes <- subset(rumen.combined, subset=labels_blue=="Myocytes")
Myocytes_ident_c <- data.frame(Myocytes$orig.ident, Myocytes$seurat_clusters)
write.table(Myocytes_ident_c, "BW_AW_Myocytes.txt", quote = F, sep = "\t")

Chondrocytes <- subset(rumen.combined, subset=labels_blue=="Chondrocytes")
Chondrocytes_ident_c <- data.frame(Chondrocytes$orig.ident, Chondrocytes$seurat_clusters)
write.table(Chondrocytes_ident_c, "BW_AW_Chondrocytes.txt", quote = F, sep = "\t")

Astrocytes <- subset(rumen.combined, subset=labels_blue=="Astrocytes")
Astrocytes_ident_c <- data.frame(Astrocytes$orig.ident, Astrocytes$seurat_clusters)
write.table(Astrocytes_ident_c, "BW_AW_Astrocytes.txt", quote = F, sep = "\t")

CLP <- subset(rumen.combined, subset=labels_blue=="CLP")
CLP_ident_c <- data.frame(CLP$orig.ident, CLP$seurat_clusters)
write.table(CLP_ident_c, "BW_AW_CLP.txt", quote = F, sep = "\t")

mv_Endothelial <- subset(rumen.combined, subset=labels_blue=="mv Endothelial cells")
mv_Endothelial_ident_c <- data.frame(mv_Endothelial$orig.ident, mv_Endothelial$seurat_clusters)
write.table(mv_Endothelial_ident_c, "BW_AW_mv_Endothelial.txt", quote = F, sep = "\t")

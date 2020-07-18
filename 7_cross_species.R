setwd("/yhgao/rumen/")

#===============#
#  Seurat step  #
#===============#
library(Seurat)
library(Matrix)

# Reading files
HAStomach1.data <- read.table("GSM4008651_Adult-Stomach1_dge.txt.gz", row.names = 1, header = T)
HAStomach2.data <- read.table("GSM4008652_Adult-Stomach2_dge.txt.gz", row.names = 1, header = T)
HAStomach3_1.data <- read.table("GSM4008653_Adult-Stomach3-1_dge.txt.gz", row.names = 1, header = T)
HAStomach3_2.data <- read.table("GSM4008654_Adult-Stomach3-2_dge.txt.gz", row.names = 1, header = T)
HAStomach3_3.data <- read.table("GSM4008655_Adult-Stomach3-3_dge.txt.gz", row.names = 1, header = T)
HFStomach1.data <- read.table("GSM4008712_Fetal-Stomach1_dge.txt.gz", row.names = 1, header = T)
HFStomach2.data <- read.table("GSM4008713_Fetal-Stomach2_dge.txt.gz", row.names = 1, header = T)

BW.data <- Read10X(data.dir = "/yhgao/rumen/BW/")
AW.data <- Read10X(data.dir = "/yhgao/rumen/AW/")
BW <- CreateSeuratObject(counts = BW.data, project = "BW", min.cells = 3, min.features = 200)
AW <- CreateSeuratObject(counts = AW.data, project = "AW", min.cells = 3, min.features = 200)

# QC
HAStomach1_500less <- HAStomach1.data[,colSums(HAStomach1.data) < 500 & colSums(HAStomach1.data) > 100]
HAStomach1_500more <- HAStomach1.data[,colSums(HAStomach1.data) >=500]
HAStomach2_500less <- HAStomach2.data[,colSums(HAStomach2.data) < 500 & colSums(HAStomach2.data) > 100]
HAStomach2_500more <- HAStomach2.data[,colSums(HAStomach2.data) >=500]
HAStomach3_1_500less <- HAStomach3_1.data[,colSums(HAStomach3_1.data) < 500 & colSums(HAStomach3_1.data) > 100]
HAStomach3_1_500more <- HAStomach3_1.data[,colSums(HAStomach3_1.data) >=500]
HAStomach3_2_500less <- HAStomach3_2.data[,colSums(HAStomach3_2.data) < 500 & colSums(HAStomach3_2.data) > 100]
HAStomach3_2_500more <- HAStomach3_2.data[,colSums(HAStomach3_2.data) >=500]
HAStomach3_3_500less <- HAStomach3_3.data[,colSums(HAStomach3_3.data) < 500 & colSums(HAStomach3_3.data) > 100]
HAStomach3_3_500more <- HAStomach3_3.data[,colSums(HAStomach3_3.data) >=500]
HFStomach1_500less <- HFStomach1.data[,colSums(HFStomach1.data) < 500 & colSums(HFStomach1.data) > 100]
HFStomach1_500more <- HFStomach1.data[,colSums(HFStomach1.data) >=500]
HFStomach2_500less <- HFStomach2.data[,colSums(HFStomach2.data) < 500 & colSums(HFStomach2.data) > 100]
HFStomach2_500more <- HFStomach2.data[,colSums(HFStomach2.data) >=500]

## Common genes
HAStomach1_genes <- as.data.frame(rownames(HAStomach1_500more))
colnames(HAStomach1_genes) <- "genes"
HAStomach2_genes <- as.data.frame(rownames(HAStomach2_500more))
colnames(HAStomach2_genes) <- "genes"
HAStomach3_1_genes <- as.data.frame(rownames(HAStomach3_1_500more))
colnames(HAStomach3_1_genes) <- "genes"
HAStomach3_2_genes <- as.data.frame(rownames(HAStomach3_2_500more))
colnames(HAStomach3_2_genes) <- "genes"
HAStomach3_3_genes <- as.data.frame(rownames(HAStomach3_3_500more))
colnames(HAStomach3_3_genes) <- "genes"
HFStomach1_genes <- as.data.frame(rownames(HFStomach1_500more))
colnames(HFStomach1_genes) <- "genes"
HFStomach2_genes <- as.data.frame(rownames(HFStomach2_500more))
colnames(HFStomach2_genes) <- "genes"

BW_umi <- as.data.frame(BW@assays$RNA@data)
BW_genes <- as.data.frame(rownames(BW_umi))
colnames(BW_genes) <- "genes"
AW_umi <- as.data.frame(AW@assays$RNA@data)
AW_genes <- as.data.frame(rownames(AW_umi))
colnames(AW_genes) <- "genes"

Stomach_rumen1 <- merge(HAStomach1_genes, HAStomach2_genes)
Stomach_rumen2 <- merge(Stomach_rumen1, HAStomach3_1_genes)
Stomach_rumen3 <- merge(Stomach_rumen2, HAStomach3_2_genes)
Stomach_rumen4 <- merge(Stomach_rumen3, HAStomach3_3_genes)
Stomach_rumen5 <- merge(Stomach_rumen4, HFStomach1_genes)
Stomach_rumen6 <- merge(Stomach_rumen5, HFStomach2_genes)
Stomach_rumen7 <- merge(Stomach_rumen6, BW_genes)
Stomach_rumen8 <- merge(Stomach_rumen7, AW_genes)


HAStomach1_500more$genes <- rownames(HAStomach1_500more)
HAStomach1_500more_common <- merge(Stomach_rumen8, HAStomach1_500more)
rownames(HAStomach1_500more_common) <- HAStomach1_500more_common$genes
HAStomach1_500more_common <- HAStomach1_500more_common[,-1]
HAStomach2_500more$genes <- rownames(HAStomach2_500more)
HAStomach2_500more_common <- merge(Stomach_rumen8, HAStomach2_500more)
rownames(HAStomach2_500more_common) <- HAStomach2_500more_common$genes
HAStomach2_500more_common <- HAStomach2_500more_common[,-1]
HAStomach3_1_500more$genes <- rownames(HAStomach3_1_500more)
HAStomach3_1_500more_common <- merge(Stomach_rumen8, HAStomach3_1_500more)
rownames(HAStomach3_1_500more_common) <- HAStomach3_1_500more_common$genes
HAStomach3_1_500more_common <- HAStomach3_1_500more_common[,-1]
HAStomach3_2_500more$genes <- rownames(HAStomach3_2_500more)
HAStomach3_2_500more_common <- merge(Stomach_rumen8, HAStomach3_2_500more)
rownames(HAStomach3_2_500more_common) <- HAStomach3_2_500more_common$genes
HAStomach3_2_500more_common <- HAStomach3_2_500more_common[,-1]
HAStomach3_3_500more$genes <- rownames(HAStomach3_3_500more)
HAStomach3_3_500more_common <- merge(Stomach_rumen8, HAStomach3_3_500more)
rownames(HAStomach3_3_500more_common) <- HAStomach3_3_500more_common$genes
HAStomach3_3_500more_common <- HAStomach3_3_500more_common[,-1]
HFStomach1_500more$genes <- rownames(HFStomach1_500more)
HFStomach1_500more_common <- merge(Stomach_rumen8, HFStomach1_500more)
rownames(HFStomach1_500more_common) <- HFStomach1_500more_common$genes
HFStomach1_500more_common <- HFStomach1_500more_common[,-1]
HFStomach2_500more$genes <- rownames(HFStomach2_500more)
HFStomach2_500more_common <- merge(Stomach_rumen8, HFStomach2_500more)
rownames(HFStomach2_500more_common) <- HFStomach2_500more_common$genes
HFStomach2_500more_common <- HFStomach2_500more_common[,-1]

BW_umi$genes <- rownames(BW_umi)
BW_umi_common <- merge(Stomach_rumen8, BW_umi)
rownames(BW_umi_common) <- BW_umi_common$genes
BW_umi_common <- BW_umi_common[,-1]
AW_umi$genes <- rownames(AW_umi)
AW_umi_common <- merge(Stomach_rumen8, AW_umi)
rownames(AW_umi_common) <- AW_umi_common$genes
AW_umi_common <- AW_umi_common[,-1]


HAStomach1 <- CreateSeuratObject(Matrix(as.matrix(HAStomach1_500more_common), sparse=T), names.delim = "\\.", project = "HAStomach1")
HAStomach2 <- CreateSeuratObject(Matrix(as.matrix(HAStomach2_500more_common), sparse=T), names.delim = "\\.", project = "HAStomach2")
HAStomach3_1 <- CreateSeuratObject(Matrix(as.matrix(HAStomach3_1_500more_common), sparse=T), names.delim = "\\.", project = "HAStomach3_1")
HAStomach3_2 <- CreateSeuratObject(Matrix(as.matrix(HAStomach3_2_500more_common), sparse=T), names.delim = "\\.", project = "HAStomach3_2")
HAStomach3_3 <- CreateSeuratObject(Matrix(as.matrix(HAStomach3_3_500more_common), sparse=T), names.delim = "\\.", project = "HAStomach3_3")
HFStomach1 <- CreateSeuratObject(Matrix(as.matrix(HFStomach1_500more_common), sparse=T), names.delim = "\\.", project = "HFStomach1")
HFStomach2 <- CreateSeuratObject(Matrix(as.matrix(HFStomach2_500more_common), sparse=T), names.delim = "\\.", project = "HFStomach2")
BW <- CreateSeuratObject(Matrix(as.matrix(BW_umi_common), sparse=T), names.delim = "\\.", project = "BW")
AW <- CreateSeuratObject(Matrix(as.matrix(AW_umi_common), sparse=T), names.delim = "\\.", project = "AW")

Stomach_rumen <- merge(HAStomach1, y =c(HAStomach2, HAStomach3_1, HAStomach3_2, HAStomach3_3,
                                        HFStomach1, HFStomach2, BW, AW),
                       add.cell.ids = c("HAStomach1","HAStomach2","HAStomach3_1","HAStomach3_2", 
                                        "HAStomach3_3","HFStomach1","HFStomach2","BW", "AW"),
                       project = "merge")

Stomach_rumen.list <- SplitObject(Stomach_rumen, split.by = "orig.ident")
for (i in 1:length(Stomach_rumen.list)) {
  Stomach_rumen.list[[i]] <- NormalizeData(Stomach_rumen.list[[i]], verbose = FALSE)
  Stomach_rumen.list[[i]] <- FindVariableFeatures(Stomach_rumen.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

Stomach_rumen.anchors <- FindIntegrationAnchors(object.list = Stomach_rumen.list, dims = 1:30)
Stomach_rumen.combined <- IntegrateData(anchorset = Stomach_rumen.anchors, dims = 1:30)
DefaultAssay(Stomach_rumen.combined) <- "integrated"

Stomach_rumen.combined <- ScaleData(Stomach_rumen.combined, verbose = FALSE)
Stomach_rumen.combined <- FindVariableFeatures(object = Stomach_rumen.combined)
Stomach_rumen.combined <- RunPCA(Stomach_rumen.combined, npcs = 30, verbose = FALSE)
Stomach_rumen.combined <- RunUMAP(Stomach_rumen.combined, reduction = "pca", dims = 1:10)
Stomach_rumen.combined <- FindNeighbors(Stomach_rumen.combined, reduction = "pca", dims = 1:10)
Stomach_rumen.combined <- FindClusters(Stomach_rumen.combined, resolution = 0.4)


pdf("BW_AW_figs5a1.pdf",width=6,height=4)
DimPlot(Stomach_rumen.combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)
dev.off()

pdf("BW_AW_figs5a2.pdf",width=7,height=6)
DimPlot(Stomach_rumen.combined, reduction = "umap", label = T, pt.size = 0.5, label.size = 8)
dev.off()

pdf("BW_AW_figs5a3.pdf",width=8,height=6)
DimPlot(Stomach_rumen.combined, reduction = "umap", split.by = "orig.ident", pt.size = 0.5, ncol = 3, label = T)
dev.off()


#================#
#  SingleR step  #
#================#
library(SingleR)
library(scRNAseq)
library(scater)
library(cowplot)

blue.se <- BlueprintEncodeData()
Stomach_rumen.combined_blue <- SingleR(test = as.SingleCellExperiment(Stomach_rumen.combined), ref = blue.se,labels = blue.se$label.fine)

Stomach_rumen.combined$SingleR.pruned.calls <- Stomach_rumen.combined_blue$pruned.labels
Stomach_rumen.combined$SingleR.calls <- Stomach_rumen.combined_blue$labels
to.remove <- pruneScores(Stomach_rumen.combined_blue)
new.pruned <- Stomach_rumen.combined_blue$labels
new.pruned[pruneScores(Stomach_rumen.combined_blue, nmads=5)] <- NA
all.markers <- metadata(Stomach_rumen.combined_blue)$de.genes
Stomach_rumen.combined$labels_blue <- Stomach_rumen.combined_blue$labels

pdf("BW_AW_figs5b1.pdf",width=15,height=12)
DimPlot(Stomach_rumen.combined, reduction = "umap", label = F, pt.size = 2, group.by = "labels_blue") +
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

pdf("BW_AW_figs5b2.pdf",width=12,height=10)
DimPlot(Stomach_rumen.combined, reduction = "umap", label = T, pt.size = 0.5, split.by = "labels_blue", 
                  ncol = 6, label.size = 4) +  NoLegend() +
  theme(legend.text = element_text(size=15),
      legend.position = "right",
      legend.key.size = unit(0.5, 'lines'),
      legend.key.width = unit(0.01,'cm'),
      plot.title = element_text(hjust = 0.5, size = 0.3),
      axis.text.x = element_text(size = 16, angle = 0, hjust = 0.6),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size=16,face="bold"),
      axis.title.y = element_text(size=13,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))
dev.off()

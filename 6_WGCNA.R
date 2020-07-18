setwd("/yhgao/rumen/")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)

#=============#
# Seurat step #
#=============#
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

rumen.anchors <- FindIntegrationAnchors(object.list = BW_AW.list, dims = 1:20)
rumen.combined <- IntegrateData(anchorset = rumen.anchors, dims = 1:20)
DefaultAssay(rumen.combined) <- "integrated"
rumen.combined <- ScaleData(rumen.combined, verbose = FALSE)
rumen.combined <- RunPCA(rumen.combined, npcs = 30, verbose = FALSE)
rumen.combined <- RunUMAP(rumen.combined, reduction = "pca", dims = 1:10)
rumen.combined <- FindNeighbors(rumen.combined, reduction = "pca", dims = 1:10)
rumen.combined <- FindClusters(rumen.combined, resolution = 0.24)
DefaultAssay(rumen.combined) <- "RNA"


#============#
# WGCNA step #
#============#
# Preprocess
datadf <- as.matrix(rumen.combined@assays$RNA@data )
idd1 <- rumen.combined@meta.data
Inter.id1 <- cbind(rownames(idd1),idd1$seurat_clusters)
rownames(Inter.id1) <- rownames(idd1)
colnames(Inter.id1) <- c("CellID","Celltype")
Inter.id1 <- as.data.frame(Inter.id1)
Inter1 <- datadf[,Inter.id1$CellID]
Inter2 <- as.matrix(Inter1)

pseudocell.size = 100 
new_ids_list1 = list()
length(levels(Inter.id1$Celltype))

for (i in 1:length(levels(Inter.id1$Celltype))) {
  cluster_id = levels(Inter.id1$Celltype)[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)     
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)    
  new_ids_list1[[i]] <- pseudo_ids      
}

new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)

new_colnames <- rownames(new_ids)
gc()
all.data <- datadf[, as.character(new_colnames)]
all.data <- t(all.data)
new.data <- aggregate(list(all.data[, 1:length(all.data[1,])]), list(name = new_ids[,1]), FUN=mean)
rownames(new.data) <- new.data$name
new.data <- new.data[, -1]
new_ids_length <- as.matrix(new_ids_length)

short <- which(new_ids_length < 100)
new_good_ids <- as.matrix(new_ids_length[-short,])
result <- t(new.data)[, rownames(new_good_ids)]

rumen.combined <- FindVariableFeatures(rumen.combined, nfeatures = 2000)
Cluster1 <- result[intersect(Seurat::VariableFeatures(rumen.combined), rownames(result)),]


# The main step of WGCNA
# Preprocess
library(WGCNA)

type <- "unsigned"
corType <- "pearson"
corFnc <- ifelse(corType == "pearson", cor, bicor)
maxPOutliers <- ifelse(corType == "pearson", 1, 0.05)
robustY <- ifelse(corType == "pearson", T, F)
dataExpr <- as.matrix(Cluster1)

m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
dim(dataExpr)

gsg <- goodSamplesGenes(dataExpr, verbose = 3)

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

# Selection of softPower
powers <- c(c(1:10), seq(from = 12, to =30, by = 2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, networkType = "signed", verbose = 5)
power <- sft$powerEstimate
softPower <- power

# Construction of network
cor <- WGCNA::cor
net <- blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                        TOMType = "unsigned", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs=TRUE, corType = corType, 
                        maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                        saveTOMFileBase = paste0("dataExpr", ".tom"),
                        verbose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

tiff(file="BW_AW_fig5a.tiff", width = 12, height = 5, units = "cm", res = 600, pointsize = 8, compression= "lzw")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Correlation between each module and each cluster
MEDissThres <- 0.25
MEs <- net$MEs
MEs_col <- MEs
library(stringr)
colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col <- orderMEs(MEs_col)
merge <- mergeCloseModules(dataExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

load(net$TOMFiles, verbose = T)
TOM <- as.matrix(TOM)
dissTOM <- 1-TOM
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
table(moduleColors)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], main = "Network heatmap plot, all genes") 

myrumen <- CreateSeuratObject(Cluster1)
mergedMEs_name <- rownames(mergedMEs)
mergedMEs_name <- substr(mergedMEs_name, 1,1)  
mergedMEs_name <- as.data.frame(mergedMEs_name)
mergedMEs <- cbind(mergedMEs_name, mergedMEs)

for (i in 1:6) {
  mergedMEs_all <- mergedMEs[mergedMEs_name==i,]
  mergedMEs_all <- mergedMEs_all[,-1]
  myrumen_all <- subset(myrumen, idents = i)
  moduleTraitCor_noFP1 <- cor(mergedMEs_all, myrumen_all@meta.data$nFeature_RNA, use = "p")
  moduleTraitPvalue_noFP1 <- corPvalueStudent(moduleTraitCor_noFP1, nSamples)
  textMatrix_noFP1 <- paste(signif(moduleTraitCor_noFP1, 2), "\n(", signif(moduleTraitPvalue_noFP1, 1), ")", sep = "")
  moduleTraitCor_noFP2 <- cor(mergedMEs_all, myrumen_all@meta.data$nCount_RNA, use = "p")
  moduleTraitPvalue_noFP2 <- corPvalueStudent(moduleTraitCor_noFP2, nSamples)
  textMatrix_noFP2 <- paste(signif(moduleTraitCor_noFP2, 2), "\n(", signif(moduleTraitPvalue_noFP2, 1), ")", sep = "")
  moduleTraitCor_noFP <- cbind(moduleTraitCor_noFP1, moduleTraitCor_noFP2)
  moduleTraitCor_noFP <- as.data.frame(moduleTraitCor_noFP)
  colnames(moduleTraitCor_noFP) <- c("nFeature_RNA", "nCount_RNA")
  textMatrix_noFP <- cbind(textMatrix_noFP1, textMatrix_noFP2)
  tiff(file=paste("BW_AW_moudle",i-1,".tiff", sep = ""), width = 10, height = 10, units = "cm", res = 600, pointsize = 8,compression= "lzw")
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor_noFP, 
                 xLabels = colnames(moduleTraitCor_noFP), 
                 yLabels = names(mergedMEs_all), 
                 ySymbols = names(mergedMEs_all), 
                 colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = textMatrix_noFP,
                 setStdMargins = FALSE, 
                 cex.text = 1, 
                 zlim = c(-1,1), 
                 main = paste("c",i-1, sep="")) 
  dev.off()
}

# Export the genes information
moduleColors_uniq <- unique(moduleColors)
for (color in moduleColors_uniq) {
  module <- color
  probes <- colnames(dataExpr)
  inModule <- (moduleColors == module)
  modProbes <- probes[inModule]
  write.table(modProbes, paste("module_gene_", module, ".txt", sep = ""), quote = F, row.names = F, col.names = F)
}


#=========================#
#   Enrichment analysis   #
#=========================#
library(clusterProfiler)
library(org.Bt.eg.db)

moduleColors_uniq

module <- "yellow"
module_gene <- read.table(paste("module_gene_", module, ".txt", sep = ""), header = F)
module_gene$V1 <- as.character(module_gene$V1)
module_gene_id <- bitr(module_gene$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Bt.eg.db")
module_gene_list <- module_gene_id$ENTREZID
module_gene_list[duplicated(module_gene_list)]
module_gene_go <- enrichGO(module_gene_list, OrgDb = org.Bt.eg.db, ont='ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')

tiff(file = paste("module_gene_", module, "_go.tiff", sep = ""), width = 40, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
dotplot(module_gene_go, showCategory=50, font.size =20) + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18))
dev.off()

setwd("/yhgao/rumen/")

library(dplyr)
library(Seurat)
library(patchwork)
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


#==================#
#   SCENIC Step    #
#==================#
# Step 1. Prepare input file
cells.use <- WhichCells(object = rumen.combined) 
exprMat <- GetAssayData(object = rumen.combined, assay= "RNA", slot = "data")[, cells.use] 
exprMat <- as(Class = 'matrix', object = exprMat)
gene <- rumen.combined$nFeature_RNA
umi <- rumen.combined$nCount_RNA
clus <- rumen.combined$seurat_clusters
cellInfo <- data.frame(gene, umi, clus)

# Step 2. Initialization
cellTypeColumn <- "clus"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
colVars <- list(CellType=c("0"="forestgreen", "1"="darkorange", "2"="magenta4", "3"="hotpink", "4"="red3", "5"="skyblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

library(SCENIC)
org="hgnc" 
dbDir="/project/dor_rna/CJL/Project_Li2/Yahui/scRNA/SCENIC/cisTarget_databases/hg19_mc9nr" # RcisTarget databases location
myDatasetTitle="SCENIC example on Human" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Step 3. Construction of gene coexpression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)

# Step 4. Construct and calculate gene regulation network
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Step 5. Convert the network activity to the format of ON/OFF
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions) # Set the threshold

nPcs <- c(5)
scenicOptions@settings$seed <- 123 # same seed for all of them
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#==================
# heatmap of Step 5
library(AUCell)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
thresholds <- loadInt(scenicOptions, "aucell_thresholds")
thresholds <- getThresholdSelected(thresholds)

regulonsCells <- setNames(lapply(names(thresholds), function(x) {
  trh <- thresholds[x]
  names(which(getAUC(regulonAUC)[x, ] > trh))
}), names(thresholds))

regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[, 1], regulonActivity[, 2]))
class(binaryRegulonActivity) <- "matrix"
regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists = "null", verbose = F)

library(NMF)
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists = "null")
cellInfo <- data.frame(cellInfo)
colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists = "null")

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}


regulon_label1 <- "all"

for (selRegs in regulon_label1) {
  if (length(regulonSelection[[selRegs]]) > 0) {
    regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in%  rownames(binaryRegulonActivity))]
    binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]], , drop = F]
    if (nrow(binaryMat) > 0) {
      fileName <- paste0(getOutName(scenicOptions, "s4_binaryActivityHeatmap"), selRegs)
      fileName <- .openDevHeatmap(fileName = fileName, devType = getSettings(scenicOptions, "devType"))
      rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
      colv <- ifelse(ncol(binaryMat) >= 2, T, NA)
      NMF::aheatmap(binaryMat, scale = "none", 
                    revC = TRUE, main = selRegs, annCol = cellInfo[colnames(binaryMat),  , drop = F], annColor = colVars, Rowv = rowv, 
                    Colv = colv, color = c("white", "black"), filename = fileName)
      if (getSettings(scenicOptions, "devType") != 
          "pdf") 
        dev.off()
    }}}

BW_AW_cluster <- as.data.frame(rumen.combined$seurat_clusters)
BW_AW_cluster$cells <- rownames(BW_AW_cluster)
colnames(BW_AW_cluster) <- c("clusters", "cells")
BW_AW_cluster <- BW_AW_cluster[order(BW_AW_cluster[,1]),]
binaryMat_col <- as.data.frame(colnames(binaryMat))
colnames(binaryMat_col) <- "cells"
binaryMat_col1 <- merge(BW_AW_cluster, binaryMat_col, sort = F)
binaryMat_col1 <- binaryMat_col1[order(binaryMat_col1[,2]),]

tiff(file="BW_AW_fig2b.tiff", width = 22, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
aheatmap(binaryMat, scale = "none", revC = T, main = selRegs, labCol = NA,
         annCol = cellInfo[binaryMat_col1$cells, ,drop = T], annColor = colVars, Rowv = rowv, 
         Colv = NA, color = c("white", "black"))
dev.off()

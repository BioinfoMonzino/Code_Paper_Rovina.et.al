rm(list = ls())
cat("\014")
graphics.off()

#load libraries
library(Seurat)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(limma)
library(DaMiRseq)
library(ggrepel)
library(plotly)

#import Data
setwd("<workdir>")
source("./custom_seurat_functions.R")

npcs_redux <- 30

set.seed(1)

dmd_mtx <- ReadMtx(mtx = "./DMD_filtered_feature_bc_matrix/matrix.mtx",
                   cells = "./DMD_filtered_feature_bc_matrix/barcodes.tsv",
                   features = "./DMD_filtered_feature_bc_matrix/genes.tsv",)

iso_mtx <- ReadMtx(mtx = "./ISO_filtered_feature_bc_matrix/matrix.mtx",
                   cells = "./ISO_filtered_feature_bc_matrix/barcodes.tsv",
                   features = "./ISO_filtered_feature_bc_matrix/genes.tsv")

dmd_seu <- CreateSeuratObject(counts = dmd_mtx,
                              min.cells = 30,
                              min.features = 500, project = "DMD")

iso_seu <- CreateSeuratObject(counts = iso_mtx,
                              min.cells = 30,
                              min.features = 500, project = "ISO")
rm(dmd_mtx)
rm(iso_mtx)

gc()

###########     QC    #######
# before cell filtering
dmd_seu[["percent.mt"]]  <- PercentageFeatureSet(dmd_seu, pattern = "^MT")
iso_seu[["percent.mt"]]  <- PercentageFeatureSet(iso_seu, pattern = "^MT")

VlnPlot(dmd_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(iso_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# afer cell filtering
dmd_seu <- subset(dmd_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
iso_seu <- subset(iso_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

VlnPlot(dmd_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(iso_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

############## Normalize separate slots (samples)
integr_list <- list()
integr_list[["dmd_seu"]] <- dmd_seu
integr_list[["iso_seu"]] <- iso_seu

integr_list <- lapply(X = integr_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
  x <- ScaleData(x, verbose = F, vars.to.regress = "percent.mt")
  x <- RunPCA(x, npcs = npcs_redux, verbose = F)
#  x <- RunUMAP(x, reduction = "pca", dims = 1:npcs_redux, verbose = F)
  x <- RunUMAP(x, reduction = "pca", dims = 1:npcs_redux, verbose = F,n.components = 3)
  #x <- RunTSNE(x, tsne.method = "Rtsne", dims = 1:npcs_redux, verbose = F)
  
})

########################### Data Integration
features <- SelectIntegrationFeatures(object.list = integr_list)
integr_anchors <- FindIntegrationAnchors(object.list = integr_list,
                                         dims = 1:npcs_redux,
integr_seurat  <- IntegrateData(anchorset = integr_anchors, dims = 1:npcs_redux)

# change Assay and integrate
DefaultAssay(integr_seurat) <- "integrated"
integr_seurat <- ScaleData(integr_seurat, verbose = F)
integr_seurat <- RunPCA(integr_seurat, npcs = npcs_redux, verbose = F)
integr_seurat <- RunUMAP(integr_seurat, reduction = "pca", dims = 1:npcs_redux, verbose = F,n.components=3)

# chose best k
ElbowPlot(integr_seurat,ndims = npcs_redux)

# best k = 20
integr_seurat <- FindNeighbors(integr_seurat, dims = 1:npcs_redux, k.param = 20, verbose = F)
integr_seurat <- FindClusters(integr_seurat, verbose = F)


DimPlot(integr_seurat,reduction = "umap",label = T,pt.size = 1,label.box = T) + NoLegend()
DimPlot(integr_seurat,label = T,pt.size = 0.5,split.by = "orig.ident") + NoLegend()


## ################### Cell ferquency by cluster ################ 
count_table <- as.data.frame(table(integr_seurat@meta.data$seurat_clusters,
                                   integr_seurat@meta.data$orig.ident))

count_tableDMD <- count_table[which(count_table$Var2 %in% "DMD"),]
count_tableISO <- count_table[which(count_table$Var2 %in% "ISO"),]

count_tableDMD$Freq_Perc <- round(100*count_tableDMD$Freq/sum(count_tableDMD$Freq),1)
count_tableISO$Freq_Perc <- round(100*count_tableISO$Freq/sum(count_tableISO$Freq),1)

count_table <- (rbind(count_tableDMD,count_tableISO))

ggplot(data=count_table, aes(x=Var1, y=Freq_Perc, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge())

# cl 0
prop.test(x = c(1120, 307), n = c(5573, 4238))
# cl 1
prop.test(x = c(848, 467), n = c(5573, 4238))
# cl 2
prop.test(x = c(424, 719), n = c(5573, 4238))
# cl 3
prop.test(x = c(245, 651), n = c(5573, 4238))
# cl 4
prop.test(x = c(173, 691), n = c(5573, 4238))
# cl 5
prop.test(x = c(506, 234), n = c(5573, 4238))
# cl 6
prop.test(x = c(669, 39), n = c(5573, 4238))
# cl 7
prop.test(x = c(316, 299), n = c(5573, 4238))
# cl 8
prop.test(x = c(512, 86), n = c(5573, 4238))
# cl 9
prop.test(x = c(29, 416), n = c(5573, 4238))
# cl 10
prop.test(x = c(306, 39), n = c(5573, 4238))
# cl 11
prop.test(x = c(241, 43), n = c(5573, 4238))
# cl 12
prop.test(x = c(146, 97), n = c(5573, 4238))
# cl 13
prop.test(x = c(38, 150), n = c(5573, 4238))

# save results
save.image("<working_dir>/res.RData")


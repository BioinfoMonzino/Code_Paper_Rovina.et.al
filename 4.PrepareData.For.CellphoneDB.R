rm(list = ls())
cat("\014")
graphics.off()

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

setwd("<working_dir>")
source("./custom_seurat_functions.R")

load("<working_dir>/res.RData")


###############################################################################
gc()

integr_seurat <- RenameIdents(integr_seurat,
                              `0` = "Cardiomyocyte", 
                              `1` = "Fibroblast", 
                              `2` = "Cardiomyocyte", 
                              `3` = "SMC",
                              `4` = "SMC",
                              `5` = "Cardiomyocyte",
                              `6` = "Cardiac_Neuron", 
                              `7` = "Fibroblast",
                              `8` = "Cardiomyocyte", 
                              `9` = "Cardiac_Neuron",
                              `10` = "Endothelial_Cell", 
                              `11` = "Endothelial_Cell", 
                              `12` = "Progenitor_Fibrobast",  
                              `13` = "Progenitor_Cell"
                              #
                              
)

integr_seurat$cell_type <- Idents(integr_seurat)
integr_seurat$celltype.group <- paste(Idents(integr_seurat), 
                                      integr_seurat$orig.ident, 
                                      sep = "_")

metadata_seu<- integr_seurat@meta.data

##############
# take raw data and normalise it
CM_FB_DMD <-which(metadata_seu[,7] %in% c("Cardiomyocyte","Fibroblast") & metadata_seu[,1] %in% "DMD")
CM_FB_ISO <-which(metadata_seu[,7] %in% c("Cardiomyocyte","Fibroblast") & metadata_seu[,1] %in% "ISO")

# Seurat Norm
seu_norm_CM_FB_DMD <- integr_seurat@assays$RNA@data[, CM_FB_DMD] 
seu_norm_CM_FB_ISO <- integr_seurat@assays$RNA@data[, CM_FB_ISO] 
gc()
####################################### Write Output
############ DMD
write.table(round(seu_norm_CM_FB_DMD,2), 
            "CM_FB.DMD.seu.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

metadata_seu_DMD <- metadata_seu[CM_FB_DMD,]

write.table(metadata_seu_DMD, 
            "CM_FB.DMD.Metadata.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

############ ISO
write.table(round(seu_norm_CM_FB_ISO,2), 
            "CM_FB.ISO.seu.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

metadata_seu_ISO <- metadata_seu[CM_FB_ISO,]

write.table(metadata_seu_ISO, 
            "CM_FB.ISO.Metadata.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

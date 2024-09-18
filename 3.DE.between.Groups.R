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
library(scCustomize)
library(viridisLite)

setwd("<working_dir>")
source("./custom_seurat_functions.R")

load("<working_dir>/res.RData")

###############################################################################

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

metadata_seu<- integr_seurat@meta.data
integr_seurat <- AddMetaData(object = integr_seurat, metadata = metadata_seu[,6,drop=F])

##########################################################################à
###### differential analysis by condition
integr_seurat$celltype.group <- paste(Idents(integr_seurat), 
                                      integr_seurat$orig.ident, 
                                      sep = "_")
integr_seurat$celltype <- Idents(integr_seurat)
Idents(integr_seurat) <- "celltype.group"

df_data_counts <- integr_seurat@assays$RNA@counts
cluster_groups <- levels(integr_seurat@meta.data$celltype)

# DE: for each Cell-type
for (j in 1:length(cluster_groups)){
  cat(paste0(cluster_groups[j],"\n"))
  gc()
  check_cluster <- cluster_groups[j]
  
  index_cl <- which(integr_seurat@meta.data$celltype %in% check_cluster)
  df_cl <- as.data.frame(t(df_data_counts[,index_cl]))
  meta_cl <- integr_seurat@meta.data[index_cl,c("celltype.group","orig.ident")]
  data_meta <- cbind(meta_cl,df_cl)
  meta_cl_rnd <- data_meta[,c(1,2)]
  data_cl_rnd <- data_meta[,-c(1,2)]
  colnames(meta_cl_rnd) <- c("group","class")
  
  # DaMiRseq
  SE<-DaMiR.makeSE(t(data_cl_rnd), meta_cl_rnd)
  data_norm <- DaMiR.normalization(SE, type = "logcpm" , minCounts=1,fSample = 0.05, hyper = "no")
  
  # DE: Limma
  design <- model.matrix(~0+class, meta_cl_rnd)
  contrasts <- makeContrasts(dmd.vs.iso=classDMD-classISO,
                             levels=design)
  fit <- lmFit(assay(data_norm), design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  res_table <- topTable(fit2,coef = 1,adjust="BH",number=10000000000000000000)
  
  # Genrate Volcano Plot
  res_table_plot <- res_table
  
  p_th <- 0.05 #10^-2
  adjp_th <- 10^-5
  Beta_th <- 0
  
  
  sign.lfc.p.adj.vect <- sign.lfc.p.vect <-sign.vect.adj.p <- sign.vect.p <- rep("NOT.SIGN", length(rownames(res_table_plot)))
  sign.vect.p[which(res_table_plot$P.Value<=p_th)] <- "SIGN"
  sign.lfc.p.vect[which(res_table_plot$P.Value<=p_th & res_table_plot$logFC >= Beta_th)] <- "SIGN_UP"
  sign.lfc.p.vect[which(res_table_plot$P.Value<=p_th & res_table_plot$logFC <= -Beta_th)] <- "SIGN_DOWN"
  sign.lfc.p.vect[which(res_table_plot$adj.P.Val<=adjp_th & res_table_plot$logFC >= Beta_th)] <- "SIGN_UPUP"
  sign.lfc.p.vect[which(res_table_plot$adj.P.Val<=adjp_th & res_table_plot$logFC <= -Beta_th)] <- "SIGN_DOWNDOWN"
  res_table_plot$sign.p <- as.factor(sign.vect.p)
  res_table_plot$sign.adj.p <- as.factor(sign.vect.adj.p)
  res_table_plot$sign.lfc.p <- as.factor(sign.lfc.p.vect)
  res_table_plot$sign.lfc.adj.p <- as.factor(sign.lfc.p.adj.vect)
  res_table_plot$name <- rownames(res_table_plot)
  
  if(length(levels(res_table_plot$sign.lfc.p)) ==3){
    gg1 <- ggplot(aes(x=logFC, y=-log10(P.Value),color=sign.lfc.p.vect),data=res_table_plot, label=name) +
      geom_point(aes(alpha=as.numeric(sign.p)),size=5)+
      scale_color_manual(values=c("gray50","lightblue","salmon")) +
      geom_hline(yintercept = -log10(p_th),color='gray50',linetype = "dashed") + 
      geom_vline(xintercept = Beta_th,color='gray50',linetype = "dashed") +
      geom_vline(xintercept = -Beta_th,color='gray50',linetype = "dashed") +
      theme(legend.position = "none") +
      xlim(-3,3) +
      ylim(0,70) +
      xlab("log2 (FC)")
    
    
  }else{
    # plot Volcano plot 
    gg1 <- ggplot(aes(x=logFC, y=-log10(P.Value),color=sign.lfc.p.vect),data=res_table_plot, label=name) +
      geom_point(aes(alpha=as.numeric(sign.p)),size=5)+
      scale_color_manual(values=c("gray50","lightblue","blue","salmon","red")) +
      geom_hline(yintercept = -log10(p_th),color='gray50',linetype = "dashed") + 
      geom_vline(xintercept = Beta_th,color='gray50',linetype = "dashed") +
      geom_vline(xintercept = -Beta_th,color='gray50',linetype = "dashed") +
      geom_hline(yintercept = -log10(res_table_plot$P.Value[max(which(res_table_plot$adj.P.Val<adjp_th))]),color='gray50',linetype = "dashed") + 
      theme(legend.position = "none") +
      xlab("log2 (FC)") #+
  }
  
  # Write Output
  pdf(file = paste0("<working_dir>/clust.",cluster_groups[j],".DMD.vs.ISO.pdf"),width = 10,height = 10)
  print(gg1)
  print(hist(res_table[,4], main="histogram of pValues", breaks=20, col="orange", xlab="p-value"))
  dev.off()
  write.table(res_table, paste0("<working_dir>/clust.",cluster_groups[j],".DMD.vs.ISO.txt"),sep="\t",quote=F,col.names = NA)
}







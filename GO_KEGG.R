############ GO/KEGG#########################################################
library(GO.db)
library(GenomeInfoDbData)
library(HDO.db)
library(clusterProfiler)
library(AnnotationHub) 
library(org.Hs.eg.db) 
library(org.Mm.eg.db) 
library(ggplot2)
library(DOSE)
library(ggplot2)
library(tidyverse)
library(dplyr)


# Run_GO_KEGG, plot_GO_bar, plot_GO_bilateral_bar
source("GO_KEGG_mathors.R")

# marker
setwd("E:\\实验室项目相关\\HCC_scRNAseq\\04_Fibroblast_Cell_Cluster")
DEGs_TAMsubsets2 <- read.csv("Fibro_cluster_Markers_Findall.csv")

table(DEGs_TAMsubsets2$cluster)
colnames(DEGs_TAMsubsets2)
dir.create('GO KEGG') 

########### cluster GO_KEGG #####################################################################
unique(DEGs_TAMsubsets2$cluster)
sel_clus = unique(DEGs_TAMsubsets2$cluster)
sel_clus
updown = T 
wz = "hsa"     # hsa；mmu,
for (i in 1:length(sel_clus)){
  clu = sel_clus[i]
  GO_KEGG_result = Run_GO_KEGG(DEGs_TAMsubsets2,  
                               wz = wz,     # hsa；mmu,
                               sel_clus = clu, 
                               updown = updown)     
  write.csv(GO_KEGG_result,file=paste0('GO KEGG/Cluster ', clu, ' GO and KEGG result.csv'))
 }

ONTOLOGYs = c("BP", "MF", "CC", "KEGG")
x_value = "Count" # "Count", "GeneRatio", "pvalue", "qvalue", or "p.adjust"
tile_names = sel_clus


all_results <- list()
for (i in 1:length(sel_clus)){
  tile_name = tile_names[i]
  clu = sel_clus[i]
  GO_KEGG_result = read.csv(paste0('GO KEGG/Cluster ', clu, ' GO and KEGG result.csv'))
  GO_KEGG_result$cluster <- clu
  all_results[[i]] <- GO_KEGG_result
  
  for (j in 1:length(ONTOLOGYs)){
    if (updown){
      p = plot_GO_bilateral_bar(GO_KEGG_result, tile_name, ONTOLOGYs[j], x_value = x_value)
      ggsave(paste0('GO KEGG/Cluster ', clu, " ", ONTOLOGYs[j], ' top20.pdf'),p,width = 9,height = 6.5)
    }else{
      p = plot_GO_bar(GO_KEGG_result, tile_name, ONTOLOGYs[j], x_value = x_value, col = "#C50B05")  # col = "#C50B05", top_n = 10
      ggsave(paste0('GO KEGG/Cluster ', clu, " ", ONTOLOGYs[j], ' top10.pdf'),p,width = 5.5,height = 4)
    }
  }
}

all_results_df <- do.call(rbind, all_results)
write.csv(all_results_df, file = "GO KEGG/All_Clusters_GO_KEGG_Result.csv", row.names = FALSE)








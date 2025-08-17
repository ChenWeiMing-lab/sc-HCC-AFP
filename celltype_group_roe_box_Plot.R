# 导入数据 -------------
rm(list = ls())
gc()
library(dplyr)
library(Matrix)
library(Seurat)
library(tidydr)
library(tidyverse)
setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")
surival_data = read.csv("./data/surival_data.csv", row.names = 1)
head(surival_data)
dim(surival_data)
table(surival_data$kmeans_group)
table(surival_data$CEA_Group)

# AFP分两组，分界点20
head(surival_data$AFP_status)
surival_data$AFP_group <- ifelse(surival_data$AFP_μg_L > 20, "High_Risk", "Low_Risk")
table(surival_data$AFP_status)

# # AFP分四组，分界点20，200，400
# head(surival_data$AFP_μg_L)
# surival_data$AFP_4group <- cut(surival_data$AFP_μg_L, breaks=c(-Inf, 20, 200, 400, Inf), 
#                                labels=c("Below_20", "20_to_200", "200_to_400", "Above_400"), 
#                                right=FALSE)
# table(surival_data$AFP_4group)

# 不同细胞类型 -------------
# Tumor  B cells  T_cells  CD4+ T cells  CD8+ T cells  Plasma cells  Mesenchymal cells
# NK cells   Macrophage   DC   Endothelial  Mast cells    Monocyte
setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")
cell_type = c("T_cells")
Hep = readRDS(paste0("/home/chenweiming/Project/HCC_scRNAseq/luo/data/HCC_",gsub(" ", "_", cell_type), ".rds"))
# Hep = readRDS("/home/chenweiming/Project/HCC_scRNAseq/luo/data/HCC_scRNA.rds")
Hep 
head(Hep@meta.data)
table(Hep@meta.data$Seu_Clusters_str)
table(Hep@meta.data$AFP_status)
# Hep@meta.data$CA199_status = paste0("statu_", Hep@meta.data$CA199_status)
Hep@meta.data$AFP_status_group = ifelse(Hep@meta.data$AFP_status == 1, "AFP_Pos", "AFP_Neg")

table(Hep@meta.data$AFP_status_group)
col_name = "AFP_status_group" #  AFP_4group CA199_status CEA_Group kmeans_group
# barcode = rownames(Hep@meta.data)
# # 合并分组信息
# Hep@meta.data <- Hep@meta.data %>%
#   left_join(surival_data[, c("Sample", col_name)], by = "Sample")
# rownames(Hep@meta.data) = barcode
# table(Hep@meta.data[[col_name]])
# table(Hep@meta.data$Seu_Clusters_str)
# sum(is.na(Hep@meta.data[[col_name]]))
# head(Hep@meta.data)

# cols = c("#db5024", "#25acda")
cols = c("#25acda", "#db5024")
# cols = c("#00b8a9","#ffdd7e","#f6416c","#1965B0")
p1 <- DimPlot(Hep, group.by = col_name, cols = cols, raster=FALSE) + # "#25acda", "#db5024"
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
  theme(panel.grid = element_blank()) +
  ggtitle(col_name) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),  # 使标题居中并加粗
    legend.text = element_text(size = 13),  # 调整图例文本大小
    legend.title = element_text(size = 15)  # 调整图例标题大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # 调整图例的整体尺寸
p1
ggsave(paste0(cell_type, "/", col_name, "_UMAP.pdf"), plot = p1, width = 6, height = 4.5)


# ROE -------
source("ROE_Heatmap_Function.R")
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:7]) # brewer.pal(9, "YlOrRd")[1:6] 

plot.data <- ROIE(table(Hep@meta.data[,c("Seu_Clusters_str", col_name)])) %>%  
  reshape2::melt() %>%  
  mutate(value = pmin(value, 3))    
head(plot.data)

plot.data <- plot.data %>%
  mutate(Var1 = factor(Var1, levels = as.character(unique(Var1[Var2 == "AFP_Pos" ][order(-value[Var2 == "AFP_Pos" ])]))))  
head(plot.data)
levels(plot.data$Var1)
levels(plot.data$Var2)

p=ggplot(plot.data, aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() + 
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_gradientn(colours = myPalette(100)) +
  scale_x_discrete(limits = c("AFP_Pos", "AFP_Neg")) +  
  labs(x = "", y = "",fill="Ro/e") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))
p
ggsave(paste0(cell_type, "/", col_name, "_RO_E_Celltype_Heatmap.pdf"), plot = p, device = "pdf",width = 8.5,height = 14,units = "cm")


col_name

table(Hep$Seu_Clusters_str)
table(Hep@meta.data[[col_name]])
table(Hep$Seu_Clusters_str, Hep@meta.data[[col_name]])  
sum(is.na(Hep@meta.data[[col_name]]))
meta <- Hep@meta.data %>% dplyr::select(orig.ident, col_name, Seu_Clusters_str)
colnames(meta) = c("orig.ident", "Group", "celltype") 
head(meta)


library(dplyr)
meta_percent <- meta %>%
  group_by(orig.ident, Group, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%          
  group_by(orig.ident) %>%
  mutate(values = count / sum(count)) %>%               
  dplyr::select(orig.ident, values, Group, celltype)           

meta_percent <- unique(meta_percent)
head(meta_percent)
sum(is.na(meta_percent))
levels(meta_percent$Group)

source("/home/chenweiming/Project/HCC_scRNAseq/luo/code/CellType_Group_BoxPlot_Function.R")
# cols = rev(c("#00b8a9","#ffdd7e","#f6416c","#1965B0"))
cols = rev(c("#327db7","#f18c8d")) #  "#9dd5cc","#327db7","#f18c8d","#ffdd7e"
p1 <- Groups_box_plot(
  data = meta_percent,
  Group_box_colors = cols,
  Group_top = "celltype", Group_top_order = levels(plot.data$Var1),
  Group_box = "Group", 
  Group_box_order = c("AFP_Pos", "AFP_Neg"),
  title_name = "Cell percentages across Samples", 
  sig_lable = T#, test_method = "t.test"  #wilcox.test 
)
p1
ggsave(paste0(cell_type, "/", col_name, "_Celltype_Boxplot.pdf"), plot = p1, device = "pdf",width = 18,height = 10,units = "cm") 





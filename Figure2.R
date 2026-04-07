# ===============================Figure2======================================

# ---- library ----
library(Seurat)
library(dplyr)
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)

# ---- load data ----
load("PancancerLEC.Rda")

# ---- Heatmap ----
markers <- FindAllMarkers(LEC,
                          logfc.threshold = 0.6,
                          min.pct = 0.25,
                          only.pos = T)
head(markers)
sig_markers <- markers %>%
  group_by(cluster)%>%
  top_n(n = 4, wt = avg_log2FC)
head(sig_markers)

genes <- unique(sig_markers$gene) 
aver_dt <- AverageExpression(LEC,
                             features = genes,
                             group.by = 'celltype',
                             slot = 'data') 
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[genes, ]
aver_dt[1:6,1:6]

gene_anno <- data.frame(gene_anno = sig_markers$cluster,
                        row.names = sig_markers$gene)

cell_anno <- data.frame(cell_anno = colnames(aver_dt),
                        row.names = colnames(aver_dt))
head(gene_anno);head(cell_anno)
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno,
         annotation_row = gene_anno) 

c4a_gui()
celltype_col <- c4a('cold', 6)
names(celltype_col) <- cell_anno$cell_anno

anno_col <- list(cell_anno = celltype_col,
                 gene_anno = celltype_col)
anno_col
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col  = cell_anno, 
         annotation_row = gene_anno,
         annotation_colors = anno_col, 
         color = mycol, 
         border_color = 'white') 
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
aver_dtt <- t(scale(t(aver_dt)))
cols <- c4a('cold', 6)
names(cols) <- cell_anno$cell_anno
cell <- data.frame(colnames(aver_dtt))
colnames(cell) <- 'cell'
col_anno <- HeatmapAnnotation(df = cell,
                              show_annotation_name = F,
                              gp = gpar(col = 'white', lwd = 2),
                              col = list(cell = cols))
Heatmap(aver_dtt,
        name = 'expression',
        col = mycol2,
        cluster_columns = F,
        cluster_rows = F,
        column_names_side = c('top'), 
        column_names_rot = 60, 
        row_names_gp = gpar(fontsize = 12, fontface = 'italic'), 
        rect_gp = gpar(col = "white", lwd = 1.5),
        top_annotation = col_anno) + row_anno
cellwidth = 1
cellheight = 0.5
cn = dim(aver_dtt)[2]
rn = dim(aver_dtt)[1]
w = cellwidth*cn
h = cellheight*rn
Heatmap(aver_dtt,
        width = unit(w, "cm"),
        height  = unit(h, "cm"),
        name = 'expression', 
        col = mycol2,
        cluster_columns = F,
        cluster_rows = F,
        column_names_side  = c('top'), 
        column_names_rot = 30, 
        row_title = 'Markers', 
        rect_gp = gpar(col = "white", lwd = 1.5),
        heatmap_legend_param = list(legend_height = unit(2.8, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        top_annotation = col_anno,
        row_split = rep(LETTERS[1:6], each = 4), #对行切片，each根据量调整
        row_gap = unit(2, "mm"),
        border = T,
        border_gp = gpar(col = "black", lwd = 1.2)) +
  row_anno



# ---- Feature Plot ----
FeaturePlot(LEC, features = c("LYVE1"),pt.size = 0.5) +
  scale_color_gradientn(colors = c("#fef9ef", "#f9dbb4","#f1b782","#e67c54","#d13f2d","#a2292d","#7c2124"))

FeaturePlot(LEC, features = c("VIM"),pt.size = 0.5) +
  scale_color_gradientn(colors = c("#f0faff", "#b7ddf2","#83bfe1","#529fd1","#2c7ac1","#135ba8","#08306b"))



# ===============================Supplemental Figure3======================================

# ---- library ----
library(Seurat)
library(data.table)
library(AnnoProbe)
library(infercnv)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(scales)
library(dplyr)
library(purrr)

# ---- load data ----
load("PancancerLEC.Rda")
load("Pancancer.Rda")

# ---- Cellchat ----
options(stringsAsFactors = FALSE)  
color.use <- c( 
  "T/NK_cell"                = "#F061DD",
  "B_cell"                   = "#c0a7e4",
  "Plasma"                   = "#A55194",
  "Epithelial"               = "#a2292d",
  "Myeloid"                  = "#17BECF",
  "Mast_cell"                = "#e67c54",
  "Fibroblast"               = "#98DF8A",
  "Endothelial"              = "#FF9896",
  "other LEC"                = "#F7B6D2",
  "LYVE-1+ VIM+ LEC"         = "#00A087",
  "Neurons"                  = "#74adc0"
)

color.use <- scPalette(nlevels(con))
names(color.use) <- levels(con)
DimPlot(con, reduction = "umap", pt.size = 0.5, label = TRUE, cols = color.use)
data.input <- GetAssayData(con, layer = "data")  
meta <- data.frame(
  labels = Idents(con),       
  samples = con$orig.ident,   
  row.names = colnames(con)  
)

unique(meta$labels)  
unique(meta$samples) 
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels",  
                           datatype = "RNA")
CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB 
showDatabaseCategory(CellChatDB)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.1, thresh.fc = 0, thresh.p = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, nboot = 100)
cellchat <- aggregateNet(cellchat)
df.net <- subsetCommunication(cellchat)
head(df.net) 
cellchat <- computeCommunProbPathway(cellchat)
df.net_P <- subsetCommunication(cellchat, slot.name = "netP")
head(df.net_P) 

#chordDiagram
celltype_colors <- c(
  "T/NK_cell"                = "#F061DD",
  "B_cell"                   = "#c0a7e4",
  "Plasma"                   = "#A55194",
  "Epithelial"               = "#a2292d",
  "Myeloid"                  = "#17BECF",
  "Mast_cell"                = "#e67c54",
  "Fibroblast"               = "#98DF8A",
  "Endothelial"              = "#FF9896",
  "other LEC"                = "#F7B6D2",
  "LYVE-1+ VIM+ LEC"         = "#00A087",
  "Neurons"                  = "#74adc0"
)

mat <- cellchat@net$count
mat_mtlec <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat_mtlec["LYVE-1+ VIM+ LEC", ] <- mat["LYVE-1+ VIM+ LEC", ]

df <- as.data.frame(as.table(mat_mtlec))
colnames(df) <- c("source", "target", "value")
df <- df[df$value > 0, ]

raw_max <- max(df$value)
upper_bound <- ceiling(raw_max / 20) * 20

new_colors <- c("#2c7ac1", "#83bfe1", "#b7ddf2",
                "#f9dbb4", "#f1b782", "#d13f2d")
breaks <- seq(0, upper_bound, length.out = length(new_colors))
col_fun <- colorRamp2(breaks, new_colors)
df$edge_col <- col_fun(df$value)

all_nodes <- union(df$source, df$target)
node_colors <- celltype_colors[all_nodes]
node_colors[is.na(node_colors)] <- "gray"

mat_plot <- with(df, tapply(value, list(source, target), sum))
mat_plot[is.na(mat_plot)] <- 0

circos.clear()
circos.par(gap.degree = 5, track.margin = c(0.02, 0.02))

chordDiagram(
  x = mat_plot,
  grid.col = node_colors,
  transparency = 0.1,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  directional = 1,
  direction.type = c("arrows"),
  diffHeight = -0.04,
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  col = df$edge_col
)
lgd_col <- Legend(
  col_fun = col_fun,
  at = c(0, upper_bound),
  labels = c("0", as.character(upper_bound)),
  title = "Interaction strength",
  direction = "vertical",
  legend_height = unit(3, "cm")
)
lgd_celltype <- Legend(
  labels = names(celltype_colors),
  legend_gp = gpar(fill = celltype_colors),
  title = "Cell Types",
  direction = "vertical",
  legend_height = unit(6, "cm"),
  nrow = 6,
  grid_height = unit(0.4, "cm"),
  grid_width = unit(0.8, "cm")
)

draw(lgd_col, x = unit(1, "npc") - unit(2, "cm"), y = unit(0.75, "npc"))
draw(lgd_celltype, x = unit(1, "npc") - unit(2, "cm"), y = unit(0.35, "npc"))

# ---- inferCNV ----
my_sub <- "Epithelial"
sub.cells <- subset(con, idents = my_sub)
sub.cells <- sub.cells |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = 1:15) |>
  FindClusters(resolution = 0.2) |>     
  RunUMAP(dims = 1:15)
DimPlot(sub.cells, reduction = "umap", label = TRUE)

epi_sampled <- sub.cells
fibro <- subset(con, subset = celltype == "Fibroblast")
plasma <- subset(con, subset = celltype == "Plasma")

set.seed(1234)
fibro_sample <- sample(colnames(fibro), min(800, ncol(fibro)))
plasma_sample <- sample(colnames(plasma), min(800, ncol(plasma)))
epiMat <- GetAssayData(epi_sampled, assay = "RNA", layer = "counts")
fibroMat <- GetAssayData(fibro[, fibro_sample], assay = "RNA", layer = "counts")
plasmaMat <- GetAssayData(plasma[, plasma_sample], assay = "RNA", layer = "counts")
common_genes <- Reduce(intersect, list(rownames(epiMat), rownames(fibroMat), rownames(plasmaMat)))
expr_mat <- cbind(epiMat[common_genes, ], fibroMat[common_genes, ], plasmaMat[common_genes, ])
write.table(as.data.frame(expr_mat), "expFile.txt", sep = "\t", quote = FALSE)

groupinfo <- data.frame(
  v1 = colnames(expr_mat),
  v2 = c(rep("epi", ncol(epiMat)), rep("ref-1", ncol(fibroMat)), rep("ref-2", ncol(plasmaMat)))
)
write.table(groupinfo, "groupFiles.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

gene_info <- annoGene(common_genes, "SYMBOL", "human") %>%
  filter(!duplicated(SYMBOL)) %>%
  arrange(chr, start) %>%
  select(SYMBOL, chr, start, end)
write.table(gene_info, "geneFile.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "expFile.txt",
  annotations_file = "groupFiles.txt",
  gene_order_file = "geneFile.txt",
  ref_group_names = c("ref-1", "ref-2"),
  delim = "\t"
)

infercnv_result <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv_epi_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  plot_steps = FALSE,
  write_expr_matrix = TRUE  
)

obs <- read.table("infercnv_epi_output/infercnv.observations.txt", header = TRUE)
ref <- read.table("infercnv_epi_output/infercnv.references.txt", header = TRUE)
expr <- cbind(obs, ref)
expr.scale <- scale(t(expr))

tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min), '-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale, 2, min)
expr_1 <- t(2 * sweep(tmp1, 2, tmp2, "/") - 1)

cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score) <- "cnv_score"
cnv_score <- rownames_to_column(cnv_score, var = "cell")
cnv_score$cell <- gsub("\\.", "-", cnv_score$cell)

meta <- epi_sampled@meta.data |>
  rownames_to_column("cell") |>
  inner_join(cnv_score, by = "cell") |>
  column_to_rownames("cell")

p <- ggboxplot(meta, x = "seurat_clusters", y = "cnv_score", fill = "seurat_clusters") +
  scale_y_continuous(limits = c(0, max(meta$cnv_score) * 1.1)) +
  xlab("Cluster") + ylab("CNV Score") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none")
p


fibro_cells <- colnames(fibroMat)
plasma_cells <- colnames(plasmaMat)
ref_score <- cnv_score %>%
  filter(cell %in% c(fibro_cells, plasma_cells)) %>%
  mutate(
    seurat_clusters = ifelse(cell %in% fibro_cells, "Ref-Fibro", "Ref-Plasma")
  )
meta_combined <- rbind(
  meta[, c("seurat_clusters", "cnv_score")] %>% 
    mutate(seurat_clusters = paste0("Epi-", seurat_clusters)),  # 可选：标记上皮聚类
  ref_score[, c("seurat_clusters", "cnv_score")]
)

p <- ggboxplot(
  meta_combined, 
  x = "seurat_clusters", 
  y = "cnv_score", 
  fill = "seurat_clusters",  
  palette = "npg",           
  outlier.shape = NA         
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # 优化Y轴范围
  labs(x = "Cell Group", y = "CNV Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # 倾斜X轴标签
    legend.position = "none"
  )
print(p)

cluster_levels <- unique(meta_combined$seurat_clusters)
n_colors <- length(cluster_levels)
base_colors <- scales::hue_pal()(n_colors)
names(base_colors) <- cluster_levels

stat.df <- meta_combined %>%
  filter(seurat_clusters != "Ref-Plasma") %>%
  group_by(seurat_clusters) %>%
  summarise(y_pos = quantile(cnv_score, 0.95), .groups = "drop") %>%
  mutate(
    p_value = map_dbl(
      seurat_clusters,
      ~ if (grepl("^Epi", .x)) {
        data_sub <- filter(meta_combined, seurat_clusters %in% c(.x, "Ref-Fibro")) %>%
          mutate(group = seurat_clusters == .x)
        wilcox.test(cnv_score ~ group, data = data_sub)$p.value
      } else {
        NA_real_
      }
    ),
    p_adj = p.adjust(p_value, method = "bonferroni"),
    signif = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_label = y_pos + max(meta_combined$cnv_score) * 0.1
  )

p_labeled <- ggplot(meta_combined, 
                    aes(x = seurat_clusters, y = cnv_score, fill = seurat_clusters)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.85) +
  geom_text(
    data = stat.df,
    aes(x = seurat_clusters, y = y_label, label = signif),
    size = 5,
    vjust = 0,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = base_colors) +  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)),
    limits = c(0, max(stat.df$y_label, na.rm = TRUE) * 1.05)
  ) +
  labs(
    x = "Cell Group", 
    y = "CNV Score",
    title = "CNV Score"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", size = 1)
  )

print(p_labeled)

malignant_clusters <- c("10", "12", "18", "20", "22")
malignant_cells_all <- meta %>%
  rownames_to_column("cell") %>%
  filter(seurat_clusters %in% malignant_clusters) %>%
  pull(cell)

sub.cells$cnv_status <- "Non-malignant epithelial"
sub.cells$cnv_status[colnames(sub.cells) %in% malignant_cells_all] <- "Tumor cell"
sub.cells$cnv_status <- factor(
  sub.cells$cnv_status,
  levels = c("Non-malignant epithelial", "Tumor cell")
)
p_sub_cnv <- DimPlot(
  sub.cells,
  reduction = "umap",
  group.by = "cnv_status",
  cols = c("#1EBDCF", "#E64B35"),
  pt.size = 0.5,
  raster = FALSE
) +
  ggtitle("CNV analysis of epithelial") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p_sub_cnv

con$cnv_label <- "Normal cell"
con$cnv_label[colnames(con) %in% malignant_cells_all] <- "Tumor cell"
con$celltype_cnv <- as.character(con$celltype)
con$celltype_cnv[colnames(con) %in% malignant_cells_all] <- "Tumor cell"
epi_cells_all <- colnames(con)[con$celltype == "Epithelial"]
non_malignant_epi_cells <- setdiff(epi_cells_all, malignant_cells_all)
con$celltype_cnv[colnames(con) %in% non_malignant_epi_cells] <- "Non-malignant epithelial"
con$celltype_cnv <- factor(con$celltype_cnv)



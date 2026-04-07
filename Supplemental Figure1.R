# ===============================Supplemental Figure1======================================
# ---- library ----
library(Seurat)
library(ggplot2)
library(patchwork)

# ---- load data ----
load("Pancancer.Rda")
load("PancancerEC.Rda")

# ---- Dotplot ----
features <- c(
  "EPCAM", "KRT17",                    # Epithelial
  "COL1A1", "COL1A2", "TAGLN",         # Fibroblast
  "PECAM1", "CLDN5",                   # Endothelial
  "CD3D", "CD3G", "NKG7", "GNLY",      # T/NK
  "MS4A1", "CD79A", "JCHAIN", "IGKC",  # B/Plasma
  "LYZ", "CD14", "C1QA" ,              # Myeloid
  "TPSAB1", "CPA3",                    # Mast cell
  "MAP2", "REFOX3"                     # Neurons
)

Idents(con) <- con$celltype
p <- DotPlot(con, features = features)
data <- p$data
data$id <- factor(data$id, levels = celltype_order)

my_gradient_colors <- c("#2c7ac1", "#529fd1", "#83bfe1", "#b7ddf2","#f0faff"
                        ,"#fef9ef","#f9dbb4", "#f1b782", "#e67c54", "#d13f2d","#a2292d"
)

p_main <- ggplot(data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colors = my_gradient_colors, limits = c(0, 3), oob = scales::squish) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    legend.position = "right"
  ) 
+coord_flip()
p_main

anno_df <- data.frame(
  cluster = celltype_order,
  group = celltype_order,
  other = "celltype"
)

p_anno <- ggplot(anno_df, aes(x = cluster, y = other, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = celltype_colors) +
  scale_y_discrete(position = "right") +
  theme_void() +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

final_plot <- p_anno / p_main + plot_layout(heights = c(0.05, 1))
final_plot


DotPlot(EC, features = c("LYVE1", "CCL21", "PROX1"),dot.scale = 8) + RotatedAxis()+scale_color_gradientn(colors = c("blue", "white", "red"))


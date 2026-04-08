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

features1 <- c(
  "ACKR1", "SELE", "TNFAIP3",          #SELE-VenEC
  "CLU","NR2F2","IL1R1",               #CLU-VenEC
  "LYVE1", "PROX1", "CCL21",           #LYVE1-LEC
  "SEMA3G", "FBLN5", "GJA5",           #SEMA3G-ArtEC
  "RGCC", "PLPP1","PLVAP",             #PLPP1-CapEC
  "SLC2A3","CD36", "CA4",              #SLC2A3-CapEC
  "ESM1","INSR", "VWA1",               #ESM1-CapEC
  "ANO2", "EMCN", "CALCRL"             #ANO2-CapEC
)

celltype_order <- c("Neurons" , "Mast_cell","Myeloid" ,"Plasma","B_cell","T/NK_cell", "Endothelial" ,"Fibroblast" ,"Epithelial")
celltype_order1 <- c("ANO2-CapEC","ESM1-CapEC","SLC2A3-CapEC","PLPP1-CapEC","SEMA3G-ArtEC","LYVE1-LEC","CLU-VenEC","SELE-VenEC")

Idents(con) <- con$celltype
p <- DotPlot(con, features = features)
data <- p$data
data$id <- factor(data$id, levels = celltype_order)

my_gradient_colors <- c("#2c7ac1", "#529fd1", "#83bfe1", "#b7ddf2","#f0faff"
                        ,"#fef9ef","#f9dbb4", "#f1b782", "#e67c54", "#d13f2d","#a2292d"
)

p_main <- ggplot(data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colors = my_gradient_colors, limits = c(0, 2), oob = scales::squish) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    legend.position = "right"
  ) 
p_main


Idents(EC) <- EC$celltype
p <- DotPlot(EC, features = features1)
data <- p$data
data$id <- factor(data$id, levels = celltype_order1)
p_main <- ggplot(data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colors = my_gradient_colors, limits = c(0, 1.5), oob = scales::squish) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    legend.position = "right"
  ) 
p_main

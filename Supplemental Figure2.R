# ===============================Supplemental Figure2======================================

# ---- library ----
library(Seurat)

# ---- load data ----
load("PancancerLEC.Rda")

# ---- Featureplot ----
FeaturePlot(LEC, features = c("S100A4"),pt.size = 0.5) +
  scale_color_gradientn(colors = c("#e5f5e0", "#a1d99b", "#31a354", "#006d2c"))
FeaturePlot(LEC, features = c("SNAI1"),pt.size = 0.5) +
  scale_color_gradientn(colors = c("#fff7bc", "#fec44f", "#fe9929", "#cc4c02"))
FeaturePlot(LEC, features = c("TAGLN2"),pt.size = 0.5) +
  scale_color_gradientn(colors = c("#EAF6FB", "#BFDDEA", "#73B3D8", "#34708B"))
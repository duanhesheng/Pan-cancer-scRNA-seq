# ===============================Figure3======================================

# ---- library ----
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)
library(scales)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
library(Seurat)
library(dplyr)
library(patchwork)
library(future)

# ---- load data ----
load("PancancerLEC.Rda")
load("Pancancer.Rda")

# ---- Slingshot ----
scale.data <- GetAssayData(LEC, layer = "scale.data", assay = "RNA")
scale.gene <- rownames(scale.data)
counts <- GetAssayData(LEC, layer = "counts", assay = "RNA")
counts <- counts[scale.gene, ]
sce <- SingleCellExperiment(assays = list(counts = counts))
umap <- Embeddings(LEC, reduction = "umap")
colnames(umap) <- c("UMAP-1", "UMAP-2")
reducedDims(sce) <- SimpleList(UMAP = umap)
colData(sce)$sampleId <- LEC$orig.ident
colData(sce)$celltype <- LEC$celltype  
plot(umap, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1)

sce <- slingshot(sce,
                 clusterLabels = 'celltype',
                 reducedDim = 'UMAP',
                 start.clus = "NLEC",
                 end.clus = NULL)
colnames(colData(sce))

plot(reducedDims(sce)$UMAP, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = brewer.pal(9, "Set1"))
legend("right",
       legend = paste0("lineage", 1:length(slingCurves(SlingshotDataSet(sce)))),
       col = brewer.pal(9, "Set1"),
       inset = 0.8,
       pch = 16)

color_vector <- c(
  "intermediate LEC" = "#F39B7F",
  "Tip-like LEC"     = "#3C5488",
  "apLEC"            = "#E64B35",
  "pEndoMT LEC"      = "#00A087",
  "NLEC"             = "#8491B4",
  "Stalk-like LEC"   = "#4DBBD5"
)

colData(sce)$celltype <- factor(colData(sce)$celltype, levels = names(color_vector))
plotcol <- color_vector[as.character(colData(sce)$celltype)]
plotcol[is.na(plotcol)] <- "lightgrey"
n_lineage <- length(slingCurves(SlingshotDataSet(sce)))
lineage_colors <- brewer.pal(n = max(3, n_lineage), name = "Set1")[1:n_lineage]

plot(reducedDims(sce)$UMAP,
     col = plotcol,
     pch = 16,
     asp = 1,
     xlab = "UMAP-1", ylab = "UMAP-2",
     main = "Slingshot Trajectory by Celltype")
for (i in seq_len(n_lineage)) {
  curve <- slingCurves(SlingshotDataSet(sce))[[i]]
  lines(curve$s[curve$ord, ], lwd = 2, col = lineage_colors[i])
  end_idx <- tail(curve$ord, 2)
  arrows(x0 = curve$s[end_idx[1], 1], y0 = curve$s[end_idx[1], 2],
         x1 = curve$s[end_idx[2], 1], y1 = curve$s[end_idx[2], 2],
         length = 0.1, col = lineage_colors[i], lwd = 2)
}
legend("bottomleft",
       legend = names(color_vector),
       col = color_vector,
       pch = 16,
       bty = "n",
       title = "Celltype")


for (i in seq_len(n_lineage)) {
  curve <- slingCurves(SlingshotDataSet(sce))[[i]]
  ord <- curve$ord
  s <- curve$s[ord, ]
  pseudotime <- curve$pseudotime[ord]
  line_colors <- colorRampPalette(c("#e0f3db", "#a8ddb5", "#43a2ca", "#0868ac"))(nrow(s) - 1)
  for (j in 1:(nrow(s) - 1)) {
    lines(s[j:(j + 1), ], col = line_colors[j], lwd = 5, lend = "round")
  }
  arrow_base <- s[nrow(s) - 1, ]
  arrow_tip <- s[nrow(s), ]
  
  arrows(x0 = arrow_base[1], y0 = arrow_base[2],
         x1 = arrow_tip[1], y1 = arrow_tip[2],
         length = 0.2, col = line_colors[nrow(s) - 1], lwd = 5)
}
alpha_val <- 0.9
celltype_vec <- as.character(colData(sce)$celltype)
plotcol <- adjustcolor(color_vector[celltype_vec], alpha.f = alpha_val)
plotcol[is.na(plotcol)] <- adjustcolor("lightgrey", alpha.f = alpha_val)
plot(reducedDims(sce)$UMAP,
     col = plotcol,
     pch = 16,
     asp = 1,
     xlab = "UMAP-1", ylab = "UMAP-2",
     main = "Slingshot Trajectory by Celltype")
for (i in seq_len(n_lineage)) {
  curve <- slingCurves(SlingshotDataSet(sce))[[i]]
  lines(curve$s[curve$ord, ], lwd = 2, col = lineage_colors[i])
  
  end_idx <- tail(curve$ord, 2)
  arrows(x0 = curve$s[end_idx[1], 1], y0 = curve$s[end_idx[1], 2],
         x1 = curve$s[end_idx[2], 1], y1 = curve$s[end_idx[2], 2],
         length = 0.1, col = lineage_colors[i], lwd = 2)
}
legend("bottomright",
       legend = names(color_vector),
       col = color_vector,  
       pch = 16,
       bty = "n",
       title = "Celltype")
lineage_colors <- list(
  c("#9F6693", "#83bfe1"),   
  c("#9F6693", "#f1b782"),   
  c("#9F6693", "#a1d99b")    
)
lwd_line <- 3
arrow_length <- 0.25  
arrow_width <- 0.08   
for (i in seq_len(n_lineage)) {
  curve <- slingCurves(SlingshotDataSet(sce))[[i]]
  ord <- curve$ord
  s <- curve$s[ord, ]
  pseudotime <- curve$pseudotime[ord]
  
  color_pair <- lineage_colors[[i]]
  n_segments <- max(nrow(s) - 1, 2)
  line_colors <- colorRampPalette(color_pair)(n_segments)
  for (j in 1:(nrow(s) - 1)) {
    lines(s[j:(j + 1), ], col = line_colors[j], lwd = lwd_line, lend = "round")
  }
  arrow_base <- s[nrow(s) - 1, ]
  arrow_tip  <- s[nrow(s), ]
  dx <- arrow_tip[1] - arrow_base[1]
  dy <- arrow_tip[2] - arrow_base[2]
  norm <- sqrt(dx^2 + dy^2)
  if (norm > 1e-3) {
    ux <- dx / norm
    uy <- dy / norm
    perp_x <- -uy
    perp_y <- ux
    tip_x <- arrow_tip[1]
    tip_y <- arrow_tip[2]
    base_x <- tip_x - arrow_length * ux
    base_y <- tip_y - arrow_length * uy
    left_x  <- base_x + arrow_width * perp_x
    left_y  <- base_y + arrow_width * perp_y
    right_x <- base_x - arrow_width * perp_x
    right_y <- base_y - arrow_width * perp_y
    polygon(x = c(tip_x, left_x, right_x),
            y = c(tip_y, left_y, right_y),
            col = line_colors[n_segments],
            border = NA)
  }
}


# ---- Cellchat ----
options(stringsAsFactors = FALSE)  
color.use <- c( 
  "T/NK_cell"                = "#F061DD",
  "B_cell"                   = "#c0a7e4",
  "Plasma"                   = "#A55194",
  "Tumor cell"               = "#a2292d",
  "Non-malignant epithelial" = "#F9DBB4",
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
  "Tumor cell"               = "#a2292d",
  "Non-malignant epithelial" = "#F9DBB4",
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


con.withln <- subset(con, subset = group_annotation == "tumorwithln")
con.withoutln <- subset(con, subset = group_annotation == "tumorwithoutln")

#tumorwithln
data.input <- GetAssayData(con.withln, layer = "data")
meta <- data.frame(
  labels = Idents(con.withln),
  samples = factor(con.withln$orig.ident),
  row.names = Cells(con.withln)
)

cellchat.withln <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "RNA")
cellchat.withln@DB <- CellChatDB.human
cellchat.withln <- subsetData(cellchat.withln)
cellchat.withln <- identifyOverExpressedGenes(cellchat.withln, thresh.pc = 0.1, thresh.fc = 0, thresh.p = 0.05)
cellchat.withln <- identifyOverExpressedInteractions(cellchat.withln)
cellchat.withln <- computeCommunProb(cellchat.withln, type = "truncatedMean", trim = 0.1, nboot = 100)
cellchat.withln <- filterCommunication(cellchat.withln, min.cells = 20)
cellchat.withln <- computeCommunProbPathway(cellchat.withln)
cellchat.withln <- aggregateNet(cellchat.withln)
save(cellchat.withln, file = "cellchat_tumorwithln_obj.Rda")

#tumorwithoutln
data.input <- GetAssayData(con.withoutln, layer = "data")
meta <- data.frame(
  labels = Idents(con.withoutln),
  samples = factor(con.withoutln$orig.ident),
  row.names = Cells(con.withoutln)
)

cellchat.withoutln <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "RNA")
cellchat.withoutln@DB <- CellChatDB.human
cellchat.withoutln <- subsetData(cellchat.withoutln)
cellchat.withoutln <- identifyOverExpressedGenes(cellchat.withoutln, thresh.pc = 0.1, thresh.fc = 0, thresh.p = 0.05)
cellchat.withoutln <- identifyOverExpressedInteractions(cellchat.withoutln)
cellchat.withoutln <- computeCommunProb(cellchat.withoutln, type = "truncatedMean", trim = 0.1, nboot = 100)
cellchat.withoutln <- filterCommunication(cellchat.withoutln, min.cells = 20)
cellchat.withoutln <- computeCommunProbPathway(cellchat.withoutln)
cellchat.withoutln <- aggregateNet(cellchat.withoutln)
save(cellchat.withoutln, file = "cellchat_tumorwithoutln_obj.Rda")

#netVisual_heatmap
netVisual_heatmap(
  merge_cellchat,
  comparison = c(1, 2),
  measure = c("count", "weight"),
  signaling = NULL,
  slot.name = c("netP", "net"),
  color.use = NULL,
  color.heatmap = c("white", "#d13f2d"),
  title.name = NULL,
  width = NULL,
  height = NULL,
  font.size = 8,
  font.size.title = 10,
  cluster.rows = FALSE,
  cluster.cols = FALSE,
  sources.use = NULL,
  targets.use = NULL,
  remove.isolate = FALSE,
  row.show = NULL,
  col.show = NULL
)

p1 <- netVisual_heatmap(cellchat.withln, measure = "count", color.heatmap = "Blues", 
                        title.name = "Tumor with LN - Number of interactions")

p2 <- netVisual_heatmap(cellchat.withoutln, measure = "count", color.heatmap = "Blues", 
                        title.name = "Tumor without LN - Number of interactions")

p1 + p2


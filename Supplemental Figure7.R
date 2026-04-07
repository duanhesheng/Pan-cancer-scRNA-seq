# ===============================Supplemental Figure7======================================
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

# ---- pySCENIC preparation ----
counts_mat <- GetAssayData(LEC, assay = "RNA", layer = "counts")
write.csv(t(as.matrix(counts_mat)), file = "for.scenic.data.csv")

# ---- pySCENIC ----
regulons_incidMat<-get_regulons(loom,column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
regulons<-regulonsToGeneLists(regulons_incidMat)
regulonAUC<-get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds<-get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings<-get_embeddings(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)
DimPlot(LEC,reduction="umap",label=T)

sub_regulonAUC<-regulonAUC[,match(colnames(LEC),colnames(regulonAUC))]
dim(sub_regulonAUC)

identical(colnames(sub_regulonAUC),colnames(LEC))

cellClusters<-data.frame(row.names=colnames(LEC),
                         seurat_clusters=as.character(LEC$celltype))
cellTypes<-data.frame(row.names=colnames(LEC),
                      celltype=LEC$celltype)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4]
save(sub_regulonAUC,cellTypes,cellClusters,LEC,
     file='for_rss_and_visual.Rdata')

selectedResolution <- "celltype" 
cellsPerGroup <- split(rownames(cellTypes),
                       cellTypes[,selectedResolution])

sub_regulonAUC<-sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled<-t(scale(t(regulonActivity_byGroup),
                                        center=T,scale=T))
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

#TF activity
Heatmap(
  regulonActivity_byGroup_Scaled,
  name ="z-score",
  col=colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11,"Spectral"))),
  show_row_names=TRUE,
  show_column_names=TRUE,
  row_names_gp=gpar(fontsize=6),
  clustering_method_rows="ward.D2",
  clustering_method_columns="ward.D2",
  row_title_rot=0,
  cluster_rows=TRUE,
  cluster_row_slices=FALSE,
  cluster_columns=FALSE)

#TF expression
# ---- TF expression ----
DefaultAssay(LEC) <- "RNA"

tf_names <- rownames(regulonActivity_byGroup_Scaled)
tf_names <- gsub("\\s*\\(.*\\)$", "", tf_names)   # 去掉如 " (123g)"
tf_names <- gsub("_extended$", "", tf_names)      # 去掉 "_extended"
tf_names <- unique(tf_names)
expr_mat <- tryCatch(
  GetAssayData(LEC, assay = "RNA", layer = "data"),
  error = function(e) NULL
)
if (is.null(expr_mat)) {
  LEC <- NormalizeData(LEC, assay = "RNA", verbose = FALSE)
  expr_mat <- GetAssayData(LEC, assay = "RNA", layer = "data")
}
tf_names <- tf_names[tf_names %in% rownames(expr_mat)]
tf_expr_byGroup <- sapply(cellsPerGroup, function(cells) {
  rowMeans(as.matrix(expr_mat[tf_names, cells, drop = FALSE]))
})

tf_expr_byGroup <- as.matrix(tf_expr_byGroup)
tf_expr_byGroup_Scaled <- t(scale(t(tf_expr_byGroup), center = TRUE, scale = TRUE))
tf_expr_byGroup_Scaled <- na.omit(tf_expr_byGroup_Scaled)

Heatmap(
  tf_expr_byGroup_Scaled,
  name ="z-score",
  col=colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11,"Spectral"))),
  show_row_names=TRUE,
  show_column_names=TRUE,
  row_names_gp=gpar(fontsize=6),
  clustering_method_rows="ward.D2",
  clustering_method_columns="ward.D2",
  row_title_rot=0,
  cluster_rows=TRUE,
  cluster_row_slices=FALSE,
  cluster_columns=FALSE
)

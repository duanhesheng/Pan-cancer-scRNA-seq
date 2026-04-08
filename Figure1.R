# ===============================Figure1======================================

# ---- library ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)
library(scales)
library(tibble)
library(ggrepel)
library(grid)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(EnhancedVolcano)

# ---- load data ----
load("Pancancer.Rda")
load("PancancerEC.Rda")
load("PancancerLEC.Rda") 

# ---- color preparation ----
color_pan <- c(
  "T/NK_cell"   = "#F061DD",
  "B_cell"      = "#c0a7e4",
  "Plasma"      = "#A55194",
  "Epithelial"  = "#a2292d",
  "Myeloid"     = "#17BECF",
  "Mast_cell"   = "#e67c54",
  "Fibroblast"  = "#98DF8A",
  "Endothelial" = "#FF9896",
  "Neurons"     = "#74adc0"
)

color_ec <- c(
 "PLPP1-CapEC"  = "#4DBBD5",
 "ESM1-CapEC"   = "#E377C2",
 "SELE-VenEC"   = "#00A087",
 "ANO2-CapEC"   = "#3C5488",
 "CLU-VenEC"    = "#F39B7F",
 "LYVE1-LEC"    = "#E64B35",
 "SEMA3G-ArtEC" = "#91D1C2",
 "SLC2A3-CapEC" = "#9467BD"
)

color_lec <- c(
  "intermediate LEC" = "#F39B7F",
  "Tip-like LEC"     = "#3C5488",
  "apLEC"            = "#E64B35",
  "pEndoMT LEC"      = "#00A087",
  "NLEC"             = "#8491B4",
  "Stalk-like LEC"   = "#4DBBD5"
)

# ---- UMAP ----
plot_umap <- function(obj, color_map, pt_size = 0.05) {
  umap_df <- as.data.frame(obj@reductions$umap@cell.embeddings)
  umap_df$celltype <- obj$celltype
  
  umap_df <- umap_df %>% 
    filter(celltype %in% names(color_map))
  
  ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = pt_size, alpha = 0.6, show.legend = TRUE) +
    scale_color_manual(values = color_map) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      text = element_text(family = "serif", size = 12),
      aspect.ratio = 1
    )
}

# All cells
p1 <- plot_umap(con, color_pan, pt_size = 0.05)
p1
# EC
p2 <- plot_umap(EC, color_ec, pt_size = 0.5)
p2
# LEC
p3 <- plot_umap(LEC, color_lec, pt_size = 2)
p3

# ---- Dotplot ----
allmarkers <- FindAllMarkers(LEC, only.pos = TRUE, min.pct = 0.6, logfc.threshold = 0.25)

features1 <- c(
  "SIM2", "DKK3", "PLAC9", "CD9", "STAB2", "PRKG1", "CD36","TFPI", "SERPINE1", "HSPA6", "CXCL8","CXCL1", "SLC27A3", "LENG8", "SPARC", "SPRY1", "COL1A1", "COL1A2",    "HLA-DPA1", "S100A4",   "VIM", "SEPN1", "FAM63B" )

features2 <- c("PECAM1", "CCL21", "LYVE1", "PROX1",  # Endothelial
              "SNAI1", "S100A4",   "VIM", "TAGLN2"  # Mesenchymal
              )
              
celltype_order <- c("NLEC" , "intermediate LEC" , "Tip-like LEC" ,"Stalk-like LEC","apLEC","pEndoMT LEC")
Idents(LEC) <- LEC$celltype
#features1
p <- DotPlot(LEC, features = features1)
#features2
#p <- DotPlot(LEC, features = features2)

data <- p$data
data$id <- factor(data$id, levels = celltype_order)
my_gradient_colors <- c("#2c7ac1", "#529fd1", "#83bfe1", "#b7ddf2","#f0faff","#fef9ef","#f9dbb4", "#f1b782", "#e67c54", "#d13f2d")
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
+coord_flip()
p_main


# ---- Volcano Plot ----
deg = FindMarkers(LEC,ident.1 =  "MT LEC",
                  ident.2 = c("apLEC", "Stalk-like LEC", "Tip-like LEC", "intermediate LEC", "NLEC"), logfc.threshold = 0, min.pct = 0.1)
head(deg[order(deg$p_val_adj),])

EnhancedVolcano(deg,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val_adj')
colnames(deg)[colnames(deg) == "avg_log2FC"] <- "logFC"
colnames(deg)[colnames(deg) == "p_val_adj"] <- "P.Adj"
deg$symbol <- rownames(deg)
nrow(deg)
logFC_t = 1
p_t = 0.05
k1 = (deg$P.Adj < p_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Adj < p_t)&(deg$logFC > logFC_t)
deg = mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)
ex = deg$P.Adj
deg$log10P.Adj = -log10(ex)
Up <- filter(deg,change=="up")
up_genes <- Up$symbol
Down<- filter(deg,change=="down")
down_genes <- Down$symbol
paste0("The number of up gene is ",length(up_genes))
paste0("The number of down gene is ",length(down_genes))
top10sig <- filter(deg,change!="stable") %>% distinct(symbol,.keep_all = T) %>% top_n(10,abs(logFC))
top10sig
up <- filter(top10sig,change=="up")
up
down <- filter(top10sig,change=="down")
down
deg$size <- case_when(!(deg$symbol %in% top10sig$symbol)~ 1,
                      deg$symbol %in% top10sig$symbol ~ 2)
head(deg)
deg <- filter(deg,size==1)
deg$change <- factor(deg$change,
                     levels = c("up","down","stable"),
                     ordered = T)

p0 <-ggplot(data=deg,aes(logFC,log10P.Adj,color=change))
p1 <- p0+geom_point(size=2)
p1
mycolor <- c("#E64B35FF","#4DBBD5FF","gray80")
p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.9))
p2
p3 <- p2+geom_point(data=up,aes(logFC,log10P.Adj),
                    color="#E64B35FF",size=2,alpha=0.9)+
  geom_point(data=down,aes(logFC,log10P.Adj),
             color="#4DBBD5FF",size=2,alpha=0.9)
p3
p4 <- p3+geom_hline(yintercept = c(-log10(0.05)),
                    size = 0.5,
                    color = "black",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 0.5,
             color = "black",
             lty = "dashed")
p4
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 10),
        axis.line = element_line(size = 0.6,colour = "gray30"),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
p4+mytheme

# ---- KEGG ----
s2e = bitr(deg$symbol, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
nrow(deg) 
deg = inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
nrow(deg)
gene_diff = deg$ENTREZID[deg$change != "stable"] 
ekk <- enrichKEGG(gene = gene_diff,organism = 'hsa')
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
ego <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,
                ont = "ALL",readable = TRUE)
class(ekk)
barplot(ego, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
dotplot(ekk,showCategory = 40)
ekk_filtered <- subset(ekk, !grepl("^hsa05", ekk@result$ID))
ekk_new <- new("enrichResult",
               result = ekk_filtered,        
               pvalueCutoff = ekk@pvalueCutoff, 
               pAdjustMethod = ekk@pAdjustMethod, 
               gene = ekk@gene,               
               universe = ekk@universe,        
               organism = ekk@organism,        
               ontology = "KEGG",           
               keytype = "kegg",               
               gene2Symbol = ekk@gene2Symbol  
)
res = ekk_new
plot_enrich_lollipop_neglogP <- function(res, top_n = 12, p_cutoff = 0.05,
                                         alpha_bar = 0.45, box_lwd = 0.8){
  df <- if (inherits(res, "enrichResult")) res@result else res
  stopifnot(all(c("Description","p.adjust","Count") %in% colnames(df)))
  
  df2 <- df %>%
    filter(!is.na(p.adjust)) %>%
    mutate(
      neglogP  = -log10(p.adjust),
      size_bin = cut(Count, breaks = c(-Inf, 7.5, 12.5, 17.5, Inf),
                     labels = c("5","10","15","20"), right = TRUE)
    ) %>%
    filter(p.adjust <= p_cutoff) %>%
    arrange(desc(neglogP), desc(Count)) %>%
    slice_head(n = top_n) %>%
    mutate(term_factor = factor(Description, levels = rev(Description)))
  
  xmax   <- max(df2$neglogP) * 1.05
  xbreak <- pretty(c(0, xmax))
  

  p_left <- ggplot(df2, aes(y = term_factor, x = neglogP)) +
    # 横线颜色映射到 neglogP（与圆点一致）
    geom_segment(aes(x = 0, xend = neglogP, y = term_factor, yend = term_factor,
                     colour = neglogP),
                 linewidth = 0.8) +
    geom_point(aes(size = size_bin, fill = neglogP), shape = 21, stroke = 0) +
    scale_size_manual(values = c(3.0, 4.5, 6.0, 7.5),
                      breaks = c("5","10","15","20"), name = "Count") +
    scale_fill_gradient(low = "#FFE3B6", high = "#E31A1C",
                        name = expression(-log[10](adj.P))) +
    scale_colour_gradient(low = "#FFE3B6", high = "#E31A1C", guide = "none") +
    scale_x_continuous(limits = c(0, xmax), breaks = xbreak,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = expression(-log[10](adjusted~P)), y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.y  = element_blank(),
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      # 仅保留 x 轴小竖杆；去掉 y 轴小横杆
      axis.ticks.x = element_line(colour = "grey40", linewidth = 0.4),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(3, "pt"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = box_lwd),
      plot.margin  = margin(t = 4, r = 6, b = 4, l = 6)
    )
  
  p_right <- ggplot(df2, aes(y = term_factor, x = neglogP, fill = neglogP)) +
    geom_col(width = 0.82, alpha = alpha_bar) +
    geom_text(aes(label = term_factor), x = 0.15, hjust = 0, size = 4.2) +
    scale_fill_gradient(low = "#FFE3B6", high = "#E31A1C",
                        name = expression(-log[10](adj.P))) +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.y  = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = box_lwd),
      plot.margin  = margin(t = 4, r = 6, b = 4, l = 0)
    )
  
  (p_left + p_right) + plot_layout(widths = c(1, 1.05), guides = "collect") &
    theme(legend.position = "right")
}

p <- plot_enrich_lollipop_neglogP(ekk_new, top_n = 12, p_cutoff = 0.05)
p

# ---- Bar Chart ----
color <- c(
  "intermediate LEC" = "#F39B7F", 
  "Tip-like LEC"     = "#3C5488",  
  "apLEC"            = "#E64B35", 
  "pEndoMT LEC"      = "#00A087", 
  "NLEC"             = "#8491B4",  
  "Stalk-like LEC"   = "#4DBBD5"   
)

LEC@meta.data %>%
  dplyr::filter(!is.na(group_annotation)) %>%
  ggplot(aes(x = group_annotation, fill = celltype)) +
  geom_bar(position = "fill", alpha = 0.6) +  
  scale_fill_manual(values = color) +
  RotatedAxis() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.text = element_text(color = 'black', size = 12)
  )

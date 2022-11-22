#pheatmap
library(pheatmap)
load("D1Result_of_4methods.RData")
source("metabolite_stat.R")
mbcluster <- read.xlsx("data/IBD/41564_2018_306_MOESM3_ESM.xlsx", sheetIndex = 1, startRow = 2)
rownames(mbcluster) <- mbcluster$Metabolomic.Feature
m <- read.table("data/IBD/IBD_meta.txt", header = T, row.names = 1, sep = "\t", comment.char = "")
mbcluster$compoundid <- gsub("-", ".", mbcluster$Metabolomic.Feature)
rownames(mbcluster) <- mbcluster$Metabolomic.Feature
d <- data.frame(group = m[rownames(mbmscale),"Diagnosis"], mbmscale[, intersect(colnames(mbmscale), colnames(mbp))])
mbmscale_dsumm <- metabolite_stat(d)
mbmscale_dsumm <- merge(mbmscale_dsumm, mbcluster, by.x = "row.names", by.y = "compoundid", all.x = T)
rownames(mbmscale_dsumm) <- mbmscale_dsumm$Row.names
mbmscale_dsumm$Row.names <- NULL

interest <- mbmscale_dsumm[mbmscale_dsumm$Putative.Chemical.Class %in% 
                                  c("Indoles and derivatives", "Bile acids, alcohols and derivatives", "Sphingolipids"), "Metabolomic.Feature"] #Sphingolipids, Indoles and derivatives, Bile acids, alcohols and derivatives
sph_mbmscale <- mbmscale[, interest]
sph_mbp <- mbp[, interest]
colnames(sph_mbmscale) <- paste0(colnames(sph_mbmscale), "_measured")
colnames(sph_mbp) <- paste0(colnames(sph_mbp), "_predicted")
identical(rownames(sph_mbmscale), rownames(sph_mbp))
sph <- cbind(sph_mbmscale, sph_mbp)
sph_wellpred <- intersect(interest, predres$wellPredicted)
annotation_col = data.frame(
  predicted = rep(c(rep("Wellpredicted", length(sph_wellpred)),
                    rep("Poorlypredicted", length(interest[! interest %in% sph_wellpred]))),
                  each = 2),
  class = rep(c(mbcluster[sph_wellpred, "Putative.Chemical.Class"], 
                mbcluster[interest[! interest %in% sph_wellpred], "Putative.Chemical.Class"]),
              each = 2)
)
rownames(annotation_col) = paste(rep(c(sph_wellpred, interest[! interest %in% sph_wellpred]),each = 2),
                                 rep(c("measured", "predicted")), sep = "_")
annotation_col$class <- factor(annotation_col$class, 
                               levels = c("Indoles and derivatives", "Bile acids, alcohols and derivatives", "Sphingolipids"))
annotation_col <- annotation_col[order(annotation_col$class), ]
annotation_row = data.frame(
  group = m[rownames(m) %in% rownames(sph), "Diagnosis"]
)
rownames(annotation_row) = rownames(m)[rownames(m) %in% rownames(sph)]
ann_colors = list(
  group = c("Control" = "#6cba56", 
            "UC" = "#ed80e6", "CD" = "#ad10ad"),
  predicted = c('Wellpredicted' = "#56c197",
                'Poorlypredicted' = "#9deadb"),
  class = c("Indoles and derivatives" = "#546e7a",
            "Bile acids, alcohols and derivatives" = "#ff8a80", 
            "Sphingolipids" = "#f8bbd0")
)
sph <- sph[rownames(annotation_row), rownames(annotation_col)]

ms_hc <- hclust(dist(sph[rownames(annotation_row), rownames(annotation_col)[grep("measured", rownames(annotation_col))]]))
sampleorder <- data.frame(sid = ms_hc$labels[ms_hc$order], 
                          order = 1:length(ms_hc$labels))
rownames(sampleorder) <- as.character(sampleorder$sid)
sampleorder$group <- m[rownames(sampleorder), "Diagnosis"]#mbmscale_dsumm[rownames(sampleorder), "Putative.Chemical.Class"]
sampleorder$group <- factor(sampleorder$group, levels = c("Control", "UC", "CD"))
sampleorder <- sampleorder[order(sampleorder$group, sampleorder$order), ]
mm_hc <- hclust(dist(t(sph[rownames(annotation_row), rownames(annotation_col)[grep("measured", rownames(annotation_col))]])))
mborder <- data.frame(compoundid = mm_hc$labels[mm_hc$order], 
                          order = 1:length(mm_hc$labels))
rownames(mborder) <- as.character(mborder$compoundid) %>% gsub("_measured", "", .)
mborder$class <- mbcluster[rownames(mborder), "Putative.Chemical.Class"]
mborder$class <- factor(mborder$class, 
                               levels = c("Indoles and derivatives", "Bile acids, alcohols and derivatives", "Sphingolipids"))
mborder$predicted <- ifelse(rownames(mborder) %in% sph_wellpred, "Wellpredicted", "Poorlypredicted")
mborder <- mborder[order(mborder$class, mborder$predicted, mborder$order), ]
sph_m <- sph[rownames(sampleorder), paste0(rownames(mborder), "_measured")] %>% t()
sph_p <- sph[rownames(sampleorder), paste0(rownames(mborder), "_predicted")] %>% t()
p1 <- pheatmap(sph_m, 
               annotation_row = annotation_col, annotation_col = annotation_row,
               cluster_cols = F, cluster_rows = F, show_rownames = F,border_color = NA,
               show_colnames = F, annotation_colors = ann_colors)
p2 <- pheatmap(sph_p, 
               annotation_row = annotation_col, annotation_col = annotation_row,
               cluster_cols = F, cluster_rows = F, show_rownames = F,border_color = NA,
               show_colnames = F, annotation_colors = ann_colors)
library(ggplotify)
library(patchwork)
as.ggplot(p1) + as.ggplot(p2) + 
  plot_annotation(tag_levels = list(c("Measured", "Predicted"))) #+ plot_layout(guides = 'collect') 
ggsave("sup4_pheatmap.pdf", width = 16, height = 5.6)

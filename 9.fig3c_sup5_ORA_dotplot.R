############################################################################
library(clusterProfiler)
library(dplyr)
library(ggplot2)
ibdclass3 <- read.table("metablites2794_summary.txt", header = T, row.names = 1, sep = "\t", fill = T, quote = "", check.names = F)
rownames(ibdclass3) <- ibdclass3$Metabolomic.Feature
# PPM: poorly-predicted metabolites, WPM: well-predicted metabolites
d <- ibdclass3[!is.na(ibdclass3$Putative.Chemical.Class) & (ibdclass3$MMINP %in% c("PPM", "WPM")), ]
tmp <- d[, c("Putative.Chemical.Class", "Metabolomic.Feature")] %>% distinct()
colnames(tmp) <- c("iterm", "sub")
tmp2 <- d[, c("class", "Putative.Chemical.Class")] %>% distinct()
tmp2 <- tmp2[tmp2$class != tmp2$Putative.Chemical.Class, ]
colnames(tmp2) <- c("iterm", "sub")
cluster2class <- rbind(tmp, tmp2)
cluster2class$annotation <- c(rep("cluster", nrow(tmp)), rep("putativeClass", nrow(tmp2)))

#putative class
a3 <- enricher(d$Metabolomic.Feature[d$MMINP_con_cd_fdr < 0.05 & d$MMINP_diff_cd_con > 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
               TERM2NAME = cluster2class[, c("iterm", "iterm")], 
               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a32 <- enricher(d$Metabolomic.Feature[d$MMINP_con_cd_fdr < 0.05 & d$MMINP_diff_cd_con < 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
                TERM2NAME = cluster2class[, c("iterm", "iterm")], 
                pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a3@result$group <- "CD"
a32@result$group <- "Control"
cdora <- rbind(a3@result, a32@result)
cdora$Ratio <- enrichplot:::parse_ratio(cdora$GeneRatio)
cdora$type <- "Predicted"
a4 <- enricher(d$Metabolomic.Feature[d$mbmscale_con_cd_fdr < 0.05 & d$mbmscale_diff_cd_con > 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
               TERM2NAME = cluster2class[, c("iterm", "iterm")], 
               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a42 <- enricher(d$Metabolomic.Feature[d$mbmscale_con_cd_fdr < 0.05 & d$mbmscale_diff_cd_con < 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
                TERM2NAME = cluster2class[, c("iterm", "iterm")], 
                pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a4@result$group <- "CD"
a42@result$group <- "Control"
cdora2 <- rbind(a4@result, a42@result)
cdora2$Ratio <- enrichplot:::parse_ratio(cdora2$GeneRatio)
cdora2$type <- "Measured"
mergedcdora <- rbind(cdora, cdora2)
mergedcdora$ID <- gsub("\"", "", mergedcdora$ID)
mergedcdora$ID <- factor(mergedcdora$ID, levels = mergedcdora[order(mergedcdora$pvalue), "ID"] %>% unique %>% rev)
ggplot(mergedcdora[mergedcdora$pvalue < 0.05, ], aes(group, ID, color = pvalue))+
  geom_point(aes(size = Ratio))+
  theme_bw()+
  labs(x = NULL, y = NULL)+
  theme(text = element_text(color = "black"))+
  facet_wrap(~type)+
  scale_color_continuous(low = "#e64a19", 
                         high = "#00acc1",
                         trans='log10',
                         guide = guide_colorbar(reverse = TRUE)) #mid = "#f9a825", low="#00acc1", high="#e64a19"
ggsave("fig3c_CDvsCon_ORA_dotplot.pdf", width = 6, height = 3.2)

#uc vs control
a3 <- enricher(d$Metabolomic.Feature[d$MMINP_con_uc_fdr < 0.05 & d$MMINP_diff_uc_con > 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
               TERM2NAME = cluster2class[, c("iterm", "iterm")], 
               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a32 <- enricher(d$Metabolomic.Feature[d$MMINP_con_uc_fdr < 0.05 & d$MMINP_diff_uc_con < 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
                TERM2NAME = cluster2class[, c("iterm", "iterm")], 
                pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a3@result$group <- "UC"
a32@result$group <- "Control"
ucora <- rbind(a3@result, a32@result)
ucora$Ratio <- enrichplot:::parse_ratio(ucora$GeneRatio)
ucora$type <- "Predicted"
a4 <- enricher(d$Metabolomic.Feature[d$mbmscale_con_uc_fdr < 0.05 & d$mbmscale_diff_uc_con > 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
               TERM2NAME = cluster2class[, c("iterm", "iterm")], 
               pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a42 <- enricher(d$Metabolomic.Feature[d$mbmscale_con_uc_fdr < 0.05 & d$mbmscale_diff_uc_con < 0], TERM2GENE = cluster2class[, c("iterm", "sub")], 
                TERM2NAME = cluster2class[, c("iterm", "iterm")], 
                pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "none"
)
a4@result$group <- "UC"
a42@result$group <- "Control"
ucora2 <- rbind(a4@result, a42@result)
ucora2$Ratio <- enrichplot:::parse_ratio(ucora2$GeneRatio)
ucora2$type <- "Measured"
mergeducora <- rbind(ucora, ucora2)
mergeducora$ID <- gsub("\"", "", mergeducora$ID)
mergeducora$ID <- factor(mergeducora$ID, levels = mergeducora[order(mergeducora$pvalue), "ID"] %>% unique %>% rev)
mergeducora$group <- factor(mergeducora$group, levels = c("UC", "Control"))
ggplot(mergeducora[mergeducora$pvalue < 0.05, ], aes(group, ID, color = pvalue))+
  geom_point(aes(size = Ratio))+
  theme_bw()+
  labs(x = NULL, y = NULL)+
  theme(text = element_text(color = "black"))+
  facet_wrap(~type)+
  scale_color_continuous(low = "#e64a19", 
                         high = "#00acc1",
                         trans='log10',
                         guide = guide_colorbar(reverse = TRUE)) #mid = "#f9a825", low="#00acc1", high="#e64a19"
ggsave("sup5_UCvsCon_ORA_dotplot.pdf", width = 6, height = 3.2)

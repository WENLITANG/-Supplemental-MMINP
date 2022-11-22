library(ggplot2)
ibdclass3 <- read.table("metablites2794_summary.txt", header = T, row.names = 1, sep = "\t", fill = T, quote = "", check.names = F)
rownames(ibdclass3) <- ibdclass3$Metabolomic.Feature
ibdclass <- ibdclass3[complete.cases(ibdclass3$Putative.Chemical.Class), ]

## diff: CD vs Control
dpred <- data.frame(group = "Predicted", ibdclass[, c("MMINP_con_cd_fdr", "Putative.Chemical.Class", "Metabolomic.Feature","MMINP_diff_cd_con")])
dmeasured <- data.frame(group = "Measured", ibdclass[, c("mbmscale_con_cd_fdr", "Putative.Chemical.Class", "Metabolomic.Feature","mbmscale_diff_cd_con")])
colnames(dpred) <- colnames(dmeasured) <- c("group", "con_cd_fdr", "Putative.Chemical.Class", "Metabolomic.Feature", "diff_cd_con")
a <- intersect(na.omit(dpred[dpred$con_cd_fdr<0.05, "Metabolomic.Feature"]), 
               na.omit(dmeasured[dmeasured$con_cd_fdr<0.05, "Metabolomic.Feature"]))
d <- rbind(dmeasured[a, ], dpred[a, ])
d2 <- d[d$Putative.Chemical.Class %in% names(which(sort(table(ibdclass[a, "Putative.Chemical.Class"]))>5)), ]
d2$Putative.Chemical.Class <- factor(d2$Putative.Chemical.Class, 
                                     levels = names(sort(tapply(d2[d2$group %in% "Predicted","diff_cd_con"], 
                                                                d2[d2$group %in% "Predicted","Putative.Chemical.Class"], median))))
d2$group <- factor(d2$group, levels = c("Predicted", "Measured"))
xaxis <- paste0(levels(d2$Putative.Chemical.Class), "(", table(d2$Putative.Chemical.Class)/2, ")")
ggplot(d2, aes(x = Putative.Chemical.Class, y = diff_cd_con, color = group))+
  geom_boxplot(outlier.color = NA)+
  geom_point(position = position_dodge(0.8), size = 0.1, alpha = 0.8)+
  coord_flip()+geom_hline(yintercept = 0)+theme_bw()+
  scale_x_discrete(labels = xaxis)+
  scale_color_manual(values = c('Predicted' = "#00acc1", #'#8dbdda', 
                                'Measured' = "#9e9d24" #'#83c202'
  ))+
  xlab("Putative Chemical Class\n(number of members)")+
  ylab("Median difference of metabolites\nbetween CD and Control")
ggsave("fig3b_CDvsCon_diffClass_member5.pdf", width = 6, height = 4)

d2 <- d[d$Putative.Chemical.Class %in% names(which((table(ibdclass[a, "Putative.Chemical.Class"])>2) & (table(ibdclass[a, "Putative.Chemical.Class"])<5))), ]
d2$Putative.Chemical.Class <- factor(d2$Putative.Chemical.Class, 
                                     levels = names(sort(tapply(d2[d2$group %in% "Predicted","diff_cd_con"], 
                                                                d2[d2$group %in% "Predicted","Putative.Chemical.Class"], median))))
d2$group <- factor(d2$group, levels = c("Predicted", "Measured"))
xaxis <- paste0(levels(d2$Putative.Chemical.Class), "(", table(d2$Putative.Chemical.Class)/2, ")")
ggplot(d2, aes(x = Putative.Chemical.Class, y = diff_cd_con, color = group))+
  geom_boxplot(outlier.color = NA)+
  geom_point(position = position_dodge(0.8), size = 0.1, alpha = 0.8)+
  coord_flip()+geom_hline(yintercept = 0)+theme_bw()+
  scale_x_discrete(labels = xaxis)+
  scale_color_manual(values = c('Predicted' = "#00acc1", #'#8dbdda', 
                                'Measured' = "#9e9d24" #'#83c202'
  ))+
  xlab("Putative Chemical Class\n(number of members)")+
  ylab("Median difference of metabolites\nbetween CD and Control")
ggsave("sup3d_CDvsCon_diffClass_member3-4.pdf", width = 6, height = 4)


## UC vs Control
dpred <- data.frame(group = "Predicted", ibdclass[, c("MMINP_con_uc_fdr", "Putative.Chemical.Class", "Metabolomic.Feature","MMINP_diff_uc_con")])
dmeasured <- data.frame(group = "Measured", ibdclass[, c("mbmscale_con_uc_fdr", "Putative.Chemical.Class", "Metabolomic.Feature","mbmscale_diff_uc_con")])
colnames(dpred) <- colnames(dmeasured) <- c("group", "con_uc_fdr", "Putative.Chemical.Class", "Metabolomic.Feature", "diff_uc_con")
b <- intersect(na.omit(dpred[dpred$con_uc_fdr<0.05, "Metabolomic.Feature"]), 
               na.omit(dmeasured[dmeasured$con_uc_fdr<0.05, "Metabolomic.Feature"]))
d <- rbind(dmeasured[b, ], dpred[b, ])
d2 <- d[d$Putative.Chemical.Class %in% names(which(sort(table(ibdclass[b, "Putative.Chemical.Class"]))>2)), ]
d2$Putative.Chemical.Class <- factor(d2$Putative.Chemical.Class, 
                                     levels = names(sort(tapply(d2[d2$group %in% "Predicted","diff_uc_con"], 
                                                                d2[d2$group %in% "Predicted","Putative.Chemical.Class"], median))))
d2$group <- factor(d2$group, levels = c("Predicted", "Measured"))
xaxis <- paste0(levels(d2$Putative.Chemical.Class), "(", table(d2$Putative.Chemical.Class)/2, ")")
ggplot(d2, aes(x = Putative.Chemical.Class, y = diff_uc_con, color = group))+
  geom_boxplot(outlier.color = NA)+
  geom_point(position = position_dodge(0.8), size = 0.1, alpha = 0.8)+
  coord_flip()+geom_hline(yintercept = 0)+theme_bw()+
  scale_x_discrete(labels = xaxis)+
  scale_color_manual(values = c('Predicted' = "#00acc1", #'#8dbdda', 
                                'Measured' = "#9e9d24")) +
  xlab("Putative Chemical Class\n(number of members)")+
  ylab("Median difference of metabolites\nbetween UC and Control")
ggsave("sup3e_UCvsCon_diffClass_member3.pdf", width = 6, height = 4)

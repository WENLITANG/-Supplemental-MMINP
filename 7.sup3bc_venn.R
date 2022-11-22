load("D1Result_of_4methods.RData")
source("metabolite_stat.R")
d <- data.frame(group = m[rownames(mbp), "Diagnosis"], mbp)
mbp_dsumm2 <- metabolite_stat(d)
d <- data.frame(group = m[rownames(mbm),"Diagnosis"], mbmscale[, intersect(colnames(mbp), colnames(mbmscale))])
mbmscale_dsumm2 <- metabolite_stat(d)

library(ggvenn)
library(patchwork)
p2 <- ggvenn(data = list(MMINP=rownames(mbp_dsumm2)[mbp_dsumm2$con_cd_fdr < 0.05],
                         Measured=rownames(mbmscale_dsumm2)[mbmscale_dsumm2$con_cd_fdr < 0.05]),
             fill_color = c("#00acc1", "#9e9d24"), #fill_alpha = 0.5, stroke_alpha = 0.5,
             show_percentage = F, text_size = 6) #4*4
p1 <- ggvenn(data = list(MMINP=rownames(mbp_dsumm2)[mbp_dsumm2$con_uc_fdr < 0.05],
                         Measured=rownames(mbmscale_dsumm2)[mbmscale_dsumm2$con_uc_fdr < 0.05]),
             fill_color = c("#00acc1", "#9e9d24"),
             show_percentage = F, text_size = 6) #4*4
p2 / p1 + 
  plot_annotation(tag_levels = list(c("B. CD vs Control", "C. UC vs Control"))) &
  theme(plot.tag.position = c(0,1),
        plot.tag = element_text(size = 20, vjust = 0.5, hjust = -0.1))
ggsave("sup3bc_venn.pdf", height = 6, width = 6)

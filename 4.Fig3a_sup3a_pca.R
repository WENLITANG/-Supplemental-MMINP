###############################################################PCA
library(vegan)
library(ggplot2)
library(patchwork)
load("D1Result_of_4methods.RData")
con <- m[rownames(mbm), ] %>% subset(Diagnosis == "Control") %>% rownames()
cd <- m[rownames(mbm), ] %>% subset(Diagnosis == "CD") %>% rownames()
uc <- m[rownames(mbm), ] %>% subset(Diagnosis == "UC") %>% rownames()
mbpdist <- vegan::vegdist(mbp[c(con, uc, cd), ], method = "euclidean") %>% as.matrix() %>% as.data.frame()
mbmscaledist <- vegan::vegdist(mbmscale[c(con, uc, cd), colnames(mbmscale) %in% colnames(mbp)], method = "euclidean") %>% as.matrix() %>% as.data.frame()

##################################################################### Fig3A CD vs Control
sap <- c(con, cd)
d1 <- mbmscale[sap, colnames(mbmscale) %in% colnames(mbp)]
d1 <- d1[, colSums(d1)>0]
a.pca <- prcomp(d1, center = F,scale = F)
plotda <- as.data.frame(a.pca$x[, 1:5])
plotda$group <- m[sap, "Diagnosis"]
eig <- a.pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:5] %>% round(digits = 2)
adon <- adonis(mbmscaledist[sap, sap]~m[sap, "Diagnosis"], permutations = 10000)$aov.tab %>% as.matrix()
p22 <- ggplot(plotda, aes(PC1, PC2, group = group, fill = group, color = group))+
  geom_point(size = 3, shape = 21, color = "black")+
  #stat_ellipse(type = "t",  level = 0.95, size = 1)+ #linetype = 2,
  #geom_encircle(s_shape = 1, expand = 0)+
  labs(x = paste0("PC1(", eig[1], "%)"),
       y = paste0("PC2(", eig[2], "%)"),
       subtitle = paste0("R2: ", round(adon[1, 5], 2), 
                         "  P: ", ifelse(round(adon[1, 6], 2) == 0, "< 0.001", round(adon[1, 6], 2))),
       title = "Measured(scaled)")+
  theme(text = element_text(size = 12, color = "black"))+
  scale_color_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  scale_fill_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_legend(title = NULL))+
  theme_bw()

a.pca <- prcomp(mbp[sap, colnames(mbp)], center = F,scale = F)
plotda <- as.data.frame(a.pca$x[, 1:5])
plotda$group <- m[sap, "Diagnosis"]
eig <- a.pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:5] %>% round(digits = 2)
adon <- adonis(mbpdist[sap, sap]~m[sap, "Diagnosis"], permutations = 10000)$aov.tab %>% as.matrix()
p32 <- ggplot(plotda, aes(PC1, PC2, group = group, fill = group, color = group))+
  geom_point(size = 3, shape = 21, color = "black")+
  #stat_ellipse(type = "t",  level = 0.95, size = 1)+ #linetype = 2,
  #geom_encircle(s_shape = 1, expand = 0)+
  labs(x = paste0("PC1(", eig[1], "%)"),
       y = paste0("PC2(", eig[2], "%)"),
       subtitle = paste0("R2: ", round(adon[1, 5], 2), 
                         "  P: ", ifelse(round(adon[1, 6], 2) == 0, "< 0.001", round(adon[1, 6], 2))),
       title = "Predicted")+
  theme(text = element_text(size = 12, color = "black"))+
  scale_color_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  scale_fill_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_legend(title = NULL))+
  theme_bw()

p22 / p32 + plot_layout(guides = "collect")
ggsave("fig3a_con_vs_cd_pca.pdf", height = 6.4, width = 4)


##################################################################### Fig3A UC vs Control
sap <- c(con, uc)
d1 <- mbmscale[sap, colnames(mbmscale) %in% colnames(mbp)]
d1 <- d1[, colSums(d1)>0]
a.pca <- prcomp(d1, center = F,scale = F)
plotda <- as.data.frame(a.pca$x[, 1:5])
plotda$group <- m[sap, "Diagnosis"]
eig <- a.pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:5] %>% round(digits = 2)
adon <- adonis(mbmscaledist[sap, sap]~m[sap, "Diagnosis"], permutations = 10000)$aov.tab %>% as.matrix()
p23 <- ggplot(plotda, aes(PC1, PC2, group = group, fill = group, color = group))+
  geom_point(size = 3, shape = 21, color = "black")+
  #stat_ellipse(type = "t",  level = 0.95, size = 1)+ #linetype = 2,
  #geom_encircle(s_shape = 1, expand = 0)+
  labs(x = paste0("PC1(", eig[1], "%)"),
       y = paste0("PC2(", eig[2], "%)"),
       subtitle = paste0("R2: ", round(adon[1, 5], 2), 
                         "  P: ", ifelse(round(adon[1, 6], 2) == 0, "< 0.001", round(adon[1, 6], 2))),
       title = "Measured(scaled)")+
  theme(text = element_text(size = 12, color = "black"))+
  scale_color_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  scale_fill_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_legend(title = NULL))+
  theme_bw()

a.pca <- prcomp(mbp[sap, colnames(mbp)], center = F,scale = F)
plotda <- as.data.frame(a.pca$x[, 1:5])
plotda$group <- m[sap, "Diagnosis"]
eig <- a.pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:5] %>% round(digits = 2)
adon <- adonis(mbpdist[sap, sap]~m[sap, "Diagnosis"], permutations = 10000)$aov.tab %>% as.matrix()
p33 <- ggplot(plotda, aes(PC1, PC2, group = group, fill = group, color = group))+
  geom_point(size = 3, shape = 21, color = "black")+
  #stat_ellipse(type = "t",  level = 0.95, size = 1)+ #linetype = 2,
  #geom_encircle(s_shape = 1, expand = 0)+
  labs(x = paste0("PC1(", eig[1], "%)"),
       y = paste0("PC2(", eig[2], "%)"),
       subtitle = paste0("R2: ", round(adon[1, 5], 2), 
                         "  P: ", ifelse(round(adon[1, 6], 2) == 0, "< 0.001", round(adon[1, 6], 2))),
       title = "Predicted")+
  theme(text = element_text(size = 12, color = "black"))+
  scale_color_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  scale_fill_manual(values = c("UC" = "#FFA500", "CD" = "#A52A2A", "Control" = "#6495ED"))+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         shape = guide_legend(title = NULL))+
  theme_bw()
p23 + p33 + plot_layout(guides = 'collect')
ggsave("sup3a_con_vs_uc_pca.pdf", height = 3, width = 6)

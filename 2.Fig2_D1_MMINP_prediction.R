library(xlsx)
library(magrittr)
load("D1Result_of_4methods.RData")
mbcluster <- read.xlsx("data/IBD/41564_2018_306_MOESM3_ESM.xlsx", sheetIndex = 1, startRow = 2)
rownames(mbcluster) <- mbcluster$Metabolomic.Feature
m <- read.table("data/IBD/IBD_meta.txt", header = T, row.names = 1, sep = "\t", comment.char = "")
mbcluster$compoundid <- gsub("-", ".", mbcluster$Metabolomic.Feature)
rownames(mbcluster) <- mbcluster$Metabolomic.Feature
mbpredict <- mbcluster[colnames(mbtscale), ]
mbpredict$assigned <- ifelse(is.na(mbpredict$Putative.Chemical.Class), "Unassigned", "Assigned")
mbpredict$MMINP <- "NIM" #metabolites removed during the modeling process 
mbpredict[rownames(model1$trainres$res), "MMINP"] <- "PPM" #poorly-predicted metabolites
mbpredict[predres$wellPredicted, "MMINP"] <- "WPM" #well-predicted metabolites
mbpredict$MMINP <- factor(mbpredict$MMINP, levels = c("WPM", "PPM", "NIM"))
mbpredict$putativeclass2 <- ifelse(is.na(mbpredict$Putative.Chemical.Class), "Unassigned",
                                   ifelse(mbpredict$Putative.Chemical.Class %in% names(which(table(mbpredict$Putative.Chemical.Class)>4)),
                                          mbpredict$Putative.Chemical.Class, "Others"))


############################################################################################ Fig2A pie
library(psych)
d <- rbind(table(mbpredict$assigned, mbpredict$MMINP), count = table(mbpredict$MMINP)) %>% 
  t() %>% as.data.frame()
d$all2 <- d$count / sum(d$count) 
d$Assigned2 <- d$Assigned / sum(d$count)
d$Unassigned2 <- d$Unassigned / sum(d$count)
d$MMINP <- rownames(d)
d2 <- reshape2::melt(d, id = c("MMINP", "count", "Assigned", "Unassigned", "all2")) %>% as.data.frame()
d2$variable <- gsub("2", "", d2$variable)
d2$variable2 <- paste(d2$MMINP, d2$variable, sep = "_")
d2$label1 <- paste0(d2$MMINP, "\n(", round(d$all2*100, 2), "%)")
d2$label2 <- paste0(d2$variable, "\n(", round(d2$value*100, 2), "%)")
d2$x <- "a"
d2$y <- "b"
d2$MMINP <- factor(d2$MMINP, levels = rev(c("WPM", "PPM", "NIM"))) #"Wellpredicted", "Predicted", "Unpredictable"
rownames(d2) <- d2$variable2
d2 <- d2[paste(rep(levels(d2$MMINP), each = 2), c("Assigned", "Unassigned"), sep = "_"), ]
d2$variable2 <- factor(d2$variable2, levels = rownames(d2))

# Build the relative position of the tag
for (i in seq(nrow(d2), 1)) {
  if (i == nrow(d2)) {
    d2$assign.y2[i] = d2$value[i] / 2
  }else{
    d2$assign.y2[i] = sum(d2$value[(i + 1):nrow(d2)]) + d2$value[i] / 2
  }
}
d2$all2.y1[5:6] <- d2$all2[5] / 2
d2$all2.y1[3:4] <- d2$all2[5] + d2$all2[3] / 2
d2$all2.y1[1:2] <- d2$all2[5] + d2$all2[3] + d2$all2[1] / 2

library(shadowtext)
library(ggplot2)
ggplot(d2) +
  geom_bar(aes(x, all2/2, fill = MMINP), stat = 'identity', width = 1.2) +
  geom_shadowtext(aes(1, as.numeric(all2.y1), label = label1),
                  size =6, bg.r = 0.1, bg.colour = "white", color = "black") +
  geom_bar(aes(y, value, fill = variable2), 
           stat = 'identity', width = 0.8, color = 'white') +
  geom_shadowtext(aes(2, as.numeric(assign.y2),label = label2),
                  size = 6, color = 'black', bg.r = 0.1, bg.color = "white") +
  scale_y_continuous(labels = scales::percent) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_manual(values = c('WPM' = "#8DBDDA", #"#CAE9F1",
                               'WPM_Assigned' = "#86c697",
                               'WPM_Unassigned' = "#FEDC85",
                               'PPM' = "#CAE9F1", #"#8DBDDA",
                               'PPM_Assigned' = "#86c697",
                               'PPM_Unassigned' = "#FEDC85",
                               'NIM' = "#84A0AB",
                               'NIM_Assigned' = "#86c697",
                               'NIM_Unassigned' = "#FEDC85"))+
  theme(legend.position = 'none') 
ggsave("fig2a_circle.pdf", width = 6, height = 6)

############################################################################################ Fig2B barplot
wpm <- mbpredict[mbpredict$MMINP %in% "Wellpredicted" & (mbpredict$assigned %in% "Assigned"), ]
toshow <- names(which(table(wpm$Putative.Chemical.Class)>4))
wpm$putativeclass2 <- ifelse(wpm$Putative.Chemical.Class %in% toshow,
                             wpm$Putative.Chemical.Class, "Others")
showlevel <- names(sort(table(wpm$putativeclass2))) %>% .[!. %in% "Others"]
assigned_MMINP <- mbpredict[mbpredict$Putative.Chemical.Class %in% toshow, ]
assigned_MMINP$Putative.Chemical.Class <- factor(assigned_MMINP$Putative.Chemical.Class, levels = showlevel)
assigned_MMINP$MMINP <- factor(assigned_MMINP$MMINP, levels = c("NIM", "PPM", "WPM"))
ggplot(assigned_MMINP, aes(Putative.Chemical.Class, fill = MMINP))+
  geom_bar(width = 0.8, color = "black", alpha = 0.8)+ 
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values = c("#84A0AB", "#CAE9F1", "#8DBDDA"))+
  xlab("Putative Chemical Class")+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0,0)) 
ggsave("fig2b_barplot.pdf", width = 10, height = 5)

######################################################################################### Fig2C procrustes
library(vegan)
distmethod <- "euclidean"
mbmdist <- vegdist(mbmscale[, names(which(colSums(mbm[, colnames(mbp)]) > 0))], method = distmethod)
mbpdist <- vegdist(mbp, method = distmethod)
mbm_pca <- cmdscale(mbmdist)
mbp_pca <- cmdscale(mbpdist)
proc <- procrustes(mbm_pca, mbp_pca)
#extract information 
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
colnames(Y) <- gsub("Dim", "PC", colnames(Y))
X <- data.frame(proc$rotation)
set.seed(123)
prot <- protest(X = mbm_pca, Y = mbp_pca, permutations = how(nperm = 999))
pvalue <- ifelse(prot$signif < 0.001, "< 0.001",round(prot$signif, 3)) 
m2 <- round(prot$ss, 3)
#add group
Y$group <- m[rownames(Y), "Diagnosis"]
#plot
ggplot(Y, aes(color = group)) +
  geom_point(aes(X1, X2, color = group), size = 3, shape = 16) +
  geom_point(aes(PC1, PC2, color = group), size = 3, shape = 1) +
  scale_color_manual(values = c("CD" = "#e64a19", #"#E5352B",
                                "Control" =  "#00acc1", #"#0081B4", 
                                "UC" =  "#f9a825"#"#EF9020"
                                ))+
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, 'cm')),
               size = 0.6, show.legend = F) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        title = element_text(size = 16), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), legend.text = element_text(size = 12),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '',
       title = distmethod, 
       subtitle = paste0("M^2: ", m2, "   p: ", pvalue))
ggsave(paste0("fig2c_", distmethod, "_procrustes.pdf"), width = 5.6, height = 5)

######################################################################################### Sup2B: mg VS mb, procrustes
library(vegan)
distmethod <- "euclidean"
#mbmdist <- vegdist(mbmscale[, names(which(colSums(mbm[, colnames(mbp)]) > 0))], method = distmethod)
mgpdist <- vegdist(mgpscale, method = distmethod)
#mbm_pca <- cmdscale(mbmdist)
mgp_pca <- cmdscale(mgpdist)
proc <- procrustes(mbm_pca, mgp_pca)
#extract information 
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
colnames(Y) <- gsub("Dim", "PC", colnames(Y))
X <- data.frame(proc$rotation)
set.seed(123)
prot <- protest(X = mbm_pca, Y = mgp_pca, permutations = how(nperm = 999))
pvalue <- ifelse(prot$signif < 0.001, "< 0.001",round(prot$signif, 3)) 
m2 <- round(prot$ss, 3)
#add group
Y$group <- m[rownames(Y), "Diagnosis"]
#plot
ggplot(Y, aes(color = group)) +
  geom_point(aes(X1, X2, color = group), size = 3, shape = 16) +
  geom_point(aes(PC1, PC2, color = group), size = 3, shape = 1) +
  scale_color_manual(values = c("CD" = "#e64a19", #"#E5352B",
                                "Control" =  "#00acc1", #"#0081B4", 
                                "UC" =  "#f9a825"#"#EF9020"
  ))+
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, 'cm')),
               size = 0.6, show.legend = F) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        title = element_text(size = 16), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), legend.text = element_text(size = 12),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '',
       title = distmethod, 
       subtitle = paste0("M^2: ", m2, "   p: ", pvalue))
ggsave(paste0("Sup2b_mgp_mbm_", distmethod, "_procrustes.pdf"), width = 5.6, height = 5)

library(ggplot2)
library(aplot)
m <- read.table("4methods_compare.txt", header = T, row.names = 1, sep = "\t")
m$index <- rownames(m)

## Number
d <- reshape2::melt(m[1:2, ])
d$variable <- factor(d$variable, levels = c("O2PLS", "MMINP", "MelonnPan", "ENVIM"))
d$index <- factor(d$index, levels = c("predicted", "TP: wellpredicted"))

p1 <- ggplot(d, aes(x = variable, y = value, group = index, label = value, color = index))+
  geom_linerange(aes(ymin = 0, ymax = value, x = variable), 
                 position = position_dodge(width = 0.3),
                 color = "lightgray", size = 2)+
  geom_point(position = position_dodge(0.3), 
    size = 6)+
  geom_text(position = position_dodge(width = 0.3), vjust = -1, color = "black")+
  theme_classic()+
  ylab("Number")+
  xlab(NULL)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14))+
  scale_color_manual(values = c("#00bfc4", "#f8766d"),
                     name = "Number",
                     breaks = c("predicted", "TP: wellpredicted"),
                     labels = c("Predicted", "Wellpredicted"))

## Ratio
dimnames(m)
d_dot <- reshape2::melt(m[10:12, ])
d_dot$variable <- factor(d_dot$variable, levels = c("O2PLS", "MMINP", "MelonnPan", "ENVIM"))
p2 <- ggplot(d_dot, aes(x = variable, y = value, color = index, group = index))+
  geom_point(size = 3, alpha = 0.8)+
  geom_line(linetype = "dashed", size = 0.6)+
  theme_minimal()+
  ylab("Ratio")+
  xlab(NULL)+
  #coord_flip()
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14))+
  scale_color_manual(values = c("#00acc1", "#b79f00", "#e64a19"), 
                     breaks = c("predicted_all_ratio", "wellpredicted_all_ratio", "wellpredicted_predicted_ratio"), 
                     labels = c("Predicted/All", "Wellpredicted/All", "Wellpredicted/Predicted"), 
                     name = "Ratio")

## Number + Ratio
p1 %>% insert_top(p2, height = 0.3)
ggsave("fig4a.pdf", height = 6, width = 7.2)


## Precision, F1score
d_dot2 <- reshape2::melt(m[6:9, ])
d_dot2$variable <- factor(d_dot2$variable, levels = c("O2PLS", "MMINP", "MelonnPan", "ENVIM"))
d_dot2$index <- factor(d_dot2$index, levels = c("precision", "accuracy", "recall", "F1score"), 
                       labels = c("Precision", "Accuracy", "Recall", "F1 score"))
ggplot(d_dot2, aes(x = variable, y = value, fill = variable))+
  geom_bar(position = position_dodge(), stat = "identity", color = "black", 
           width = 0.8, alpha = 0.8, show.legend = F)+
  theme_bw()+
  facet_wrap(~index)+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values = c("#b79f00", "#e64a19", "#00acc1", "#619cff"), name = NULL)+
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 10, color = "black"))
ggsave("fig4b.pdf", height = 6, width = 4)


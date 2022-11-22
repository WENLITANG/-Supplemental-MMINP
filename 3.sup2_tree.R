library(dplyr)
library(ggtree)
library(do)
library(tidytree)
library(ggplot2)
library(ggtreeExtra)
library(ggnewscale)
colors = c("#80B1D3","#B3DE69", #"#FFFFB3",
           "#8DD3C7","#4daf4a",
           "#377eb8","#BEBADA","#FB8072","#FDB462","#FCCDE5",
           "#BC80BD","#CCEBC5","#FFED6F","#CD4F39","#BC41A4",
           "#4F94CD","#E41A1C","#00CD66","#CD3278","#CD8A96",
           "#00C5CD","#CDCD00","#CD85CD","#CD853F","#8B5A2B",
           "#5CACEE","#EE5C42","#00EE76","#EE4A8C","#EED8AE",
           "#00E5EE", "#cccc00", #"#EEEE00",
           "#EED2EE","#EE9A49","#E41A1C",
           "#377EB8","#FF6A6A","#87CEFA","#6E8B3D", "#ffc267",#"#FFEBCD",
           "#B2DFEE")

ibdclass3 <- read.table("metablites2794_summary.txt", header = T, sep = "\t", quote = "")
# PPM: poorly-predicted metabolites, WPM: well-predicted metabolites
d1 <- ibdclass3[ibdclass3$MMINP %in% c("PPM", "WPM"), c("kingdom", "super_class", "class", "Putative.Chemical.Class")]
d1 <- d1[complete.cases(d1$Putative.Chemical.Class), ]
d1$super_class <- paste0("SC_", d1$super_class)
d1$class <- paste0("C_", d1$class)
d1$Putative.Chemical.Class <- paste0("PC_", d1$Putative.Chemical.Class)
tmp <- d1[, c("kingdom", "super_class")] %>% distinct()
colnames(tmp) <- c("from", "to")
tmp1 <- d1[, c("super_class", "class")] %>% distinct()
colnames(tmp1) <- c("from", "to")
d2 <- rbind(tmp, tmp1)
tmp <- d1[, c("class", "Putative.Chemical.Class")] %>% distinct()
colnames(tmp) <- c("from", "to")
d2 <- rbind(d2, tmp)
label <- unlist(d2) %>% unique() %>% data.frame()
colnames(label) <- "class"
d2 <- Replace(data = d2, from = c(" ", ",", "\"", "\\(", "\\)", "-"), to = "_", pattern = c("\": "))
label$label <- Replace(data = label$class, from = c(" ", ",", "\"", "\\(", "\\)", "-"), to = "_", pattern = c("\": "))
label$type <- lapply(label$label, function(x) strsplit(x, "_") %>% unlist %>% .[1]) %>% unlist()
n <- c(table(d1$kingdom), table(d1$super_class), 
       table(d1$class), table(d1$Putative.Chemical.Class)) %>% as.data.frame()
colnames(n) <- "Members"
label <- merge(label, n , by.x = "class", by.y = "row.names")
label$class2 <- Replace(data = label$class, from = c("^SC_", "^C_", "^PC_"), to = "", pattern = c("\": "))
label$class3 <- Replace(data = label$class, from = c("^SC_", "^C_", "^PC_"), to = "")
label$colorlabel <- ifelse(label$Members >= 5, label$class2, "Others")
label$colorlabel2 <- ifelse(label$Members >= 10, label$class2, "Others")
label$colorlabel3 <- ifelse(label$Members >= 10 & label$type %in% "PC", label$class2, "Others")
label$colorlabel_class <- ifelse(label$Members >= 10 & label$type %in% "C", label$class2, NA)
ratio <- table(ibdclass3$Putative.Chemical.Class, ibdclass3$MMINP) %>% as.data.frame.array()
ratio$Predictable <- ratio$PPM + ratio$WPM
ratio$Ratio <- ratio$WPM / ratio$Predictable
label2 <- merge(label, ratio, by.x = "class3", by.y = "row.names", all.x = T)


d3 <- as.treedata(d2)
d32 <- as.phylo(d3)
d4 <- as_tibble(d3)
d5 <- merge(d4, label2, by = "label")
d5$label2 <- paste(d5$node, d5$class2, sep = ":")
mycolors <- c(rev(colors[1:length(unique(d5$colorlabel2))-1]))
names(mycolors) <- c(unique(d5$colorlabel3)[!unique(d5$colorlabel3) %in% "Others"], 
                     unique(d5$colorlabel2)[!unique(d5$colorlabel2) %in% unique(d5$colorlabel3)])
d5$colorlabel_class <- factor(d5$colorlabel_class, levels = unique(d5$colorlabel_class))
d5$offset <- ifelse(d5$type == "PC", 0.1,
                    ifelse(d5$type == "C", 0.3,
                           ifelse(d5$type == "SC", 0.5, NA)))
d5$extend <- ifelse(d5$type == "PC", 0.3,
                    ifelse(d5$type == "C", 0.5,
                           ifelse(d5$type == "SC", 0.7, NA)))
d5 <- arrange(d5, node)
p <- ggtree(d32, layout = "circular") %<+% d5 + 
  geom_tree(aes(size = Members)) + #color=colorlabel3,
  geom_tiplab(aes(subset=(colorlabel2 != "Others" & type == 'PC'), label = node), size = 3, hjust = 0, offset = 0.05)+
  geom_point(aes(subset=(colorlabel2 != "Others"  & type == 'PC'), color=colorlabel3), size=3)+
  labs(color = 'color', size = 'size')+
  scale_color_manual(values = mycolors, breaks = d5$colorlabel2, labels = d5$label2)#+
  guides(color = "none")

dta <- d5[d5$type %in% "PC", ] %>% arrange(desc(node))
dta 

p <- ggtree(d32, layout = "circular")
p <- flip(p, 28, 31)
p <- flip(p, 29, 34)
p <- flip(p, 30, 37)
p <- flip(p, 31, 40)
p <- flip(p, 32, 43)
p <- flip(p, 34, 46)
p <- flip(p, 57, 61)
p <- flip(p, 67, 72)
p <- flip(p, 209, 229)
p <- flip(p, 298, 307)
 p <- flip(p, 310, 312)
p %<+% d5 + geom_hilight(#data=td_filter(type=='C' & Members>=10), 
  mapping=aes(subset=(Members >= 10  & type %in% c('PC', 'C', 'SC')), extend = extend,
              node=node, fill=colorlabel2)#, extend=.3
)+
  scale_fill_manual(values = mycolors, breaks = d5$colorlabel2, labels = d5$label2)+
  geom_cladelab(
    mapping = aes(subset = (type %in% c('PC', 'C', 'SC') & Members >= 10), offset = offset,
                  node = node, label = node),
    geom = 'shadowtext',
    angle = "auto",
    horizontal = FALSE,
    size = 1.5,
    hjust = 0.5,
    #offset = 0.1,
    fontsize = 3,
    barsize = 0,
    barcolor = NA,
    bg.colour = 'white'
  )+
  geom_point2(aes(subset=(colorlabel2 != "Others"  & type == 'PC'), color=colorlabel3), size=2)+
  guides(color = "none")+
  new_scale_fill() +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=label, x= Ratio, fill = Predictable), 
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
             offset = 0.3)+
  scale_fill_continuous(low = "#b0d5df", 
                        high = "#0f95b0", trans = "log")+
  labs(fill = 'log(Predicted)')
ggsave("sup2_tree.pdf", height = 20, width = 20)
write.table(d5, "sup2_treedata.txt", row.names = F, sep = "\t", quote = F)

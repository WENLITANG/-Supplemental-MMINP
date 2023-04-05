# library(OmicsPLS)
# source("MMINP_predict.R")
# source("MMINP_preprocess.R")
# source("MMINP_train.R")
# metab <- openxlsx::read.xlsx("data/IBD/41564_2018_306_MOESM4_ESM.xlsx", sheet = 1, startRow = 2, rowNames = T)
# metab2 <- as.data.frame(t(apply(metab[8:nrow(metab), ], 2, as.numeric))) 
# colnames(metab2) <- rownames(metab)[8:nrow(metab)]
# metag <- openxlsx::read.xlsx("data/IBD/41564_2018_306_MOESM8_ESM.xlsx", sheet = 1, startRow = 2, rowNames = T)
# metag2 <- as.data.frame(t(apply(metag[9:nrow(metag), ], 2, as.numeric))) 
# colnames(metag2) <- rownames(metag)[9:nrow(metag)]
# metag2_pro <- as.data.frame(t(apply(metag2, 1, function(x) x/sum(x))))
# metab2_pro <- as.data.frame(t(apply(metab2, 1, function(x) x/sum(x))))

load("pre.RData")
library(OmicsPLS)
library(magrittr)
library(forecast)
trainS <- rownames(metag2)[grep("PRISM", rownames(metag2))]
for(k in 1:50){
  all_res <- matrix(nrow = 26, ncol = 12, 
                    dimnames = list(paste0("r", 1:26), 
                                    c("r", paste("fit0",c("train","tSp_n", "tSp_y", "pS_n", "pS_y"), sep = "."),
                                      paste("fit1",c("numb", "train", "tSp_n", "tSp_y", "pS_n", "pS_y"), sep = "."))))
  trainSt <- sample(trainS, 155*2/3)
  trainSp <- trainS[!trainS %in% trainSt]
  predS <- rownames(metag2)[grep("Validation", rownames(metag2))]
  #train
  mbt <- metab2_pro[trainSt, ]
  mgt <- metag2_pro[trainSt, ]
  mbtrm <- mbt[, apply(mbt, 2, function(x) length(which(x>0.0001))>nrow(mbt)*0.1)] #0.0001 * 1000000
  mgtrm <- mgt[, apply(mgt, 2, function(x) (length(which(x>0.0001))>nrow(mgt)*0.1))]
  mbtscale <- MMINP.preprocess(mbtrm, normalized = F, transformed = T, scaled = T)
  mgtscale <- MMINP.preprocess(mgtrm, normalized = F, transformed = T, scaled = T)
  components <- get_Components(mgtscale, mbtscale, compmethod = "cvo2m", n = 3:10, nx = 0:10, ny = 0:10,
                                nr_folds = 5, nr_cores = 10)
  fit0 <- o2m(mgtscale, mbtscale, n = as.numeric(components$n) , nx = as.numeric(components$nx) , ny = as.numeric(components$ny))
  #spearman r
  rcutoff <- c(0.1, 0.2, seq(0.3, 0.5, 0.01), seq(0.6, 0.8, 0.1))
  for (i in 1:length(rcutoff)) {
    r <- rcutoff[i]
    all_res[i, 1] <- r
    pred <- MMINP.predict(fit0, mgtscale)
    trainres <- compareFeatures(mbtscale, pred, rsignif = r)
    all_res[i, 2] <- length(trainres$wellPredicted)
    if(all_res[i, 2] < 5) next
    #validation
    mbm <- metab2_pro[predS, ]
    mgp <- metag2_pro[predS, ]
    mgpscale <- MMINP.preprocess(mgp[, rownames(fit0$W.)], normalized = F, transformed = T, scaled = T)
    mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = T, scaled = T)
    mbp <- MMINP.predict(fit0, mgpscale)
    predres <- compareFeatures(mbp, mbm)
    table(predres$res$signif)
    a <- predres$res[predres$res$compound %in% trainres$wellPredicted, ]
    table(a$signif)
    all_res[i, 5] <- length(grep("no", a$signif))
    all_res[i, 6] <- length(grep("yes", a$signif))
    #trainSp
    mbm <- metab2_pro[trainSp, ]
    mgp <- metag2_pro[trainSp, ]
    mgpscale <- MMINP.preprocess(mgp[, rownames(fit0$W.)], normalized = F, transformed = T, scaled = T)
    mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = T, scaled = T)
    mbp <- MMINP.predict(fit0, mgpscale)
    predres <- compareFeatures(mbp, mbm)
    table(predres$res$signif)
    a <- predres$res[predres$res$compound %in% trainres$wellPredicted, ]
    table(a$signif)
    all_res[i, 3] <- length(grep("no", a$signif))
    all_res[i, 4] <- length(grep("yes", a$signif))
    
    #iteration
    trainnumb <- 1
    while(length(trainres$wellPredicted) < ncol(pred)){
      metab_well <- mbtscale[, trainres$wellPredicted]
      # components <- get_Components(mgtscale, metab_well, compmethod = "cvo2m", n = 3:8, nx = 0:5, ny = 0:5,
      #                              nr_folds = 3, nr_cores = 10)
      if(length(trainres$wellPredicted)<5) break
      fit1 <- tryCatch({o2m(mgtscale, metab_well, n = as.numeric(components$n) , nx = as.numeric(components$nx) , ny = as.numeric(components$ny))},
               error = function(e) {
                 0})
      if(is.numeric(fit1)) break
      pred <- MMINP.predict(fit1, mgtscale)
      trainnumb <- trainnumb + 1
      print(trainnumb)
      trainres <- compareFeatures(metab_well, pred, rsignif = r)
    }
    all_res[i, 7] <- trainnumb - 1
    if(is.numeric(fit1)) next
    all_res[i, 8] <- length(trainres$wellPredicted)
    #trainSp
    mbm <- metab2_pro[trainSp, ]
    mgp <- metag2_pro[trainSp, ]
    mgpscale <- MMINP.preprocess(mgp[, rownames(fit1$W.)], normalized = F, transformed = T, scaled = T)
    mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = T, scaled = T)
    mbp <- MMINP.predict(fit1, mgpscale)
    predres <- compareFeatures(mbp, mbm)
    table(predres$res$signif)
    a <- predres$res[predres$res$compound %in% trainres$wellPredicted, ]
    table(a$signif)
    all_res[i, 9] <- length(grep("no", a$signif))
    all_res[i, 10] <- length(grep("yes", a$signif))
    #validation
    mbm <- metab2_pro[predS, ]
    mgp <- metag2_pro[predS, ]
    mgpscale <- MMINP.preprocess(mgp[, rownames(fit1$W.)], normalized = F, transformed = T, scaled = T)
    mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = T, scaled = T)
    mbp <- MMINP.predict(fit1, mgpscale)
    predres <- compareFeatures(mbp, mbm)
    table(predres$res$signif)
    a <- predres$res[predres$res$compound %in% trainres$wellPredicted, ]
    table(a$signif)
    all_res[i, 11] <- length(grep("no", a$signif))
    all_res[i, 12] <- length(grep("yes", a$signif))
  }
  all_res <- as.data.frame(all_res)
  all_res$fit0.train_all <- all_res$fit0.train / ncol(mbtscale)
  all_res$fit0.tSp_y_train <- all_res$fit0.tSp_y / all_res$fit0.train
  all_res$fit0.pS_y_train <- all_res$fit0.pS_y / all_res$fit0.train
  all_res$fit1.train_all <- all_res$fit1.train / ncol(mbtscale)
  all_res$fit1.tSp_y_train <- all_res$fit1.tSp_y / all_res$fit1.train
  all_res$fit1.pS_y_train <- all_res$fit1.pS_y / all_res$fit1.train
  write.table(all_res, paste0("rcutoff_selected/", k, "mbt", ncol(mbtscale), "_", paste(components, collapse = "_"), ".txt"), sep = "\t", quote = F, row.names = F)
}

####################### merge
filelist <- list.files(path = "rcutoff_selected/", pattern = ".txt", full.names = T)
res <- matrix(nrow = length(filelist), ncol = 19)
i <- 1
for (f in filelist) {
  m <- read.table(f, header = T, sep = "\t")
  m$absdiff <- abs(m$fit1.tSp_y_train - m$fit1.train_all)
  res[i, ] <- c(unlist(m[which.min(m$absdiff), ]))
  i <- i+1
}
res <- as.data.frame(res)
colnames(res) <- colnames(m)
res <- res[order(res$r), ]
write.table(res, "merged_rcutoff50.txt", row.names = F, sep = "\t", quote = F)

####################################################################### plot
library(ggplot2)
# sup1a
m <- read.table(filelist[27], header = T, sep = "\t") #"34mbt2817_1.84933269097928_8_6_5.txt"
m$absdiff <- abs(m$fit1.tSp_y_train - m$fit1.train_all)
d <- reshape2::melt(m[, c("r", "fit1.train_all", "fit1.tSp_y_train")], c("r"))
ggplot(d, aes(x = r, y = value, color = variable))+
  geom_point(size = 2, alpha = 0.8)+
  geom_line(size = 0.75)+
  theme_bw()+
  theme(text = element_text(size = 20))+ylab("Ration")+
  xlab("Thread of well-fitted metabolites")+
  scale_color_manual(values = c("#00acc1", "#e64a19"), 
                     breaks = c("fit1.train_all", "fit1.tSp_y_train"), 
                     labels = c("Well-fitted/All", "Well-predicted/Predicted"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.3, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave("sup1a_rcutoff_example.pdf", height = 5.4, width = 6) 

# sup1c
res$group <- cut(res$r, breaks = c(0.25, 0.35, 0.45, 0.55))
ggplot(res, aes(fit1.tSp_y_train, fit1.train_all, fill = group))+
  geom_point(size = 3, alpha = 0.9, show.legend = T, shape = 21)+
  scale_fill_manual(values = c("#e6cecc", "#AD615B", "#B2B6BF"))+
  xlab("Ratio of well-predicted to predicted metabolites")+
  ylab("Ratio of well-fitted to all metabolites")+
  # xlim(0, 1)+
  # ylim(0, 1)+
  theme_bw()
ggsave("sup1c_rcutoff_ratio.pdf", height = 5.4, width = 6) 

# sup1b
library(dplyr)
library(shadowtext)
summary(res$r)
numb <- cut(res$r, breaks = c(0.25, 0.35, 0.45, 0.55)) %>% table() %>% as.data.frame()
numb$x <- c(0.3, 0.4, 0.5)
ggplot(data=res,aes(x=r)) + 
  geom_histogram(binwidth=0.1,fill=c("#e6cecc", "#AD615B", "#B2B6BF"), 
                 color="#e9ecef", alpha=0.9)+
  geom_shadowtext(data = numb, aes(x = x, y = Freq, label= Freq), 
                  size =5, bg.r = 0.1, bg.colour = "white", color = "black", vjust = 1.2)+
  theme_bw()+
  xlab("Threshold of well-fitted metabolites")+
  ylab("Frequence")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, colour = "black"))
ggsave("sup1b_rcutoff_histogram.pdf", height = 5.4, width = 4.2)

library(MMINP)
library(OmicsPLS)
library(melonnpan)
source("melonnpan_predict_modified.R")
source("ENVIM-main/ENVIM.R")
load("data/D1'_D2_D3_data.RData")
smoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  return(x)
}
mg_pro <- as.data.frame(t(apply(f_ibdg, 1, function(x) x/sum(x))))
mb_pro <- as.data.frame(t(apply(f_ibdb, 1, function(x) x/sum(x))))
m <- read.table("data/IBD_meta.txt", header = T, row.names = 1, sep = "\t")
metag2_pro <- mg_pro[rownames(m), ]
metab2_pro <- mb_pro[rownames(m), ]
metag2_pro <- metag2_pro[, colSums(metag2_pro)>0]
metab2_pro <- metab2_pro[, colSums(metab2_pro)>0]
set.seed(1234)
predS <- sample(rownames(metag2_pro), 65)
trainSs <- rownames(metag2_pro)[!rownames(metag2_pro) %in% predS] 

dir.create("sample_size")
dir.create("sample_size/Melon")
dir.create("sample_size/ENVIM")
nlist <- seq(20, 150, 10)
a <- c("trainnumb", "iteration", "prednumb", "samples*genes*metabolites",
       "components", paste(rep(c("O2PLS", "MMINP", "Melon", "ENVIM"), each = 4), 
                           c("welltrain", "wellpred", "welltrain_ratio", "wellpred_ratio"), sep = "_"))
a <- as.data.frame(matrix(a, nrow = 1))
write.table(a, "sample_size/results.txt", sep = "\t", quote = F, row.names = F, col.names = F)

res <- matrix(ncol = 21, nrow = length(nlist)*50)
k <- 1
for(ii in c(1, 50)){
  for(i in 1:length(nlist)){
    n <- nlist[i]
    set.seed(ii)
    trainS <- sample(trainSs, n)
    mbt <- metab2_pro[trainS, ]
    mgt <- metag2_pro[trainS, ]
    mbtrm <- mbt[, apply(mbt, 2, function(x) length(which(x>0.0001))>nrow(mbt)*0.1)] #0.0001 * 1000000
    mgtrm <- mgt[, apply(mgt, 2, function(x) (length(which(x>0.0001))>nrow(mgt)*0.1))]
    mbtscale <- MMINP.preprocess(mbtrm, normalized = F, transformed = 'boxcox', scaled = T)
    mgtscale <- MMINP.preprocess(mgtrm, normalized = F, transformed = 'boxcox', scaled = T)
    dim(mbtscale)
    dim(mgtscale)
    mbm <- metab2_pro[predS, ]
    mgp <- metag2_pro[predS, ]
    mgp <- mgp[, colSums(mgp)>0]
    dim(mgp)
    mgpscale <- MMINP.preprocess(mgp, normalized = F, transformed = 'boxcox', scaled = T)
    mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = 'boxcox', scaled = T)
    
    ## O2PLS, MMINP
    if(n < 60){
      components <- get_Components(mgtscale, mbtscale, compmethod = "cvo2m", n = 3:10, nx = 0:10, ny = 0:10,
                                   nr_folds = 2, nr_cores = 5) 
    }else if(n <= 100){
      components <- get_Components(mgtscale, mbtscale, compmethod = "cvo2m", n = 3:10, nx = 0:10, ny = 0:10,
                                   nr_folds = 3, nr_cores = 5) 
    }else{
      components <- get_Components(mgtscale, mbtscale, compmethod = "cvo2m", n = 3:10, nx = 0:10, ny = 0:10,
                                   nr_folds = 5, nr_cores = 5) 
    }
    components
    res[k, 1] <- n
    res[k, 2] <- ii
    res[k, 3] <- length(predS)
    res[k, 4] <- paste(c(dim(mgtscale), ncol(mbtscale)), collapse = "*")
    res[k, 5] <- paste(components, collapse = " ")
    fit0 <- o2m(mgtscale, mbtscale, n = as.numeric(components$n) , nx = as.numeric(components$nx) , ny = as.numeric(components$ny))
    ### 1.O2PLS
    mbp0 <- MMINP.predict(fit0, mgpscale)
    predres0 <- compareFeatures(mbp0, mbmscale)
    res[k, 6] <- nrow(predres0$res)
    res[k, 7] <- length(predres0$wellPredicted) #paste(table(predres$res$signif), collapse = " vs ")
    res[k, 8] <- nrow(predres0$res) / ncol(mbtscale)
    res[k, 9] <- length(predres0$wellPredicted) / nrow(predres0$res)
    
    ### 2.MMINP
    #spearman r
    r <- 0.4 
    pred <- predict(fit0, mgtscale, XorY = c("X"))
    trainres <- compareFeatures(mbtscale, pred, rsignif = r)
    trainnumb <- 1
    while(length(trainres$wellPredicted) < ncol(pred)){
      metab_well <- mbtscale[, trainres$wellPredicted]
      # components <- get_Components(mgtscale, metab_well, compmethod = "cvo2m", n = 3:8, nx = 0:5, ny = 0:5,
      #                              nr_folds = 3, nr_cores = 10)
      fit1 <- o2m(mgtscale, metab_well, n = as.numeric(components$n) , nx = as.numeric(components$nx) , ny = as.numeric(components$ny))
      pred <- predict(fit1, mgtscale, XorY = c("X"))
      trainnumb <- trainnumb + 1
      print(trainnumb)
      trainres <- compareFeatures(metab_well, pred, rsignif = r)
    }
    length(trainres$wellPredicted)
    mbp <- MMINP.predict(fit1, mgpscale[, colnames(mgpscale) %in% rownames(fit1$W.)])
    predres <- compareFeatures(mbp, mbmscale)
    table(predres$res$signif)
    res[k, 10] <- length(trainres$wellPredicted)
    res[k, 11] <- length(predres$wellPredicted) #paste(table(predres$res$signif), collapse = " vs ")
    res[k, 12] <- length(trainres$wellPredicted) / ncol(mbtscale)
    res[k, 13] <- length(predres$wellPredicted) / length(trainres$wellPredicted)

    ### 3.melonnpan
    dir.create(melondir <- paste0("sample_size/Melon/", n, "_", ii, "/"))
    melonnpan.train(metab = mbtrm, metag = mgtrm, outputDirectory = melondir, cores = 5)
    weight <- read.table(paste0(melondir, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
    colnames(weight) <- gsub("\\.", "-", colnames(weight)) 
    dim(weight)
    res[k, 14] <- ncol(weight)
    metabpredict <- melonnpan.predict.modified(metag = mgp, weight.matrix = weight, train.metag = mgtrm,
                                               corr.method = "spearman", output = paste0(melondir, "pred"))
    melonn_pred <- metabpredict$pred
    rownames(melonn_pred) <- melonn_pred$ID
    melonn_pred$ID <- NULL
    melonn_res <- compareFeatures(melonn_pred, mbm[rownames(melonn_pred), colnames(mbm) %in% colnames(melonn_pred)])
    table(melonn_res$res$signif) 
    res[k, 15] <- length(melonn_res$wellPredicted)
    res[k, 16] <- ncol(weight) / ncol(mbtscale)
    res[k, 17] <- length(melonn_res$wellPredicted) / ncol(weight)
    
    ### 4.ENVIM
    metab.train <- apply(mbtrm, 2, smoothZero)
    metab.test <- mbm[, colnames(mbtrm)]
    gid <- intersect(colnames(mgtrm), colnames(mgp))
    bid <- intersect(colnames(metab.train), colnames(metab.test))
    dir.create(output <- paste0("sample_size/ENVIM/", n, "_", ii, "/"))
    ENVIM(microbio.train = mgtrm[, gid],
          microbio.test = mgp[, gid],
          metab.train = metab.train[, bid],
          metab.test = metab.test[, bid],
          seed = 1234,
          outputdirectory = output,
          fold_rf = 10,
          fold_ENVIM = 10)
    envimres <- read.csv("summary.csv", header = T)
    res[k, 18] <- sum(envimres$Train_spearman_cor>0.3, na.rm = T)
    res[k, 19] <- sum((envimres$Train_spearman_cor>0.3) & (envimres$Test_spearman_cor > 0.3), na.rm = T)
    res[k, 20] <- sum(envimres$Train_spearman_cor>0.3, na.rm = T) / ncol(mbtscale)
    res[k, 21] <- sum((envimres$Train_spearman_cor>0.3) & (envimres$Test_spearman_cor > 0.3), na.rm = T) / sum(envimres$Train_spearman_cor>0.3, na.rm = T)
    setwd("../../../")
    
    ## this loop, done
    write.table(as.data.frame(matrix(res[k, ], nrow = 1)), "sample_size/results.txt", sep = "\t", quote = F, row.names = F, col.names = F, append = T)
    
    ## next
    k <- k+1
    print(ii)
    print(n)
  }
}
res <- as.data.frame(res)
colnames(res) <- a
write.table(res, "sample_size/all_results.txt", sep = "\t", quote = F, row.names = F)

############################################################################################################### plot
library(ggplot2)
library(patchwork)
m <- read.table("sample_size/all_results.txt", header = T, sep = "\t")
n1 <- c("trainnumb", paste(c("O2PLS", "MMINP", "Melon", "ENVIM"), "welltrain_ratio", sep = "_"))
m1 <- reshape2::melt(m, id.vars = "trainnumb", measure.vars = n1[-1])
m1$trainnumb <- factor(m1$trainnumb, levels = seq(20, 150, 10))
p1 <- ggplot(m1, aes(x = trainnumb, y = value, color = variable))+
  geom_boxplot(outlier.size = 1, position = position_dodge(0.9))+
  theme_bw()+
  #scale_color_nejm()+
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"))+
  ylab("Predicted/All ratio")+
  xlab("Size of training samples")+
  scale_color_manual(values = c("#b79f00", "#e64a19", "#00d9f4", "#004ac7"))+#ylim(0,1)
  geom_smooth(aes(x = as.numeric(m1[,"trainnumb"])), position = position_dodge(0.9))
#ggsave("welltrainRatio_result.pdf", width = 8, height = 5)

n2 <- c("trainnumb", paste(c("O2PLS", "MMINP", "Melon", "ENVIM"), "wellpred_ratio", sep = "_"))
m2 <- reshape2::melt(m, id.vars = "trainnumb", measure.vars = n2[-1])
m2$trainnumb <- factor(m2$trainnumb, levels = seq(20, 150, 10))
p2 <- ggplot(m2, aes(x = trainnumb, y = value, color = variable))+
  geom_boxplot(outlier.size = 1, position = position_dodge(0.9))+
  theme_bw()+
  #scale_color_nejm()+
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"))+
  ylab("Wellpredicted/Predicted ratio")+
  xlab("Size of training samples")+
  scale_color_manual(values = c("#b79f00", "#e64a19", "#00d9f4", "#004ac7"))+#ylim(0,1)
  geom_smooth(aes(x = as.numeric(m2[,"trainnumb"])), position = position_dodge(0.9))

m$all <- m$O2PLS_welltrain
m$O2PLS_wellpred_All <- m$O2PLS_wellpred / m$all
m$MMINP_wellpred_All <- m$MMINP_wellpred / m$all
m$Melon_wellpred_All <- m$Melon_wellpred / m$all
m$ENVIM_wellpred_All <- m$ENVIM_wellpred / m$all
#write.table(m, "samplesize.txt", row.names = F, sep = "\t", quote = F)
n3 <- c("trainnumb", paste(c("O2PLS", "MMINP", "Melon", "ENVIM"), "wellpred_All", sep = "_"))
m3 <- reshape2::melt(m, id.vars = "trainnumb", measure.vars = n3[-1])
m3$trainnumb <- factor(m3$trainnumb, levels = seq(20, 150, 10))
p3 <- ggplot(m3, aes(x = trainnumb, y = value, color = variable))+
  geom_boxplot(outlier.size = 1, position = position_dodge(0.9))+
  theme_bw()+
  #scale_color_nejm()+
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"))+
  ylab("Wellpredicted/All ratio")+
  xlab("Size of training samples")+
  scale_color_manual(values = c("#b79f00", "#e64a19", "#00d9f4", "#004ac7"))+#ylim(0,1)
  geom_smooth(aes(x = as.numeric(m3[,"trainnumb"])), position = position_dodge(0.9))

p1+p3+p2+plot_layout(guides = 'collect')
ggsave("fig5_samplesize.pdf", width = 13.2, height = 4)

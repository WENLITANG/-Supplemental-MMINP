library(MMINP)
library(OmicsPLS)
library(melonnpan)
source("melonnpan_predict_modified.R")
source("ENVIM-main/ENVIM.R")
source("ENVIM-main/ENVIM_predict.R")
load("data/D1'_D2_D3_data.RData")
smoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  return(x)
}

mg_pro <- as.data.frame(t(apply(f_ibdg, 1, function(x) x/sum(x))))
mb_pro <- as.data.frame(t(apply(f_ibdb, 1, function(x) x/sum(x))))
m <- read.table("script_and_meta/IBD_meta.txt", header = T, row.names = 1, sep = "\t")
metag2_pro <- mg_pro[rownames(m), ]
metab2_pro <- mb_pro[rownames(m), ]
metag2_pro <- metag2_pro[, colSums(metag2_pro)>0]
metab2_pro <- metab2_pro[, colSums(metab2_pro)>0]
control <- rownames(m)[m$Diagnosis %in% "Control"]
uc <- rownames(m)[m$Diagnosis %in% "UC"]
cd <- rownames(m)[m$Diagnosis %in% "CD"]
conbt <- metab2_pro[control, ]
conbt_filter <- conbt[, apply(conbt, 2, function(x) length(which(x>0.0001))>nrow(conbt)*0.1)]
congt <- metag2_pro[control, ]
congt_filter <- congt[, apply(congt, 2, function(x) length(which(x>0.0001))>nrow(congt)*0.1)]
ucbt <- metab2_pro[uc, ]
ucbt_filter <- ucbt[, apply(ucbt, 2, function(x) length(which(x>0.0001))>nrow(ucbt)*0.1)]
ucgt <- metag2_pro[uc, ]
ucgt_filter <- ucgt[, apply(ucgt, 2, function(x) length(which(x>0.0001))>nrow(ucgt)*0.1)]
cdbt <- metab2_pro[cd, ]
cdbt_filter <- cdbt[, apply(cdbt, 2, function(x) length(which(x>0.0001))>nrow(cdbt)*0.1)]
cdgt <- metag2_pro[cd, ]
cdgt_filter <- cdgt[, apply(cdgt, 2, function(x) length(which(x>0.0001))>nrow(cdgt)*0.1)]

##################################################################
########################################## MMINP
#con
conbtscale <- MMINP.preprocess(conbt_filter, normalized = F, transformed = "boxcox")
congtscale <- MMINP.preprocess(data = congt_filter, normalized = F, transformed = 'boxcox', scaled = T)
mcon <- MMINP.train(metag = congtscale, metab = conbtscale,
                    n = 3:10, nx = 0:10, ny = 0:10,
                    nr_folds = 3, nr_cores = 10, seed = 1234)
mcon 

#uc
ucbtscale <- MMINP.preprocess(ucbt_filter, normalized = F, transformed = 'boxcox')
ucgtscale <- MMINP.preprocess(data = ucgt_filter, normalized = F, transformed = 'boxcox', scaled = T)
muc <- MMINP.train(metag = ucgtscale, metab = ucbtscale,
                    n = 3:10, nx = 0:10, ny = 0:10,
                    nr_folds = 3, nr_cores = 10, seed = 1234)
muc 

#cd
cdbtscale <- MMINP.preprocess(cdbt_filter, normalized = F, transformed = 'boxcox')
cdgtscale <- MMINP.preprocess(data = cdgt_filter, normalized = F, transformed = 'boxcox', scaled = T)
mcd <- MMINP.train(metag = cdgtscale, metab = cdbtscale,
                   n = 3:10, nx = 0:10, ny = 0:10, 
                   nr_folds = 3, nr_cores = 10, seed = 1234)
mcd 
#-----------------------------predict other data set
#conmodel
mcon_cdp <- MMINP.predict(mcon, newdata = cdgtscale, minGeneSize = 0.3)
mcon_rescd <- compareFeatures(predicted = mcon_cdp, measured = cdbtscale)
length(mcon_rescd$wellPredicted) 

mcon_ucp <- MMINP.predict(mcon, newdata = ucgtscale)
mcon_resuc <- compareFeatures(predicted = mcon_ucp, measured = ucbtscale)
length(mcon_resuc$wellPredicted) 

#cdmodel
mcd_conp <- MMINP.predict(mcd, newdata = congtscale, minGeneSize = 0.2)
mcd_rescon <- compareFeatures(predicted = mcd_conp, measured = conbtscale)
length(mcd_rescon$wellPredicted) 

mcd_ucp <- MMINP.predict(mcd, newdata = ucgtscale, minGeneSize = 0.3)
mcd_resuc <- compareFeatures(predicted = mcd_ucp, measured = ucbtscale)
length(mcd_resuc$wellPredicted) 

#ucmodel
muc_conp <- MMINP.predict(muc, newdata = congtscale, minGeneSize = 0.2)
muc_rescon <- compareFeatures(predicted = muc_conp, measured = conbtscale)
length(muc_rescon$wellPredicted) 

muc_cdp <- MMINP.predict(muc, newdata = cdgtscale)
muc_rescd <- compareFeatures(predicted = muc_cdp, measured = cdbtscale)
length(muc_rescd$wellPredicted) 


########################################## O2PLS
mcon_o2pls <- o2m(X = congtscale, Y = conbtscale, 
                  n = mcon$components$n, nx = mcon$components$nx, ny = mcon$components$ny)

muc_o2pls <- o2m(X = ucgtscale, Y = ucbtscale, 
                  n = muc$components$n, nx = muc$components$nx, ny = muc$components$ny)

mcd_o2pls <- o2m(X = cdgtscale, Y = cdbtscale, 
                 n = mcd$components$n, nx = mcd$components$nx, ny = mcd$components$ny)
#-----------------------------predict other data set
#conmodel
mcon_o2pls_cdp <- MMINP.predict(mcon_o2pls, newdata = cdgtscale)
mcon_o2pls_rescd <- compareFeatures(predicted = mcon_o2pls_cdp, measured = cdbtscale)
length(mcon_o2pls_rescd$wellPredicted) 

mcon_o2pls_ucp <- MMINP.predict(mcon_o2pls, newdata = ucgtscale)
mcon_o2pls_resuc <- compareFeatures(predicted = mcon_o2pls_ucp, measured = ucbtscale)
length(mcon_o2pls_resuc$wellPredicted) 

#cdmodel
mcd_o2pls_conp <- MMINP.predict(mcd_o2pls, newdata = congtscale, minGeneSize = 0.2)
mcd_o2pls_rescon <- compareFeatures(predicted = mcd_o2pls_conp, measured = conbtscale)
length(mcd_o2pls_rescon$wellPredicted) 

mcd_o2pls_ucp <- MMINP.predict(mcd_o2pls, newdata = ucgtscale)
mcd_o2pls_resuc <- compareFeatures(predicted = mcd_o2pls_ucp, measured = ucbtscale)
length(mcd_o2pls_resuc$wellPredicted) 

#ucmodel
muc_o2pls_conp <- MMINP.predict(muc_o2pls, newdata = congtscale, minGeneSize = 0.2)
muc_o2pls_rescon <- compareFeatures(predicted = muc_o2pls_conp, measured = conbtscale)
length(muc_o2pls_rescon$wellPredicted) 

muc_o2pls_cdp <- MMINP.predict(muc_o2pls, newdata = cdgtscale)
muc_o2pls_rescd <- compareFeatures(predicted = muc_o2pls_cdp, measured = cdbtscale)
length(muc_o2pls_rescd$wellPredicted) 


########################################## melonnpan
dir.create("sample_status/")
dir.create(dir_con <- "sample_status/con/")
melonnpan.train(metab = conbt_filter, metag = congt_filter, outputDirectory = dir_con, cores = 10, nfolds = 3)
weight_con <- read.table(paste0(dir_con, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_con)

dir.create(dir_uc <- "sample_status/uc/")
melonnpan.train(metab = ucbt_filter, metag = ucgt_filter, outputDirectory = dir_uc, cores = 10, nfolds = 3)
weight_uc <- read.table(paste0(dir_uc, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_uc) 

dir.create(dir_cd <- "sample_status/cd/")
melonnpan.train(metab = cdbt_filter, metag = cdgt_filter, outputDirectory = dir_cd, cores = 10, nfolds = 3)
weight_cd <- read.table(paste0(dir_cd, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_cd) 
#--------------------------------------predict other data set
#conmodel
meloncon_cd <- melonnpan.predict.modified(metag = cdgt, weight.matrix = weight_con, train.metag = congt_filter,
                                          corr.method = "spearman", output = paste0(dir_con, "cd"))
meloncon_cdp <- meloncon_cd$pred
rownames(meloncon_cdp) <- meloncon_cdp$ID
meloncon_cdp$ID <- NULL
melonconres_cdp <- compareFeatures(predicted = meloncon_cdp, measured = cdbt)
length(melonconres_cdp$wellPredicted)

meloncon_uc <- melonnpan.predict.modified(metag = ucgt, weight.matrix = weight_con, train.metag = congt_filter,
                                           corr.method = "spearman", output = paste0(dir_con, "uc"))
meloncon_ucp <- meloncon_uc$pred
rownames(meloncon_ucp) <- meloncon_ucp$ID
meloncon_ucp$ID <- NULL
melonconres_ucp <- compareFeatures(predicted = meloncon_ucp, measured = ucbt)
length(melonconres_ucp$wellPredicted) 

#cdmodel
meloncd_con <- melonnpan.predict.modified(metag = congt, weight.matrix = weight_cd, train.metag = cdgt_filter,
                                          corr.method = "spearman", output = paste0(dir_cd, "con"))
meloncd_conp <- meloncd_con$pred
rownames(meloncd_conp) <- meloncd_conp$ID
meloncd_conp$ID <- NULL
meloncdres_conp <- compareFeatures(predicted = meloncd_conp, measured = conbt)
length(meloncdres_conp$wellPredicted) 

meloncd_uc <- melonnpan.predict.modified(metag = ucgt, weight.matrix = weight_cd, train.metag = cdgt_filter,
                                          corr.method = "spearman", output = paste0(dir_cd, "uc"))
meloncd_ucp <- meloncd_uc$pred
rownames(meloncd_ucp) <- meloncd_ucp$ID
meloncd_ucp$ID <- NULL
meloncdres_ucp <- compareFeatures(predicted = meloncd_ucp, measured = ucbt)
length(meloncdres_ucp$wellPredicted) 

#ucmodel
melonuc_con <- melonnpan.predict.modified(metag = congt, weight.matrix = weight_uc, train.metag = ucgt_filter,
                                           corr.method = "spearman", output = paste0(dir_uc, "con"))
melonuc_conp <- melonuc_con$pred
rownames(melonuc_conp) <- melonuc_conp$ID
melonuc_conp$ID <- NULL
melonucres_conp <- compareFeatures(predicted = melonuc_conp, measured = conbt)
length(melonucres_conp$wellPredicted)

melonuc_cd <- melonnpan.predict.modified(metag = cdgt, weight.matrix = weight_uc, train.metag = ucgt_filter,
                                          corr.method = "spearman", output = paste0(dir_uc, "cd"))
melonuc_cdp <- melonuc_cd$pred
rownames(melonuc_cdp) <- melonuc_cdp$ID
melonuc_cdp$ID <- NULL
melonucres_cdp <- compareFeatures(predicted = melonuc_cdp, measured = cdbt)
length(melonucres_cdp$wellPredicted) 

########################################## ENVIM
#-----------------------------predict other data set
conb.train <- apply(conbt_filter, 2, smoothZero)
dir.create(output <- paste0("sample_status/ENVIM_con_uc/"))
ENVIM_predict(microbio.train = congt_filter,
              microbio.test = ucgt[, apply(ucgt, 2, function(x) sum(x>0)>2)],
              metab.train = conb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_con_uc <- read.csv("sample_status/ENVIM_con_uc/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimconuc <- compareFeatures(envim_con_uc, ucbt)
length(res_envimconuc$wellPredicted)

dir.create(output <- paste0("sample_status/ENVIM_con_cd/"))
ENVIM_predict(microbio.train = congt_filter,
              microbio.test = cdgt[, apply(cdgt, 2, function(x) sum(x>0)>2)],
              metab.train = conb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_con_cd <- read.csv("sample_status/ENVIM_con_cd/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimconcd <- compareFeatures(envim_con_cd, cdbt)
length(res_envimconcd$wellPredicted)


ucb.train <- apply(ucbt_filter, 2, smoothZero)
dir.create(output <- paste0("sample_status/ENVIM_uc_cd/"))
ENVIM_predict(microbio.train = ucgt_filter,
              microbio.test = cdgt[, apply(cdgt, 2, function(x) sum(x>0)>2)],
              metab.train = ucb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_uc_cd <- read.csv("sample_status/ENVIM_uc_cd/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimucdc <- compareFeatures(envim_uc_cd, cdbt)
length(res_envimucdc$wellPredicted)

dir.create(output <- paste0("sample_status/ENVIM_uc_con/"))
ENVIM_predict(microbio.train = ucgt_filter,
              microbio.test = congt[, apply(congt, 2, function(x) sum(x>0)>2)],
              metab.train = ucb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_uc_con <- read.csv("sample_status/ENVIM_uc_con/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimucdon <- compareFeatures(envim_uc_con, conbt)
length(res_envimucdon$wellPredicted)

cdb.train <- apply(cdbt_filter, 2, smoothZero)
dir.create(output <- paste0("sample_status/ENVIM_cd_con/"))
ENVIM_predict(microbio.train = cdgt_filter,
              microbio.test = congt[, apply(congt, 2, function(x) sum(x>0)>2)],
              metab.train = cdb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_cd_con <- read.csv("sample_status/ENVIM_cd_con/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimcdcon <- compareFeatures(envim_cd_con, conbt)
length(res_envimcdcon$wellPredicted)

dir.create(output <- paste0("sample_status/ENVIM_cd_uc/"))
ENVIM_predict(microbio.train = cdgt_filter,
              microbio.test = ucgt[, apply(ucgt, 2, function(x) sum(x>0)>2)],
              metab.train = cdb.train,
              #metab.test = conbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_cd_uc <- read.csv("sample_status/ENVIM_cd_uc/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimcduc <- compareFeatures(envim_cd_uc, ucbt)
length(res_envimcduc$wellPredicted)

##########################################
save.image(file = "sample_status/all.RData")

############################################################################################# plot
#MMINP
mcon
muc
mcd
length(mcon_resuc$wellPredicted)
length(mcon_rescd$wellPredicted)
length(muc_rescon$wellPredicted)
length(muc_rescd$wellPredicted)
length(mcd_rescon$wellPredicted)
length(mcd_resuc$wellPredicted)

#o2pls
length(mcon_o2pls_resuc$wellPredicted)
length(mcon_o2pls_rescd$wellPredicted)
length(muc_o2pls_rescon$wellPredicted)
length(muc_o2pls_rescd$wellPredicted)
length(mcd_o2pls_rescon$wellPredicted)
length(mcd_o2pls_resuc$wellPredicted)

#melonnpan
length(weight_con)
length(weight_uc)
length(weight_cd)
length(melonconres_ucp$wellPredicted)
length(melonconres_cdp$wellPredicted)
length(melonucres_conp$wellPredicted)
length(melonucres_cdp$wellPredicted)
length(meloncdres_conp$wellPredicted)
length(meloncdres_ucp$wellPredicted)

#ENVIM
length(res_envimconuc$wellPredicted)
length(res_envimconcd$wellPredicted)
length(res_envimucdon$wellPredicted)
length(res_envimucdc$wellPredicted)
length(res_envimcdcon$wellPredicted)
length(res_envimcduc$wellPredicted)


library(ggplot2)
mm <- read.table("sample_conditions.txt", header = T, sep = "\t")
mm$model <- factor(mm$model, levels = c("Control", "UC", "CD"))
mm$predicted <- factor(mm$predicted, levels = c("numb", "Control", "UC", "CD"))
d <- reshape2::melt(mm)
ggplot(d, aes(model, value, fill = predicted))+
  geom_bar(stat = 'identity', position = position_dodge(0.9), width = 0.8, alpha = 0.8)+
  scale_fill_manual(values = c("numb" = "#546e7a", 
                               "CD" = "#e64a19", #"#E5352B",
                               "Control" =  "#00acc1", #"#0081B4", 
                               "UC" =  "#f9a825"))+
  theme_bw()+
  theme(axis.title = element_text(size = 18, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  facet_wrap(~variable)+
  xlab("Model")+
  ylab("Number")
ggsave("fig6_sample_status.pdf", width = 8.4, height = 8)

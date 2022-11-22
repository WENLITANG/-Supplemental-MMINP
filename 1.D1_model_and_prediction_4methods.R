library(MMINP)
library(OmicsPLS)
library(forecast)
m <- read.table("data/IBD/IBD_meta.txt", header = T, row.names = 1, sep = "\t", comment.char = "")
metab <- openxlsx::read.xlsx("data/IBD/41564_2018_306_MOESM4_ESM.xlsx", sheet = 1, startRow = 2, rowNames = T)
metab2 <- as.data.frame(t(apply(metab[8:nrow(metab), ], 2, as.numeric))) 
colnames(metab2) <- rownames(metab)[8:nrow(metab)]
metag <- openxlsx::read.xlsx("data/IBD/41564_2018_306_MOESM8_ESM.xlsx", sheet = 1, startRow = 2, rowNames = T)
metag2 <- as.data.frame(t(apply(metag[9:nrow(metag), ], 2, as.numeric))) 
colnames(metag2) <- rownames(metag)[9:nrow(metag)]
metag2_pro <- as.data.frame(t(apply(metag2, 1, function(x) x/sum(x))))
metab2_pro <- as.data.frame(t(apply(metab2, 1, function(x) x/sum(x))))

#samples
trainS <- rownames(metag2)[grep("PRISM", rownames(metag2))]
predS <- rownames(metag2)[grep("Validation", rownames(metag2))]

#train model
mbt <- metab2_pro[trainS, ]
mgt <- metag2_pro[trainS, ]
mbtrm <- mbt[, apply(mbt, 2, function(x) length(which(x>0.0001))>nrow(mbt)*0.1)] #0.0001 * 1000000
mgtrm <- mgt[, apply(mgt, 2, function(x) (length(which(x>0.0001))>nrow(mgt)*0.1))]
mbtscale <- MMINP.preprocess(mbtrm, normalized = F, transformed = 'boxcox', scaled = T)
mgtscale <- MMINP.preprocess(mgtrm, normalized = F, transformed = 'boxcox', scaled = T)

#MMINP
##train
components <- get_Components(mgtscale, mbtscale, compmethod = "cvo2m", n = 3:10, nx = 0:10, ny = 0:10,
                             nr_folds = 5, nr_cores = 10)
components
# minCV nx ny n
# 1 1.84129  6  4 6
model1 <- MMINP.train(mgtscale, mbtscale, n = as.numeric(components$n), 
                      nx = as.numeric(components$nx), ny = as.numeric(components$ny), 
                      nr_folds = 3, nr_cores = 3)
model1 #2015

##predict
mbm <- metab2_pro[predS, ]
mgp <- metag2_pro[predS, ]
mgpscale <- MMINP.preprocess(mgp, normalized = F, transformed = 'boxcox', scaled = T)
mbmscale <- MMINP.preprocess(mbm, normalized = F, transformed = 'boxcox', scaled = T)
mbp <- MMINP.predict(model1, mgpscale)
predres <- compareFeatures(mbp, mbm)
length(predres$wellPredicted) #1233
predres2 <- compareFeatures(mbp, mbmscale)
length(predres2$wellPredicted) #1233

#vipres <- O2PLSvip(mgtscale, mbtscale, model = model1)

#########################O2PLS
fit0 <- o2m(mgtscale, mbtscale, n = as.numeric(components$n),
            nx = as.numeric(components$nx), ny = as.numeric(components$ny))
vipres0 <- O2PLSvip(mgtscale, mbtscale, model = fit0)
o2pls_pred <- MMINP.predict(fit0, mgpscale)
o2pls_res <- compareFeatures(o2pls_pred, mbmscale)
length(o2pls_res$wellPredicted)

###### melonnpan
library(melonnpan)
source("melonnpan_predict_modified.R")
dir.create(melondir <- "IBD_melon/")
melonnpan.train(metab = mbtrm, metag = mgtrm, outputDirectory = melondir, cores = 10)
weight <- read.table(paste0(melondir, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
colnames(weight) <- gsub("\\.", "-", colnames(weight))
dim(weight) #1664
metabpredict <- melonnpan.predict.modified(metag = mgp, weight.matrix = weight, train.metag = mgtrm,
                                           corr.method = "spearman", output = paste0(melondir, "predict5"))
melonn_pred <- metabpredict$pred
rownames(melonn_pred) <- melonn_pred$ID
melonn_pred$ID <- NULL
melonn_res <- compareFeatures(predicted = melonn_pred, measured = mbm[rownames(melonn_pred), colnames(melonn_pred)])
length(melonn_res$wellPredicted) #976

###### ENVIM
# Run on the Terminal:
# nohup Rscript.exe splitD1MetabENVIM.R -n 0 > 0.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 1 > 1.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 2 > 2.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 3 > 3.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 4 > 4.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 5 > 5.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 6 > 6.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 7 > 7.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 8 > 8.log 1>&1 &
# nohup Rscript.exe splitD1MetabENVIM.R -n 9 > 9.log 1>&1 &
# 
# R:
# filelist <- list.files("IBD_ENVIM/", pattern = "^summary.csv$", recursive = T, full.names = T)
# envim_res <- c()
# for(f in filelist){
#   tmp <- read.csv(f)
#   envim_res <- rbind(envim_res, tmp)
# }
# write.table(envim_res, "IBD_ENVIM/ENVIM_merged_summary.csv", row.names = F, sep = "\t", quote = F)
envim_res <- read.csv("IBD_ENVIM/ENVIM_merged_summary.csv", header = T)
rownames(envim_res) <- envim_res$Compound_name

###################upset
library(UpSetR)
pdf("sup7_upset.pdf", width = 6, height = 6)
upset(fromList(list(O2PLS = gsub("-", "_", o2pls_res$wellPredicted), 
                    MMINP = gsub("-", "_", predres2$wellPredicted),
                    MelonnPan = gsub("-", "_", melonn_res$wellPredicted),
                    ENVIM = rownames(envim_res)[envim_res$Train_spearman_cor > 0.3 & envim_res$Test_spearman_cor > 0.3])),
      keep.order = TRUE, text.scale = 1.5)
dev.off()

####################
save.image(file = "D1Result_of_4methods.RData")


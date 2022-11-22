library(MMINP)
library(melonnpan)
source("melonnpan_predict_modified.R")
source("ENVIM-main/ENVIM.R")
source("ENVIM-main/ENVIM_predict.R")

load("data/D1'_D2_D3_data.RData")
set.seed(1234)
p_ibd <- sample(rownames(f_ibdb), nrow(f_ibdb)/3)
t_ibd <- rownames(f_ibdb)[!rownames(f_ibdb) %in% p_ibd]
set.seed(1234)
p_bio <- sample(rownames(f_biomb), nrow(f_biomb)/3)
t_bio <- rownames(f_biomb)[!rownames(f_biomb) %in% p_bio]
set.seed(1234)
p_cc <- sample(rownames(f_ccmb), nrow(f_ccmb)/3)
t_cc <- rownames(f_ccmb)[!rownames(f_ccmb) %in% p_cc]

###################################################################################### MMINP
# D1': ibd
ibdbt <- f_ibdb[t_ibd, ]
ibdgt <- f_ibdg[t_ibd, ]
ibdbtscale <- MMINP.preprocess(ibdbt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
ibdgtscale <- MMINP.preprocess(data = ibdgt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
ibdgp <- f_ibdg[p_ibd, ]
ibdgpscale <- MMINP.preprocess(ibdgp, transformed = "boxcox")
ibdbm <- f_ibdb[p_ibd, ]
ibdbm_pro <- MMINP.preprocess(ibdbm, normalized = T, transformed = "none", scaled = F)
ibdbmscale <- MMINP.preprocess(ibdbm_pro, normalized = F, transformed = "boxcox", scaled = T)
mibd <- MMINP.train(metag = ibdgtscale, metab = ibdbtscale, 
                    n = 3:10, nx = 0:10, ny = 0:10, 
                    nr_folds = 7, nr_cores = 10, seed = 1234)
mibd 
ibdp <- MMINP.predict(mibd, newdata = ibdgpscale)
resibd <- compareFeatures(predicted = ibdp, measured = ibdbmscale)
length(resibd$wellPredicted) 

# D3: bio
biobt <- f_biomb[t_bio, ]
biogt <- f_biomg[t_bio, ]
biobtscale <- MMINP.preprocess(biobt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
biogtscale <- MMINP.preprocess(data = biogt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
biogp <- f_biomg[p_bio, ]
biogpscale <- MMINP.preprocess(biogp, transformed = "boxcox")
biobm <- f_biomb[p_bio, ]
biobm_pro <- MMINP.preprocess(biobm, normalized = T, transformed = 'none', scaled = F)
biobmscale <- MMINP.preprocess(biobm_pro, normalized = F, transformed = "boxcox", scaled = T)
mbio <- MMINP.train(metag = biogtscale, metab = biobtscale, 
                    n = 3:10, nx = 0:10, ny = 0:10, 
                    nr_folds = 5, nr_cores = 10, seed = 1234)
mbio 
biop <- MMINP.predict(mbio, newdata = biogpscale)
resbio <- compareFeatures(predicted = biop, measured = biobmscale)
length(resbio$wellPredicted) 


# D2: cc
ccbt <- f_ccmb[t_cc, ]
ccgt <- f_ccmg[t_cc, ]
ccbtscale <- MMINP.preprocess(ccbt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
ccgtscale <- MMINP.preprocess(data = ccgt, prev = 0.1, abund = 0.0001, transformed = "boxcox")
ccgp <- f_ccmg[p_cc, ]
ccgpscale <- MMINP.preprocess(ccgp, transformed = "boxcox")
ccbm <- f_ccmb[p_cc, ]
ccbm_pro <- MMINP.preprocess(ccbm, normalized = T, transformed = 'none', scaled = F)
ccbmscale <- MMINP.preprocess(ccbm_pro, normalized = F, transformed = "boxcox", scaled = T)
mcc <- MMINP.train(metag = ccgtscale, metab = ccbtscale, 
                   n = 3:10, nx = 0:10, ny = 0:10, 
                   nr_folds = 10, nr_cores = 10, seed = 1234)
mcc 
ccp <- MMINP.predict(mcc, newdata = ccgpscale)
rescc <- compareFeatures(predicted = ccp, measured = ccbmscale)
length(rescc$wellPredicted) 


#--------------------------------------predict other data set
#ibdmodel
mibd_ccp <- MMINP.predict(mibd, newdata = ccgpscale)
mibd_rescc <- compareFeatures(predicted = mibd_ccp, measured = ccbmscale)
length(mibd_rescc$wellPredicted) 

mibd_biop <- MMINP.predict(mibd, newdata = biogpscale)
mibd_resbio <- compareFeatures(predicted = mibd_biop, measured = biobmscale)
length(mibd_resbio$wellPredicted) 

#ccmodel
mcc_ibdp <- MMINP.predict(mcc, newdata = ibdgpscale, minGeneSize = 0.2)
mcc_resibd <- compareFeatures(predicted = mcc_ibdp, measured = ibdbmscale)
length(mcc_resibd$wellPredicted) 

mcc_biop <- MMINP.predict(mcc, newdata = biogpscale)
mcc_resbio <- compareFeatures(predicted = mcc_biop, measured = biobmscale)
length(mcc_resbio$wellPredicted) 

#biomodel
mbio_ibdp <- MMINP.predict(mbio, newdata = ibdgpscale, minGeneSize = 0.2)
mbio_resibd <- compareFeatures(predicted = mbio_ibdp, measured = ibdbmscale)
length(mbio_resibd$wellPredicted) 

mbio_ccp <- MMINP.predict(mbio, newdata = ccgpscale)
mbio_rescc <- compareFeatures(predicted = mbio_ccp, measured = ccbmscale)
length(mbio_rescc$wellPredicted) 

###################################################################################### O2PLS
mibd_o2pls <- o2m(X = ibdgtscale, Y = ibdbtscale, 
                  n = mibd$components$n, nx = mibd$components$nx, ny = mibd$components$ny)
mibd_o2pls 
ibdp_o2pls <- MMINP.predict(mibd_o2pls, newdata = ibdgpscale)
resibd_o2pls <- compareFeatures(predicted = ibdp_o2pls, measured = ibdbmscale)
length(resibd_o2pls$wellPredicted) 

mbio_o2pls <- o2m(X = biogtscale, Y = biobtscale, 
                    n = mbio$components$n, nx = mbio$components$nx, ny = mbio$components$ny)
mbio_o2pls 
biop_o2pls <- MMINP.predict(mbio_o2pls, newdata = biogpscale)
resbio_o2pls <- compareFeatures(predicted = biop_o2pls, measured = biobmscale)
length(resbio_o2pls$wellPredicted) 

mcc_o2pls <- o2m(X = ccgtscale, Y = ccbtscale, 
                  n = mcc$components$n, nx = mcc$components$nx, ny = mcc$components$ny)
mcc_o2pls 
ccp_o2pls <- MMINP.predict(mcc_o2pls, newdata = ccgpscale)
rescc_o2pls <- compareFeatures(predicted = ccp_o2pls, measured = ccbmscale)
length(rescc_o2pls$wellPredicted) 

#--------------------------------------predict other data set
#ibdmodel
mibd_o2pls_ccp <- MMINP.predict(mibd_o2pls, newdata = ccgpscale)
mibd_o2pls_rescc <- compareFeatures(predicted = mibd_o2pls_ccp, measured = ccbmscale)
length(mibd_o2pls_rescc$wellPredicted) 

mibd_o2pls_biop <- MMINP.predict(mibd_o2pls, newdata = biogpscale)
mibd_o2pls_resbio <- compareFeatures(predicted = mibd_o2pls_biop, measured = biobmscale)
length(mibd_o2pls_resbio$wellPredicted) 

#ccmodel
mcc_o2pls_ibdp <- MMINP.predict(mcc_o2pls, newdata = ibdgpscale, minGeneSize = 0.2)
mcc_o2pls_resibd <- compareFeatures(predicted = mcc_o2pls_ibdp, measured = ibdbmscale)
length(mcc_o2pls_resibd$wellPredicted) 

mcc_o2pls_biop <- MMINP.predict(mcc_o2pls, newdata = biogpscale)
mcc_o2pls_resbio <- compareFeatures(predicted = mcc_o2pls_biop, measured = biobmscale)
length(mcc_o2pls_resbio$wellPredicted) 

#biomodel
mbio_o2pls_ibdp <- MMINP.predict(mbio_o2pls, newdata = ibdgpscale, minGeneSize = 0.2)
mbio_o2pls_resibd <- compareFeatures(predicted = mbio_o2pls_ibdp, measured = ibdbmscale)
length(mbio_o2pls_resibd$wellPredicted) 

mbio_o2pls_ccp <- MMINP.predict(mbio_o2pls, newdata = ccgpscale)
mbio_o2pls_rescc <- compareFeatures(predicted = mbio_o2pls_ccp, measured = ccbmscale)
length(mbio_o2pls_rescc$wellPredicted) 

###################################################################################### melonnpan
dir.create("individual")
dir.create(dir_ibd <- "individual/ibd/")
ibdbt_pro <- MMINP.preprocess(ibdbt, normalized = T, prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
ibdgt_pro <- MMINP.preprocess(ibdgt, normalized = T, prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
ibdgp_pro <- MMINP.preprocess(ibdgp, normalized = T, #prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
melonnpan.train(metab = ibdbt_pro, metag = ibdgt_pro, outputDirectory = dir_ibd, cores = 10)
weight_ibd <- read.table(paste0(dir_ibd, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_ibd) 
melon_ibd <- melonnpan.predict.modified(metag = ibdgp_pro, weight.matrix = weight_ibd, train.metag = ibdgt_pro,
                                        corr.method = "spearman", output = paste0(dir_ibd, "ibd"))
melon_ibdp <- melon_ibd$pred
rownames(melon_ibdp) <- melon_ibdp$ID
melon_ibdp$ID <- NULL
melonres_ibdp <- compareFeatures(predicted = melon_ibdp, measured = ibdbm_pro)
length(melonres_ibdp$wellPredicted) 

dir.create(dir_bio <- "individual/bio/")
biobt_pro <- MMINP.preprocess(biobt, normalized = T, prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
biogt_pro <- MMINP.preprocess(biogt, normalized = T, prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
biogp_pro <- MMINP.preprocess(biogp, normalized = T, #prev = 0.1, abund = 0.0001, 
                              transformed = 'none', scaled = F)
melonnpan.train(metab = biobt_pro, metag = biogt_pro, outputDirectory = dir_bio, cores = 10)
weight_bio <- read.table(paste0(dir_bio, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_bio) 
melon_bio <- melonnpan.predict.modified(metag = biogp_pro, weight.matrix = weight_bio, train.metag = biogt_pro,
                                        corr.method = "spearman", output = paste0(dir_bio, "bio"))
melon_biop <- melon_bio$pred
rownames(melon_biop) <- melon_biop$ID
melon_biop$ID <- NULL
melonres_biop <- compareFeatures(predicted = melon_biop, measured = biobm_pro)
length(melonres_biop$wellPredicted) 


dir.create(dir_cc <- "individual/cc/")
ccbt_pro <- MMINP.preprocess(ccbt, normalized = T, prev = 0.1, abund = 0.0001, 
                             transformed = 'none', scaled = F)
ccgt_pro <- MMINP.preprocess(ccgt, normalized = T, prev = 0.1, abund = 0.0001, 
                             transformed = 'none', scaled = F)
ccgp_pro <- MMINP.preprocess(ccgp, normalized = T, #prev = 0.1, abund = 0.0001, 
                             transformed = 'none', scaled = F)
melonnpan.train(metab = ccbt_pro, metag = ccgt_pro, outputDirectory = dir_cc, cores = 10)
weight_cc <- read.table(paste0(dir_cc, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_cc) 
melon_cc <- melonnpan.predict.modified(metag = ccgp_pro, weight.matrix = weight_cc, train.metag = ccgt_pro,
                                       corr.method = "spearman", output = paste0(dir_cc, "cc"))
melon_ccp <- melon_cc$pred
rownames(melon_ccp) <- melon_ccp$ID
melon_ccp$ID <- NULL
melonres_ccp <- compareFeatures(predicted = melon_ccp, measured = ccbm_pro)
length(melonres_ccp$wellPredicted) 

#--------------------------------------predict other data set
#ibdmodel
melonibd_cc <- melonnpan.predict.modified(metag = ccgp_pro, weight.matrix = weight_ibd, train.metag = ibdgt_pro,
                                          corr.method = "spearman", output = paste0(dir_ibd, "cc"))
melonibd_ccp <- melonibd_cc$pred
rownames(melonibd_ccp) <- melonibd_ccp$ID
melonibd_ccp$ID <- NULL
melonibdres_ccp <- compareFeatures(predicted = melonibd_ccp, measured = ccbm_pro)
length(melonibdres_ccp$wellPredicted) 

#biogp_pro <- MMINP.preprocess(biogp_pro, normalized = F, prev = 3/nrow(biogp_pro), transformed = 'none', scaled = F)
melonibd_bio <- melonnpan.predict.modified(metag = biogp_pro, weight.matrix = weight_ibd, train.metag = ibdgt_pro,
                                           corr.method = "spearman", output = paste0(dir_ibd, "bio"))
melonibd_biop <- melonibd_bio$pred
rownames(melonibd_biop) <- melonibd_biop$ID
melonibd_biop$ID <- NULL
melonibdres_biop <- compareFeatures(predicted = melonibd_biop, measured = biobm_pro)
length(melonibdres_biop$wellPredicted) 

#ccmodel
meloncc_ibd <- melonnpan.predict.modified(metag = ibdgp_pro, weight.matrix = weight_cc, train.metag = ccgt_pro,
                                          corr.method = "spearman", output = paste0(dir_cc, "ibd"))
meloncc_ibdp <- meloncc_ibd$pred
rownames(meloncc_ibdp) <- meloncc_ibdp$ID
meloncc_ibdp$ID <- NULL
melonccres_ibdp <- compareFeatures(predicted = meloncc_ibdp, measured = ibdbm_pro)
length(melonccres_ibdp$wellPredicted) 

meloncc_bio <- melonnpan.predict.modified(metag = biogp_pro, weight.matrix = weight_cc, train.metag = ccgt_pro,
                                          corr.method = "spearman", output = paste0(dir_cc, "bio"))
meloncc_biop <- meloncc_bio$pred
rownames(meloncc_biop) <- meloncc_biop$ID
meloncc_biop$ID <- NULL
melonccres_biop <- compareFeatures(predicted = meloncc_biop, measured = biobm_pro)
length(melonccres_biop$wellPredicted) 

#biomodel
melonbio_ibd <- melonnpan.predict.modified(metag = ibdgp_pro, weight.matrix = weight_bio, train.metag = biogt_pro,
                                           corr.method = "spearman", output = paste0(dir_bio, "ibd"))
melonbio_ibdp <- melonbio_ibd$pred
rownames(melonbio_ibdp) <- melonbio_ibdp$ID
melonbio_ibdp$ID <- NULL
melonbiores_ibdp <- compareFeatures(predicted = melonbio_ibdp, measured = ibdbm_pro)
length(melonbiores_ibdp$wellPredicted) 

melonbio_cc <- melonnpan.predict.modified(metag = ccgp_pro, weight.matrix = weight_bio, train.metag = biogt_pro,
                                          corr.method = "spearman", output = paste0(dir_bio, "cc"))
melonbio_ccp <- melonbio_cc$pred
rownames(melonbio_ccp) <- melonbio_ccp$ID
melonbio_ccp$ID <- NULL
melonbiores_ccp <- compareFeatures(predicted = melonbio_ccp, measured = ccbm_pro)
length(melonbiores_ccp$wellPredicted) 

###################################################################################### ENVIM
smoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  return(x)
}
ibdb.train <- apply(ibdbt_pro, 2, smoothZero)
ibdb.test <- ibdbm_pro[, colnames(ibdbt_pro)]
gid_ibd <- intersect(colnames(ibdgp_pro), colnames(ibdgt_pro))
bid_ibd <- intersect(colnames(ibdb.train), colnames(ibdb.test))
dir.create(output <- paste0("individual/ENVIM_ibd/"))
ENVIM(microbio.train = ibdgt_pro[, gid_ibd],
      microbio.test = ibdgp_pro[, gid_ibd],
      metab.train = ibdb.train[, bid_ibd],
      metab.test = ibdb.test[, bid_ibd],
      seed = 1234,
      outputdirectory = output,
      fold_rf = 10,
      fold_ENVIM = 10)
setwd("../../")

biob.train <- apply(biobt_pro, 2, smoothZero)
biob.test <- biobm_pro[, colnames(biobt_pro)]
gid_bio <- intersect(colnames(biogp_pro), colnames(biogt_pro))
bid_bio <- intersect(colnames(biob.train), colnames(biob.test))
dir.create(output <- paste0("individual/ENVIM_bio/"))
ENVIM(microbio.train = biogt_pro[, gid_bio],
      microbio.test = biogp_pro[, gid_bio],
      metab.train = biob.train[, bid_bio],
      metab.test = biob.test[, bid_bio],
      seed = 1234,
      outputdirectory = output,
      fold_rf = 10,
      fold_ENVIM = 10)
setwd("../../")

ccb.train <- apply(ccbt_pro, 2, smoothZero)
ccb.test <- ccbm_pro[, colnames(ccbt_pro)]
gid_cc <- intersect(colnames(ccgp_pro), colnames(ccgt_pro))
bid_cc <- intersect(colnames(ccb.train), colnames(ccb.test))
dir.create(output <- paste0("individual/ENVIM_cc/"))
ENVIM(microbio.train = ccgt_pro[, gid_cc],
      microbio.test = ccgp_pro[, gid_cc],
      metab.train = ccb.train[, bid_cc],
      metab.test = ccb.test[, bid_cc],
      seed = 1234,
      outputdirectory = output,
      fold_rf = 10,
      fold_ENVIM = 10)
setwd("../../")

envim_ibd <- read.csv("individual/ENVIM_ibd/summary.csv", header = T)
sum(envim_ibd$Train_spearman_cor > 0.3 & envim_ibd$Test_spearman_cor > 0.3)
envim_bio <- read.csv("individual/ENVIM_bio/summary.csv", header = T)
sum(envim_bio$Train_spearman_cor > 0.3 & envim_bio$Test_spearman_cor > 0.3)
envim_cc <- read.csv("individual/ENVIM_cc/summary.csv", header = T)
sum(envim_cc$Train_spearman_cor > 0.3 & envim_cc$Test_spearman_cor > 0.3, na.rm = T)

#--------------------------------------predict other data set
dir.create(output <- paste0("individual/ENVIM_ibd_bio2/"))
ENVIM_predict(microbio.train = ibdgt_pro,
              microbio.test = biogp_pro[, apply(biogp_pro, 2, function(x) length(unique(x))>2)],
              metab.train = ibdb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 10,
              fold_ENVIM = 10)
setwd("../../")
envim_ibd_bio <- read.csv("individual/ENVIM_ibd_bio2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimibdbio <- compareFeatures(envim_ibd_bio, biobm_pro)
length(res_envimibdbio$wellPredicted)

dir.create(output <- paste0("individual/ENVIM_ibd_cc2/"))
ENVIM_predict(microbio.train = ibdgt_pro,
              microbio.test = ccgp_pro[, apply(ccgp_pro, 2, function(x) length(unique(x))>2)],
              metab.train = ibdb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 10,
              fold_ENVIM = 10)
setwd("../../")
envim_ibd_cc <- read.csv("individual/ENVIM_ibd_cc2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimibdcc <- compareFeatures(envim_ibd_cc, ccbm_pro)
length(res_envimibdcc$wellPredicted)

dir.create(output <- paste0("individual/ENVIM_bio_ibd2/"))
ENVIM_predict(microbio.train = biogt_pro,
      microbio.test = ibdgp_pro[, apply(ibdgp_pro, 2, function(x) length(unique(x))>2)],
      metab.train = biob.train,
      #metab.test = ibdbm_pro[, bid],
      seed = 1234,
      outputdirectory = output,
      fold_rf = 10,
      fold_ENVIM = 10)
setwd("../../")
envim_bio_ibd <- read.csv("individual/ENVIM_bio_ibd2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimbioibd <- compareFeatures(envim_bio_ibd, ibdbm_pro)
length(res_envimbioibd$wellPredicted)

dir.create(output <- paste0("individual/ENVIM_bio_cc2/"))
ENVIM_predict(microbio.train = biogt_pro,
              microbio.test = ccgp_pro[, apply(ccgp_pro, 2, function(x) length(unique(x))>2)],
              metab.train = biob.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 10,
              fold_ENVIM = 10)
setwd("../../")
envim_bio_cc <- read.csv("individual/ENVIM_bio_cc2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimbiocc <- compareFeatures(envim_bio_cc, ccbm_pro)
length(res_envimbiocc$wellPredicted)

dir.create(output <- paste0("individual/ENVIM_cc_ibd2/"))
ENVIM_predict(microbio.train = ccgt_pro,
              microbio.test = ibdgp_pro[, apply(ibdgp_pro, 2, function(x) length(unique(x))>2)],
              metab.train = ccb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 10,
              fold_ENVIM = 10)
setwd("../../")
envim_cc_ibd <- read.csv("individual/ENVIM_cc_ibd2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimccibd <- compareFeatures(envim_cc_ibd, ibdbm_pro)
length(res_envimccibd$wellPredicted)

dir.create(output <- paste0("individual/ENVIM_cc_bio2/"))
ENVIM_predict(microbio.train = ccgt_pro,
              microbio.test = biogp_pro[, apply(biogp_pro, 2, function(x) length(unique(x))>2)],
              metab.train = ccb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 10,
              fold_ENVIM = 10)
setwd("../../")
envim_cc_bio <- read.csv("individual/ENVIM_cc_bio2/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimccbio <- compareFeatures(envim_cc_bio, biobm_pro)
length(res_envimccbio$wellPredicted)


###################################################################################### overlap
#O2PLS
intersect(rownames(mibd_o2pls$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mibd_o2pls$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mibd_o2pls$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mibd_o2pls$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mibd_o2pls$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mibd_o2pls$C.), colnames(ccbmscale)) %>% length()

intersect(rownames(mbio_o2pls$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mbio_o2pls$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mbio_o2pls$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mbio_o2pls$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mbio_o2pls$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mbio_o2pls$C.), colnames(ccbmscale)) %>% length()

intersect(rownames(mcc_o2pls$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mcc_o2pls$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mcc_o2pls$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mcc_o2pls$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mcc_o2pls$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mcc_o2pls$C.), colnames(ccbmscale)) %>% length()

#MMINP
intersect(rownames(mibd$model$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mibd$model$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mibd$model$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mibd$model$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mibd$model$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mibd$model$C.), colnames(ccbmscale)) %>% length()

intersect(rownames(mbio$model$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mbio$model$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mbio$model$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mbio$model$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mbio$model$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mbio$model$C.), colnames(ccbmscale)) %>% length()

intersect(rownames(mcc$model$W.), colnames(ibdgp_pro)) %>% length()
intersect(rownames(mcc$model$W.), colnames(biogp_pro)) %>% length()
intersect(rownames(mcc$model$W.), colnames(ccgp_pro)) %>% length()
intersect(rownames(mcc$model$C.), colnames(ibdbmscale)) %>% length()
intersect(rownames(mcc$model$C.), colnames(biobmscale)) %>% length()
intersect(rownames(mcc$model$C.), colnames(ccbmscale)) %>% length()

#MelonnPan
intersect(rownames(weight_ibd), colnames(ibdgp_pro)) %>% length()
intersect(rownames(weight_ibd), colnames(biogp_pro)) %>% length()
intersect(rownames(weight_ibd), colnames(ccgp_pro)) %>% length()
intersect(colnames(weight_ibd), colnames(ibdbm_pro)) %>% length()
intersect(colnames(weight_ibd), colnames(biobm_pro)) %>% length()
intersect(colnames(weight_ibd), colnames(ccbm_pro)) %>% length()

intersect(rownames(weight_bio), colnames(ibdgp_pro)) %>% length()
intersect(rownames(weight_bio), colnames(biogp_pro)) %>% length()
intersect(rownames(weight_bio), colnames(ccgp_pro)) %>% length()
intersect(colnames(weight_bio), colnames(ibdbm_pro)) %>% length()
intersect(colnames(weight_bio), colnames(biobm_pro)) %>% length()
intersect(colnames(weight_bio), colnames(ccbm_pro)) %>% length()

intersect(rownames(weight_cc), colnames(ibdgp_pro)) %>% length()
intersect(rownames(weight_cc), colnames(biogp_pro)) %>% length()
intersect(rownames(weight_cc), colnames(ccgp_pro)) %>% length()
intersect(colnames(weight_cc), colnames(ibdbm_pro)) %>% length()
intersect(colnames(weight_cc), colnames(biobm_pro)) %>% length()
intersect(colnames(weight_cc), colnames(ccbm_pro)) %>% length()

#ENVIM
intersect(colnames(envim_ibd_bio), colnames(ibdbm_pro)) %>% length()
intersect(colnames(envim_ibd_bio), colnames(biobm_pro)) %>% length()
intersect(colnames(envim_ibd_bio), colnames(ccbm_pro)) %>% length()

intersect(colnames(envim_bio_cc), colnames(ibdbm_pro)) %>% length()
intersect(colnames(envim_bio_cc), colnames(biobm_pro)) %>% length()
intersect(colnames(envim_bio_cc), colnames(ccbm_pro)) %>% length()

intersect(colnames(envim_cc_bio), colnames(ibdbm_pro)) %>% length()
intersect(colnames(envim_cc_bio), colnames(biobm_pro)) %>% length()
intersect(colnames(envim_cc_bio), colnames(ccbm_pro)) %>% length()

######################################################################################
save.image(file = "individual/result_of_4methods.RData")

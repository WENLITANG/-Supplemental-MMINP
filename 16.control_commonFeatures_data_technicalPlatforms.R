library(MMINP)
library(melonnpan)
source("melonnpan_predict_modified.R")
source("ENVIM-main/ENVIM.R")
source("ENVIM-main/ENVIM_predict.R")
load("data/D1'_D2_D3_data.RData")

meta <- read.table("data/merged_hc.meta.txt", header = T, row.names = 1, sep = "\t", check.names = F)
# D1': ibd
# D2: cc
# D3: bio
hc_ibdb <- MMINP.preprocess(f_ibdb[rownames(f_ibdb) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 
hc_biomb <- MMINP.preprocess(f_biomb[rownames(f_biomb) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 
hc_ccmb <- MMINP.preprocess(f_ccmb[rownames(f_ccmb) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 
hc_ibdg <- MMINP.preprocess(f_ibdg[rownames(f_ibdg) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 
hc_biomg <- MMINP.preprocess(f_biomg[rownames(f_biomg) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 
hc_ccmg <- MMINP.preprocess(f_ccmg[rownames(f_ccmg) %in% rownames(meta), ], normalized = F, transformed = "none", scaled = F) 

set.seed(1234)
id_bio <- sample(rownames(hc_biomb), 56)
id_biop <- rownames(hc_biomb)[!rownames(hc_biomb) %in% id_bio]

set.seed(1234)
id_cc <- sample(rownames(hc_ccmb), 56)
id_ccp <- rownames(hc_ccmb)[!rownames(hc_ccmb) %in% id_cc]

gid <- intersect(colnames(hc_ibdg)[apply(hc_ibdg, 2, function(x) sum(x>0)>2)], 
                 colnames(hc_biomg)[apply(hc_biomg[id_bio, ], 2, function(x) sum(x>0)>2)]) %>% 
  intersect(colnames(hc_ccmg)[apply(hc_ccmg[id_cc, ], 2, function(x) sum(x>0)>2)])
bid <- intersect(colnames(hc_ibdb)[apply(hc_ibdb, 2, function(x) sum(x>0)>2)], 
                 colnames(hc_biomb)[apply(hc_biomb[id_bio, ], 2, function(x) sum(x>0)>2)]) %>% 
  intersect(colnames(hc_ccmb)[apply(hc_ccmb[id_cc, ], 2, function(x) sum(x>0)>2)])
hc_ibdb <- hc_ibdb[, bid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)
hc_ibdg <- hc_ibdg[, gid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)
hc_biomb <- hc_biomb[, bid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)
hc_biomg <- hc_biomg[, gid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)
hc_ccmb <- hc_ccmb[, bid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)
hc_ccmg <- hc_ccmg[, gid] %>% MMINP.preprocess(normalized = T, transformed = "none", scaled = F)


########################################################## MMINP
#ibd
ibdbtscale <- MMINP.preprocess(hc_ibdb, normalized = F, transformed = "boxcox")
ibdgtscale <- MMINP.preprocess(data = hc_ibdg, normalized = F, transformed = 'boxcox', scaled = T)
mibd <- MMINP.train(metag = ibdgtscale, metab = ibdbtscale,
                    n = 3:10, nx = 0:10, ny = 0:10,
                    nr_folds = 3, nr_cores = 10, seed = 1234)
mibd 

#bio
biobtscale <- MMINP.preprocess(hc_biomb[id_bio, ], normalized = F, transformed = 'boxcox')
biogtscale <- MMINP.preprocess(data = hc_biomg[id_bio, ], normalized = F, transformed = 'boxcox', scaled = T)
biogpscale <- MMINP.preprocess(data = hc_biomg[id_biop, ], normalized = F, transformed = 'boxcox', scaled = T)
biobpscale <- MMINP.preprocess(data = hc_biomb[id_biop, ], normalized = F, transformed = 'boxcox', scaled = T)
mbio <- MMINP.train(metag = biogtscale, metab = biobtscale,
                    n = 3:10, nx = 0:10, ny = 0:10,
                    nr_folds = 3, nr_cores = 10, seed = 1234)
mbio 
mbio_biop <- MMINP.predict(mbio, newdata = biogpscale, minGeneSize = 0.2)
mbio_resbio <- compareFeatures(predicted = mbio_biop, measured = biobpscale)
length(mbio_resbio$wellPredicted) 

#cc
ccbtscale <- MMINP.preprocess(hc_ccmb[id_cc, ], normalized = F, transformed = 'boxcox')
ccgtscale <- MMINP.preprocess(data = hc_ccmg[id_cc, ], normalized = F, transformed = 'boxcox', scaled = T)
ccgpscale <- MMINP.preprocess(data = hc_ccmg[id_ccp, ], normalized = F, transformed = 'boxcox', scaled = T)
ccbpscale <- MMINP.preprocess(data = hc_ccmb[id_ccp, ], normalized = F, transformed = 'boxcox', scaled = T)
mcc <- MMINP.train(metag = ccgtscale, metab = ccbtscale,
                   n = 3:10, nx = 0:10, ny = 0:10, 
                   nr_folds = 3, nr_cores = 10, seed = 1234)
mcc 
mcc_ccp <- MMINP.predict(mcc, newdata = ccgpscale, minGeneSize = 0.2)
mcc_rescc <- compareFeatures(predicted = mcc_ccp, measured = ccbpscale)
length(mcc_rescc$wellPredicted) 
#-----------------------------------predict other data set
#ibdmodel
mibd_ccp <- MMINP.predict(mibd, newdata = ccgpscale, minGeneSize = 0.3)
mibd_rescc <- compareFeatures(predicted = mibd_ccp, measured = ccbpscale)
length(mibd_rescc$wellPredicted) 

mibd_biop <- MMINP.predict(mibd, newdata = biogpscale)
mibd_resbio <- compareFeatures(predicted = mibd_biop, measured = biobpscale)
length(mibd_resbio$wellPredicted) 

#ccmodel
mcc_ibdp <- MMINP.predict(mcc, newdata = ibdgtscale, minGeneSize = 0.2)
mcc_resibd <- compareFeatures(predicted = mcc_ibdp, measured = ibdbtscale)
length(mcc_resibd$wellPredicted) 

mcc_biop <- MMINP.predict(mcc, newdata = biogpscale, minGeneSize = 0.3)
mcc_resbio <- compareFeatures(predicted = mcc_biop, measured = biobpscale)
length(mcc_resbio$wellPredicted) 

#biomodel
mbio_ibdp <- MMINP.predict(mbio, newdata = ibdgtscale, minGeneSize = 0.2)
mbio_resibd <- compareFeatures(predicted = mbio_ibdp, measured = ibdbtscale)
length(mbio_resibd$wellPredicted) 

mbio_ccp <- MMINP.predict(mbio, newdata = ccgpscale)
mbio_rescc <- compareFeatures(predicted = mbio_ccp, measured = ccbpscale)
length(mbio_rescc$wellPredicted) 

########################################################## O2PLS
mibd_o2pls <- o2m(X = ibdgtscale, Y = ibdbtscale, 
                  n = mibd$components$n, nx = mibd$components$nx, ny = mibd$components$ny)
mibd_o2pls 

mbio_o2pls <- o2m(X = biogtscale, Y = biobtscale, 
                  n = mbio$components$n, nx = mbio$components$nx, ny = mbio$components$ny)
mbio_o2pls 
biop_o2pls <- MMINP.predict(mbio_o2pls, newdata = biogpscale)
resbio_o2pls <- compareFeatures(predicted = biop_o2pls, measured = biobpscale)
length(resbio_o2pls$wellPredicted) 

mcc_o2pls <- o2m(X = ccgtscale, Y = ccbtscale, 
                 n = mcc$components$n, nx = mcc$components$nx, ny = mcc$components$ny)
mcc_o2pls 
ccp_o2pls <- MMINP.predict(mcc_o2pls, newdata = ccgpscale)
rescc_o2pls <- compareFeatures(predicted = ccp_o2pls, measured = ccbpscale)
length(rescc_o2pls$wellPredicted) 

#-----------------------------------predict other data set
#ibdmodel
mibd_o2pls_ccp <- MMINP.predict(mibd_o2pls, newdata = ccgpscale)
mibd_o2pls_rescc <- compareFeatures(predicted = mibd_o2pls_ccp, measured = ccbpscale)
length(mibd_o2pls_rescc$wellPredicted) 

mibd_o2pls_biop <- MMINP.predict(mibd_o2pls, newdata = biogpscale)
mibd_o2pls_resbio <- compareFeatures(predicted = mibd_o2pls_biop, measured = biobpscale)
length(mibd_o2pls_resbio$wellPredicted) 

#ccmodel
mcc_o2pls_ibdp <- MMINP.predict(mcc_o2pls, newdata = ibdgtscale, minGeneSize = 0.2)
mcc_o2pls_resibd <- compareFeatures(predicted = mcc_o2pls_ibdp, measured = ibdbtscale)
length(mcc_o2pls_resibd$wellPredicted) 

mcc_o2pls_biop <- MMINP.predict(mcc_o2pls, newdata = biogpscale)
mcc_o2pls_resbio <- compareFeatures(predicted = mcc_o2pls_biop, measured = biobpscale)
length(mcc_o2pls_resbio$wellPredicted) 

#biomodel
mbio_o2pls_ibdp <- MMINP.predict(mbio_o2pls, newdata = ibdgtscale, minGeneSize = 0.2)
mbio_o2pls_resibd <- compareFeatures(predicted = mbio_o2pls_ibdp, measured = ibdbtscale)
length(mbio_o2pls_resibd$wellPredicted) 

mbio_o2pls_ccp <- MMINP.predict(mbio_o2pls, newdata = ccgpscale)
mbio_o2pls_rescc <- compareFeatures(predicted = mbio_o2pls_ccp, measured = ccbpscale)
length(mbio_o2pls_rescc$wellPredicted) 


########################################################## melonnpan
dir.create("control_commonFeatures/")
dir.create(dir_ibd <- "control_commonFeatures/ibd/")
melonnpan.train(metab = hc_ibdb, metag = hc_ibdg, outputDirectory = dir_ibd, cores = 10, nfolds = 3)
weight_ibd <- read.table(paste0(dir_ibd, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_ibd) 

dir.create(dir_bio <- "control_commonFeatures/bio/")
melonnpan.train(metab = hc_biomb[id_bio, ], metag = hc_biomg[id_bio, ], outputDirectory = dir_bio, cores = 10, nfolds = 3)
weight_bio <- read.table(paste0(dir_bio, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_bio) 
melonbio_bio <- melonnpan.predict.modified(metag = hc_biomg[id_biop, ], weight.matrix = weight_bio, train.metag = hc_biomg[id_bio, ],
                                           corr.method = "spearman", output = paste0(dir_bio, "bio"))
melonbio_biop <- melonbio_bio$pred
rownames(melonbio_biop) <- melonbio_biop$ID
melonbio_biop$ID <- NULL
melonbiores_biop <- compareFeatures(predicted = melonbio_biop, measured = hc_biomb[id_biop, ])
length(melonbiores_biop$wellPredicted) 


dir.create(dir_cc <- "control_commonFeatures/cc/")
melonnpan.train(metab = hc_ccmb[id_cc, ], metag = hc_ccmg[id_cc, ], outputDirectory = dir_cc, cores = 10, nfolds = 3)
weight_cc <- read.table(paste0(dir_cc, "MelonnPan_Trained_Weights.txt"), header = T, row.names = 1, sep = "\t", quote = "")
length(weight_cc) 
meloncc_cc <- melonnpan.predict.modified(metag = hc_ccmg[id_ccp, ], weight.matrix = weight_cc, train.metag = hc_ccmg[id_cc, ],
                                         corr.method = "spearman", output = paste0(dir_cc, "cc"))
meloncc_ccp <- meloncc_cc$pred
rownames(meloncc_ccp) <- meloncc_ccp$ID
meloncc_ccp$ID <- NULL
melonccres_ccp <- compareFeatures(predicted = meloncc_ccp, measured = hc_ccmb[id_ccp, ])
length(melonccres_ccp$wellPredicted) 

#-----------------------------------predict other data set
#ibdmodel
melonibd_cc <- melonnpan.predict.modified(metag = hc_ccmg[id_ccp, ], weight.matrix = weight_ibd, train.metag = hc_ibdg,
                                          corr.method = "spearman", output = paste0(dir_ibd, "cc"))
melonibd_ccp <- melonibd_cc$pred
rownames(melonibd_ccp) <- melonibd_ccp$ID
melonibd_ccp$ID <- NULL
melonibdres_ccp <- compareFeatures(predicted = melonibd_ccp, measured = hc_ccmb[id_ccp, ])
length(melonibdres_ccp$wellPredicted) 

melonibd_bio <- melonnpan.predict.modified(metag = hc_biomg[id_biop, ], weight.matrix = weight_ibd, train.metag = hc_ibdg,
                                           corr.method = "spearman", output = paste0(dir_ibd, "bio"))
melonibd_biop <- melonibd_bio$pred
rownames(melonibd_biop) <- melonibd_biop$ID
melonibd_biop$ID <- NULL
melonibdres_biop <- compareFeatures(predicted = melonibd_biop, measured = hc_biomb[id_biop, ])
length(melonibdres_biop$wellPredicted) 

#ccmodel
meloncc_ibd <- melonnpan.predict.modified(metag = hc_ibdg, weight.matrix = weight_cc, train.metag = hc_ccmg[id_cc, ],
                                          corr.method = "spearman", output = paste0(dir_cc, "ibd"))
meloncc_ibdp <- meloncc_ibd$pred
rownames(meloncc_ibdp) <- meloncc_ibdp$ID
meloncc_ibdp$ID <- NULL
melonccres_ibdp <- compareFeatures(predicted = meloncc_ibdp, measured = hc_ibdb)
length(melonccres_ibdp$wellPredicted) 

meloncc_bio <- melonnpan.predict.modified(metag = hc_biomg[id_biop, ], weight.matrix = weight_cc, train.metag = hc_ccmg[id_cc, ],
                                          corr.method = "spearman", output = paste0(dir_cc, "bio"))
meloncc_biop <- meloncc_bio$pred
rownames(meloncc_biop) <- meloncc_biop$ID
meloncc_biop$ID <- NULL
melonccres_biop <- compareFeatures(predicted = meloncc_biop, measured = hc_biomb[id_biop, ])
length(melonccres_biop$wellPredicted) 

#biomodel
melonbio_ibd <- melonnpan.predict.modified(metag = hc_ibdg, weight.matrix = weight_bio, train.metag = hc_biomg[id_bio, ],
                                           corr.method = "spearman", output = paste0(dir_bio, "ibd"))
melonbio_ibdp <- melonbio_ibd$pred
rownames(melonbio_ibdp) <- melonbio_ibdp$ID
melonbio_ibdp$ID <- NULL
melonbiores_ibdp <- compareFeatures(predicted = melonbio_ibdp, measured = hc_ibdb)
length(melonbiores_ibdp$wellPredicted) 

melonbio_cc <- melonnpan.predict.modified(metag = hc_ccmg[id_ccp, ], weight.matrix = weight_bio, train.metag = hc_biomg[id_bio, ],
                                          corr.method = "spearman", output = paste0(dir_bio, "cc"))
melonbio_ccp <- melonbio_cc$pred
rownames(melonbio_ccp) <- melonbio_ccp$ID
melonbio_ccp$ID <- NULL
melonbiores_ccp <- compareFeatures(predicted = melonbio_ccp, measured = hc_ccmb[id_ccp, ])
length(melonbiores_ccp$wellPredicted) 

########################################################## ENVIM
#-----------------------------------predict other data set
smoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  return(x)
}
ibdb.train <- apply(hc_ibdb, 2, smoothZero)
dir.create(output <- paste0("control_commonFeatures/ENVIM_ibd_bio/"))
ENVIM_predict(microbio.train = hc_ibdg,
              microbio.test = hc_biomg[id_biop, ],
              metab.train = ibdb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_ibd_bio <- read.csv("control_commonFeatures/ENVIM_ibd_bio/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimibdbio <- compareFeatures(envim_ibd_bio, hc_biomb[id_biop, ])
length(res_envimibdbio$wellPredicted)

dir.create(output <- paste0("control_commonFeatures/ENVIM_ibd_cc/"))
ENVIM_predict(microbio.train = hc_ibdg,
              microbio.test = hc_ccmg[id_ccp, ],
              metab.train = ibdb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_ibd_cc <- read.csv("control_commonFeatures/ENVIM_ibd_cc/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimibdcc <- compareFeatures(envim_ibd_cc, hc_ccmb[id_ccp, ])
length(res_envimibdcc$wellPredicted)


biob.train <- apply(hc_biomb[id_bio, ], 2, smoothZero)
dir.create(output <- paste0("control_commonFeatures/ENVIM_bio/"))
ENVIM_predict(microbio.train = hc_biomg[id_bio, ],
              microbio.test = hc_biomg[id_biop, ],
              metab.train = biob.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_bio_bio <- read.csv("control_commonFeatures/ENVIM_bio/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimbiobio <- compareFeatures(envim_bio_bio, hc_biomb[id_biop, ])
length(res_envimbiobio$wellPredicted)

dir.create(output <- paste0("control_commonFeatures/ENVIM_bio_cc/"))
ENVIM_predict(microbio.train = hc_biomg[id_bio, ],
              microbio.test = hc_ccmg[id_ccp, ],
              metab.train = biob.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_bio_cc <- read.csv("control_commonFeatures/ENVIM_bio_cc/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimbiocc <- compareFeatures(envim_bio_cc, hc_ccmb[id_ccp, ])
length(res_envimbiocc$wellPredicted)

dir.create(output <- paste0("control_commonFeatures/ENVIM_bio_ibd/"))
ENVIM_predict(microbio.train = hc_biomg[id_bio, ],
              microbio.test = hc_ibdg,
              metab.train = biob.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_bio_ibd <- read.csv("control_commonFeatures/ENVIM_bio_ibd/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimbioibd <- compareFeatures(envim_bio_ibd, hc_ibdb)
length(res_envimbioibd$wellPredicted)

ccb.train <- apply(hc_ccmb[id_cc, ], 2, smoothZero)
dir.create(output <- paste0("control_commonFeatures/ENVIM_cc/"))
ENVIM_predict(microbio.train = hc_ccmg[id_cc, ],
              microbio.test = hc_ccmg[id_ccp, ],
              metab.train = ccb.train,
              #metab.test = hc_ccmb[id_ccp, ],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_cc <- read.csv("control_commonFeatures/ENVIM_cc/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimcc <- compareFeatures(envim_cc, hc_ccmb[id_ccp, ])
length(res_envimcc$wellPredicted)

dir.create(output <- paste0("control_commonFeatures/ENVIM_cc_ibd/"))
ENVIM_predict(microbio.train = hc_ccmg[id_cc, ],
              microbio.test = hc_ibdg,
              metab.train = ccb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_cc_ibd <- read.csv("control_commonFeatures/ENVIM_cc_ibd/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimccibd <- compareFeatures(envim_cc_ibd, hc_ibdb)
length(res_envimccibd$wellPredicted)

dir.create(output <- paste0("control_commonFeatures/ENVIM_cc_bio/"))
ENVIM_predict(microbio.train = hc_ccmg[id_cc, ],
              microbio.test = hc_biomg[id_biop, ],
              metab.train = ccb.train,
              #metab.test = ibdbm_pro[, bid],
              seed = 1234,
              outputdirectory = output,
              fold_rf = 3,
              fold_ENVIM = 10)
setwd("../../")
envim_cc_bio <- read.csv("control_commonFeatures/ENVIM_cc_bio/Pred_test.csv", header = T, row.names = 1, check.names = F)
res_envimccbio <- compareFeatures(envim_cc_bio, hc_biomb[id_biop, ])
length(res_envimccbio$wellPredicted)


##########################################################
save.image(file = "control_commonFeatures/all.RData")


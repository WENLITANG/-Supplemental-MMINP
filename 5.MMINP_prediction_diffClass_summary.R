###############################################################################
#add: annotation
load("D1Result_of_4methods.RData")
source("metabolite_stat.R")
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

ibdclass <- read.table("data/IBD_putativeclass_HMDB.txt", header = T, sep = "\t", quote = "")
tmp <- merge(mbpredict, predres$res, by.x = "Metabolomic.Feature", by.y = "compound", all.x = T)
tmp <- merge(tmp, ibdclass, by.x = "Putative.Chemical.Class", by.y = "putative", all.x = T)
ibdclass <- merge(tmp, melonn_res$res, by.x = "Metabolomic.Feature", by.y = "compound", all.x = T)
rownames(ibdclass) <- ibdclass$Metabolomic.Feature

#add: statistical result(relative abundance)
tmp <- metab2_pro[trainS, rownames(ibdclass)]
tmp$group <- m[rownames(tmp) ,"Diagnosis"]
tmp <- tmp[, c("group", rownames(ibdclass))]
trainS_stat <- metabolite_stat(tmp)
trainS_stat$mean <- apply(tmp[, -1], 2, mean)
trainS_stat$median <- apply(tmp[, -1], 2, median)
colnames(trainS_stat) <- paste0("trainS_", colnames(trainS_stat))

tmp <- metab2_pro[predS, rownames(ibdclass)]
tmp$group <- m[rownames(tmp) ,"Diagnosis"]
tmp <- tmp[, c("group", rownames(ibdclass))]
predS_stat <- metabolite_stat(tmp)
predS_stat$mean <- apply(tmp[, -1], 2, mean)
predS_stat$median <- apply(tmp[, -1], 2, median)
colnames(predS_stat) <- paste0("predS_", colnames(predS_stat))

measured_stat <- merge(trainS_stat, predS_stat, by = "row.names")
rownames(measured_stat) <- measured_stat$Row.names
measured_stat$Row.names <- NULL
ibdclass2 <- merge(ibdclass, measured_stat, by = "row.names")
ibdclass2$sub_class <- ifelse(ibdclass2$type != "class", ibdclass2$sub_class, "unknown")
rownames(ibdclass2) <- ibdclass2$Row.names
ibdclass2$Row.names <- NULL

#add: prediction
tmp <- model1$trainres$res
colnames(tmp) <- paste0("MMINPmodel_", colnames(tmp))
ibdclass2 <- merge(ibdclass2, tmp, by = "row.names", all.x = T)
rownames(ibdclass2) <- ibdclass2$Row.names
ibdclass2$Row.names <- NULL

#add: statistical result(preprocessed)
d <- data.frame(group = m[rownames(mbp), "Diagnosis"], mbp)
mbp_dsumm2 <- metabolite_stat(d)
colnames(mbp_dsumm2) <- paste("MMINP", colnames(mbp_dsumm2), sep = "_")

d <- data.frame(group = m[rownames(mbm),"Diagnosis"], mbmscale[, intersect(colnames(mbp), colnames(mbmscale))])
mbmscale_dsumm2 <- metabolite_stat(d)
colnames(mbmscale_dsumm2) <- paste("mbmscale", colnames(mbmscale_dsumm2), sep = "_")
tmp1 <- merge(mbp_dsumm2, mbmscale_dsumm2, by = "row.names", all = T)
rownames(tmp1) <- tmp1$Row.names
tmp1$Row.names <- NULL

ibdclass3 <- merge(ibdclass2, tmp1, by.x = "compoundid", by.y = "row.names", all = T)
write.table(ibdclass3, "metablites2794_summary.txt", row.names = F, sep = "\t", quote = F)

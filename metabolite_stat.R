metabolite_stat <- function(d){
  dstat <- apply(d[, -1], 2, function(x) unlist(tapply(x, d$group, function(y) 
    c(mean=mean(y), median = median(y)))))
  kw <- apply(d[, -1], 2, function(x) 
    tryCatch({kruskal.test(x, d$group)$p.value}, error = function(e) {NA}))
  d1 <- subset(d, group != "CD")
  con_uc <- apply(d1[, -1], 2, function(x) 
    tryCatch({wilcox.test(x~d1$group)$p.value}, error = function(e) {NA}))
  d2 <- subset(d, group != "UC")
  con_cd <- apply(d2[, -1], 2, function(x) 
    tryCatch({wilcox.test(x~d2$group)$p.value}, error = function(e) {NA}))
  d3 <- subset(d, group != "Control")
  uc_cd <- apply(d3[, -1], 2, function(x) 
    tryCatch({wilcox.test(x~d3$group)$p.value}, error = function(e) {NA}))
  dsumm <- data.frame(t(dstat), kw, con_uc, con_cd, uc_cd)
  dsumm$kw_fdr <- p.adjust(dsumm$kw, method = "fdr")
  dsumm$con_uc_fdr <- p.adjust(dsumm$con_uc, method = "fdr")
  dsumm$con_cd_fdr <- p.adjust(dsumm$con_cd, method = "fdr")
  dsumm$uc_cd_fdr <- p.adjust(dsumm$uc_cd, method = "fdr")
  dsumm$diff_cd_con <- dsumm$CD.median - dsumm$Control.median
  dsumm$diff_cd_uc <- dsumm$CD.median - dsumm$UC.median
  dsumm$diff_uc_con <- dsumm$UC.median - dsumm$Control.median
  dsumm$fc_cd_con <- (dsumm$CD.median - dsumm$Control.median) / abs(dsumm$Control.median)
  dsumm$fc_cd_uc <- (dsumm$CD.median - dsumm$UC.median) / abs(dsumm$UC.median)
  dsumm$fc_uc_con <- (dsumm$UC.median - dsumm$Control.median) / abs(dsumm$Control.median)
  dsumm <- dsumm[order(dsumm$kw), ]
  return(dsumm)
}
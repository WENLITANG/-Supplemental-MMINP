#-------------------------------
#!usr/bin/R
#--------------------------------------------
#doInstall<-"FALSE" # Change to FALSE if you don't want packages installed.
#toInstall<-c("optparse","stats")
#if(doInstall){install.packages(toInstall, repos = "http://cran.us.r-project.org")}
#toLib<-toInstall
#lapply(toLib, library, character.only = TRUE)
library(optparse)
# parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-n", "--number"), type="character", help="split number. [required].")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$number)) stop('Please supply a number.')
#if(is.null(opts$outfile)) opts$outfile = "addPathwayinfo.txt"
i <- as.numeric(opts$number)
print(i)
########################################################################################
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

############################### ENVIM
source("ENVIM-main/ENVIM.R")
library(do)
mbm <- metab2_pro[predS, ]
mgp <- metag2_pro[predS, ]
microbio.train = mgtrm
microbio.test = mgp[, colnames(mgtrm)]
colnames(microbio.train) <- Replace(colnames(microbio.train), from = c(":", " ", "-", "\\[", "\\]", "\\(", "\\)", ","), to = "_")
colnames(microbio.test) <- Replace(colnames(microbio.test), from = c(":", " ", "-", "\\[", "\\]", "\\(", "\\)", ","), to = "_")
colnames(microbio.test) <- paste0("EC", colnames(microbio.test))
colnames(microbio.train) <- paste0("EC", colnames(microbio.train))

smoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  return(x)
}
metab.train <- apply(mbtrm, 2, smoothZero)
metab.test <- mbm[, colnames(mbtrm)]
colnames(metab.train) <- Replace(colnames(metab.train), from = "-", to = "_")
colnames(metab.test) <- Replace(colnames(metab.test), from = "-", to = "_")

dir.create(output <- paste0("IBD_ENVIM/split", i))
b <- i*300+1
e <- b+299
e <- ifelse(e>2794, 2794, e)
#e <- b+1
ENVIM(microbio.train = microbio.train,
      microbio.test = microbio.test,
      metab.train = metab.train[, b:e],
      metab.test = metab.test[, b:e],
      seed = 1234,
      outputdirectory = output,
      fold_rf = 10,
      fold_ENVIM = 10)
setwd("../..")


################################################################# diff
library(forecast)
library(dplyr)
load("D1Result_of_4methods.RData")
source("metabolite_stat.R")
#pathway
e <- read.table("data/20220924_kegg/enzyme_pathway.txt", header = F, sep = "\t")
colnames(e) <- c("enzyme", "pathway")
e <- apply(e, 2, function(x) 
  lapply(x, function(d) 
    strsplit(d, ":") %>% unlist %>% .[2]) %>% unlist) %>% as.data.frame()
e <- e[grep("map", e$pathway), ]
interest <- paste0("map", c("00120", "00121", "00380", "00400", "00600", "04071",
                            "01062", "00061", "00062", "00071", "01040", "01212"))
e2 <- e[e$pathway %in% interest, ]
pathway <- matrix(c("map00120", "Primary bile acid biosynthesis",
                    "map00121", "Secondary bile acid biosynthesis",
                    "map00380", "Tryptophan metabolism",
                    "map00400", "Phenylalanine, tyrosine and tryptophan biosynthesis",
                    "map00600", "Sphingolipid metabolism",
                    "map04071", "Sphingolipid signaling pathway",
                    "map01062", "Biosynthesis of terpenoids and steroids",
                    "map00061", "Fatty acid biosynthesis",
                    "map00062", "Fatty acid elongation",
                    "map00071", "Fatty acid degradation",
                    "map01040", "Biosynthesis of unsaturated fatty acids",
                    "map01212", "Fatty acid metabolism"), ncol = 2, byrow = T) %>% as.data.frame()
colnames(pathway) <- c("ID", "annotation")
pathway$name <- paste(pathway$ID, pathway$annotation, sep = ":")
pathway$type <- "pathway"
#trans to KO
ec_ko <- read.table("data/20220924_kegg/enzyme_ko.txt", header = F, sep = "\t")
colnames(ec_ko) <- c("enzyme", "ko")
ec_ko <- apply(ec_ko, 2, function(x) 
  lapply(x, function(d) 
    strsplit(d, ":") %>% unlist %>% .[2]) %>% unlist) %>% as.data.frame()
#kegg, bacteria list
taxon <- read.table("data/20220924_kegg/bacteria_list.txt", header = F, sep = "\t")
colnames(taxon) <- c("taxid", "shortname", "taxonomy")
taxon$genus <- lapply(taxon$taxonomy, function(x) strsplit(x, " ") %>% unlist %>% .[1]) %>% unlist()
taxon$species <- lapply(taxon$taxonomy, function(x) {a <- strsplit(x, " ") %>% unlist
                                                     b <- paste(a[1], a[2])}) %>% unlist()
#kegg, species:ec
#bac120_taxonomy_r207.tsv: download from https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_taxonomy_r207.tsv
gtdb <- read.table("data/bac120_taxonomy_r207.tsv", header = F, sep = "\t")
gtdb_taxonomy <- lapply(gtdb$V2, function(x) strsplit(x, ";") %>% unlist) %>% 
  unlist() %>% matrix(ncol = 7, byrow = T) %>% as.data.frame() %>% dplyr::distinct()
colnames(gtdb_taxonomy) <- c("domain", "phylum", "class", "object", "family", "genus", "species")
gtdb_taxonomy <- do::Replace(gtdb_taxonomy, from = ".__", to = "")
gtdb_taxonomy[nrow(gtdb_taxonomy)+1, ] <- c("Bacteria", "Firmicutes", "Clostridia", "Eubacteriales", "Oscillospiraceae", "Subdoligranulum", "manual")
gtdb_taxonomy[nrow(gtdb_taxonomy)+1, ] <- c("Archaea", "Euryarchaeota", "Methanobacteria", "Methanobacteriales", "Methanobacteriaceae", "Methanobrevibacter", "manual")

# diff
dta <- data.frame(group = m[rownames(mgpscale), "Diagnosis"], mgpscale, check.names = F)
mgp_dsumm <- metabolite_stat(dta)
#write.table(mgp_dsumm, "diff_mgpscale.txt", col.names = NA, sep = "\t", quote = F)
rownames(mgp_dsumm)[mgp_dsumm$con_uc < 0.05] %>% na.omit() %>%
  lapply(., function(x) strsplit(x, ":") %>% unlist %>% .[1]) %>% unlist()
ec_uccon <- data.frame(name = na.omit(rownames(mgp_dsumm)[mgp_dsumm$con_uc < 0.05]))
ec_uccon$ec <- ec_uccon$name %>% na.omit() %>% 
  lapply(., function(x) strsplit(x, ":") %>% unlist %>% .[1]) %>% unlist()
intec <- e2[e2$enzyme %in% ec_uccon$ec, ]
pathway_ec <- data.frame(intec[, c("pathway", "enzyme")], Value = 1)
colnames(pathway_ec) <- c("Source", "Target", "Value")
#4.2.1.20, RIX, RIM, Roseburia intestinalis, indole

## diff species
species <- openxlsx::read.xlsx("data/IBD/41564_2018_306_MOESM6_ESM.xlsx", sheet = 1, startRow = 2, rowNames = T)
species2 <- as.data.frame(t(apply(species[9:nrow(species), ], 2, as.numeric))) 
colnames(species2) <- rownames(species)[9:nrow(species)]
species2_pro <- as.data.frame(t(apply(species2, 1, function(x) x/sum(x))))
all(rownames(mgpscale) %in% rownames(species2_pro))
sp <- species2_pro[rownames(mgpscale), ]
sptrm <- sp[, apply(sp, 2, function(x) length(which(x>0.0001))>nrow(sp)*0.1)] #0.0001 * 1000000
spdta <- data.frame(group = m[rownames(mgpscale), "Diagnosis"], sptrm, check.names = F)
sp_dsumm <- metabolite_stat(spdta)
s_uccon <- rownames(sp_dsumm)[sp_dsumm$con_uc < 0.05]

## diff enzyme, diff species, correlation
sp_mgp_uccon <- psych::corr.test(sp[, s_uccon], mgpscale[, ec_uccon$name], method = "spearman", adjust = "fdr")
dim(sp_mgp_uccon$p)
sp_mgp <- cbind(reshape2::melt(sp_mgp_uccon$r), reshape2::melt(sp_mgp_uccon$p))
colnames(sp_mgp) <- c("species", "ec", "spearman_r", "sp", "e", "spearman_p")
identical(sp_mgp$species, sp_mgp$sp) & identical(sp_mgp$ec, sp_mgp$e)
sp_mgp$ecnumb <- sp_mgp$ec %>% as.character() %>%
  lapply(., function(x) strsplit(x, ":") %>% unlist %>% .[1]) %>% unlist()
sp_mgp2 <- sp_mgp[(sp_mgp$species %in% s_uccon) & (sp_mgp$ecnumb %in% intec$enzyme), ]
sp_mgp2$p_adjust <- p.adjust(sp_mgp2$spearman_p, method = "fdr")
ec_species <- sp_mgp2[(sp_mgp2$species %in% s_uccon) & (sp_mgp2$ecnumb %in% intec$enzyme) & (sp_mgp2$p_adjust < 0.05), 
                      c("ecnumb", "species", "spearman_r")]
colnames(ec_species) <- c("Source", "Target", "Value")
ec_info <- dplyr::distinct(sp_mgp[sp_mgp$ecnumb %in% ec_species$Source, c("ecnumb", "ec")])
colnames(ec_info) <- c("ID", "annotation")
ec_info$name <- ec_info$annotation
ec_info$type <- "enzyme"
species_info <- data.frame(ID = unique(ec_species$Target), annotation = unique(ec_species$Target),
                           name = unique(ec_species$Target), type = "species")
species_info2 <- data.frame(ID = species_info$ID, taxonomy = gsub("_", " ", species_info$ID))
species_info2$genus <- lapply(species_info2$taxonomy, function(x) strsplit(x, " ") %>% unlist %>% .[1]) %>% unlist()
gsub("_", " ", species_info$ID) %in% taxon$taxonomy
species_info2$species <- lapply(species_info2$taxonomy, function(x) {
  a <- strsplit(x, " ") %>% unlist
  b <- paste(a[1], a[2])}) %>% unlist()
species_info2$shortname <- lapply(species_info2$taxonomy, function(x) 
  ifelse(x %in% taxon$taxonomy, taxon$shortname[taxon$taxonomy %in% x], NA)) %>% unlist()
species_info2$shortname2 <- lapply(species_info2$species, function(x) 
  ifelse(x %in% taxon$species, taxon$shortname[taxon$species %in% x], NA)) %>% unlist()
taxon2 <- taxon[taxon$genus %in% species_info2$genus, ]
# write.table(taxon2, "diff_uc_con_species_kegg.txt", row.names = F, sep = "\t", quote = F)
# write.table(species_info2, "diff_uc_con_species.txt", row.names = F, sep = "\t", quote = F)

## download specis:ec (kegg)  
#https://rest.kegg.jp/link/fpla/enzyme
#https://rest.kegg.jp/list/genome
dir.create("kegg_uc")
for (i in na.omit(species_info2$shortname2)) {
  download.file(paste0("https://rest.kegg.jp/link/", i, "/enzyme"), destfile = paste0("kegg_uc/", i, "_ec.txt"))
}
species_kegg_ec <- c()
for (f in list.files("kegg_uc/", pattern = "*_ec.txt", full.names = T)) {
  tmp <- read.table(f, header = F, sep = "\t")
  tmp$V1 <- gsub("ec:", "", tmp$V1)
  species_kegg_ec <- rbind(species_kegg_ec, tmp[tmp$V1 %in% ec_info$ID, ])
}
species_kegg_ec$shortname <- lapply(species_kegg_ec$V2, function(x) 
  strsplit(x, ":") %>% unlist %>% .[1]) %>% unlist()
species_kegg_ec2 <- merge(species_kegg_ec, species_info2, by.x = "shortname", by.y = "shortname2", all.y = T)
species_kegg_ec2$value <- 1
library(reshape2)
species_kegg_ec3 <- dcast(species_kegg_ec2[, c("V1", "ID", "value")], ID~V1, fill = 0)
rownames(species_kegg_ec3) <- species_kegg_ec3$ID
species_kegg_ec3$ID <- NULL


d_heatmap <- dcast(ec_species, Target~Source, fill = 0)
rownames(d_heatmap) <- d_heatmap$Target
d_heatmap$Target <- NULL
library(pheatmap)
pheatmap(d_heatmap)
# Generate annotations for rows and columns
table(pathway_ec$Target)
mgp_dsumm2 <- mgp_dsumm[ec_uccon$name[ec_uccon$ec %in% ec_info$ID], ]
mgp_dsumm2$names <- rownames(mgp_dsumm2)
rownames(mgp_dsumm2) <- lapply(rownames(mgp_dsumm2), function(x) 
  strsplit(x, ":") %>% unlist %>% .[1]) %>% unlist()
pathway_ec2 <- pathway_ec[!pathway_ec$Target %in% "3.5.1.24", ]
rownames(pathway_ec2) <- pathway_ec2$Target
pathway_ec2["3.5.1.24", ] <- c("mix", "3.5.1.24", 1)
annotation_col = data.frame(
  diff_uc_con = mgp_dsumm2[colnames(d_heatmap), "diff_uc_con"], 
  pathway = pathway_ec2[colnames(d_heatmap), "Source"]
)
annotation_col$group <- ifelse(annotation_col$diff_uc_con > 0, "UC", "Control")
rownames(annotation_col) = colnames(d_heatmap)
annotation_col$ec <- colnames(d_heatmap)
annotation_col <- merge(annotation_col, pathway[, c("ID", "name")], by.x = "pathway", by.y = "ID", all.x = T)
annotation_col <- merge(annotation_col, ec_info[, c("ID", "name")], by.x = "ec", by.y = "ID", all.x = T)
rownames(annotation_col) <- annotation_col$ec
annotation_col["3.5.1.24", "name.x"] <- "mix" #map00120, map00121
annotation_col <- annotation_col[order(annotation_col$group, annotation_col$pathway), ]


## add specis:ec (kegg)
species_kegg_ec3[, setdiff(colnames(d_heatmap), colnames(species_kegg_ec3))] <- 0
species_kegg_ec3 <- species_kegg_ec3[rownames(d_heatmap), colnames(d_heatmap)]
species_kegg_ec3[species_kegg_ec3 == 0] <- ""
species_kegg_ec3[species_kegg_ec3 > 0] <- "+"

## add taxonomy
species_info3 <- merge(species_info2, distinct(gtdb_taxonomy[, 1:6]), by = "genus", all.x = T)
species_info3[is.na(species_info3$domain), ]
rownames(species_info3) <- species_info3$ID

sp_dsumm2 <- sp_dsumm[rownames(d_heatmap), ]
annotation_row = data.frame(
  diff_uc_con = sp_dsumm2[rownames(d_heatmap), "diff_uc_con"], 
  phylum = species_info3[rownames(d_heatmap), "phylum"],
  phylum = species_info3[rownames(d_heatmap), "class"],
  object = species_info3[rownames(d_heatmap), "object"]
)
rownames(annotation_row) = rownames(d_heatmap)
annotation_row$group <- ifelse(annotation_row$diff_uc_con > 0, "UC", "Control")
library(RColorBrewer)
ann_colors = list(
  group = c(UC = "#ed80e6", Control = "#6cba56"), 
  name.x = c('map00061:Fatty acid biosynthesis' = "#66C2A5", 
             'map00071:Fatty acid degradation' = "#E5C494",
             #'map00380:Tryptophan metabolism' = "#FC8D62", 
             'map00400:Phenylalanine, tyrosine and tryptophan biosynthesis' = "#d149f2", 
             'map00121:Secondary bile acid biosynthesis' = "#E78AC3", 
             'map00600:Sphingolipid metabolism' = "#A6D854", 
             'mix' = "#FFD92F")#brewer.pal(7, "Set2")
)

annotation_col$group <- ifelse(annotation_col$diff_uc_con > 0, "UC", 
                               ifelse(annotation_col$diff_uc_con < 0, "Control", NA))
annotation_row$group <- ifelse(annotation_row$diff_uc_con > 0, "UC", 
                               ifelse(annotation_row$diff_uc_con < 0,
                                      "Control", NA))
annotation_row2 <- annotation_row[complete.cases(annotation_row$group), ]
annotation_col2 <- annotation_col[complete.cases(annotation_col$group), ]
annotation_row2 <- annotation_row2[order(annotation_row2$group), ]
annotation_col2 <- annotation_col2[order(annotation_col2$group, annotation_col2$pathway), ]
pheatmap(d_heatmap[rownames(annotation_row2), rownames(annotation_col2)], annotation_col = annotation_col[, c("name.x", "group")], 
         annotation_row = annotation_row[, "group", drop = F], fontsize_number = 10, cluster_rows = F,
         labels_col = annotation_col$name.y, cluster_cols = F, annotation_colors = ann_colors,
         display_numbers = species_kegg_ec3[rownames(annotation_row2), rownames(annotation_col2)],
         filename = "sup6b_UCvsControl_diff_species_genes_correlation.pdf", width = 10.4, height = 4.2)
dev.off()

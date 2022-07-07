#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for GLM based clustering of ASVs
# Removing ASVs of importance and obtaining the most contributing taxa
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "RColorBrewer", "parallel", "cowplot")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))
load(paste0("./data/asvs_CAS_DE_asv_E.RData"))
load(paste0("./data/asvs_CAS_DE_asv_RP.RData"))
load(paste0("./data/asvs_CAS11_DE_asv_E.RData"))
load(paste0("./data/asvs_CAS11_DE_asv_RP.RData"))
load(paste0("./data/asvs_CAS13_DE_asv_E.RData"))
load(paste0("./data/asvs_CAS13_DE_asv_RP.RData"))
# load(paste0("./data/GH_DE_asv_RP.RData"))
# load(paste0("./data/GH_DE_asv_E.RData"))

# REad the taxonomy table from the ATSPhereCC
taxa_at_cc <- read.table("./output/taxonomy_syncom.txt", sep = "\t", header = TRUE, as.is = TRUE)

pro <- c("gh")

# pro <- names(meta$datasets)
figs.out <- paste0(figs, paste0("/gh/differential_analysis/uspiked/asv/", pro))
stats.out <- paste0(stats, paste0("/gh/differential_analysis/uspiked/asv/", pro))
out.out <- paste0(out, paste0("/gh/differential_analysis/uspiked/asv/", pro))

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

temp.meta <- rbind.data.frame(meta$datasets[["CAS"]]) %>%
filter(Compartment != "rhizosphere", !SampleID %in% outlier)
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$Replicate))
temp.meta$Experiment <- factor(as.character(temp.meta$Experiment))
temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# Combine the LFC data
idx <- unique(c(
  row.names(CAS_DE_asv_RP$logFC_P), row.names(CAS_DE_asv_E$logFC_P), 
  row.names(CAS11_DE_asv_RP$logFC_P), row.names(CAS11_DE_asv_E$logFC_P),
  row.names(CAS13_DE_asv_RP$logFC_P), row.names(CAS13_DE_asv_E$logFC_P)
))

lfc <- data.frame(tag = idx) %>% cbind.data.frame(., 
  CAS_DE_asv_RP$logFC_P[.$tag,],
  CAS_DE_asv_E$logFC_P[.$tag,],
  CAS11_DE_asv_RP$logFC_P[.$tag,],
  CAS11_DE_asv_E$logFC_P[.$tag,],
  CAS13_DE_asv_RP$logFC_P[.$tag,],
  CAS13_DE_asv_E$logFC_P[.$tag,]
)

# Use short terms for factoring
temp.meta$Genotype <- genotype$short[match(temp.meta$Genotype, genotype$names)]
temp.meta$Compartment <- compartment$short[match(temp.meta$Compartment, compartment$names)]
temp.meta$group <- paste(temp.meta$Genotype, temp.meta$Compartment, temp.meta$Soil_Batch, sep = "_")

# Sort the LFC matrix
idx <- grep(x = colnames(lfc), pattern = "_logFC")
lfc_i <- lfc[, idx]
colnames(lfc_i) <- str_replace(colnames(lfc_i), "_logFC", "")

# Sorting by LFC asv INDEX
d <- dist((lfc_i))
clust_obj <- hclust(as.dist(d), method = "ward.D2")
lfc_sorted <- rownames(lfc_i)[clust_obj$order]


# Overall mean relative abundance
# taxa.mat <- all_exp_cut$taxa
# taxa_level <- "lineage_family"

# taxa.mat_temp <- temp.mat %>% 
# add_column(tag = taxa.mat[match(row.names(.), taxa.mat$asvb), taxa_level], .before = 1) %>%
# gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
# group_by(tag, SampleID) %>%
# summarise(cRA = sum(RA)) %>%
# spread(key = SampleID, value = cRA, fill = 0, convert = FALSE) %>%
# data.frame(., stringsAsFactors = FALSE)

# row.names(taxa.mat_temp) <- taxa.mat_temp$tag
# asv.ra <- apply(taxa.mat_temp[,-1], 2, function(x) x/sum(x))
asv.ra <- apply(temp.mat[,match(colnames(temp.mat), temp.meta$SampleID)], 2, function(x) x/sum(x))
# idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))

ra <- as.data.frame(asv.ra) %>% 
add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = temp.meta$group[match(.$SampleID, temp.meta$SampleID)]) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(ra)[-1] <- paste0("meanra_", colnames(ra[,-1]), "_RA")
row.names(ra) <- ra$tag

# Integrate clusters, taxonomy with the LFC data and the ASV data
lfc_annotated <- cbind.data.frame(lfc, 
  lineage_family = all_exp_cut$taxa$lineage_family[match(row.names(lfc), all_exp_cut$taxa$asvb)],
  taxa_atsphere_strain_lineage = taxa_at_cc$taxa_atsphere_strain[match(row.names(lfc), taxa_at_cc$asvb)])

lfc_annotated$taxa_atsphere_strain_lineage <- str_replace_na(lfc_annotated$taxa_atsphere_strain_lineage, "Unclassified")

# Relative abundance
lfc_annotated.ra <- cbind.data.frame(lfc_annotated, 
  ra[row.names(lfc_annotated), which(colnames(ra) != "tag")])

# Make list of essentials
asv_clustered <- list(

  lfc_annotated = lfc_annotated,
  lfc_annotated.ra = lfc_annotated.ra,
  ra.table = ra,
  lfc_sorted = lfc_sorted

) # Use these objects to spill the canvas in a new sheet

# Make a composite heatmap
save(list = "asv_clustered", file = "./data/GLM_clustered_asv.RData")

# END OF SCRIPT
sessionInfo()
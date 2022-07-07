#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "reshape2", "sva")

lapply(pkgs, require, character.only = T)

# Loading feature data
asv_tab <- read.table(paste0(asv.16s, "ASV_table.txt"), sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
asv_lib <- read.table(paste0(asv.16s, "ASV_lib_stat.txt"), sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
# asv_blast <- read.csv("./output/asv_blast_output.csv", sep = ",", header = F, stringsAsFactors = F)
# asv_blast.cdhit <- read.csv("./output/asv_blast_output.cdhit.csv", sep = "\t", header = F, stringsAsFactors = F)

# Filter the spike-ins
spike <- asv_tab[spike_16s[,2],]
stats_spike <- as.data.frame(t(spike))
spike_detects <- row.names(stats_spike)[stats_spike[,1] != 0]
spike_id <- row.names(spike)

p <- data.frame(total_asvs = colSums(asv_tab[row.names(asv_tab) %in% spike_id,]), SampleIDs = as.factor(names(colSums(asv_tab)))) %>%
filter(total_asvs > 0) %>%
ggplot(aes(x = total_asvs, y = SampleIDs)) +
# geom_vline(xintercept = b, lwd = 1, lty = "solid", colour = c_red) +
geom_bar(stat = "identity") +
# scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black)) +
# scale_fill_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black)) +
theme_RTN +
theme(axis.text.y = element_text(size = 2))
ggsave(p, filename = "./figures/diagnosis/ASVS_across_samples_spike_in.png",
  units = img,
  device = "png",
  width = 3,
  height = 7,
  bg = "white")

# Remove the sipike reads
asv_tab <- asv_tab[!row.names(asv_tab) %in% spike_id, ]

# Metadata 
taxa_s <- read.table(paste0(asv.16s, "ASV_taxonomy_silva.txt"), sep = "\t", stringsAsFactors = F, fill = T, header = T)
taxa_gg <- read.table(paste0(asv.16s, "ASV_taxonomy_gg.txt"), sep = "\t", stringsAsFactors = F, fill = T, header = T)
taxa_rdp <- read.table(paste0(asv.16s, "ASV_taxonomy_rdp.txt"), sep = "\t", stringsAsFactors = F, fill = T, header = T)
taxa_syncom <- read.table(paste0(asv.16s, "ASV_robust_syncom.txt"), sep = "\t", stringsAsFactors = F, fill = T, header = T)
map_syncom <- read.table(map_16s_cc, sep = "\t", stringsAsFactors = F, fill = T, header = T)
#map <- read.table(paste0(asv.16s, "asvb_map.txt"), sep = "\t", stringsAsFactors = F, fill = T, header = F)

nb_s <- row.names(taxa_s)[which(is.na(taxa_s$Kingdom))]
nb_gg <- row.names(taxa_gg)[which(is.na(taxa_gg$Kingdom))]
nb_rdp <- row.names(taxa_rdp)[which(is.na(taxa_rdp$Kingdom))]

# ASVs that are predicted to be in Syncom
asv_hit_syncom <- row.names(taxa_syncom)[!is.na(taxa_syncom$Species)]
asv_tab_sync <- asv_tab[asv_hit_syncom,]

# Remove noise
noise1 <- unique(

  apply(taxa_s, 2, function(x){

    temp <- row.names(taxa_s)[grep(x, pattern = "mitochondria|chloroplast|archaea|eukaryota", ignore.case = TRUE)]
    return(temp)
  }) %>% unlist(.)
)

noise2 <- unique(

  apply(taxa_gg, 2, function(x){

    temp <- row.names(taxa_gg)[grep(x, pattern = "mitochondria|chloroplast|archaea|eukaryota", ignore.case = TRUE)]
    return(temp)
  }) %>% unlist(.)
)

noise3 <- unique(

  apply(taxa_rdp, 2, function(x){

    temp <- row.names(taxa_rdp)[grep(x, pattern = "mitochondria|chloroplast|archaea|eukaryota", ignore.case = TRUE)]
    
    return(temp)
  }) %>% unlist(.)
)

noise <- unique(c(noise1, noise2, noise3, nb_rdp, nb_gg, nb_s))
taxa <- taxa_s[which(!row.names(taxa_s) %in% noise), ]

# Write noise
taxa[which(row.names(taxa) %in% noise), ] %>% write.table(., paste(out, "noise_dataset.txt", sep = ""), sep = "\t", quote = F)

# Remove noise elements
asv_tab <- asv_tab[-which(row.names(asv_tab) %in% noise),]

# ASV libs stats summary
b = 500


p <- data.frame(total_asvs = colSums(asv_tab), SampleIDs = as.factor(names(colSums(asv_tab)))) %>%
mutate(cut = log10(total_asvs+1) < log10(b)) %>%
ggplot(aes(x = log10(total_asvs+1), y = SampleIDs)) +
geom_vline(xintercept = log10(b), lwd = 1, lty = "solid", colour = c_red) +
geom_bar(stat = "identity", aes(colour = cut, fill = cut)) +
scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black)) +
scale_fill_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black)) +
theme_RTN +
theme(axis.text.y = element_text(size = 2))
ggsave(p, filename = "./figures/diagnosis/ASVS_across_samples.png",
  units = img,
  device = "png",
  width = 3,
  height = 7,
  bg = "white")

parallel::mclapply(colnames(asv_tab), function(x){


  p <- qplot(log10(asv_tab[,x]+1), main = x, geom = "histogram", bins = 100)
  ggsave(p, filename = paste0("./figures/diagnosis/ASV_histogram_", x,".png"),
  units = img,
  device = "png",
  width = 4,
  height = 2,
  bg = "white")

}, mc.cores = 24)


# Optimisation of sample size
id <- row.names(asv_lib)[asv_lib$Abundance >= b]
message(paste("No. of samples removed at cutoff ", b, " : ", nrow(asv_lib) - length(id), " out of ", nrow(asv_lib), "", sep = ""))
idx <- which(design$SampleID %in% id)
asv_tab <- asv_tab[,id]
metadata <- design[idx,]

colnames(metadata) <- str_replace(colnames(metadata), "Study", "Experiment")
colnames(metadata) <- str_replace(colnames(metadata), "Host_", "")

# Experiment filters
CAS_pro <- metadata %>% dplyr::filter(Experiment %in% c("pro01", "pro3", "pro20"), 
!Genotype %in% c("Capsella_rubella", "nai2-2", "bhlh04/05/06"))

LRO_gh_pro <- metadata %>% dplyr::filter(Experiment %in% c("pro8"), 
!Genotype %in% c("Capsella_rubella", "nai2-2", "bhlh04/05/06"))

LRO_all <- metadata %>% dplyr::filter(Experiment %in% c("pro8", "pro15"), 
!Genotype %in% c("Capsella_rubella", "nai2-2", "bhlh04/05/06"))

LRO_pro <- metadata %>% dplyr::filter(Experiment %in% c("pro15"))

ITA_pro <- metadata %>% dplyr::filter(Experiment %in% c("pro9"))

strex_pro <- metadata[str_detect(metadata$Experiment, "^STREX"),]

pro_all <- metadata %>% dplyr::filter(!Genotype %in% c("Capsella_rubella", "nai2-2", "bhlh04/05/06"))

# Note the samples with lower reads in Run 1 and filtered out and the new ID samples are replaced
design_strex$genotype <- str_replace_all(design_strex$genotype, "^pyk$", "pyk10")
design_strex$genotype <- factor(design_strex$genotype, levels = genotype$short[genotype$short %in% design_strex$genotype])
design_strex$experiment <- as.factor(design_strex$experiment)
design_strex$tech_rep <- as.factor(design_strex$tech_rep)
design_strex$well_id <- as.factor(design_strex$well_id)
design_strex$replicate <- as.factor(design_strex$replicate)
design_strex$Treatment <- as.factor(design_strex$Treatment)
design_strex$random <- as.factor(paste(design_strex$experiment, design_strex$tech_rep, design_strex$replicate, design_strex$well_id, sep = "_"))

levels <- sapply(unique(design_strex$Treatment), function(x) {
  paste(x, (genotype$short[which(genotype$short %in% design_strex$genotype)]), sep = "_")
}, simplify = TRUE) %>% as.vector(.)
design_strex$group <- as.factor(paste(design_strex$Treatment, design_strex$genotype, sep = "_"))
levels(design_strex$group) <- levels[which(levels %in% as.character(design_strex$group))]

# Output of the metadata operation
meta <- list(
  datasets = list(
    CAS = CAS_pro, 
    LRO = LRO_gh_pro, 
    LRO_all = LRO_all,
    LRO_field = LRO_pro, 
    ITA = ITA_pro, 
    STREX = strex_pro, 
    ALL = pro_all
    ),
  only_spiked = list(stats_spike, spike_detects),
  design_strex = design_strex
)

#===== asv - Realtive abundance =====
# Rarefy -- Remove low abundant amplicons
ra <- apply(asv_tab, 2, function(x) x/sum(x))
# idx <- rowSums(ra * 100 > threshold) >= 1
# asv_cut <- asv_tab[idx,]
asv_cut <- asv_tab

# Absolute abundance calculation
idy <- as.character(strex_pro$SampleID)
aa <- sweep(as.matrix(asv_tab)[names(idx), idy], 2, 1+as.matrix(spike[,idy]), "/")
idx <- which(rowSums(abs(aa)) > 0)
idy <- row.names(t(spike))[t(spike) > 1]
aa <- as.data.frame(aa)[names(idx),(idy)]

# Make a relative abundance table
ra_cut <- apply(asv_cut, 2, function(x) x/sum(x))


#===== asv syncom - Realtive abundance =====
# Rarefy -- Remove low abundant amplicons
ra_sync <- apply(asv_tab_sync, 2, function(x) x/sum(x))
idx <- names(rowSums(ra_sync * 100 > threshold) >= 1)
asv_cut_sync <- asv_tab_sync[idx,]
ra_cut_sync <- apply(asv_cut_sync, 2, function(x) x/sum(x))

# Absolute abundance calculation
idy <- as.character(strex_pro$SampleID)
asv_spike_sync <- asv_tab_sync[,spike_detects]
aa_sync <- sweep(as.matrix(asv_tab_sync)[names(idx),idy], 2, 1+as.matrix(spike[,idy]), "/")
idx <- rowSums(aa_sync * 100 > threshold) >= 1
aa_sync <- aa_sync[idx,]

# Make a relative abundance table

# Write table fors of output

# ASV
# Unspiked
write.table(file = "./output/asv.ra.txt", x = ra, sep = "\t", quote = F)
write.table(file = "./output/asv.ra.cut.txt", x = ra_cut, sep = "\t", quote = F)

write.table(file = "./output/asv.aa.txt", x = asv_tab, sep = "\t", quote = F)
write.table(file = "./output/asv.aa.cut.txt", x = asv_cut, sep = "\t", quote = F)

# spiked
write.table(file = "./output/asv.spiked.aa.txt", x = aa, sep = "\t", quote = F)

# ASV that are mapped to AtSphere CC using assignSpecies
# Unspiked
write.table(file = "./output/asv.aa.sync.txt", x = asv_tab_sync, sep = "\t", quote = F)
write.table(file = "./output/asv.aa.sync.cut.txt", x = asv_cut_sync, sep = "\t", quote = F)

write.table(file = "./output/asv.ra.sync.txt", x = ra_sync, sep = "\t", quote = F)
write.table(file = "./output/asv.ra.sync.cut.txt", x = ra_cut_sync, sep = "\t", quote = F)

# spiked
write.table(file = "./output/asv.spiked.aa.sync.txt", x = aa_sync, sep = "\t", quote = F)

# Make a seperate function or R script for COMBAT normalization
# send the fixed group and the known batch to the function to normalise
# NOTE:: USE COMBAT ONLY WHILE INTEGRATING THE HORIZONTAL DATASET
# Subset data for separate analysis

# # ========= COMBAT
# idx <- intersect(metadata$SampleID, colnames(asv_tab))
# asv_cut <- asv_cut[,idx]

# ra <- apply(asv_cut, 2, function(x) x/sum(x))
# idx <- rowSums(ra * 100 > threshold) >= 1

# asv_cut <- asv_cut[idx,]
# # ra_cut <- ra_cut[idx,]


# # Remove the low abundant ASVs within filtered group

# # Make ComBat normalisation and correct for the batch
# levels <- lapply(genotype$names, function(x) paste(x, soilbatch$names, sep = "_")) %>% unlist(., use.names = F)
# metadata$group <- paste(metadata$Genotype, metadata$Soil_Batch, sep = "_")
# metadata$group <- factor(metadata$group, levels = levels[levels %in% metadata$group])
# metadata$Soil_Batch <- as.factor(metadata$Soil_Batch)

# levels <- lapply(genotype$names, function(x) paste(x, compartment$names, sep = "_")) %>% unlist(., use.names = F)
# metadata$fixed_group <- paste(metadata$Genotype, metadata$Compartment, sep = "_")
# metadata$fixed_group <- factor(metadata$fixed_group, levels = levels[levels %in% metadata$fixed_group])
# metadata$batch <- as.factor(paste(metadata$Soil_Batch, metadata$Experiment, sep = "_"))


# # Normalisng batch effect
# mod_fixed <- model.matrix(~fixed_group, metadata)
# # combat_obj <- sva::ComBat(dat = asv_cca,
# #   batch = dat$batch, 
# #   mod = mod_fixed, 
# #   par.prior = TRUE, 
# #   prior.plots = FALSE)
# combat_obj <- sva::ComBat_seq(counts = as.matrix(asv_cut),
#   covar_mod = mod_fixed, 
#   group = NULL,
#   batch = metadata$batch)

# asv_cut <- combat_obj
# ra_cut <- apply(asv_cut, 2, function(x) x/sum(x))


# write.table(file = "./output/asv.ra.cut.txt", x = ra_cut, sep = "\t", quote = F)
# write.table(file = "./output/asv.cut.txt", x = asv_cut, sep = "\t", quote = F)

# # COMBAT =========

#===== Taxonomy table =====
taxa <- taxa %>% add_column(asvb = row.names(taxa), .before = 1) %>%
  filter(!Kingdom %in% c("Eukaryota", "Archaea"))
row.names(taxa) <- taxa$asvb

taxa_cut <- taxa[row.names(taxa) %in% row.names(asv_cut),]
asv_sorted <- taxa %>% group_by(Kingdom, Phylum, Class, Order, Family, asvb) %>% 
                            summarise(lineage_asvb = 
                                      paste(Kingdom, Phylum, Class, Order, Family, asvb, 
                                            sep = ";")) %>% arrange(lineage_asvb) %>%
                            data.frame(.)

asv_sorted_cut <- asv_sorted[asv_sorted$asvb %in% row.names(asv_cut),]

family_sorted <- taxa %>% group_by(Kingdom, Phylum, Class, Order, Family) %>% 
                            mutate(lineage_family = 
                                      paste(Kingdom, Phylum, Class, Order, Family, 
                                            sep = ";")) %>% arrange(lineage_family) %>%
                            data.frame(.)

family_sorted_cut <- family_sorted[family_sorted$asvb %in% row.names(asv_cut),]

taxa$lineage_asvb <- asv_sorted$lineage_asv[match(taxa$asvb, asv_sorted$asvb)]
taxa$lineage_family <- family_sorted$lineage_family[match(taxa$asvb, family_sorted$asvb)]
taxa <- taxa %>% mutate(tag_family = paste(Class, Order, Family, asvb, sep = ";"))

taxa_cut$lineage_asvb <- asv_sorted_cut$lineage_asv[match(taxa_cut$asvb, asv_sorted_cut$asvb)]
taxa_cut$lineage_family <- family_sorted_cut$lineage_family[match(taxa_cut$asvb, family_sorted_cut$asvb)]

taxa_cut <- taxa_cut %>% mutate(tag_family = paste(Class, Order, Family, asvb, sep = ";"))

# Spiked
taxa_cut_aa <- taxa[row.names(taxa) %in% row.names(aa),]
asv_sorted <- taxa %>% group_by(Kingdom, Phylum, Class, Order, Family, asvb) %>% 
                            summarise(lineage_asvb = 
                                      paste(Kingdom, Phylum, Class, Order, Family, asvb, 
                                            sep = ";")) %>% arrange(lineage_asvb) %>%
                            data.frame(.)

asv_sorted_cut <- asv_sorted[asv_sorted$asvb %in% row.names(aa),]

family_sorted <- taxa %>% group_by(Kingdom, Phylum, Class, Order, Family) %>% 
                            mutate(lineage_family = 
                                      paste(Kingdom, Phylum, Class, Order, Family, 
                                            sep = ";")) %>% arrange(lineage_family) %>%
                            data.frame(.)

family_sorted_cut <- family_sorted[family_sorted$asvb %in% row.names(aa),]

taxa$lineage_asvb <- asv_sorted$lineage_asv[match(taxa$asvb, asv_sorted$asvb)]
taxa$lineage_family <- family_sorted$lineage_family[match(taxa$asvb, family_sorted$asvb)]
taxa <- taxa %>% mutate(tag_family = paste(Class, Order, Family, asvb, sep = ";"))

taxa_cut_aa$lineage_asvb <- asv_sorted_cut$lineage_asv[match(taxa_cut_aa$asvb, asv_sorted_cut$asvb)]
taxa_cut_aa$lineage_family <- family_sorted_cut$lineage_family[match(taxa_cut_aa$asvb, family_sorted_cut$asvb)]

taxa_cut_aa <- taxa_cut_aa %>% mutate(tag_family = paste(Class, Order, Family, asvb, sep = ";"))

# AtSphere CC

# Unspiked
taxa_atcc <- taxa_syncom[asv_hit_syncom, ] %>%
add_column(asvb = row.names(.), .before = 1) %>%
select(asvb, Species) %>%
mutate(taxa_db = taxa_cut$lineage_asvb[match(.$asvb, taxa_cut$asvb)],
  taxa_db_lineage = taxa_cut$lineage[match(.$asvb, taxa_cut$asvb)]
  ) %>%
cbind.data.frame(., map_syncom[match(.$Species, map_syncom$ID), c("kingdom","phylum","class","order", "family")]) %>%
mutate(taxa_atsphere = paste(kingdom, phylum, class, order, family, asvb, sep = ";"),
  taxa_atsphere_lineage = paste(kingdom, phylum, class, order, family, sep = ";"),
  taxa_atsphere_strain = paste(kingdom, phylum, class, order, family, Species, asvb, sep = ";"),
  taxa_atsphere_strain_lineage = paste(kingdom, phylum, class, order, family, Species, sep = ";")
)
row.names(taxa_atcc) <- taxa_atcc$asvb

# Spiked
taxa_atcc_spike <- taxa_atcc[row.names(taxa_atcc) %in% row.names(aa_sync),]

# Save the taxonomy in table format
write.table(file = "./output/taxonomy.txt", taxa, sep = "\t", quote = F)
write.table(file = "./output/taxonomy_cut.txt", taxa_cut, sep = "\t", quote = F)
write.table(file = "./output/taxonomy_spiked_cut.txt", taxa_cut_aa, sep = "\t", quote = F)
write.table(file = "./output/taxonomy_syncom.txt", taxa_atcc, sep = "\t", quote = F)
write.table(file = "./output/taxonomy_syncom_spiked.txt", taxa_atcc_spike, sep = "\t", quote = F)

#===== RA-Table =====

# Unspiked
ra_tab <- melt(ra)
colnames(ra_tab) <- c("asvb", "SampleID", "RA")
write.table(file = "./output/ra_table.txt", x = ra_tab, sep = "\t", quote = F)

ra_tab_cut <- melt(ra_cut)
colnames(ra_tab_cut) <- c("asvb", "SampleID", "RA")
write.table(file = "./output/ra_table.cut.txt", x = ra_tab_cut, sep = "\t", quote = F)

# Spiked
aa_tab <- melt(aa)
colnames(aa_tab) <- c("asvb", "SampleID", "AA")
write.table(file = "./output/aa_table.txt", x = aa_tab, sep = "\t", quote = F)

# Unspiked
ra_tab <- melt(ra_sync)
colnames(ra_tab) <- c("asvb", "SampleID", "RA")
write.table(file = "./output/ra_table_sync.txt", x = ra_tab, sep = "\t", quote = F)

ra_tab_cut <- melt(ra_cut_sync)
colnames(ra_tab_cut) <- c("asvb", "SampleID", "RA")
write.table(file = "./output/ra_table.cut_sync.txt", x = ra_tab_cut, sep = "\t", quote = F)

# Spiked
aa_tab <- melt(aa_sync)
colnames(aa_tab) <- c("asvb", "SampleID", "AA")
write.table(file = "./output/aa_table_sync.txt", x = aa_tab, sep = "\t", quote = F)


# ===== ASV Taxa lineage =====

# Unspiked
asv_taxa <- as.data.frame(asv_cut) %>% add_column(asvb = row.names(.), .before = 1)
asv_taxa <- melt(asv_taxa)
colnames(asv_taxa) <- c("asvb", "SampleID", "counts")
asv_taxa$lineage_family <- taxa_cut$lineage_family[match(asv_taxa$asvb, taxa_cut$asvb)]
asv_taxa <- asv_taxa %>% group_by(lineage_family, SampleID) %>% summarise(counts = sum(counts)) %>% spread(key = "SampleID", value = "counts", convert = FALSE) %>% data.frame(., stringsAsFactors = F)
asv_taxa$lineage_family[is.na(asv_taxa$lineage)] <- "Undetermined"
rownames(asv_taxa) <- asv_taxa$lineage_family
asv_taxa <- asv_taxa[,-which(colnames(asv_taxa) == "lineage_family")]

# ra_taxa <- apply(asv_taxa, 2, function(x) x/sum(x))
# idx <- rowSums(ra_taxa * 100 > threshold) >= 1
# asv_taxa_cut <- asv_taxa[idx,]
asv_taxa_cut <- asv_taxa
ra_taxa_cut <- apply(as.matrix(asv_taxa_cut), 2, function(x) x/sum(x))

ra_tab_taxa_cut <- as.data.frame(ra_taxa_cut) %>% add_column(asvb = row.names(asv_taxa_cut), .before = 1) %>% data.frame(.)
ra_tab_taxa_cut <- melt(ra_tab_taxa_cut)
colnames(ra_tab_taxa_cut) <- c("asvb", "SampleID", "RA")

write.table(file = "./output/ra.taxa.cut.txt", x = ra_taxa_cut, sep = "\t", quote = F)
write.table(file = "./output/ra_table.taxa.cut.txt", x = ra_tab_taxa_cut, sep = "\t", quote = F)
write.table(file = "./output/asv.taxa.cut.txt", x = asv_taxa, sep = "\t", quote = F)

# Absolute
asv_taxa <- as.data.frame(aa) %>% add_column(asvb = row.names(.), .before = 1)
asv_taxa <- melt(asv_taxa)
colnames(asv_taxa) <- c("asvb", "SampleID", "counts")
asv_taxa$lineage_family <- taxa_cut_aa$lineage_family[match(asv_taxa$asvb, taxa_cut_aa$asvb)]
asv_taxa <- asv_taxa %>% group_by(lineage_family, SampleID) %>% 
summarise(counts = sum(counts)) %>% 
spread(key = "SampleID", value = "counts", convert = FALSE) %>% 
data.frame(., stringsAsFactors = F)
asv_taxa$lineage_family[is.na(asv_taxa$lineage)] <- "Undetermined"
rownames(asv_taxa) <- asv_taxa$lineage_family
asv_taxa <- asv_taxa[,-which(colnames(asv_taxa) == "lineage_family")]

# aa_taxa <- apply(asv_taxa, 2, function(x) sum(x))
# idx <- rowSums(asv_taxa) >= 1
# asv_taxa_cut <- asv_taxa[idx,]
asv_taxa_cut <- asv_taxa
aa_taxa_cut <- apply(as.matrix(asv_taxa_cut), 2, function(x) x/sum(x))

aa_tab_taxa_cut <- as.data.frame(aa_taxa_cut) %>% add_column(asvb = row.names(asv_taxa_cut), .before = 1) %>% data.frame(.)
aa_tab_taxa_cut <- melt(aa_tab_taxa_cut)
colnames(aa_tab_taxa_cut) <- c("asvb", "SampleID", "AA")

write.table(file = "./output/aa.spiked.taxa.cut.txt", x = aa_taxa_cut, sep = "\t", quote = F)
write.table(file = "./output/aa.spiked_table.taxa.cut.txt", x = aa_tab_taxa_cut, sep = "\t", quote = F)
write.table(file = "./output/asv.spiked.taxa.cut.txt", x = asv_taxa, sep = "\t", quote = F)

# Syncom
asv_taxa <- as.data.frame(ra_cut_sync) %>% add_column(asvb = row.names(.), .before = 1)
asv_taxa <- melt(asv_taxa)
colnames(asv_taxa) <- c("asvb", "SampleID", "counts")
asv_taxa$lineage_family <- taxa_atcc$taxa_atsphere_strain_lineage[match(asv_taxa$asvb, taxa_atcc$asvb)]
asv_taxa <- asv_taxa %>% group_by(lineage_family, SampleID) %>% 
summarise(counts = sum(counts)) %>% 
spread(key = "SampleID", value = "counts", convert = FALSE) %>% 
data.frame(., stringsAsFactors = F)
asv_taxa$lineage_family[is.na(asv_taxa$lineage)] <- "Unclassified"
rownames(asv_taxa) <- asv_taxa$lineage_family
asv_taxa <- asv_taxa[,-which(colnames(asv_taxa) == "lineage_family")]

# aa_taxa <- apply(asv_taxa, 2, function(x) x/sum(x))
# idx <- rowSums(aa_taxa * 100 > threshold) >= 1
# asv_taxa_cut <- asv_taxa[idx,]
asv_taxa_cut <- asv_taxa
ra_taxa_cut_syn <- apply(as.matrix(asv_taxa_cut), 2, function(x) x/sum(x))

ra_tab_taxa_cut_syn <- as.data.frame(ra_taxa_cut_syn) %>% 
add_column(asvb = row.names(asv_taxa_cut), .before = 1) %>% data.frame(.)
ra_tab_taxa_cut_syn <- melt(ra_tab_taxa_cut_syn)
colnames(ra_tab_taxa_cut_syn) <- c("asvb", "SampleID", "RA")

write.table(file = "./output/ra.syncom.taxa.cut.txt", x = ra_taxa_cut_syn, sep = "\t", quote = F)
write.table(file = "./output/ra.syncom_table.taxa.cut.txt", x = ra_tab_taxa_cut_syn, sep = "\t", quote = F)
write.table(file = "./output/asv.syncom.taxa.cut.txt", x = asv_taxa, sep = "\t", quote = F)

# Syncom -- Absolute
asv_taxa <- as.data.frame(aa_sync) %>% add_column(asvb = row.names(.), .before = 1)
asv_taxa <- melt(asv_taxa)
colnames(asv_taxa) <- c("asvb", "SampleID", "counts")
asv_taxa$lineage_family <- taxa_atcc_spike$taxa_atsphere_strain_lineage[match(asv_taxa$asvb, taxa_atcc_spike$asvb)]
asv_taxa <- asv_taxa %>% group_by(lineage_family, SampleID) %>% 
summarise(counts = sum(counts)) %>% 
spread(key = "SampleID", value = "counts", convert = FALSE) %>% 
data.frame(., stringsAsFactors = F)
asv_taxa$lineage_family[is.na(asv_taxa$lineage)] <- "Unclassified"
rownames(asv_taxa) <- asv_taxa$lineage_family
asv_taxa <- asv_taxa[,-which(colnames(asv_taxa) == "lineage_family")]

# aa_taxa <- apply(asv_taxa, 2, function(x) sum(x))
# idx <- rowSums(asv_taxa) >= 1
# asv_taxa_cut <- asv_taxa[idx,]
asv_taxa_cut <- asv_taxa
aa_taxa_cut_syn <- apply(as.matrix(asv_taxa_cut), 2, function(x) x/sum(x))

aa_tab_taxa_cut_syn <- as.data.frame(aa_taxa_cut_syn) %>% 
add_column(asvb = row.names(asv_taxa_cut), .before = 1) %>% data.frame(.)
aa_tab_taxa_cut_syn <- melt(aa_tab_taxa_cut_syn)
colnames(aa_tab_taxa_cut_syn) <- c("asvb", "SampleID", "AA")

write.table(file = "./output/aa.spiked.syncom.taxa.cut.txt", x = aa_taxa_cut_syn, sep = "\t", quote = F)
write.table(file = "./output/aa.spiked.syncom_table.taxa.cut.txt", x = aa_tab_taxa_cut_syn, sep = "\t", quote = F)
write.table(file = "./output/asv.spiked.syncom.taxa.cut.txt", x = asv_taxa, sep = "\t", quote = F)

# ===== Sort-samples ====

# Unspiked
d <- 1 - cor(ra_cut)
hc <- hclust(as.dist(d), method = "ward.D")
samples_sorted <- colnames(ra_cut)[hc$order]

# ===== Sort-samples-taxa_prevelence ====
d <- 1 - cor(ra_taxa_cut)
hc <- hclust(as.dist(d), method = "ward.D")
samples_sorted_taxa <- colnames(ra_taxa_cut)[hc$order]

# # ===== Sort-ASVs =====
# d <- dist((ra_cut))
# hc <- hclust(as.dist(d), method = "ward.D")
# asv_cut_clusters <- row.names(ra_cut)[hc$order]

# # ===== Sort-ASVs-taxa_prevelence =====
# d <- dist((ra_taxa_cut))
# hc <- hclust(as.dist(d), method = "ward.D")
# asv_cut_clusters_taxa <- row.names(ra_taxa_cut)[hc$order]

clusters_unspiked <- list(samples_sorted = samples_sorted, 
  samples_sorted_taxa = samples_sorted_taxa
  # asv_cut_clusters = asv_cut_clusters,
  # asv_cut_clusters_taxa = asv_cut_clusters_taxa
)

# Make dashboard of experiments as an objects for the horizontal integration
all_exp <- list(

  metadata = meta,
  asv_ra = ra,
  asv = asv_tab,
  taxa = taxa

)

# Filtered dataset
all_exp_cut <- list(

  metadata = meta,
  asv_ra = ra_cut,
  asv = asv_cut,
  taxa = taxa_cut,
  samples_sorted = clusters_unspiked$samples_sorted,
  asv_sorted = clusters_unspiked$asv_cut_clusters
)

# Filtered dataset
all_exp_cut_sync <- list(

  metadata = meta,
  asv_ra = ra_cut_sync,
  asv = asv_cut_sync,
  taxa = aa_taxa_cut_syn
  # samples_sorted = clusters_unspiked$samples_sorted,
  # asv_sorted = clusters_unspiked$asv_cut_clusters
)

# Save objects of clusters
save(list = c("clusters_unspiked", "clusters_spiked"), 
  file = paste0("./data/asv_clusters_16s.Rdata"))

# Save objects of experiments
save(list = "all_exp",
  file = paste0("./data/all_experiments.RData"))

save(list = "gh_exp",
  file = paste0("./data/gh_experiments.RData"))

save(list = "all_exp_cut",
  file = paste0("./data/all_experiments_filtered_asvs.RData"))

save(list = "all_exp_cut_sync",
  file = paste0("./data/all_experiments_sync_filtered_asvs.RData"))

save(list = "meta",
  file = paste0("./data/metadata_asvs.RData"))

# !END OF SCRIPT
sessionInfo()
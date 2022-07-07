# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/parameters.R")
setwd(analysis.syncom)

# Loading required packages
pkgs <- c("tidyverse", "reshape2")

lapply(pkgs, require, character.only = T)

# Load potential contaminated samples
contaminated <- read.table("./data/contaminated_samples.txt", sep = "\t", header = TRUE, as.is = TRUE) %>% 
select(1) %>% 
unlist(., use.names = FALSE)

# Loading np matrix -> count data
mat <- read.table("./data/merged_strain_table.txt", sep = "\t", header = TRUE, as.is = TRUE, row.names = 1, comment.char = "")

# Filter low abundant strains
idy <- which(colnames(mat) %in% contaminated)
mat.ra <- mat[-which(row.names(mat) == "DH5alpha"), -idy]
mat.ra <- apply(mat.ra, 2, function(x) x/sum(x))
idx <- which(rowSums(mat.ra)*100 >= 0.01)

# Filtered and relative abundance is calculated
mat.filtered <- mat[names(idx), -idy]
mat.ra.filtered <- apply(mat.filtered, 2, function(x) x/sum(x))

# Compute absolute abundance from spike in
spike <- mat["DH5alpha", -idy]
mat.aa <- sweep(as.matrix(mat[names(idx),-idy]), 2, 1+as.matrix(spike), "/")
idx <- rowSums(abs(mat.aa)) <= 0
idy <- as.matrix(t(spike)) <= 1
mat.aa <- mat.aa[!idx, !idy]
mat.aa <- mat.aa[-which(row.names(mat) == "DH5alpha"),] # Remove the spike after normalisation

lib.size <- colSums(mat)

# Load metadata 
meta <- read.table("./data/metadata_table.txt", sep = "\t", header = TRUE, as.is = TRUE) %>%
mutate(SampleID = paste0("Sample_", SampleID)) %>%
filter(!SampleID %in% contaminated, !biol_rep %in% c(1, 4))

# Sample summary
meta_sum <- meta %>%
group_by(biol_rep, time, genotype, exudates_rep, dose) %>%
summarise(N = n()) %>%
data.frame(., stringsAsFactors = FALSE)

# Make factors

# Make a metadata
meta$genotype <- str_replace_all(meta$genotype, "pyk", "pyk10")
meta$genotype <- str_replace_all(meta$genotype, "inoculum", "inoc")
meta$lib <- "lib1"

meta <- meta %>% 
mutate(lib = as.factor(ifelse(Fwd_barcode %in% c(13:24), "lib2", lib))) %>%
mutate(
	random = as.factor(paste(lib, biol_rep, exudates_rep, sep = "_")),
	fixed = as.factor(paste(genotype, dose, time, sep = "_")),
	lib = as.factor(lib)
	) %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short),
	exudates_rep = as.factor(exudates_rep),
	biol_rep = as.factor(biol_rep),
	time = as.factor(as.character(time)),
	dose = factor(as.character(dose), levels = dosage.syncom$names)
	) %>%
data.frame(., stringsAsFactors = FALSE)

row.names(meta) <- meta$SampleID

# Load taxonomy data
taxa <- read.table("./data/taxonomy_filtered.txt", sep = "\t", header = TRUE, as.is = TRUE, row.names = 1)

idx <- which(row.names(taxa) %in% row.names(mat.filtered))
taxa.filtered <- taxa[idx,]

# Clustering samples and strains

# Based on absolute abundance
# Samples
mat.aa <- mat.aa[, which(colnames(mat.aa) %in% meta$SampleID)]
mat.filtered <- mat.filtered[, which(colnames(mat.filtered) %in% meta$SampleID)]

d <- 1 - cor(mat.aa)
hc <- hclust(as.dist(d), clust_method)
aa.sample.clusters <- colnames(mat.aa)[hc$order]

# Strains
d <- 1 - cor(t(mat.aa))
hc <- hclust(as.dist(d), clust_method)
aa.strains.clusters <- row.names(mat.aa)[hc$order]

# Based on relative abundance
# Samples
d <- 1 - cor(mat.filtered)
hc <- hclust(as.dist(d), clust_method)
ra.sample.clusters <- colnames(mat.aa)[hc$order]

# Strains
d <- dist(mat.filtered)
hc <- hclust(as.dist(d), clust_method)
ra.strains.clusters <- row.names(mat.aa)[hc$order]

# Make an object for storing the data

# Non-filtered data
sync.dat <- list(cData = mat, 
	metadata = meta, 
	taxa = taxa)

# Filtered data

# non-spiked
sync.filtered.dat <- list(cData = mat.filtered, 
	ra = mat.ra.filtered,
	metadata = meta[which(meta$SampleID %in% colnames(mat.filtered)),], 
	taxa = taxa.filtered)

# spiked
sync.spiked.dat <- list(cData = mat.aa, 
	log.aa = log2(mat.aa+0.0001),
	metadata = meta[which(meta$SampleID %in% colnames(mat.aa)),], 
	taxa = taxa[which(row.names(taxa) %in% row.names(mat.aa)),])

# Save the objects
save(list = "sync.dat", file = "./data/sync_data.Rdata")
save(list = "sync.filtered.dat", file = "./data/sync_filtered_data.Rdata")
save(list = "sync.spiked.dat", file = "./data/sync_filtered_spike_data.Rdata")

# !END OF SCRIPT
sessionInfo()
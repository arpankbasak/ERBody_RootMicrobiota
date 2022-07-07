# Script for differential analysis
# @Arpan Kumar Basak inherited from R.T Nakano

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
pkgs <- c("tidyverse", "reshape2", "vegan", "DEseq2", "edgeR")

lapply(pkgs, require, character.only = T)

source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/parameters.R")
setwd(analysis.syncom)

load("./data/sync_filtered_spike_data.Rdata")
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/functions.R")

mat <- sync.spiked.dat$cData
meta <- sync.spiked.dat$metadata %>% filter(biol_rep %in% c(2, 3))

figs <- paste0(figs, "/differential_analysis/spiked/")
stats <- paste0(stats, "/differential_analysis/spiked/")
out <- paste0(out, "/differential_analysis/spiked/")

# Grouping factors
levels <- c(sapply(dosage.syncom$names, paste, c("0", "24", "72"), sep = "_"))
levels <- c(sapply(genotype.syncom$short, paste, unique(levels), sep = "_"))
group <- factor(meta$fixed, levels = levels[levels %in% unique(as.character(meta$fixed))])
exp_rep <- as.factor(meta$random)

# Model design
model <- model.matrix(~ 0 + group + exp_rep) # Dropping intercept

colnames(model) <- str_replace(colnames(model), "group", "")

idx <- colSums(model) == 0 # separating terms with no interaction
model <- model[,!idx]
set.seed(seeder)
mat <- mat[, meta$SampleID]

# Contrast
contrasts.mat <- limma::makeContrasts(
    
    # Estimate the difference from incocula
    Col_0.05_24_competency = (Col_0.05_24 - inoc_0.05_0),
    Col_0.005_24_competency = (Col_0.005_24 - inoc_0.005_0),
    Col_0.05_72_competency = (Col_0.05_72 - inoc_0.05_0),
    Col_0.005_72_competency = (Col_0.005_72 - inoc_0.005_0),

    # Estimate the difference from incocula
    pyk10_0.05_24_competency = (pyk10_0.05_24 - inoc_0.05_0),
    pyk10_0.005_24_competency = (pyk10_0.005_24 - inoc_0.005_0),
    pyk10_0.05_72_competency = (pyk10_0.05_72 - inoc_0.05_0),
    pyk10_0.005_72_competency = (pyk10_0.005_72 - inoc_0.005_0),

    # Estimate the difference from incocula
    cyp_0.05_24_competency = (cyp_0.05_24 - inoc_0.05_0),
    cyp_0.005_24_competency = (cyp_0.005_24 - inoc_0.005_0),
    cyp_0.05_72_competencey = (cyp_0.05_72 - inoc_0.05_0),
    cyp_0.005_72_competency = (cyp_0.005_72 - inoc_0.005_0),

    # Estimate the difference for the mutants within 24h
    pyk10_0.005_24_genotype          = (  pyk10_0.005_24   - Col_0.005_24 ),
    cyp_0.005_24_genotype         = ( cyp_0.005_24   - Col_0.005_24 ),
    pyk10_0.05_24_genotype          = (  pyk10_0.05_24   - Col_0.05_24 ),
    cyp_0.05_24_genotype         = ( cyp_0.05_24   - Col_0.05_24 ),

    # Extimate the difference for the mutants within 72h
    pyk10_0.005_72_genotype          = (  pyk10_0.005_72   - Col_0.005_72 ),
    cyp_0.005_72_genotype         = ( cyp_0.005_72   - Col_0.005_72 ),
    pyk10_0.05_72_genotype          = (  pyk10_0.05_72   - Col_0.05_72 ),
    cyp_0.05_72_genotype         = ( cyp_0.05_72   - Col_0.05_72 ),
    
    
levels=model)

# Create contrast matrix - Genotype effect on soiltype
contrasts.names <- attr(contrasts.mat, "dimnames")$Contrasts
n <- length(contrasts.names)

rownames(meta) <- meta$SampleID

# Create DGEList object for edgeR
delist <- DGEList(counts = mat, group = group)
# delist <- calcNormFactors(delist)

# Estimate common and tag-wise disperse
de <- estimateGLMCommonDisp(delist, model)
de <- estimateGLMTagwiseDisp(de, model)

mat.names <- row.names(de$counts)
fit <- glmFit(de, model)

# LRT for each contrasts
LRT.list <- lapply(contrasts.names, function(x) glmLRT(fit, contrast=contrasts.mat[, which(contrasts.names == x)]))
names(LRT.list) <- contrasts.names

# LogFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
  table <- LRT.list[[x]]$table[,c(1,4)]
  table <- table %>% dplyr::mutate(PValue = p.adjust(table[,2], method=p.adj.method))
  colnames(table) <- paste(contrasts.names[x], colnames(table), sep="_")
  return(table)
})

logFC_P <- do.call(data.frame, logFC_P.list)
row.names(logFC_P) <- mat.names

# Significance picking for each tested model
contrasts_idx <- intersect(names(LRT.list), contrasts.names)
DE.list <- sapply(contrasts_idx, function(x) decideTestsDGE(LRT.list[[x]], 
                                                              adjust.method=p.adj.method, 
                                                              p.value=alpha))

temp <- as.data.frame(DE.list, row.names = mat.names)

# Number of significant differentially abundant OTUs
count.mat <- data.frame(total = sapply(temp, function(x) sum(abs(x))),
                      enriched = sapply(temp, function(x) sum(x ==  1)),
                      depleted = sapply(temp, function(x) sum(x == -1)))


# Significance table
DE <- temp

# Counts per million
log2cpm <- cpm(de, prior.count = 2, log = T)

# Cluster LFC
d <- 1 - cor(t(logFC_P))
hc <- hclust(as.dist(d), clust_method)
strain_sorted <- row.names(logFC_P)[hc$order]

# Make object for easier computation later
de_analysis <- list(DE_strains = count.mat, 
      lfc_p = logFC_P,
      strain_clusters = strain_sorted,
      cpm_mat = log2cpm
)

# Save DE object for further analysis
save(list = "de_analysis", file = "./data/GLM_analysis_spiked.RData")

# Write table for statistics
write.table(logFC_P, paste0(stats, "/stats_logFC.P.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)

write.table(count.mat, paste0(stats, "/stats_number_of_enriched_depleted.txt"), 
            quote = F, row.names = T, col.names = NA, sep="\t")

write.table(DE, paste0(stats, "/stats_enrichment_test.significance_table.txt"), sep="\t", 
            quote = T, row.names = T, col.names = NA)

write.table(log2cpm, paste0(stats, "/stats_log2cpm.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)

# !END of Script
sessionInfo()

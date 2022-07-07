#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.itsf)

# Loading required packages
pkgs <- c("tidyverse", "vegan", "RColorBrewer", "parallel", "DESeq2", "edgeR", "limma")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))
# load(paste0("./data/soilbatch_CAS_DE_asv.RData"))
# common_in_CAS <- soilbatch_DE_asv$asv_list$similar_asvs

# What are we anlysing
pro <- c("STREX", "CAS")
# pro <- names(meta$datasets)

figs.out <- paste0(figs, "/gh/beta_diversity/uspiked/cpcoa/")
stats.out <- paste0(stats, "/gh/beta_diversity/uspiked/cpcoa/")
out.out <- paste0(out, "/gh/beta_diversity/uspiked/cpcoa/")

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

p = "st"

plot_path <- paste0(figs, "/st/beta_diversity/uspiked/cpcoa/")
stat_path <- paste0(stats, "/st/beta_diversity/uspiked/cpcoa/")
out_path <- paste0(out, "/st/beta_diversity/uspiked/cpcoa/")


mclapply(c(plot_path,stat_path,out_path), function(x){

    exp_path <- paste0(x, p)
    if(!dir.exists(paths = exp_path)){

        message(paste0("Directory created ", exp_path))
        dir.create(exp_path, recursive = FALSE)

    } else{
      
      message(paste0("Directory exists ", exp_path))

    }

}, mc.cores = 4)

plot_path <- paste0(plot_path, p)
stat_path <- paste0(stat_path, p)
out_path <- paste0(out_path, p)

# Prepare the ASV and metadata
temp.meta <- (meta$datasets[["STREX"]])
temp.meta <- data.frame(meta$design_strex[match(temp.meta$SampleID, meta$design_strex$SampleID),]) %>% filter(Treatment != "Extract")
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
# temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$replicate))
temp.meta$tech_rep <- as.factor(as.character(temp.meta$tech_rep))
temp.meta$Experiment <- factor(as.character(temp.meta$experiment))
# temp.meta$Treatment <- factor(as.character(temp.meta$Treatment), levels = c("Bulk soil", "Exudate", "Extract"))
# temp.meta$Run <- as.factor(as.character(temp.meta$Run))
asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv <- temp.mat[idx,]

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$short %in% as.character(temp.meta$genotype))
# idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
# idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

temp.meta$Genotype <- factor(as.character(temp.meta$genotype), levels = genotype$short[idg])

# replace genotype and compartment names with simpler ones
idx <- match(temp.meta$Genotype, genotype$names)
idy <- match(temp.meta$Compartment, compartment$names)

df <- temp.meta %>% 
  mutate(exp_rep = paste(Experiment, tech_rep, sep = "_"))

# Group those groups that are available
idx <- which(df$SampleID %in% colnames(asv))

df <- df[order(idx),]
asv <- asv[,idx]

# Grouping factors
exp_rep <- as.factor(df$exp_rep)
group <- df$Genotype

# Model design
model <- model.matrix(~ 0 + group + exp_rep) # Dropping intercept

colnames(model) <- str_replace(colnames(model), "group", "")

idx <- colSums(model) == 0 # separating terms with no interaction
model <- model[,!idx]
set.seed(1)

# Contrast
contrasts.mat <- limma::makeContrasts(
    
    strex_pyk10_exudate_genotype = pyk10 - Col,
    strex_cyp_exudate_genotype = cyp - Col,
    strex_MS_exudate_genotype = MS - Col,
    
    strex_Col_exudate_soil = Col - Soil,
    strex_pyk10_exudate_soil = pyk10 - Soil,
    strex_cyp_exudate_soil = cyp - Soil,
    strex_MS_exudate_soil = MS - Soil,
    
levels=model)

# Create contrast matrix - Genotype effect on soiltype
contrasts.names <- attr(contrasts.mat, "dimnames")$Contrasts
n <- length(contrasts.names)

rownames(df) <- df$SampleID

# Create DGEList object for edgeR
delist <- DGEList(counts = asv, group = group)
delist <- calcNormFactors(delist, method = "TMM")

# delist$samples %>% 
# ggplot(aes(x = norm.factors, fill = as.factor(group))) +
# geom_density(alpha = 0.5)

# Estimate common and tag-wise disperse
de <- estimateGLMCommonDisp(delist, model)
de <- estimateGLMTagwiseDisp(de, model)
#asv.norm <- t(t(de$pseudo.counts) * (de$samples$norm.factors))
asv.names <- row.names(de$counts)
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
row.names(logFC_P) <- asv.names

head(logFC_P)

write.table(logFC_P, paste0(stats.out, "/ST_stats_logFC.P_combat_Exudate.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)


# Significance picking for each tested model
contrasts_idx <- intersect(names(LRT.list), contrasts.names)
DE.list <- sapply(contrasts_idx, function(x) decideTestsDGE(LRT.list[[x]], 
                                                              adjust.method=p.adj.method, 
                                                              p.value=alpha))

temp <- as.data.frame(DE.list, row.names = asv.names)


# Number of significant differentially abundant OTUs
count.mat <- data.frame(total = sapply(temp, function(x) sum(abs(x))),
                      enriched = sapply(temp, function(x) sum(x ==  1)),
                      depleted = sapply(temp, function(x) sum(x == -1)))

write.table(count.mat, paste0(stats.out, "/ST_stats_number_of_enriched_depleted_combat_Exudate.txt"), 
            quote = F, row.names = T, col.names = NA, sep="\t")


# Significance table
DE <- temp

write.table(DE, paste0(stats.out, "/ST_stats_enrichment_test.significance_table_combat_Exudate.txt"), sep="\t", 
            quote = T, row.names = T, col.names = NA)

log2cpm <- cpm(de, prior.count = 2, log = T)

write.table(log2cpm, paste0(stats.out, "/ST_stats_log2cpm_combat_Exudate.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)


# Make object an RObject as output
ST_DE_asv_Exudate <- list(
  logFC_P = logFC_P, 
  DE_table = count.mat,
  cpm = log2cpm
)

save(list = "ST_DE_asv_Exudate", file = paste0("./data/ST_DE_asv_Exudate.RData"))

# !END of Script
sessionInfo()

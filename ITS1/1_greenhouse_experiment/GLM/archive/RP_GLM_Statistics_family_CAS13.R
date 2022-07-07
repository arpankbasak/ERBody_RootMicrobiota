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

# What are we anlysing
pro <- c("gh")

# pro <- names(meta$datasets)
figs.out <- paste0(figs, paste0("/gh/differential_analysis/uspiked/family/", pro))
stats.out <- paste0(stats, paste0("/gh/differential_analysis/uspiked/family/", pro))
out.out <- paste0(out, paste0("/gh/differential_analysis/uspiked/family/", pro))

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

temp.meta <- rbind.data.frame(meta$datasets[["CAS"]])
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$Replicate))
temp.meta$Experiment <- factor(as.character(temp.meta$Experiment))
temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$name %in% as.character(temp.meta$Genotype))
idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

# Make factors
temp.meta <- temp.meta %>%
filter(!Compartment %in% c("rhizosphere", "endosphere"), Soil_Batch == "CAS13") %>%
mutate(
    
    Genotype = factor(Genotype, levels = genotype$names[idg]),
    Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
    Compartment = factor(Compartment, levels = compartment$names[idc])
) %>% filter(!SampleID %in% outlier)

# Fetch taxa
taxa.mat <- all_exp_cut$taxa
taxa_level <- "lineage_family"
taxa.mat_temp <- temp.mat %>% 
add_column(tag = taxa.mat[match(row.names(.), taxa.mat$asvf), taxa_level], .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
group_by(tag, SampleID) %>%
summarise(cRA = sum(RA)) %>%
spread(key = SampleID, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

row.names(taxa.mat_temp) <- taxa.mat_temp$tag

# Remove low abudant ASVs
asv.ra <- apply(taxa.mat_temp[,-1], 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv <- taxa.mat_temp[idx,-1]

# replace genotype and compartment names with simpler ones
idx <- match(temp.meta$Genotype, genotype$names)
idy <- match(temp.meta$Compartment, compartment$names)

df <- temp.meta %>%
  add_column(
    genotype_short = genotype$short[idx], 
    compartment_short = compartment$short[idy]) %>% 
  mutate(exp_rep = paste(Experiment, Replicate, sep = "_"),
                      group = paste(genotype_short, compartment_short, sep = "_"))

# Group those groups that are available
idx <- which(df$SampleID %in% colnames(asv))

# df <- df[order(idx),]
# asv <- asv[,idx]
asv <- asv[,df$SampleID]

# Random factors
#exp_rep <- droplevels(factor(df$replicate, levels = 1:5))

# Grouping factors
levels <- c(sapply(genotype$short, paste, unique(paste(compartment$short, sep="_")), sep = "_"))
group <- factor(df$group, levels = levels[levels %in% unique(df$group)])
exp_rep <- as.factor(df$exp_rep)

# Model design
model <- model.matrix(~ 0 + group + exp_rep) # Dropping intercept

colnames(model) <- str_replace(colnames(model), "group", "")

idx <- colSums(model) == 0 # separating terms with no interaction
model <- model[,!idx]
set.seed(1)

# Contrast
contrasts.mat <- limma::makeContrasts(
    
    # compartment_Col_E = (Col_E - soil_BS),
    compartment_Col_RP_cas13 = (Col_RP - soil_BS),

    # compartment_nai1_E = (nai1_E - soil_BS),
    compartment_nai1_RP_cas13 = (nai1_RP - soil_BS),

    # compartment_pyk10_E = (pyk10_E - soil_BS),
    compartment_pyk10_RP_cas13 = (pyk10_RP - soil_BS),

    # compartment_cyp_E = (cyp_E - soil_BS),
    compartment_cyp_RP_cas13 = (cyp_RP - soil_BS),

    # compartment_myb_E = (myb_E - soil_BS),
    compartment_myb_RP_cas13 = (myb_RP - soil_BS),

    # genotype_nai1_E          = (  nai1_E   - Col_E ),
    # #nai2_E          = (  nai2_E   - Col_E ),
    # genotype_pyk10_E         = ( pyk10_E   - Col_E ),
    # genotype_cyp_E           = (   cyp_E   - Col_E ),
    # genotype_myb_E           = (   myb_E   - Col_E ),
    #bhlh_E = (   bhlh_E   - Col_E ),
    
    
    genotype_nai1_RP_cas13     = (  nai1_RP -   Col_RP),
    # nai2_RP     = (  nai2_RP -   Col_RP),
    genotype_pyk10_RP_cas13    = ( pyk10_RP -   Col_RP),
    genotype_cyp_RP_cas13      = (   cyp_RP -   Col_RP  ),
    genotype_myb_RP_cas13      = (   myb_RP -    Col_RP),
    # bhlh_RP = (   bhlh_RP   - Col_RP ),

    # nai1_RS     = (  nai1_RS -   Col_RS),
    # # nai2_RP     = (  nai2_RP -   Col_RP),
    # pyk10_RS    = ( pyk10_RS -   Col_RS),
    # cyp_RS     = (   cyp_RS -   Col_RS  ),
    # myb_RS      = (   myb_RS -    Col_RS),
    # bhlh_RP = (   bhlh_RP   - Col_RP ),
    
    
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

write.table(logFC_P, paste0(stats.out, "/CAS13_stats_logFC.P_combat_RP.txt"), sep="\t", 
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

write.table(count.mat, paste0(stats.out, "/CAS13_stats_number_of_enriched_depleted_combat_RP.txt"), 
            quote = F, row.names = T, col.names = NA, sep="\t")


# Significance table
DE <- temp

write.table(DE, paste0(stats.out, "/CAS13_stats_enrichment_test.significance_table_combat_RP.txt"), sep="\t", 
            quote = T, row.names = T, col.names = NA)

log2cpm <- cpm(de, prior.count = 2, log = T)

write.table(log2cpm, paste0(stats.out, "/CAS13_stats_log2cpm_combat_RP.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)


# Make object an RObject as output
CAS13_DE_asv_RP <- list(
  logFC_P = logFC_P, 
  DE_table = count.mat,
  cpm = log2cpm
)

save(list = "CAS13_DE_asv_RP", file = paste0("./data/CAS13_DE_asv_RP.RData"))

# !END of Script
sessionInfo()

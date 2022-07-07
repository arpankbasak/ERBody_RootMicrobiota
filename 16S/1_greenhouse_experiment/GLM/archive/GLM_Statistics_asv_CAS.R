#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "vegan", "RColorBrewer", "parallel", "DESeq2", "edgeR", "limma")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))

# What are we anlysing
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

temp.meta <- data.frame(meta$datasets[["CAS"]])
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
# idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,colnames(temp.mat)]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$Replicate))
temp.meta$Experiment <- factor(as.character(temp.meta$Experiment))
temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$name %in% as.character(temp.meta$Genotype))
idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

term <- c("CAS13")

# Make factors
temp.meta <- temp.meta %>%
filter(!Compartment %in% c("rhizosphere"), Soil_Batch %in% term) %>%
mutate(
    
    Genotype = factor(Genotype, levels = genotype$names[idg]),
    Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
    Compartment = factor(Compartment, levels = compartment$names[idc])
) %>% filter(!SampleID %in% outlier)

# Fetch taxa
# taxa.mat <- all_exp_cut$taxa
# taxa_level <- "lineage_family"
# taxa.mat_temp <- temp.mat %>% 
# add_column(tag = taxa.mat[match(row.names(.), taxa.mat$asvf), taxa_level], .before = 1) %>%
# gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
# group_by(tag, SampleID) %>%
# summarise(cRA = sum(RA)) %>%
# spread(key = SampleID, value = cRA, fill = 0, convert = FALSE) %>%
# data.frame(., stringsAsFactors = FALSE)

# row.names(taxa.mat_temp) <- taxa.mat_temp$tag

# Remove low abudant ASVs
asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv <- temp.mat[idx,]+1

# replace genotype and compartment names with simpler ones
idx <- match(temp.meta$Genotype, genotype$names)
idy <- match(temp.meta$Compartment, compartment$names)

levels <- c(sapply(genotype$short, paste, unique(paste(compartment$short, sep="_")), sep = "_"))

df <- temp.meta %>%
  add_column(
    genotype_short = genotype$short[idx], 
    compartment_short = compartment$short[idy]) %>% 
  mutate(
    # exp_rep = paste(Experiment, Replicate, Soil_Batch, Run, sep = "_"),
    exp_rep = paste(Experiment, Replicate, sep = "_"),
    group = paste(genotype_short, compartment_short, sep = "_")
    ) %>%
  mutate(
    group = factor(group, levels = levels[levels %in% unique(group)]),
    exp_rep = as.factor(exp_rep)
    )

# Group those groups that are available
idx <- which(df$SampleID %in% colnames(asv))

# df <- df[order(idx),]
# asv <- asv[,idx]
asv <- asv[,df$SampleID]
# asv.norm <- apply(asv, 2, function(x) log2(1+ (x/sum(x))))

# mod_fixed <- model.matrix(~group + exp_rep, df)
# combat_obj <- ComBat_seq(
#   counts = as.matrix(asv),
#   # covar_mod = mod_fixed, 
#   group = df$group,
#   batch = df$exp_rep,
#   full_mod = TRUE
# )


# Random factors
# exp_rep <- droplevels(factor(df$replicate, levels = 1:5))

# Grouping factors
group <- factor(df$group, levels = levels[levels %in% unique(df$group)])
exp_rep <- as.factor(df$exp_rep)

# limma::removeBatchEffect(assay(vsd), vsd$Batch)

# Model design
model <- model.matrix(~ 0 + group + exp_rep) # Dropping intercept
idx <- colSums(model) == 0 # separating terms with no interaction
model <- model[,!idx]
set.seed(1)
colnames(model) <- str_replace(colnames(model), "group", "")

# model_null <- model.matrix(~1, df)
# svseq <- svaseq(as.matrix(asv), model, model_null, n.sv=1)$sv

# Contrast
contrasts.mat <- limma::makeContrasts(
    
    compartment_Col_E_cas = (Col_E - soil_BS),
    compartment_Col_RP_cas = (Col_RP - soil_BS),

    compartment_nai1_E_cas = (nai1_E - soil_BS),
    compartment_nai1_RP_cas = (nai1_RP - soil_BS),

    compartment_pyk10_E_cas = (pyk10_E - soil_BS),
    compartment_pyk10_RP_cas = (pyk10_RP - soil_BS),

    compartment_cyp_E_cas = (cyp_E - soil_BS),
    compartment_cyp_RP_cas = (cyp_RP - soil_BS),

    compartment_myb_E_cas = (myb_E - soil_BS),
    compartment_myb_RP_cas = (myb_RP - soil_BS),

    genotype_nai1_E_cas          = (  nai1_E   - Col_E ),
    #nai2_E          = (  nai2_E   - Col_E ),
    genotype_pyk10_E_cas         = ( pyk10_E   - Col_E ),
    genotype_cyp_E_cas           = (   cyp_E   - Col_E ),
    genotype_myb_E_cas           = (   myb_E   - Col_E ),
    # bhlh_E = (   bhlh_E   - Col_E ),
    
    
    genotype_nai1_RP_cas     = (  nai1_RP -   Col_RP),
    # # nai2_RP     = (  nai2_RP -   Col_RP),
    genotype_pyk10_RP_cas    = ( pyk10_RP -   Col_RP),
    genotype_cyp_RP_cas      = (   cyp_RP -   Col_RP  ),
    genotype_myb_RP_cas      = (   myb_RP -    Col_RP),
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
delist <- DGEList(counts = asv, group = group, samples = df)
delist <- calcNormFactors(delist, method = "RLE")

# delist$samples %>% 
# ggplot(aes(x = norm.factors, fill = as.factor(group))) +
# geom_density(alpha = 0.5)

# Estimate common and tag-wise disperse
de <- estimateGLMCommonDisp(delist, model)
de <- estimateGLMTrendedDisp(de, model)
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

write.table(logFC_P, paste0(stats.out, "/stats_logFC.P_combat.txt"), sep="\t", 
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

write.table(count.mat, paste0(stats.out, "/stats_number_of_enriched_depleted_combat.txt"), 
            quote = F, row.names = T, col.names = NA, sep="\t")


# Significance table
DE <- temp

write.table(DE, paste0(stats.out, "/stats_enrichment_test.significance_table_combat.txt"), sep="\t", 
            quote = T, row.names = T, col.names = NA)

log2cpm <- cpm(de, prior.count = 2, log = T)

write.table(log2cpm, paste0(stats.out, "/stats_log2cpm_combat.txt"), sep="\t", 
            quote = F, col.names = NA, row.names = T)

# Temporary visualisation
col_data <- data.frame(id = colnames(logFC_P)) %>%
separate(id, remove = FALSE, into = c("hypothesis", "genotype", "compartment", "x", "value")) %>%
select(-x, -value) %>%
data.frame(., row.names = 1)

mat <- cbind(
  as.matrix(logFC_P[str_detect(colnames(logFC_P), "logFC$")]), 
  as.matrix(-log10(logFC_P[str_detect(colnames(logFC_P), "PValue$")]))
  )

lfc_mat <- as.matrix(logFC_P[str_detect(colnames(logFC_P), "logFC$") & str_detect(colnames(logFC_P), "^genotype")])
id_lfc <- (rowSums(abs(lfc_mat)) > 1)

pmap <- pheatmap::pheatmap(lfc_mat[id_lfc,], 
  annotation_col = col_data,
  color = colorRampPalette((RColorBrewer::brewer.pal(n = 3, name = "PiYG")))(5))

pdf(paste0(figs.out,"/GLM_heatmap_", term,".pdf"), width=5, height=12)
grid::grid.newpage()
grid::grid.draw(pmap$gtable)
dev.off()

GGally::ggpairs(as.data.frame(lfc_mat)) %>%
ggsave(
  filename = paste0(figs.out, "/GLM_comparison_pairs_", term,".png"),
  dpi = 72, 
  device = "png",
  bg = "white",
  units = "in", 
  width = 20, 
  height = 20, 
  limitsize = FALSE)


# Make object an RObject as output
CAS_DE_asv <- list(
  logFC_P = logFC_P, 
  DE_table = count.mat,
  cpm = log2cpm
)

save(list = "CAS_DE_asv", file = paste0("./data/asv_CAS.RData"))

# !END of Script
sessionInfo()

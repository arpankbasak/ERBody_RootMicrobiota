#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for detrending the data and obtaining the varition due to genotype
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "vegan", "RColorBrewer", "parallel")
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

plot_path <- paste0(figs, "/gh/beta_diversity/uspiked/cpcoa/")
stat_path <- paste0(stats, "/gh/beta_diversity/uspiked/cpcoa/")
out_path <- paste0(out, "/gh/beta_diversity/uspiked/cpcoa/")


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
outliers <- c("strex.75", "strex.25", "strex.49", "strex.27", "strex.48")
temp.meta <- (meta$datasets[["STREX"]]) %>% filter(!SampleID %in% outliers)
temp.meta <- data.frame(meta$design_strex[match(temp.meta$SampleID, meta$design_strex$SampleID),])
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
# temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$replicate))
temp.meta$tech_rep <- as.factor(as.character(temp.meta$tech_rep))
temp.meta$Experiment <- factor(as.character(temp.meta$experiment))
temp.meta$Treatment <- factor(as.character(temp.meta$Treatment), levels = c("Bulk soil", "Exudate", "Extract"))
# temp.meta$Run <- as.factor(as.character(temp.meta$Run))
asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
temp.mat <- temp.mat[idx,]

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$short %in% as.character(temp.meta$genotype))
# idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
# idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

temp.meta$Genotype <- factor(as.character(temp.meta$genotype), levels = genotype$short[idg])

# Compartment wise
pw_stat_list <- mclapply(c("Exudate", "Extract"), function(x){

	temp.df <- temp.meta %>% filter(Treatment == x)

	idg <- which(genotype$short %in% as.character(temp.df$genotype))
    # idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))
    temp.df$Genotype <- factor(as.character(temp.df$genotype), levels = genotype$short[idg])

    row.names(temp.df) <- temp.df$SampleID
    id <- as.character(temp.df$SampleID)
    asv_cpcoa <- apply(temp.mat[, id], 2, function(x) x/sum(x))

	gen <- unique(as.character(temp.df$Genotype))[unique(as.character(temp.df$Genotype)) != "Col"]

      stat_list <- lapply(gen, function(i){

  		ix <- i
        temp.metadata <- temp.meta %>%
          filter(Treatment == x, Genotype %in% c("Col", ix))

        temp.metadata <- temp.metadata %>%
	    mutate(
	    	group = Genotype,
	    	random = as.factor(paste(Replicate, Experiment, tech_rep, sep = "_")))
	    
	    # Overall beta-diversity
	    row.names(temp.metadata) <- temp.metadata$SampleID
      id <- as.character(temp.metadata$SampleID)
	    asv_cpcoa.ra <- asv_cpcoa[, id]

	      f <- formula(t(asv_cpcoa.ra) ~ group + random)

	      ad_obj <- adonis(f,
	      sqrt.dist = TRUE, 
	      method = "bray",
	      data = temp.metadata %>% select(group, random),
	      strata = temp.metadata$group,
	      permutations = 1000)

	      stat_sum <- broom::tidy(ad_obj$aov.tab)
	      
	      stat <- stat_sum %>% 
	      dplyr::filter(term == "group") %>%
	        mutate(
	          con.var.explained = .$SumsOfSqs/stat_sum$SumsOfSqs[which(stat_sum$term == "Total")],
	          Genotype = ix,
	          term = (ix),
	          Treatment = x) %>%
	        data.frame(., stringsAsFactors = FALSE)

          return(stat)

      }) 

    stat_list <- do.call(rbind.data.frame, stat_list)
    return(stat_list)

}, mc.cores = 4)

# Make data frame
stat_df <- do.call(rbind.data.frame, pw_stat_list)

write.table(stat_df, file= paste0(stats.out, "/pw_PERMANOVA.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# stat_df$effect <- "Genotype"
# stat_df$effect[str_detect(stat_df$term, "soilbatch")] <- "Genotype*Soilbatch"
# stat_df$effect[str_detect(stat_df$term, "CAS")] <- "Genotype*Soilbatch:CAS"

p <- stat_df %>%
mutate(
	FDR = p.adjust(p.value, "bonferroni"), 
  con.var.explained = con.var.explained*100) %>%
mutate(
  significance = p.value <= alpha,
  Genotype = factor(Genotype, levels = genotype$short[genotype$short %in% .$Genotype]),
  Treatment = factor(Treatment, levels = c("Exudate", "Extract")),
  tag = paste0(as.character(round(con.var.explained, 2)), "%,\n", as.character(round(p.value, 4)))
  ) %>%
ggplot(aes(x = Genotype, y = con.var.explained)) +
geom_bar(stat = "identity", aes(colour = significance, fill = Genotype, alpha = significance)) +
geom_text(aes(label = tag, alpha = significance), nudge_y = 1, size = 1) +
facet_grid(Treatment ~., space = "free", switch = "y", scale = "free_y") +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.8, `TRUE` = 1), guide = FALSE) +
scale_fill_manual(values = genotype$colour[which(genotype$short %in% stat_df$Genotype)], guide = FALSE) +
coord_flip() +
theme_RTN +
theme() +
labs(y = "Variance explained %", x = "", colour = "", fill = "pvalue", alpha = "")
ggsave(p, file = paste0(figs.out, "ST_pwCPCOA.png"), 
              dpi = 600, device = "png", units = "in", 
              bg = "transparent", 
              width = 3, height = 5, limitsize = F)



# END OF SCRIPT
sessionInfo()
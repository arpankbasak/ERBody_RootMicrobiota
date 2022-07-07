#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for Overall Beta diversity measures
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.itsf)

# Loading required packages
pkgs <- c("tidyverse", "vegan", "RColorBrewer", "parallel")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))
# load(paste0("./data/soilbatch_CAS_DE_asv.RData"))
# common_in_CAS <- soilbatch_DE_asv$asv_list$similar_asvs

# What are we anlysin
# pro <- names(meta$datasets)

figs.out <- paste0(figs, "/st/beta_diversity/uspiked/pcoa/")
stats.out <- paste0(stats, "/st/beta_diversity/uspiked/pcoa/")
out.out <- paste0(out, "/st/beta_diversity/uspiked/pcoa/")

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

p = "st"

plot_path <- paste0(figs, "/st/beta_diversity/uspiked/pcoa/")
stat_path <- paste0(stats, "/st/beta_diversity/uspiked/pcoa/")
out_path <- paste0(out, "/st/beta_diversity/uspiked/pcoa/")


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

temp.meta <- (meta$datasets[["STREX"]])
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

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$short %in% as.character(temp.meta$genotype))
# idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
idc <- which(treatment$name %in% as.character(temp.meta$Treatment))

temp.meta$Genotype <- factor(as.character(temp.meta$genotype), levels = genotype$short[idg])

# Overall beta-diversity
id <- colnames(temp.mat) %in% as.character(temp.meta$SampleID)
asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv_pcoa <- temp.mat[,id]

temp.mds <- beta_diversity(metadata = temp.meta, asv_matrix = asv_pcoa)


# MDS 12
(plot_obj_1 <-  temp.mds$mds_df %>%
ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
ggtitle(paste0("All")) +
geom_point(alpha = 0.6, size = 3, 
    aes(colour = Treatment, 
        fill = Genotype, 
        shape = Experiment), # Soiltype
    show.legend = T) +
scale_fill_manual(values = genotype$colours[idg]) + 
scale_colour_manual(values = treatment$colours[idc]) + 
scale_shape_manual(values = c(22, 23, 24, 25)) +
theme_RTN_MDS +
labs(x = paste0("PCoA - 1: ",temp.mds$variance[1],"%"), 
     y = paste0("PCoA - 2: ",temp.mds$variance[2],"%"))) +
ggsave(filename = paste0(plot_path, "/overall_beta_diversity_all.png"), 
           dpi = 600, 
           device = "png", 
           bg = "transparent", 
           units = img, 
           width = cpcoa.box[1]+2.5, 
           height = cpcoa.box[2]+2.5,
           limitsize = FALSE)

# MDS 23
(plot_obj_2 <-  temp.mds$mds_df %>%
ggplot(aes(x = saturate(MDS3), y = saturate(MDS2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
ggtitle(paste0("All")) +
geom_point(alpha = 0.6, size = 3, 
    aes(colour = Treatment, 
        fill = Genotype, 
        shape = Experiment), # Soiltype
    show.legend = T) +
scale_fill_manual(values = genotype$colours[idg]) + 
scale_colour_manual(values = treatment$colours[idc]) + 
scale_shape_manual(values = c(22, 23, 24, 25)) +
theme_RTN_MDS +
labs(x = paste0("PCoA - 3: ",temp.mds$variance[3],"%"), 
     y = paste0("PCoA - 2: ",temp.mds$variance[2],"%"))) +
ggsave(filename = paste0(plot_path, "/overall_beta_diversity_all_23.png"), 
           dpi = 600, 
           device = "png", 
           bg = "transparent", 
           units = img, 
           width = cpcoa.box[1]+2.5, 
           height = cpcoa.box[2]+2.5, 
           limitsize = FALSE)

# Treatment wise
mclapply(c("Bulk soil", "Exudate", "Extract"), function(x){

        temp.metadata <- temp.meta %>%
        filter(Treatment == x)
        
        idg <- which(genotype$short %in% as.character(temp.metadata$genotype))
        # idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))
        temp.metadata$Genotype <- factor(as.character(temp.metadata$genotype), levels = genotype$short[idg])
        

        # Overall beta-diversity
        id <- colnames(temp.mat) %in% as.character(temp.metadata$SampleID)
        asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
        idx <- rowSums(asv.ra * 100 > threshold) >= 1
        asv_pcoa <- temp.mat[,id]

        temp.mds <- beta_diversity(metadata = temp.metadata, asv_matrix = asv_pcoa)


        # MDS 12
        temp.mds$mds_df %>%
        ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
        geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        ggtitle(paste0(x)) +
        geom_point(alpha = 0.6, size = 3, colour = c_black,
        aes( 
            fill = Genotype, 
            shape = Replicate), # Soiltype
        show.legend = T) +
        scale_fill_manual(values = genotype$colours[idg]) + 
        # scale_colour_manual(values = treatment$colours[idc]) + 
        scale_shape_manual(values = c(22, 23, 24, 25)) +
        theme_RTN_MDS +
        labs(x = paste0("PCoA - 1: ",temp.mds$variance[1],"%"), 
             y = paste0("PCoA - 2: ",temp.mds$variance[2],"%")) +
        ggsave(filename = paste0(plot_path, "/overall_beta_diversity", x,"_treatment_12.png"), 
                   dpi = 600, 
                   device = "png", 
                   bg = "transparent", 
                   units = img, 
                   width = cpcoa.box[1], 
                   height = cpcoa.box[2],
                   limitsize = FALSE)

        # MDS 23
        temp.mds$mds_df %>%
        ggplot(aes(x = saturate(MDS3), y = saturate(MDS2))) +
        geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        ggtitle(paste0(x)) +
        geom_point(alpha = 0.6, size = 3, colour = c_black,
        aes( 
            fill = Genotype, 
            shape = Replicate), # Soiltype
        show.legend = T) +
        scale_fill_manual(values = genotype$colours[idg]) + 
        # scale_colour_manual(values = treatment$colours[idc]) + 
        scale_shape_manual(values = c(22, 23, 24, 25)) +
        theme_RTN_MDS +
        labs(x = paste0("PCoA - 3: ",temp.mds$variance[3],"%"), 
             y = paste0("PCoA - 2: ",temp.mds$variance[2],"%")) +
        ggsave(filename = paste0(plot_path, "/overall_beta_diversity", x,"_treatment_23.png"), 
                   dpi = 600, 
                   device = "png", 
                   bg = "transparent", 
                   units = img, 
                   width = cpcoa.box[1], 
                   height = cpcoa.box[2], 
                   limitsize = FALSE)


}, mc.cores = 4)


# END OF SCRIPT
sessionInfo()
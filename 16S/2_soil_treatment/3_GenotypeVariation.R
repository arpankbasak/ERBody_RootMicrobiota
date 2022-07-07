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

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$short %in% as.character(temp.meta$genotype))
# idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
# idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

temp.meta$Genotype <- factor(as.character(temp.meta$genotype), levels = genotype$short[idg])
asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv_cpcoa <- apply(temp.mat[idx,], 2, function(x) x/sum(x))

levs <- sapply(genotype$short[idg], function(x) paste(x, as.character(unique(temp.meta$Replicate)), sep = "_")) %>% c(.)

# Compartment wise
lapply(c("Exudate", "Extract"), function(x){

        temp.metadata <- temp.meta %>%
        filter(Treatment == x)
        
        idg <- which(genotype$short %in% as.character(temp.metadata$genotype))
        # idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))
        temp.metadata$Genotype <- factor(as.character(temp.metadata$genotype), levels = genotype$short[idg])
        # temp.metadata$group <- factor(paste(as.character(temp.metadata$genotype), as.character(temp.metadata$Replicate), sep = "_"),
        #     levels = levs)
        # temp.metadata$group <- droplevels(temp.metadata$group)
        

        # Overall beta-diversity
        row.names(temp.metadata) <- as.character(temp.metadata$SampleID)
        id <- as.character(temp.metadata$SampleID)
        # asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
        # idx <- rowSums(asv.ra * 100 > threshold) >= 1
        asv_cpcoa <- apply(temp.mat[idx,id], 2, function(x) x/sum(x))

        f <- formula(t(asv_cpcoa) ~ Genotype + Condition(Experiment + Replicate + tech_rep))
        # cpcoa_obj <- constrained_beta_diversity(
        #   metadata = temp.metadata, 
        #   asv_matrix = asv_cpcoa, 
        #   formula = f)

        cca_obj <- capscale(formula = f, 
                        data = temp.metadata %>% select(Genotype, Experiment, Replicate, tech_rep), 
                        distance = "bray", 
                        sqrt.dist = T, 
                        add = F)
    
        # Conduct PERMANOVA
        set.seed(1)
        p.nova <- anova.cca(cca_obj)

        # Extract variables from stat objects
        eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
        pval <- p.nova$`Pr(>F)`[1]
        chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
        variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                               row.names = c("total", "constrianed", "unconstrained"))

        # Format P-Value for representation
        variance <- format(100 * variable$proportion[2], digits = 4)
        tf <- paste(as.character(f)[2], as.character(f)[1], as.character(f)[3])
        ti <- paste0("[", tf ,"] :  "," = ",variance,"% p = ",pval, 
                     ifelse(pval < 0.05, "*", " (NS)"))
       
        idx <- match(temp.metadata$SampleID, row.names(cca_obj$CCA$wa))
        cpcoa_df <- cbind.data.frame(temp.metadata,
            CCA1 = (as.numeric(cca_obj$CCA$wa[idx,1])),
            CCA2 = (as.numeric(cca_obj$CCA$wa[idx,2]))
        )

        # Use the biplot
        markers <- as.data.frame(cca_obj$CCA$centroids) %>%
        mutate(group = str_replace_all(row.names(.), "Genotype", "")) %>%
        # separate(group, into = c("A", "B"), sep = ":", convert = FALSE, remove = FALSE) %>%
        # mutate(A = factor(A)) %>%
        select(group, CAP1, CAP2)

        # CPCOA 12
        p <- cpcoa_df %>%
        # mutate(direction = (sign(CCA3) == 1)) %>%
        cbind.data.frame(., markers[match(as.character(.$Genotype), as.character(markers$group)),c("CAP1", "CAP2")]) %>%
        mutate(CAP1 = replace_na(CAP1, 0), CAP2 = replace_na(CAP2, 0)) %>%
        ggplot(aes(x = saturate(CCA1), y = saturate(CCA2))) +
        geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        ggtitle(paste0(ti)) +
        geom_segment(aes(x = saturate(CAP1), y = saturate(CAP2),
          xend = saturate(CCA1), yend = saturate(CCA2), 
          colour = Genotype), 
        lwd = 0.3, lty = "solid") +
        scale_color_manual(values = genotype$colours[idg]) + 
        geom_point(alpha = 0.8, size = 2, colour = c_black,
            aes(fill = Genotype, 
                shape = Replicate), # Soiltype
            show.legend = T) +
        scale_fill_manual(values = genotype$colours[idg]) + 
        scale_shape_manual(values = c(22, 23, 24, 25)) +
        # scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = FALSE) +
        theme_RTN_MDS +
        labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
             y = paste0("cPCoA - 2: ",eig[2],"%"))
        ggsave(p, filename = paste0(plot_path, "/cpcoa_", x,"_treatment_12.png"), 
                   dpi = 600, 
                   device = "png", 
                   bg = "transparent", 
                   units = img, 
                   width = cpcoa.box[1], 
                   height = cpcoa.box[2],
                   limitsize = FALSE)


})


# END OF SCRIPT
sessionInfo()
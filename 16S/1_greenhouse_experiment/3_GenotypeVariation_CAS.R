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

p = "gh"

plot_path <- paste0(figs, "/gh/beta_diversity/uspiked/cpcoa/asv")
stat_path <- paste0(stats, "/gh/beta_diversity/uspiked/cpcoa/asv")
out_path <- paste0(out, "/gh/beta_diversity/uspiked/cpcoa/asv")

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
temp.meta <- rbind.data.frame(meta$datasets[["CAS"]]) %>% filter(!SampleID %in% outlier)
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
mutate(
    
    Genotype = factor(Genotype, levels = genotype$names[idg]),
    Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
    Compartment = factor(Compartment, levels = compartment$names[idc])
) %>% 
filter(!SampleID %in% outlier) %>%
data.frame(., row.names = .$SampleID)


# Compartment wise
mclapply(c("rhizoplane", "endosphere"), function(x){

        temp.metadata <- temp.meta %>%
        filter(Compartment == x)
        
        idg <- which(genotype$name %in% as.character(temp.metadata$Genotype))
        idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))

        temp.metadata <- temp.metadata %>%
        mutate(
            
            Genotype = factor(Genotype, levels = genotype$names[idg]),
            Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
            # group = as.factor(paste(Genotype, Soil_Batch, sep = "_")),
            group = factor(Genotype, levels = genotype$names[idg])
        )

        # Overall beta-diversity
        id <- as.character(temp.metadata$SampleID)
        asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
        idx <- rowSums(asv.ra * 100 > threshold) >= 1
        asv_cpcoa <- apply(temp.mat[idx,id], 2, function(x) x/sum(x))

        f <- formula(t(asv_cpcoa) ~ group + Condition(Soil_Batch + Experiment + Replicate + Run))
        # cpcoa_obj <- constrained_beta_diversity(
        #   metadata = temp.metadata, 
        #   asv_matrix = asv_cpcoa, 
        #   formula = f)

        cca_obj <- capscale(
                        formula = f, 
                        data = (temp.metadata %>% select(group, Soil_Batch, Experiment, Replicate, Run)), 
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
            CCA2 = (as.numeric(cca_obj$CCA$wa[idx,2])),
            CCA3 = (as.numeric(cca_obj$CCA$wa[idx,3]))
        )

        # Use the biplot
        markers <- as.data.frame(cca_obj$CCA$centroids) %>%
        mutate(group = str_replace_all(row.names(.), "group", "")) %>%
        # separate(group, into = c("A", "B"), sep = ":", convert = FALSE, remove = FALSE) %>%
        # mutate(A = factor(A)) %>%
        select(group, CAP1, CAP2, CAP3)

        # Make ASVs that are explaining the cumulative varance >1 at family level
        ftable <- as.data.frame(cca_obj$CCA$v) %>%
        add_column(asvb = row.names(.), .before = 1) %>%
        separate(asvb, into = c("Kingdom", "Phylum", "Class", "Order", "Family", sep = ";")) %>%
        select(Family, CAP1, CAP2, CAP3)

        taxa.mat.temp <- ftable %>% 
        # separate(asvb, into = c("Kingdom", "Phylum", "Class", "Order", "Family", sep = ";")) %>%
        mutate(Family = str_replace_na(Family, "Unclassified")) %>%
        group_by(Family) %>%
        summarise(cCAP1 = sum(CAP1), cCAP2 = sum(CAP2), cCAP3 = sum(CAP3)) %>%
        filter(
          cCAP1 >= quantile((cCAP1), .95) | cCAP1 <= quantile((cCAP1), .05), 
          cCAP2 >= quantile((cCAP2), .95) | cCAP2 <= quantile((cCAP2), .05), 
          cCAP3 >= quantile((cCAP3), .95) | cCAP3 <= quantile((cCAP3), .05)
          ) %>%
         mutate(
              cCAP1 = ifelse(cCAP1 > max(cpcoa_df$CCA1), max(cpcoa_df$CCA1)+0.001, cCAP1),
              cCAP2 = ifelse(cCAP2 > max(cpcoa_df$CCA2), max(cpcoa_df$CCA2)+0.001, cCAP2),
              cCAP3 = ifelse(cCAP3 > max(cpcoa_df$CCA3), max(cpcoa_df$CCA3)+0.001, cCAP3)
            ) %>%
          mutate(
              cCAP1 = ifelse(cCAP1 < min(cpcoa_df$CCA1), min(cpcoa_df$CCA1)-0.001, cCAP1),
              cCAP2 = ifelse(cCAP2 < min(cpcoa_df$CCA2), min(cpcoa_df$CCA2)-0.001, cCAP2),
              cCAP3 = ifelse(cCAP3 < min(cpcoa_df$CCA3), min(cpcoa_df$CCA3)-0.001, cCAP3)
            ) %>%
        data.frame(., stringsAsFactors = FALSE)

        # CPCOA 12
        p <- cpcoa_df %>%
        mutate(direction = (sign(CCA3) == 1)) %>%
        cbind.data.frame(., markers[match(as.character(.$group), as.character(markers$group)),c("CAP1", "CAP2", "CAP3")]) %>%
        ggplot(aes(x = (CCA1), y = (CCA2))) +
        # geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        # geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        ggtitle(paste0(ti)) +
        geom_segment(aes(x = (CAP1), y = (CAP2),
          xend = (CCA1), yend = (CCA2), 
          colour = Genotype, 
          lty  = direction), 
        lwd = 0.3) +
        # geom_text(data = taxa.mat.temp, aes(x = saturate(cCAP1), y = saturate(cCAP2),  position = "dodge", label = Family), colour = c_black, size = 1) +
        geom_point(alpha = 0.8, size = 2, colour = c_black,
            aes(fill = Genotype, 
                shape = Soil_Batch), # Soiltype
            show.legend = T) +
        scale_color_manual(values = genotype$colours[idg]) + 
        scale_fill_manual(values = genotype$colours[idg]) + 
        scale_shape_manual(values = c(22, 23, 24, 25)) +
        scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = FALSE) +
        theme_RTN_MDS +
        labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
             y = paste0("cPCoA - 2: ",eig[2],"%"))
        ggsave(p, filename = paste0(plot_path, "/cpcoa_", x,"_compartemnt_12.png"), 
                   dpi = 600, 
                   device = "png", 
                   bg = "transparent", 
                   units = img, 
                   width = cpcoa.box[1], 
                   height = cpcoa.box[2],
                   limitsize = FALSE)

        # CPCOA 23
        p <- cpcoa_df %>%
        mutate(direction = (sign(CCA1) == 1)) %>%
        cbind.data.frame(., markers[match(as.character(.$group), as.character(markers$group)), c("CAP1", "CAP2", "CAP3")]) %>%
        ggplot(aes(x = (CCA3), y = (CCA2))) +
        # geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        # geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
        ggtitle(paste0(ti)) +
        geom_segment(aes(x = (CAP3), y = (CAP2), xend = (CCA3), yend = (CCA2), 
          colour = Genotype, 
          lty = direction), 
        lwd = 0.3) +
        # geom_text(data = taxa.mat.temp, aes(x = saturate(cCAP1), y = saturate(cCAP2),  position = "dodge", label = Family), colour = c_black, size = 1) +
        geom_point(alpha = 0.8, size = 2, colour = c_black,
            aes(fill = Genotype, 
                shape = Soil_Batch), # Soiltype
            show.legend = T) +
        scale_color_manual(values = genotype$colours[idg]) + 
        scale_fill_manual(values = genotype$colours[idg]) + 
        scale_shape_manual(values = c(22, 23, 24, 25)) +
        scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = FALSE) +
        theme_RTN_MDS +
        labs(x = paste0("cPCoA - 3: ",eig[3],"%"), 
             y = paste0("cPCoA - 2: ",eig[2],"%"))
        ggsave(p, filename = paste0(plot_path, "/cpcoa_", x,"_compartemnt_23.png"), 
                   dpi = 600, 
                   device = "png", 
                   bg = "transparent", 
                   units = img, 
                   width = cpcoa.box[1], 
                   height = cpcoa.box[2], 
                   limitsize = FALSE)

        if(x != "bulk soil"){  
          # Subset the variation within the CAS soil
          temp.metadata <- temp.meta %>%
          filter(Compartment == x, soil == "CAS")
          
          idg <- which(genotype$name %in% as.character(temp.metadata$Genotype))
          idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))

          temp.metadata <- temp.metadata %>%
          mutate(
              
              Genotype = factor(Genotype, levels = genotype$names[idg]),
              Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb])
          )

          # Overall beta-diversity
          id <- as.character(temp.metadata$SampleID)
          asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
          idx <- rowSums(asv.ra * 100 > threshold) >= 1
          asv_cpcoa <- apply(temp.mat[idx,id], 2, function(x) x/sum(x))

          f <- formula(t(asv_cpcoa) ~ Genotype + Condition(Soil_Batch + Experiment + Replicate + Run))
          # cpcoa_obj <- constrained_beta_diversity(
          #   metadata = temp.metadata, 
          #   asv_matrix = asv_cpcoa, 
          #   formula = f)

          cca_obj <- capscale(formula = f, 
                          data = (temp.metadata %>% select(Genotype, Soil_Batch, Experiment, Replicate, Run)), 
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
              CCA2 = (as.numeric(cca_obj$CCA$wa[idx,2])),
              CCA3 = (as.numeric(cca_obj$CCA$wa[idx,3]))
          )

          # Use the biplot
          markers <- as.data.frame(cca_obj$CCA$centroids) %>%
          mutate(group = factor(str_replace_all(row.names(.), "Genotype|Soil_Batch", ""), levels = genotype$name[idg])) %>%
          select(group, CAP1, CAP2, CAP3)

          # Make ASVs that are explaining the cumulative varance >1 at family level
          ftable <- as.data.frame(cca_obj$CCA$v) %>%
          add_column(asvb = row.names(.), .before = 1) %>%
          separate(asvb, into = c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";") %>%
          select(Family, CAP1, CAP2, CAP3)

          taxa.mat.temp <- ftable %>% 
          mutate(Family = str_replace_na(Family, "Unclassified")) %>%
          group_by(Family) %>%
          summarise(cCAP1 = sum(CAP1), cCAP2 = sum(CAP2), cCAP3 = sum(CAP3)) %>%
          filter(
            cCAP1 >= quantile((cCAP1), .95) | cCAP1 <= quantile((cCAP1), .05), 
            cCAP2 >= quantile((cCAP2), .95) | cCAP2 <= quantile((cCAP2), .05), 
            cCAP3 >= quantile((cCAP3), .95) | cCAP3 <= quantile((cCAP3), .05)
            ) %>%
           mutate(
                cCAP1 = ifelse(cCAP1 > max(cpcoa_df$CCA1), max(cpcoa_df$CCA1)+0.001, cCAP1),
                cCAP2 = ifelse(cCAP2 > max(cpcoa_df$CCA2), max(cpcoa_df$CCA2)+0.001, cCAP2),
                cCAP3 = ifelse(cCAP3 > max(cpcoa_df$CCA3), max(cpcoa_df$CCA3)+0.001, cCAP3)
              ) %>%
            mutate(
                cCAP1 = ifelse(cCAP1 < min(cpcoa_df$CCA1), min(cpcoa_df$CCA1)-0.001, cCAP1),
                cCAP2 = ifelse(cCAP2 < min(cpcoa_df$CCA2), min(cpcoa_df$CCA2)-0.001, cCAP2),
                cCAP3 = ifelse(cCAP3 < min(cpcoa_df$CCA3), min(cpcoa_df$CCA3)-0.001, cCAP3)
              ) %>%
          data.frame(., stringsAsFactors = FALSE)

          # CPCOA 12
          p <- cpcoa_df %>%
          cbind.data.frame(., markers[match(.$Genotype, as.character(markers$group)),]) %>%
          mutate(direction = (sign(CCA3)==1)) %>%
          ggplot(aes(x = (CCA1), y = (CCA2))) +
          # geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
          # geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
          ggtitle(paste0(ti)) +
          # geom_point(alpha = 0.8, size = 5, colour = c_black, shape = 22,
          #     aes(fill = Genotype), # Soiltype
          #     show.legend = T) +
          geom_segment(aes(x = CAP1, y = CAP2, xend = (CCA1), yend = (CCA2), colour = Genotype), lty = "solid", lwd = 0.5) +
          # geom_text(data = taxa.mat.temp, aes(x = saturate(cCAP1), y = saturate(cCAP2),  position = "dodge", label = Family), colour = c_black, size = 1) +
          geom_point(alpha = 0.8, size = 3, colour = c_black, 
              aes(fill = Genotype, 
                  shape = Experiment), # Soiltype
              show.legend = T) +
          scale_color_manual(values = genotype$colours[idg]) + 
          scale_fill_manual(values = genotype$colours[idg]) + 
          scale_shape_manual(values = c(22, 23, 24, 25)) +
          scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = FALSE) +
          theme_RTN_MDS +
          labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
               y = paste0("cPCoA - 2: ",eig[2],"%"))
          ggsave(p, filename = paste0(plot_path, "/cpcoa_CAS_", x,"_compartemnt_12.png"), 
                     dpi = 600, 
                     device = "png", 
                     bg = "transparent", 
                     units = img, 
                     width = cpcoa.box[1], 
                     height = cpcoa.box[2],
                     limitsize = FALSE)

          # CPCOA 23
          p <- cpcoa_df %>%
          cbind.data.frame(., markers[match(.$Genotype, as.character(markers$group)),]) %>%
          mutate(direction = (sign(CCA1)==1)) %>%
          ggplot(aes(x = (CCA3), y = (CCA2))) +
          # geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
          # geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
          ggtitle(paste0(ti)) +
          # geom_point(alpha = 0.8, size = 5, colour = c_black, shape = 22,
          #     aes(fill = Genotype), # Soiltype
          #     show.legend = T) +
          geom_segment(aes(x = CAP3, y = CAP2, xend = (CCA3), yend = (CCA2), colour = Genotype), lty = "solid", lwd = 0.5) +
          # geom_text(data = taxa.mat.temp, aes(x = saturate(cCAP1), y = saturate(cCAP2),  position = "dodge", label = Family), colour = c_black, size = 1) +
          geom_point(alpha = 0.8, size = 3, colour = c_black, 
              aes(fill = Genotype, 
                  shape = Experiment), # Soiltype
              show.legend = T) +
          scale_colour_manual(values = genotype$colours[idg]) + 
          scale_fill_manual(values = genotype$colours[idg]) + 
          scale_shape_manual(values = c(22, 23, 24, 25)) +
          scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = FALSE) +
          theme_RTN_MDS +
          labs(x = paste0("cPCoA - 3: ",eig[3],"%"), 
               y = paste0("cPCoA - 2: ",eig[2],"%"))
          ggsave(p, filename = paste0(plot_path, "/cpcoa_CAS_", x,"_compartemnt_23.png"), 
                     dpi = 600, 
                     device = "png", 
                     bg = "transparent", 
                     units = img, 
                     width = cpcoa.box[1], 
                     height = cpcoa.box[2], 
                     limitsize = FALSE)
        }

}, mc.cores = 4)


# END OF SCRIPT
sessionInfo()
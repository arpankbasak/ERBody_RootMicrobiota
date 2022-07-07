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

p = "gh"

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
temp.meta <- rbind.data.frame(meta$datasets[["CAS"]]) #%>% filter(!SampleID %in% outlier)
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
pw_stat_list <- mclapply(c("rhizoplane", "endosphere"), function(x){

    stat_list <- lapply(c("CAS11","CAS13"), function(k){

      temp.metadata <- temp.meta %>%
          filter(Compartment == x, Soil_Batch == k)

      # Overall beta-diversity
      id <- as.character(temp.metadata$SampleID)
      asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
      idx <- rowSums(asv.ra * 100 > threshold) >= 1
      asv_cpcoa.compartment <- apply(temp.mat[idx,id], 2, function(x) x/sum(x))
  
      stat_list <- lapply(c("nai1", "pyk10", "cyp", "myb"), function(i){


          ix <- genotype$names[genotype$short == i]
          temp.metadata <- temp.meta %>%
          filter(Compartment == x, Genotype %in% c("Col-0", ix), Soil_Batch == k)
          
          # idg <- which(genotype$name %in% as.character(temp.metadata$Genotype))
          # idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))

          # temp.metadata <- temp.metadata %>%
          # mutate(
              
          #     Genotype = factor(Genotype, levels = genotype$names[idg]),
          #     Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
          #     group = as.factor(paste(Genotype, Soil_Batch, sep = "_")),
          #     random = as.factor(paste(Experiment, Replicate, Run, sep = "_"))
          # )

          # # Overall beta-diversity
          # id <- colnames(temp.mat) %in% as.character(temp.metadata$SampleID)
          # asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
          # idx <- rowSums(asv.ra * 100 > threshold) >= 1
          # asv_cpcoa <- temp.mat[,id]

          # f <- formula(t(asv_cpcoa) ~ group + random)

          # ad_obj <- adonis2(f,
          # sqrt.dist = TRUE, 
          # method = "bray",
          # data = temp.metadata,
          # strata = temp.metadata$group,
          # permutations = 1000)

          # stat_sum <- broom::tidy(ad_obj)
          
          # stat <- stat_sum %>% 
          # dplyr::filter(term == "group") %>%
          #   mutate(
          #     con.var.explained = .$SumOfSqs/stat_sum$SumOfSqs[which(stat_sum$term == "Total")],
          #     Genotype = ix,
          #     term = paste(ix, "soilbatch", sep = ":")) %>%
          #   data.frame(., stringsAsFactors = FALSE)

          # # Soilbatch as random
          # ix <- genotype$names[genotype$short == i]
          # temp.metadata <- temp.meta %>%
          # filter(Compartment == x, Genotype %in% c("Col-0", ix))
          
          idg <- which(genotype$name %in% as.character(temp.metadata$Genotype))
          idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))

          temp.metadata <- temp.metadata %>%
          mutate(
              
              Genotype = factor(Genotype, levels = genotype$names[idg]),
              Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
              group = as.factor(paste(Genotype, Soil_Batch, sep = "_")),
              random = as.factor(paste(Experiment, Replicate, Run, sep = "_"))
          )

          set.seed(seeder)
          
          asv_cpcoa <- asv_cpcoa.compartment[,temp.metadata$SampleID]

          f <- formula(t(asv_cpcoa) ~ Genotype + random)
          perm <- how(nperm = 1000)
          # setBlocks(perm) <- with(temp.metadata, random)

          # ad_obj <- adonis2(f,
          # sqrt.dist = TRUE, 
          # by = "term",
          # method = "bray",
          # data = temp.metadata %>% select(Genotype, random),
          # # strata = temp.metadata$Genotype,
          # perm = perm
          # )

          ad_obj <- adonis(f,
          sqrt.dist = TRUE, 
          by = "term",
          method = "bray",
          data = temp.metadata %>% select(Genotype, random),
          # strata = temp.metadata$Genotype,
          permutations = 1000
          )

          # stat_sum <- broom::tidy(ad_obj)
          stat_sum <- broom::tidy(ad_obj$aov.tab)
          
          stat_temp <- stat_sum %>% 
          dplyr::filter(term == "Genotype") %>%
            mutate(
              con.var.explained = .$SumsOfSqs/stat_sum$SumsOfSqs[which(stat_sum$term == "Total")],
              Genotype = ix, Compartment = x, Soil_Batch = k,
              term = paste(ix, x, k, sep = "_")) %>%
            data.frame(., stringsAsFactors = FALSE)

          stat <- stat_temp

          # Within CAS soil
          # if(x != "bulk soil"){  
              
          #     # Subset the variation within the CAS soil
          #     temp.metadata <- temp.meta %>%
          #     filter(Compartment == x, Genotype %in% c("Col-0", ix), soil == "CAS", Soil_Batch == k)
              
          #     idg <- which(genotype$name %in% as.character(temp.metadata$Genotype))
          #     idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))

          #     temp.metadata <- temp.metadata %>%
          #     mutate(
                
          #       Genotype = factor(Genotype, levels = genotype$names[idg]),
          #       Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[idsb]),
          #       group = as.factor(paste(Genotype, Soil_Batch, sep = "_")),
          #       random = as.factor(paste(Experiment, Replicate, Run, sep = "_"))
          #     )

          #     set.seed(seeder)
          #     # Overall beta-diversity
          #     id <- colnames(temp.mat) %in% as.character(temp.metadata$SampleID)
          #     asv.ra <- apply(temp.mat[,id], 2, function(x) x/sum(x))
          #     idx <- rowSums(asv.ra * 100 > threshold) >= 1
          #     asv_cpcoa <- temp.mat[,id]

          #     f <- formula(t(asv_cpcoa) ~ group + random)

          #     ad_obj <- adonis2(f,
          #     sqrt.dist = TRUE, 
          #     method = "bray",
          #     data = temp.metadata,
          #     strata = temp.metadata$group,
          #     permutations = 1000)

          #     stat_sum <- broom::tidy(ad_obj)
              
          #     stat_temp <- stat_sum %>% 
          #     dplyr::filter(term == "group") %>%
          #       mutate(
          #         con.var.explained = .$SumOfSqs/stat_sum$SumOfSqs[which(stat_sum$term == "Total")],
          #         Genotype = ix,
          #         term = paste(ix, k, "_CAS", sep = ":")) %>%
          #       data.frame(., stringsAsFactors = FALSE)


          #   stat <- rbind.data.frame(stat, stat_temp) %>% mutate(Compartment = x)
          # }
          
          return(stat)

      }) 
    stat_list <- do.call(rbind.data.frame, stat_list)
    return(stat_list)
  })
  stat_list <- do.call(rbind.data.frame, stat_list)

}, mc.cores = 4)

# Make data frame
stat_df <- do.call(rbind.data.frame, pw_stat_list)

p <- stat_df %>%
mutate(FDR = p.adjust(p.value, "bonferroni"), 
  con.var.explained = con.var.explained*100) %>%
mutate(
  significance = p.value <= alpha,
  Genotype = factor(Genotype, levels = genotype$names[genotype$names %in% .$Genotype]),
  Compartment = factor(Compartment, levels = compartment$names[compartment$names %in% .$Compartment]),
  Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[soilbatch$names %in% .$Soil_Batch]),
  tag = paste0(as.character(round(con.var.explained, 2)), "%,\n", as.character(round(p.value, 4)))
  ) %>%
ggplot(aes(x = Genotype, y = con.var.explained)) +
geom_bar(stat = "identity", aes(colour = significance, fill = Genotype, alpha = significance)) +
geom_text(aes(label = tag, alpha = significance), nudge_y = 2, size = 1) +
facet_grid(Compartment + Soil_Batch ~., space = "free", switch = "y", scale = "free_y") +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.8, `TRUE` = 1), guide = FALSE) +
scale_fill_manual(values = genotype$colour[which(genotype$names %in% stat_df$Genotype)], guide = FALSE) +
coord_flip() +
theme_RTN +
theme() +
labs(y = "Variance explained %", x = "", colour = "", fill = "pvalue", alpha = "")
ggsave(p,file = paste0(figs.out, "CAS_pwCPCOA_grouped_asv.png"), 
              dpi = 600, device = "png", units = "in", 
              bg = "transparent", 
              width = 3, height = 5, limitsize = F)


# Print the statistics
write.table(stat_df, file = paste0(stats.out, "/pw_PERMANOVA_Table.txt"), 
  quote = FALSE, 
  row.names = FALSE, 
  sep  = "\t")

# END OF SCRIPT
sessionInfo()
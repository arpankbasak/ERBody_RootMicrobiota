#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for GLM based clustering of ASVs
# Removing ASVs of importance and obtaining the most contributing taxa
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.itsf)

# Loading required packages
pkgs <- c("tidyverse", "RColorBrewer", "parallel", "cowplot")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
# load(paste0("./data/all_experiments.RData"))
load(paste0("./data/GLM_clustered_family.RData"))

# REad the taxonomy table from the ATSPhereCC
# taxa_at_cc <- read.table("./output/taxonomy_syncom.txt", sep = "\t", header = TRUE, as.is = TRUE)

# What are we anlysing
pro <- c("st")

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

# Filter criteria
cf <- -log10(alpha)
ra_cf <- c(log2(0.05/100), log2(0.05/100), log2(0.05/100))
ra_cf10 <- c(log10(0.05/100), log10(0.5/100), log10(5/100))
lfc_cf <- 2

# Group out the cluster informations
cluster_taxa_reps <- family_clusterd$lfc_annotated.ra %>%
select(lineage_family) %>%
separate(lineage_family, into = c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", remove = TRUE) %>%
select(Order) %>%
mutate(Order = str_replace_all(.$Order, "^NA$", "Unclassified")) %>%
group_by(Order) %>%
summarise(nreps = n()) %>%
mutate(nrep_prec = 100* nreps/sum(nreps)) %>%
mutate(new_Order = ifelse(nrep_prec < 1, "Rare_taxa", Order)) %>%
arrange(desc(nrep_prec)) %>%
data.frame(., stringsAsFactors = FALSE)

# Make colours
set.seed(seeder)
col_rep <- cluster_taxa_reps %>% 
group_by(new_Order) %>% 
summarise(rep = sum(nreps)) %>% 
arrange(desc(rep)) %>% 
.$new_Order

names(col_rep) <- col_rep
col.pal <- colorRampPalette(colors=RColorBrewer::brewer.pal(n=5, "Set2"))(length(unique(col_rep)))
names(col.pal) <- unique(col_rep)
col.pal[which(names(col.pal) == "Unclassified")] <- "#A9A9A9" # Darkgrey
col.pal[which(names(col.pal) == "Rare_taxa")] <- "#D3D3D3" # Lightgrey
sliva_pal <- col.pal

# Plot the mean relative abundance of of the 
temp <- family_clusterd$lfc_annotated.ra
idx <- colnames(temp)[str_detect(colnames(temp), "meanRA$")]

top20 <- data.frame(top_tag = rowSums(temp[,idx])) %>%
mutate(tags = row.names(.)) %>%
arrange(desc(top_tag)) %>%
.$tags %>%
head(., n=20)

# Plot Canvas
set.seed(seeder)
mclapply(X = c("Extract", "Exudate"), FUN = function(x){

  # x = "Extract"
  temp <- family_clusterd$lfc_annotated.ra
  # idx <- str_detect(colnames(temp), paste("Col", x, "mean", sep = "_"))
  idy <- str_detect(colnames(temp), paste0(x, "_genotype_")) & str_detect(colnames(temp), paste0("logFC$|", "PValue$"))
  temp <- temp[,idy] %>% cbind.data.frame(., lineage_asvf = temp[,c("tag")])
  colnames(temp) <- str_replace_all(colnames(temp), "genotype_", "")

    i <- "cyp"
    j <- "pyk10"
    
    message(paste0("Spilling for compartment ", x, " Comparing ", i, " and ", j))
    
    temp.geno <- temp
    idx <- str_detect(colnames(temp.geno), paste(i, x, sep = "_"))
    colnames(temp.geno)[idx] <- c("X_lfc", "X_p")
    
    idy <- str_detect(colnames(temp.geno), paste(j, x, sep = "_"))
    colnames(temp.geno)[idy] <- c("Y_lfc", "Y_p")
    
    sig_asvs <- names(which(rowSums(temp.geno[,c("X_p", "Y_p")] < 0.05) == 1))
    temp.geno <- temp.geno %>%
      mutate(
        significant = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf & X_p <= alpha & Y_p <= alpha),
        FDR = row.names(.) %in% sig_asvs
        # siz = cut(mean_df[match(row.names(.), names(mean_df))], 
          # breaks = c(0,0.01,0.1,1,10, 100), 
          # labels = c("0.01", "0.1", "1", "10", "100"))
      ) %>%
      # separate(Atsphere_strain, into = c("At_Kingdom", "At_Phylum", "At_Class", "At_Order", "At_Family", "At_strain", "X1"), sep = ";", remove = TRUE) %>%
      separate(lineage_asvf, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "X2"), sep = ";", remove = FALSE) %>%
      select(-X2) %>%
      mutate(
        Order = str_replace_all(.$Order, "^NA$", "Unclassified")
        # tag = str_replace_all(paste(At_Order, At_strain, sep = "_"), "^NA_NA$", "Unclassified")
      ) %>%
      mutate(Order = sliva_pal[match(.$Order, names(sliva_pal))]) %>%
      data.frame(., stringsAsFactors = FALSE)

    stat <- cor.test(temp.geno$X_lfc, temp.geno$Y_lfc) 
    til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", stat$p.value)

    # SLIVA  
    p <- temp.geno %>%
    mutate(
      Family_txt = as.character(ifelse(significant == TRUE, Order, "")), 
      Order = factor(Order, levels = names(col.pal)[names(col.pal) %in% .$Order]), 
      # tag = as.character(ifelse(cumRA >= ra_cf10[1], tag, "")), 
      alp = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf),
      prev = lineage_asvf %in% top20
      # siz = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf)
    ) %>%
    ggplot(aes(x = X_lfc, y = Y_lfc)) +
    ggtitle(paste0(paste(x,i,j, sep = "_"), til)) +
    # geom_hline(yintercept = c(0), 
    #         lty = "solid", 
    #         colour = c_grey, 
    #         lwd = 1) +
    # geom_vline(xintercept = c(0), 
    #         lty = "solid", 
    #         colour = c_grey, 
    #         lwd = 1) +
    geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
    geom_point(aes(
          # size = siz,
          colour = prev), 
          shape = 21, fill = soilbatch$colours[soilbatch$names == "CAS15"], 
          alpha = 0.8, size = 3,
          na.rm = FALSE) +
    # ggrepel::geom_text_repel(aes(label = Family_txt), 
    #         colour = c_black, 
    #         size = 2,  max.overlaps = 15,
    #         box.padding = unit(0.5, "lines"),
    #         point.padding = unit(0.3, "lines")
    #         ) +
    # scale_fill_manual(values = col.pal) +
    # scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.8), guide = FALSE) +
    # scale_size_manual(values = c(0.3, .5, .7, .9, 1.2), guide = FALSE) +
    scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
    theme_RTN_MDS +
    theme(panel.spacing = unit(0.5, "lines"), 
                              axis.text.x = element_text(size = 12),
                              axis.text.y = element_text(size = 12),
                              # axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
                              strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
                              legend.text = element_text(size = 2),
                              legend.title = element_blank(),
                              axis.line.y = element_line(size = 1),
                              axis.line.x = element_line(size = 1),
                              axis.ticks.y = element_line(size = 1),
                              axis.ticks.x = element_line(size = 1),
                              plot.title = element_text(size=4)
                              ) +
    labs(x = i, 
      y = j)
    ggsave(plot = p, file = paste0(figs.out, "/Parwise_",x,"_", j,"_", i,"_logFC.png"), 
                dpi = 600, 
                units = img, 
                device = "png", 
                bg = "transparent", 
                width = 4, 
                height = 3, 
                limitsize = T)

    # AtSphere
    # temp.geno %>%
    # mutate(
    #   Family_txt = as.character(ifelse(significant == TRUE, At_Order, "")), 
    #   tag = factor(tag, levels = names(col.pal.strain)[names(col.pal.strain) %in% .$tag]), 
    #   # tag = as.character(ifelse(cumRA >= ra_cf10[1], tag, "")), 
    #   alp = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf),
    #   siz = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf)
    #   ) %>%
    # ggplot(aes(x = X_lfc, y = Y_lfc)) +
    # ggtitle(paste0(paste(x,i,j, sep = "_"), til)) +
    # # geom_hline(yintercept = c(0), 
    # #         lty = "solid", 
    # #         colour = c_grey, 
    # #         lwd = 1) +
    # # geom_vline(xintercept = c(0), 
    # #         lty = "solid", 
    # #         colour = c_grey, 
    # #         lwd = 1) +
    # geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
    # geom_point(aes(
    #       alpha = alp, 
    #       size = siz,
    #       colour = FDR,
    #       fill = tag), 
    #       shape = 21, 
    #       na.rm = FALSE) +
    # ggrepel::geom_text_repel(aes(label = Family_txt), 
    #         colour = c_black, 
    #         size = 2,  max.overlaps = 15,
    #         box.padding = unit(0.5, "lines"),
    #         point.padding = unit(0.3, "lines")
    #         ) +
    # scale_fill_manual(values = col.pal.strain, guide = FALSE) +
    # scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.8), guide = FALSE) +
    # scale_size_manual(values = c(0.3, .5, .7, .9, 1.2), guide = FALSE) +
    # scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
    # theme_RTN +
    # theme(panel.spacing = unit(0.5, "lines"), 
    #                           axis.text.x = element_text(size = 12),
    #                           axis.text.y = element_text(size = 12),
    #                           # axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
    #                           strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
    #                           legend.text = element_text(size = 2),
    #                           legend.title = element_blank(),
    #                           axis.line.y = element_line(size = 1),
    #                           axis.line.x = element_line(size = 1),
    #                           axis.ticks.y = element_line(size = 1),
    #                           axis.ticks.x = element_line(size = 1),
    #                           plot.title = element_text(size=8)
    #                           ) +
    # labs(x = i, 
    #   y = j) +
    # ggsave(file = paste0(figs.out, "/atcc_Parwise_",x,"_", j,"_", i,"_logFC.png"), 
    #             dpi = 600, 
    #             units = img, 
    #             device = "png", 
    #             bg = "transparent", 
    #             width = 4, 
    #             height = 4, 
    #             limitsize = T)

    message(paste0("DONE!! compartment ", x, " Comparing ", i, " and ", j))
      
}, mc.cores = 12)

# within genotype and across treatment
mclapply(X = c("cyp", "pyk10"), FUN = function(x){

  # x = "cyp"
  temp <- family_clusterd$lfc_annotated.ra
  # idx <- str_detect(colnames(temp), paste("Col", x, "mean", sep = "_"))
  idy <- (str_detect(colnames(temp), paste0(x, "_Exudate_genotype_")) | str_detect(colnames(temp), paste0(x, "_Extract_genotype_"))) & str_detect(colnames(temp), paste0("logFC$|", "PValue$"))
  temp <- temp[,idy] %>% cbind.data.frame(., lineage_asvb = temp[,c("tag")])
  colnames(temp) <- str_replace_all(colnames(temp), "genotype_", "")

    i <- "Extract"
    j <- "Exudate"
    
    message(paste0("Spilling for compartment ", x, " Comparing ", i, " and ", j))
    
    temp.geno <- temp
    idx <- str_detect(colnames(temp.geno), paste(x, i, sep = "_"))
    colnames(temp.geno)[idx] <- c("X_lfc", "X_p")
    
    idy <- str_detect(colnames(temp.geno), paste(x, j, sep = "_"))
    colnames(temp.geno)[idy] <- c("Y_lfc", "Y_p")
    
    sig_asvs <- names(which(rowSums(temp.geno[,c("X_p", "Y_p")] < 0.05) == 1))
    temp.geno <- temp.geno %>%
      mutate(
        significant = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf & X_p <= alpha & Y_p <= alpha),
        FDR = row.names(.) %in% sig_asvs
        # siz = cut(mean_df[match(row.names(.), names(mean_df))], 
          # breaks = c(0,0.01,0.1,1,10, 100), 
          # labels = c("0.01", "0.1", "1", "10", "100"))
      ) %>%
      # separate(Atsphere_strain, into = c("At_Kingdom", "At_Phylum", "At_Class", "At_Order", "At_Family", "At_strain", "X1"), sep = ";", remove = TRUE) %>%
      separate(lineage_asvb, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "X2"), sep = ";", remove = FALSE) %>%
      select(-X2) %>%
      mutate(
        Order = str_replace_all(.$Order, "^NA$", "Unclassified")
        # tag = str_replace_all(paste(At_Order, At_strain, sep = "_"), "^NA_NA$", "Unclassified")
      ) %>%
      mutate(Order = sliva_pal[match(.$Order, names(sliva_pal))]) %>%
      data.frame(., stringsAsFactors = FALSE)

    stat <- cor.test(temp.geno$X_lfc, temp.geno$Y_lfc) 
    til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", stat$p.value)

    # SLIVA  
    p <- temp.geno %>%
    mutate(
      Family_txt = as.character(ifelse(significant == TRUE, Order, "")), 
      Order = factor(Order, levels = names(col.pal)[names(col.pal) %in% .$Order]), 
      # tag = as.character(ifelse(cumRA >= ra_cf10[1], tag, "")), 
      alp = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf),
      prev = lineage_asvb %in% top20
      # siz = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf)
    ) %>%
    ggplot(aes(x = X_lfc, y = Y_lfc)) +
    ggtitle(paste0(paste(x,i,j, sep = "_"), til)) +
    # geom_hline(yintercept = c(0), 
    #         lty = "solid", 
    #         colour = c_grey, 
    #         lwd = 1) +
    # geom_vline(xintercept = c(0), 
    #         lty = "solid", 
    #         colour = c_grey, 
    #         lwd = 1) +
    geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
    geom_point(aes(
          # size = siz,
          colour = prev), 
          shape = 21, fill = soilbatch$colours[soilbatch$names == "CAS15"], 
          alpha = 0.8,
          na.rm = FALSE, size = 3) +
    # ggrepel::geom_text_repel(aes(label = Family_txt), 
    #         colour = c_black, 
    #         size = 2,  max.overlaps = 15,
    #         box.padding = unit(0.5, "lines"),
    #         point.padding = unit(0.3, "lines")
    #         ) +
    # scale_fill_manual(values = col.pal) +
    # scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.8), guide = FALSE) +
    # scale_size_manual(values = c(0.3, .5, .7, .9, 1.2), guide = FALSE) +
    scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
    theme_RTN_MDS +
    theme(panel.spacing = unit(0.5, "lines"), 
                              axis.text.x = element_text(size = 12),
                              axis.text.y = element_text(size = 12),
                              # axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
                              strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
                              legend.text = element_text(size = 2),
                              legend.title = element_blank(),
                              axis.line.y = element_line(size = 1),
                              axis.line.x = element_line(size = 1),
                              axis.ticks.y = element_line(size = 1),
                              axis.ticks.x = element_line(size = 1),
                              plot.title = element_text(size=4)
                              ) +
    labs(x = i, 
      y = j)
    ggsave(p, file = paste0(figs.out, "/Parwise_",x,"_", j,"_", i,"_logFC.png"), 
                dpi = 600, 
                units = img, 
                device = "png", 
                bg = "transparent", 
                width = 4, 
                height = 3, 
                limitsize = T)

    # AtSphere
    # temp.geno %>%
    # mutate(
    #   Family_txt = as.character(ifelse(significant == TRUE, At_Order, "")), 
    #   tag = factor(tag, levels = names(col.pal.strain)[names(col.pal.strain) %in% .$tag]), 
    #   # tag = as.character(ifelse(cumRA >= ra_cf10[1], tag, "")), 
    #   alp = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf),
    #   siz = (abs(X_lfc) > lfc_cf & abs(Y_lfc) > lfc_cf)
    #   ) %>%
    # ggplot(aes(x = X_lfc, y = Y_lfc)) +
    # ggtitle(paste0(paste(x,i,j, sep = "_"), til)) +
    # # geom_hline(yintercept = c(0), 
    # #         lty = "solid", 
    # #         colour = c_grey, 
    # #         lwd = 1) +
    # # geom_vline(xintercept = c(0), 
    # #         lty = "solid", 
    # #         colour = c_grey, 
    # #         lwd = 1) +
    # geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
    # geom_point(aes(
    #       alpha = alp, 
    #       size = siz,
    #       colour = FDR,
    #       fill = tag), 
    #       shape = 21, 
    #       na.rm = FALSE) +
    # ggrepel::geom_text_repel(aes(label = Family_txt), 
    #         colour = c_black, 
    #         size = 2,  max.overlaps = 15,
    #         box.padding = unit(0.5, "lines"),
    #         point.padding = unit(0.3, "lines")
    #         ) +
    # scale_fill_manual(values = col.pal.strain, guide = FALSE) +
    # scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.8), guide = FALSE) +
    # scale_size_manual(values = c(0.3, .5, .7, .9, 1.2), guide = FALSE) +
    # scale_colour_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
    # theme_RTN +
    # theme(panel.spacing = unit(0.5, "lines"), 
    #                           axis.text.x = element_text(size = 12),
    #                           axis.text.y = element_text(size = 12),
    #                           # axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
    #                           strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
    #                           legend.text = element_text(size = 2),
    #                           legend.title = element_blank(),
    #                           axis.line.y = element_line(size = 1),
    #                           axis.line.x = element_line(size = 1),
    #                           axis.ticks.y = element_line(size = 1),
    #                           axis.ticks.x = element_line(size = 1),
    #                           plot.title = element_text(size=8)
    #                           ) +
    # labs(x = i, 
    #   y = j) +
    # ggsave(file = paste0(figs.out, "/atcc_Parwise_",x,"_", j,"_", i,"_logFC.png"), 
    #             dpi = 600, 
    #             units = img, 
    #             device = "png", 
    #             bg = "transparent", 
    #             width = 4, 
    #             height = 4, 
    #             limitsize = T)

    message(paste0("DONE!! compartment ", x, " Comparing ", i, " and ", j))
      
}, mc.cores = 12)



p <- temp.ra <- temp[,idx] %>%
add_column(tag = row.names(.), .before = 1) %>%
gather(key = "key", value = "vals", -tag) %>%
separate(key, into = c("x", "Genotype", "Treatment", "x1", "x2"), sep = "_") %>%
filter(tag %in% top20) %>%
mutate(
  Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
  Treatment = factor(Treatment, levels = c("Bulk.soil", "Exudate", "Extract")),
  tag = factor(tag, levels = top20))

temp.ra %>%
ggplot(aes(x = tag, y = 100*vals)) +
geom_bar(stat = "identity", position = "dodge", aes(fill = Genotype), colour = c_dark_grey) +
facet_grid(.~Treatment, space = "free_x", scale = "free_x", switch = "x") +
scale_fill_manual(values = genotype$colours[which(genotype$short %in% temp.ra$Genotype)]) +
theme_RTN +
theme(
  panel.spacing = unit(0.5, "lines"), 
  axis.text.x = element_text(size = 6, angle = 90, hjust = 0.9, vjust = 0.5),
  axis.text.y = element_text(size = 10), 
  axis.title = element_text(size = 10)) +
labs(x = "Top 20 abundant families", y = "mean RA (%)")
ggsave(p, file = paste0(figs.out, "/rank_abundance_top20_family_level.png"), 
            dpi = 600, 
            units = img, 
            device = "png", 
            bg = "transparent", 
            width = 7.5, 
            height = 6, 
            limitsize = T)


# Correlation between the genotypes
temp <- family_clusterd$lfc_annotated.ra
idx <- colnames(temp)[str_detect(colnames(temp), "logFC$")]
d <- 1-cor(apply(temp[,idx], 2, function(x) replace_na(x, -6)))
mds <- cmdscale(as.dist(d), eig = TRUE, k = 3)
eigen <- round(100*mds$eig/sum(mds$eig), 3)

mds_df <- as.data.frame(mds$points) %>%
mutate(group = row.names(.)) %>%
separate(group, into = c("X", "Genotype", "Treatment", "comparison", "X1"), sep = "_") %>%
select(-X, -X1) %>%
mutate(
	comparison = as.factor(comparison), 
  Genotype = factor(Genotype, levels = genotype$short[genotype$short %in% .$Genotype]),
  Treatment = factor(Treatment, levels = c("Exudate", "Extract")),
  comp = paste(Genotype, Treatment, sep = "_"))

mds_df_tests <-  mds_df %>%
group_by(Genotype) %>%
summarise(cV1 = mean(V1), cV2 = mean(V2), cV3 = mean(V3)) %>%
data.frame(.)

mds_df <- mds_df %>% 
cbind.data.frame(., mds_df_tests[match(mds_df$Genotype, mds_df_tests$Genotype), c("cV1", "cV2", "cV3")])

# MDS for the correlation
p <- mds_df %>%
ggplot(aes(x = saturate(V1), y = saturate(V2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
# ggtitle(paste0(ti)) +
geom_segment(aes(
  x = saturate(cV1), 
  y = saturate(cV2),
  xend = saturate(V1), 
  yend = saturate(V2), 
  colour = Genotype, 
  lty  = comparison), 
lwd = 0.3) +
geom_point(alpha = 0.8, size = 3, colour = c_black,
    aes(fill = Genotype, 
        shape = Treatment), # Soiltype
    show.legend = T) +
scale_color_manual(values = genotype$colours[genotype$short %in% mds_df$Genotype]) + 
scale_fill_manual(values = genotype$colours[genotype$short %in% mds_df$Genotype]) + 
scale_shape_manual(values = c(22, 23, 24, 25)) +
scale_linetype_manual(values = c(`genotype` = "solid", `soil` = "dashed"), guide = FALSE) +
theme_RTN_MDS +
labs(x = paste0("Axis - 1: ",eigen[1],"%"), 
     y = paste0("Axis - 2: ",eigen[2],"%")) +
ggsave(p, filename = paste0(figs.out,"/MDS_coorealtion_12_family.png"), 
           dpi = 600, 
           device = "png", 
           bg = "transparent", 
           units = img, 
           width = cpcoa.box[1], 
           height = cpcoa.box[2],
           limitsize = FALSE)

p <- mds_df %>%
ggplot(aes(x = saturate(V3), y = saturate(V2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
# ggtitle(paste0(ti)) +
geom_segment(aes(
  x = saturate(cV3), 
  y = saturate(cV2),
  xend = saturate(V3), 
  yend = saturate(V2), 
  colour = Genotype, 
  lty  = comparison), 
lwd = 0.3) +
geom_point(alpha = 0.8, size = 3, colour = c_black,
    aes(fill = Genotype, 
        shape = Treatment), # Soiltype
    show.legend = T) +
scale_color_manual(values = genotype$colours[genotype$short %in% mds_df$Genotype]) + 
scale_fill_manual(values = genotype$colours[genotype$short %in% mds_df$Genotype]) + 
scale_shape_manual(values = c(22, 23, 24, 25)) +
scale_linetype_manual(values = c(`genotype` = "solid", `soil` = "dashed"), guide = FALSE) +
theme_RTN_MDS +
labs(x = paste0("Axis - 3: ",eigen[3],"%"), 
     y = paste0("Axis - 2: ",eigen[2],"%"))
ggsave(p, filename = paste0(figs.out,"/MDS_coorealtion_23_family.png"), 
           dpi = 600, 
           device = "png", 
           bg = "transparent", 
           units = img, 
           width = cpcoa.box[1], 
           height = cpcoa.box[2],
           limitsize = FALSE)

# Heatmaps -- make another by just filtering those asvs that have logFC > lfc_cutoff
idx <- colnames(temp)[!str_detect(colnames(temp), "meanRA$")]
temp.mat_lfc <- temp[, idx] %>%
    # add_column(asvf = row.names(.), .before = 1) %>%
    gather(key = "key", value = "val", convert = FALSE, -tag, -lineage_family) %>%
    separate(key, into = c("X", "Genotype", "Treatment", "comparison", "vals"), sep = "_") %>%
    select(-X) %>%
    spread(key = vals, value = val, fill = NA, convert = FALSE) %>%
    mutate(
      group = as.factor(comparison),
      significant = (abs(logFC) > 0.5 & PValue <= alpha),
      FDR = (PValue <= alpha),
      Genotype = factor(.$Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
      Treatment = factor(Treatment, levels = c("Exudate", "Extract")),
      # lfc_sat = saturate(logFC),
      # RA = ifelse(!is.finite(log10(RA)), log10(RA+0.0001), log10(RA)),
      # cluster = factor(cluster, levels = family_clusterd$cluster_sorted),
      tag = factor(tag, levels = family_clusterd$lfc_sorted)) %>%
    # separate(Atsphere_strain, into = c("At_Kingdom", "At_Phylum", "At_Class", "At_Order", "At_Family", "At_strain", "X1"), sep = ";", remove = TRUE) %>%
    separate(lineage_family, into = c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", remove = TRUE) %>%
    mutate(
      Order = str_replace_all(.$Order, "^NA$", "Unclassified")
      # tag = str_replace_all(paste(At_Family, At_strain, sep = "_"), "^NA_NA$", "Unclassified"),
      # cumRA = ifelse(RA < -5, -5, RA)
    ) %>%
    mutate(Order = factor(cluster_taxa_reps$new_Order[match(.$Order, cluster_taxa_reps$Order)], levels = names(col.pal))) %>%
    data.frame(., stringsAsFactors = FALSE)

idx <- colnames(temp)[str_detect(colnames(temp), "meanRA$")]
temp.mat_RA <- temp[, c("tag", "lineage_family", idx)] %>%
    # add_column(asvf = row.names(.), .before = 1) %>%
    gather(key = "key", value = "val", convert = FALSE, -tag, -lineage_family) %>%
    separate(key, into = c("X", "Genotype", "Treatment", "comparison", "vals"), sep = "_") %>%
    select(-X) %>%
    spread(key = vals, value = val, fill = NA, convert = FALSE) %>%
    mutate(
      group = as.factor(comparison),
      # significant = (abs(logFC) > 0.5 & PValue <= alpha),
      # FDR = (PValue <= alpha),
      Genotype = factor(.$Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
      Treatment = factor(Treatment, levels = c("Bulk.soil", "Exudate", "Extract")),
      # lfc_sat = saturate(logFC),
      RA = ifelse(!is.finite(log10(meanRA)), log10(meanRA+0.0001), log10(meanRA)),
      # cluster = factor(cluster, levels = family_clusterd$cluster_sorted),
      tag = factor(tag, levels = family_clusterd$lfc_sorted)) %>%
    # separate(Atsphere_strain, into = c("At_Kingdom", "At_Phylum", "At_Class", "At_Order", "At_Family", "At_strain", "X1"), sep = ";", remove = TRUE) %>%
    separate(lineage_family, into = c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", remove = TRUE) %>%
    mutate(
      Order = str_replace_all(.$Order, "^NA$", "Unclassified"),
      # tag = str_replace_all(paste(At_Family, At_strain, sep = "_"), "^NA_NA$", "Unclassified"),
      cumRA = ifelse(RA < -5, -5, RA)
    ) %>%
    mutate(Order = factor(cluster_taxa_reps$new_Order[match(.$Order, cluster_taxa_reps$Order)], levels = names(col.pal))) %>%
    data.frame(., stringsAsFactors = FALSE)

# Heatmap -- RA
(hmap_ra <- 
temp.mat_RA %>%
# filter(group == "meanra") %>% 
ggplot(aes(x = Genotype, y = tag, fill = cumRA)) +
geom_raster(alpha=1) +
    scale_fill_gradient2(
        high=c_red, 
        mid = c_yellow,
        low = c_black, 
        na.value = c_black,
        midpoint = ra_cf10[2], 
        breaks = c(-4, -3, -2, -1), 
        limits = c(-5,0), 
        labels = paste0(c(".01", ".1", "1", "10"))) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="% meanRA"))
    ggsave(hmap_ra, filename = paste0(figs.out,"/RA_asv_heatmap_lfc_sorted_family.png"),
           bg="transparent", 
           width=3, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

# Heatmap -- LFC
(hmap_lfc <- 
temp.mat_lfc %>%
filter(group == "genotype") %>%
mutate(logFC = ifelse(abs(logFC) > 6, sign(logFC) * 6, logFC )) %>%
ggplot(aes(x = Genotype, y = tag, fill = logFC)) +
geom_raster(alpha=1) +
# geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = FDR)) +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-6, -3, 0, 3, 6), 
        limits = c(-6,6), 
        labels = paste0(c("-6", "-3", "0", "3", "6"))) +
# scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = c_black), guide = FALSE) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="lfc(vCol-0)")) +
    ggsave(hmap_lfc, filename = paste0(figs.out,"/LFC_asv_heatmap_lfc_km_sorted_family.png"),
           bg="transparent", 
           width=3, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

(hmap_lfc_vSoil <- 
temp.mat_lfc %>%
filter(group == "soil") %>%
mutate(logFC = ifelse(abs(logFC) > 6, sign(logFC) * 6, logFC )) %>%
ggplot(aes(x = Genotype, y = tag, fill = logFC)) +
geom_raster(alpha=1) +
# geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = FDR)) +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-6, -3, 0, 3, 6), 
        limits = c(-6,6), 
        labels = paste0(c("-6", "-3", "0", "3", "6"))) +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = c_black), guide = FALSE) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="lfc(vs Soil)"))
    ggsave(hmap_lfc_vSoil, filename = paste0(figs.out,"/LFC_asv_heatmap_lfc_km_sorted_vs_soil_family.png"),
           bg="transparent", 
           width=4, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)


ra_group <- temp.mat_RA %>%
mutate(group = paste(tag, Genotype, Treatment, sep = "_")) %>%
select(group, meanRA)

# Heatmap for the GH_prevelant orders
p <- temp.mat_RA %>%
filter(tag %in% top20) %>%
# filter(Order != "Mortierellales") %>%
# mutate(
#   behaviour = as.factor(ifelse(!is.na(logFC), sign(logFC), -1)),
#   # ra = log10(ra),
#   logFC = ifelse(abs(logFC) > 2, sign(logFC)*2, logFC),
#   group = paste(tag, Genotype, Treatment, sep = "_")
#   ) %>%
# add_column(meanRA = ra_group$meanRA[match(.$group, ra_group$group)]) %>%
mutate(logRA = log10(meanRA)) %>%
ggplot(aes(y = Genotype, x = tag)) +
geom_point(alpha=0.8, aes(fill = logRA), shape = 21, size = 3, colour = c_black, na.rm = FALSE) +
# geom_tile(alpha=0.5, aes(fill = ra), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(Treatment ~ Class, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_red, 
        mid = c_yellow,
        low = c_black, 
        na.value = c_black,
        midpoint = -2, 
        breaks = c(-4, -3, -2, -1, 0), 
        limits = c(-4,0), 
        labels = paste0(c("<0.001%", "0.01%", "1%", "10%", "100%"))) +
# scale_shape_manual(values = c(`1` = 24, `-1` = 25)) +
scale_size_continuous(range = c(1,4)) +
scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_dark_grey)) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.y=element_text(angle = 0, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text.x=element_text(size = 4, angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/prevelant_heatmap_t20family.png"),
       bg="transparent", 
       width=5, 
       height=7,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

p <- temp.mat_lfc %>%
filter(tag %in% top20) %>%
filter(group == "genotype", !Genotype %in% c("MS", "PBS")) %>%
# filter(Order != "Mortierellales") %>%
mutate(
  behaviour = as.factor(ifelse(!is.na(logFC), sign(logFC), -1)),
  # ra = log10(ra),
  logFC = ifelse(abs(logFC) >= 0.5, sign(logFC)*0.5, logFC),
  group = paste(tag, Genotype, Treatment, sep = "_")
  ) %>%
# add_column(meanRA = ra_group$meanRA[match(.$group, ra_group$group)]) %>%
# mutate(logRA = log10(meanRA)) %>%
ggplot(aes(y = Genotype, x = tag)) +
geom_point(alpha=0.8, aes(fill = logFC), shape = 21, size = 3, colour = c_black, na.rm = FALSE) +
# geom_tile(alpha=0.5, aes(fill = ra), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(Treatment ~ Class, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-0.5, 0, 0.5), 
        limits = c(-0.5,0.5), 
        labels = paste0(c("<-0.5", "0", "0.5"))) +
# scale_shape_manual(values = c(`1` = 24, `-1` = 25)) +
scale_size_continuous(range = c(1,4)) +
scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_dark_grey)) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.y=element_text(angle = 0, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text.x=element_text(size = 4, angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/prevelant_heatmap_LFC_t20family.png"),
       bg="transparent", 
       width=5, 
       height=6,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)


# BArplot
p <- temp.mat_RA %>%
filter(tag %in% top20) %>%
group_by(Class, tag, Treatment) %>%
summarise(RA = (sum(meanRA)), meanra = log10(mean(meanRA)+0.0001)) %>%
ggplot(aes(y = RA, x = tag)) +
geom_bar(position = "dodge", stat = "identity", aes(fill = Treatment), colour = c_black) +
# geom_line(aes(y = meanra, x = tag, lty = Treatment), colour = c_red) +
facet_grid(. ~ Class, switch = "both", scales="free", space="free") +
scale_fill_manual(values = c(`Extract` = "#B4B4B4CC", `Exudate` = "#000000CC", `Bulk.soil` = "#654321CC")) +
scale_linetype_manual(values = c(`Extract` = "dashed", `Exudate` = "solid", `Bulk.soil` = "dotted")) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.y=element_text(angle = 0, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text.x=element_text(size = 4, angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/barplot_t20family.png"),
       bg="transparent", 
       width=5, 
       height=5,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

# boxplot of lfc
bxp <- temp.mat_lfc %>%
filter(group == "genotype", !Genotype %in% c("PBS", "MS")) %>%
mutate(
  logFC = ifelse(abs(logFC) > 6, sign(logFC) * 6, logFC ),
  direction = as.factor(ifelse(abs(logFC) >= 0.5, sign(logFC)*1, 0))
  ) %>%
ggplot(aes(x = Genotype, y = abs(logFC))) +
geom_point(
  alpha = 0.8, 
  shape = 21, 
  size = 0.5,
  fill = c_white, 
  colour = c_black) +
geom_boxplot(fill = c_white, outlier.alpha = 0, alpha = 0.2, size = 0.5, colour = c_black) +
scale_fill_gradient2(low = c_cudo_magenta, high = c_dark_green, midpoint = 0, mid = c_white) +
facet_grid(. ~ Treatment, switch = "both", scales="free", space="free") +
theme_RTN +
theme() +
labs(y="", x="", fill="lfc(vs Col-0)")
ggsave(bxp, filename = paste0(figs.out,"/LFC_family_boxplot_DA.png"),
           bg="transparent", 
           width=3, 
           height=3,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)


temp.temp <- temp.mat_lfc %>%
filter(group == "genotype") %>%
select(tag, logFC, Genotype, Treatment) %>%
na.omit(.) %>%
group_by(tag,Treatment) %>%
spread(key = Genotype, value = logFC, fill = 0) %>%
ggplot(aes(x = cyp, y = pyk10)) +
geom_hline(yintercept = 0, lty = "solid", colour = c_dark_grey) +
geom_vline(xintercept = 0, lty = "solid", colour = c_dark_grey) +
geom_point(aes(fill = Treatment),
  alpha = 0.8, 
  shape = 21, 
  size = 0.8, 
  colour = c_black
  ) +
geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
scale_fill_manual(values = c(`Extract` = c_white, `Exudate` = c_black)) +
theme_RTN_MDS +
theme() +
labs(y="", x="", fill="")
ggsave(temp.temp, filename = paste0(figs.out,"/pyk10_vs_cyp_family_.png"),
           bg="transparent", 
           width=5, 
           height=5,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600
           )


temp.temp <- temp.mat_lfc %>%
filter(group == "genotype") %>%
select(tag, logFC, Genotype, Treatment) %>%
na.omit(.) %>%
group_by(tag,Treatment) %>%
spread(key = Genotype, value = logFC, fill = 0) %>%
ggplot(aes(x = cyp, y = pyk10)) +
geom_hline(yintercept = 0, lty = "solid", colour = c_dark_grey) +
geom_vline(xintercept = 0, lty = "solid", colour = c_dark_grey) +
geom_point(aes(fill = Treatment),
  alpha = 0.8, 
  shape = 21, 
  size = 0.8, 
  colour = c_black
  ) +
geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
scale_fill_manual(values = c(`Extract` = c_white, `Exudate` = c_black)) +
theme_RTN_MDS +
theme() +
labs(y="", x="", fill="")
ggsave(temp.temp, filename = paste0(figs.out,"/pyk10_vs_cyp_family_.png"),
           bg="transparent", 
           width=5, 
           height=5,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600
           )

# Strip -- Taxonomy -- Order
# In Silva
# (strip_unite <- 
# temp.mat_RA %>%
# mutate(Order = factor(Order, levels = names(col.pal)[which(names(col.pal) %in% .$Order)])) %>%
# ggplot(aes(x = "", y = tag, fill = Order)) +
# geom_raster(alpha=1) +
# scale_fill_manual(values = col.pal, guide = FALSE) +
#     facet_grid(Order ~ ., switch = "both", scales="free", space="free") +
#     theme_RTN +
#     theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
#           axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
#           strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
#           axis.text.y=element_text(size = 2),
#           strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
#     labs(y="", x="", fill="")) +
#     ggsave(filename = paste0(figs.out,"/UNITE_annotated_lfc_sorted.png"),
#            bg="transparent", 
#            width=2, 
#            height=8,
#            units = img, 
#            limitsize=F, 
#            device = "png", 
#            dpi = 600)

# col.pal.strain[1] <- "white"
# (strip_atcc <- 
# temp %>%
# mutate(At_Family = factor(tag, levels = names(col.pal.strain)[which(names(col.pal.strain) %in% .$tag)])) %>%
# ggplot(aes(x = "", y = asvf, fill = At_Family)) +
# geom_raster(alpha=1) +
# scale_fill_manual(values = col.pal.strain) +
#     facet_grid(cluster ~ ., switch = "both", scales="free", space="free") +
#     theme_RTN +
#     theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
#           axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
#           strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
#           axis.text.y=element_text(size = 10),
#           strip.text.y = element_text(size = 15, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
#     labs(y="", x="", fill="")) +
#     ggsave(filename = paste0(figs.out,"/atcc_annotated_asv_heatmap_lfc_km_sorted.png"),
#            bg="transparent", 
#            width=1.75, 
#            height=8,
#            units = img, 
#            limitsize=F, 
#            device = "png", 
#            dpi = 600)

# Spill Canvas in one aisle
# cowplot::plot_grid(
#   # theme_clean(strip_clusters),
#   theme_clean(hmap_ra), 
#   # theme_clean(hmap_lfc),
#   # theme_clean(hmap_lfc_vSoil),   
#   theme_clean(strip_unite),
#   # theme_clean(strip_atcc), 
#   ncol=2, 
#   align='hv', 
#   axis = 'tblr',
#   labels=NULL, 
#   rel_widths = c(0.5, 2.5),
#   rel_heights = c(1, 1)
#   ) +
# ggsave(filename = paste0(figs.out,"/summary_glm_family_heatmap_lfc_sorted.png"),
#            bg="transparent", 
#            width=8, 
#            height=7,
#            units = img, 
#            limitsize=F, 
#            device = "png", 
#            dpi = 600)

idx <- unique(as.character(temp.mat_lfc$tag[temp.mat_lfc$significant == TRUE]))
ids <- unique(as.character(temp.mat_lfc$tag[temp.mat_lfc$significant == TRUE & 
	temp.mat_lfc$group == "genotype"]))
temp.mat <- temp.mat_RA %>% 
filter(tag %in% idx) %>%
mutate(target = tag %in% ids)

# Heatmap -- RA
(hmap_ra <- 
temp.mat %>%
# filter(group == "meanra") %>%
ggplot(aes(x = Genotype, y = tag, fill = cumRA)) +
geom_raster(alpha=1) +
scale_fill_gradient2(
        high=c_red, 
        mid = c_yellow,
        low = c_black, 
        na.value = c_black,
        midpoint = ra_cf10[2], 
        breaks = c(-4, -3, -2, -1), 
        limits = c(-5,0), 
        labels = paste0(c(".01", ".1", "1", "10"))) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="% meanRA"))
    ggsave(hmap_ra, filename = paste0(figs.out,"/RA_asv_heatmap_lfc_sorted_sig_family.png"),
           bg="transparent", 
           width=3, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)


# Strip -- Taxonomy -- Order
# In Silva
(strip_unite <- 
temp.mat %>%
mutate(Order = factor(Order, levels = names(col.pal)[which(names(col.pal) %in% .$Order)])) %>%
ggplot(aes(x = "", y = tag, fill = Order)) +
geom_raster(alpha=1) +
geom_tile(fill = NA, width = .95, size = .5, height = .95, aes(colour = target)) +
scale_fill_manual(values = col.pal, guide = FALSE) +
scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = NA)) +
    facet_grid(Order ~ ., switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, angle = 90, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill=""))
    ggsave(strip_unite, filename = paste0(figs.out,"/UNITE_annotated_lfc_sorted_sig_family.png"),
           bg="transparent", 
           width=2, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

temp.mat <- temp.mat_lfc %>% filter(tag %in% idx)

# Heatmap -- LFC
(hmap_lfc <- 
temp.mat %>%
filter(group == "genotype") %>%
mutate(logFC = ifelse(abs(logFC) > 6, sign(logFC) * 6, logFC )) %>%
ggplot(aes(x = Genotype, y = tag, fill = logFC)) +
geom_raster(alpha=1) +
# geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = FDR)) +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-6, -3, 0, 3, 6), 
        limits = c(-6,6), 
        labels = paste0(c("-6", "-3", "0", "3", "6"))) +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = c_black), guide = FALSE) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="lfc(vCol-0)"))
    ggsave(hmap_lfc, filename = paste0(figs.out,"/LFC_asv_heatmap_lfc_km_sorted_sig_family.png"),
           bg="transparent", 
           width=3, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

(hmap_lfc_vSoil <- 
temp.mat %>%
filter(group == "soil") %>%
mutate(logFC = ifelse(abs(logFC) > 6, sign(logFC) * 6, logFC )) %>%
ggplot(aes(x = Genotype, y = tag, fill = logFC)) +
geom_raster(alpha=1) +
# geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = FDR)) +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-6, -3, 0, 3, 6), 
        limits = c(-6,6), 
        labels = paste0(c("-6", "-3", "0", "3", "6"))) +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = c_black), guide = FALSE) +
    facet_grid(Order ~ Treatment, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 2, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 2, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 2),
          strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="lfc(vs Soil)"))
    ggsave(hmap_lfc_vSoil, filename = paste0(figs.out,"/LFC_asv_heatmap_lfc_km_sorted_vs_soil_sig_family.png"),
           bg="transparent", 
           width=3, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

# col.pal.strain[1] <- "white"
# (strip_atcc <- 
# temp %>%
# mutate(At_Family = factor(tag, levels = names(col.pal.strain)[which(names(col.pal.strain) %in% .$tag)])) %>%
# ggplot(aes(x = "", y = asvf, fill = At_Family)) +
# geom_raster(alpha=1) +
# scale_fill_manual(values = col.pal.strain) +
#     facet_grid(cluster ~ ., switch = "both", scales="free", space="free") +
#     theme_RTN +
#     theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
#           axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
#           strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
#           axis.text.y=element_text(size = 10),
#           strip.text.y = element_text(size = 15, hjust = 0.5, vjust = 0.5),
#           legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
#     labs(y="", x="", fill="")) +
#     ggsave(filename = paste0(figs.out,"/atcc_annotated_asv_heatmap_lfc_km_sorted.png"),
#            bg="transparent", 
#            width=1.75, 
#            height=8,
#            units = img, 
#            limitsize=F, 
#            device = "png", 
#            dpi = 600)

# Spill Canvas in one aisle
# cowplot::plot_grid(
#   # theme_clean(strip_clusters),
#   theme_clean(hmap_ra), 
#   # theme_clean(hmap_lfc),
#   # theme_clean(hmap_lfc_vSoil),   
#   theme_clean(strip_unite),
#   # theme_clean(strip_atcc), 
#   ncol=2, 
#   align='hv', 
#   axis = 'tblr',
#   labels=NULL, 
#   rel_widths = c(.5, 2.5),
#   rel_heights = c(1, 1),
#   greedy = TRUE
#   ) +
# ggsave(filename = paste0(figs.out,"/summary_glm_family_heatmap_lfc_sorted_sig.png"),
#            bg="transparent", 
#            width=7, 
#            height=8,
#            units = img, 
#            limitsize=F, 
#            device = "png", 
#            dpi = 600)


# Individual plots
load(paste0("./data/all_experiments_filtered_asvs.RData"))

temp.meta <- (meta$datasets[["STREX"]])
temp.meta <- data.frame(meta$design_strex[match(temp.meta$SampleID, meta$design_strex$SampleID),]) %>% filter(replicate != 3)
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

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$name %in% as.character(temp.meta$Genotype))
idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

# Make factors
temp.meta$Genotype <- genotype$short[match(temp.meta$genotype, genotype$short)]

# Fetch taxa
taxa.mat <- all_exp_cut$taxa
taxa_level <- "lineage_family"

asv.ra <- apply(temp.mat, 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv <- apply(temp.mat[idx,-1], 2, function(x) x/sum(x))

idt <- taxa.mat$lineage_family[str_detect(taxa.mat$lineage_family, ";Dothideo|;Sordario|;Leotio")]
taxa.mat_temp <- as.data.frame(asv) %>% 
add_column(asv = row.names(.), .before = 1) %>%
add_column(tag = taxa.mat[match(row.names(.), taxa.mat$asvf), taxa_level], .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag, -asv) %>%
# filter(tag %in% ids) %>%
filter(tag %in% idt) %>%
data.frame(., stringsAsFactors = FALSE)

taxa.df <- taxa.mat_temp %>% 
cbind.data.frame(., temp.meta[match(.$SampleID, temp.meta$SampleID), c("Genotype", "Treatment", "Experiment", "Replicate", "tech_rep")])  %>%
cbind.data.frame(., taxa.mat[match(.$tag, taxa.mat$lineage_family), c("Phylum", "Class", "Order")]) %>%
mutate(new_Order = replace_na(Order, "Unclassified"))

taxa.df$new_Order <- ifelse(taxa.df$new_Order %in% names(col.pal), taxa.df$new_Order, "Rare_taxa")
taxa.df$new_Order <- factor(taxa.df$new_Order, levels = names(col.pal))

taxa.df %>%
  mutate(RA = ifelse(!is.finite(log10(RA)), log10(RA+0.0001), log10(RA))) %>%
  ggplot(aes(x = Genotype, y = RA)) +
  # ggtitle(paste0(x, "; no. ASVs = ", n_asvs)) +
  # geom_hline(yintercept = c(mark, mark_exudate, mark_extract), 
  #   colour = c_dark_grey, 
  #   lty = "solid", 
  #   lwd = 1) +
  geom_point(
    size = 0.1, 
    alpha = 0.2, 
    aes(shape = Replicate, colour = Genotype),
    position = position_jitterdodge(jitter.width=1)) +
  geom_boxplot(
    size = 0.8, 
    alpha = 0.8, 
    outlier.alpha = 0,
    aes(colour = Genotype),
    fill = c_white) +
  facet_wrap(new_Order ~ Treatment, switch = "x", scales = "free") +
  scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(taxa.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_blank(), 
          strip.text.x = element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="log10(RA)", x="", colour="", shape = "") +
    ggsave(filename = paste0(figs.out,"/Targetted_boxplot.png"),
           bg="transparent", 
           width=7, 
           height=8,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

# # Plot individual orders
mclapply(unique(as.character(taxa.df$new_Order)), function(x){

	# Subset the data
  temp.df <- taxa.df %>% filter(new_Order == x) %>%
  mutate(logRA = ifelse(!is.finite(log10(RA)), log10(RA+0.0001), log10(RA)))

  # Plot the distibution
  mark_extract <- temp.df %>% filter(Genotype == "Col", Treatment == "Extract") %>%
  summarise(sum = mean(RA)) %>%
  .$sum

  mark_exudate <- temp.df %>% filter(Genotype == "Col", Treatment == "Exudate") %>%
  summarise(sum = mean(RA)) %>%
  .$sum

  mark <- temp.df %>% filter(Genotype == "Soil") %>%
  summarise(sum = mean(RA)) %>%
  .$sum

  x.fig.out <- paste0(figs.out, "/", x)
  n_asvs <- length(unique(temp.df$asv))

  if(!dir.exists(paths = x.fig.out)){

        message(paste0("Directory created ", x.fig.out))
      dir.create(x.fig.out, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x.fig.out))

    }
    temp.df$Genotype <- factor(temp.df$Genotype, levels = genotype$short[which(genotype$short %in% temp.df$Genotype)])

  # Plot Canvas
  temp.df %>%
  ggplot(aes(x = Genotype, y = logRA)) +
  ggtitle(paste0(x, "; no. ASVs = ", n_asvs)) +
  geom_hline(yintercept = c(mark, mark_exudate, mark_extract), 
    colour = c_dark_grey, 
    lty = "solid", 
    lwd = 1) +
  geom_point(
    size = 0.1, 
    alpha = 0.2, 
    aes(shape = Replicate, colour = Genotype),
    position = position_jitterdodge(jitter.width=1)) +
  geom_boxplot(
    size = 0.8, 
    alpha = 0.8, 
    outlier.alpha = 0,
    aes(colour = Genotype),
    fill = c_white) +
  facet_grid(.~ Treatment, switch = "x", scales = "free_x", space = "free_x") +
  scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(temp.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          # panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 4, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="log10(RA)", x="", colour="", shape = "") +
    ggsave(filename = paste0(x.fig.out,"/", x, "_", n_asvs, "_boxplot.png"),
           bg="transparent", 
           width=5, 
           height=3.5,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

    # Plot Canvas
  temp.df %>%
  group_by(Genotype, Treatment, Experiment, SampleID, tag) %>%
  summarise(RA = sum(RA)) %>%
  ggplot(aes(x = SampleID, y = RA*100)) +
  geom_bar(
    stat = "identity", 
    position = "stack",
    aes(fill = tag),
    colour = c_grey) +
  facet_grid(.~ Treatment + Genotype, switch = "x", scales = "free_x", space = "free_x") +
  # scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(temp.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  scale_fill_brewer(palette = "Spectral") +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          # panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_blank(), 
          legend.position = "right",
          strip.text.x = element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="RA%", x="", colour="", shape = "") +
    ggsave(filename = paste0(x.fig.out,"/", x, "_", "_barplot.png"),
           bg="transparent", 
           width=8, 
           height=4,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

  # Plot Canvas
  temp.df %>%
  group_by(Genotype, Treatment, Experiment, SampleID, tag) %>%
  summarise(RA = sum(RA)) %>%
  ungroup(.) %>%
  group_by(Genotype, Treatment, tag) %>%
  summarise(RA = mean(RA)) %>%
  ggplot(aes(x = Genotype, y = RA*100)) +
  geom_bar(
    stat = "identity", 
    position = "stack",
    aes(fill = tag),
    colour = c_grey) +
  facet_grid(.~ Treatment, switch = "x", scales = "free_x", space = "free_x") +
  # scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(temp.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  scale_fill_brewer(palette = "Spectral") +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          # panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5), 
          legend.position = "right",
          strip.text.x = element_text(angle = 90, size = 6, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="RA%", x="", colour="", shape = "", fill = "") +
    ggsave(filename = paste0(x.fig.out,"/", x, "_", "_mean_RA_barplot.png"),
           bg="transparent", 
           width=5, 
           height=3,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

}, mc.cores = 4)

# Ternary plots
tern.df <- taxa.df %>%
filter(Genotype %in% c("Col", "cyp", "pyk10")) %>%
group_by(tag, asv, Treatment, Experiment, Replicate, tech_rep, Genotype) %>%
summarise(RA = (sum(RA))) %>%
# mutate(RA = ifelse(!is.finite(RA),log10(0.0001), RA)) %>%
spread(key = Genotype, value = RA, fill = NA) %>%
data.frame(., stringsAsFactors = FALSE)

# Plot
require(ggtern)
tern.df  %>% 
filter((Col+cyp+pyk10)!=0) %>%
mutate(behaviour = as.factor(ifelse((pyk10 > Col) & (cyp > Col), 1, 2))) %>%
        ggtern::ggtern(aes(x = cyp, y = Col, z = pyk10), na.rm = FALSE) +
        # ggtern::tern_limits(T=0.1, L=10, R=4) +
        geom_point(aes(size = Col, fill = behaviour), colour = c_black, alpha = 0.4, shape = 21,
        # position = position_jitter_tern(), 
        show.legend = TRUE, na.rm = FALSE) +
        facet_wrap(.~Treatment, scale = "free") +
        scale_size_continuous(range = c(0.3, 3)) +
        # scale_fill_manual(values = c(`1` = c_blue, `2` = c_cudo_magenta, `3` = c_red, `4` = c_dark_green), guide = FALSE) +
        scale_fill_manual(values = c(`1` = c_dark_green, `2` = c_cudo_magenta), guide = FALSE) +
        # scale_alpha_manual(values = c(as.numeric(unique(temp$alpha))), guide = FALSE) +
        labs(size = "") +
        theme_bw() +
        theme(legend.position      = c(0,1),
              legend.justification = c(0, 1)) +
        ggtern::ggsave(filename = paste0(figs.out,"/ternary_diagram_ra_targetted.png"),
               width=8, height=4, units = "in", 
               dpi = 300, 
               device = "png", limitsize = FALSE)



# END OF SCRIPT
sessionInfo()

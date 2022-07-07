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
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "RColorBrewer", "parallel", "cowplot")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))
load(paste0("./data/CAS_DE_family_E.RData"))

# REad the taxonomy table from the ATSPhereCC
# taxa_at_cc <- read.table("./output/taxonomy_syncom.txt", sep = "\t", header = TRUE, as.is = TRUE)

# What are we anlysing
pro <- c("gh_silva")

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

# Filter criteria
cf <- -log10(alpha)
ra_cf <- c(log2(0.05/100), log2(0.05/100), log2(0.05/100))
ra_cf10 <- c(log10(0.05/100), log10(0.5/100), log10(5/100))
lfc_cf <- 2

# Group out the cluster informations

# Plot Canvas
set.seed(seeder)

temp.meta <- rbind.data.frame(meta$datasets[["CAS"]]) %>%
filter(Compartment != "rhizosphere", !SampleID %in% outlier)
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$Replicate))
temp.meta$Experiment <- factor(as.character(temp.meta$Experiment))
temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# Use short terms for factoring
temp.meta$Genotype <- genotype$short[match(temp.meta$Genotype, genotype$names)]
temp.meta$Compartment <- compartment$short[match(temp.meta$Compartment, compartment$names)]
temp.meta$group <- paste(temp.meta$Genotype, temp.meta$Compartment, temp.meta$Soil_Batch, sep = "_")

# LFC matrix
lfc <- CAS_DE_family$logFC_P[str_detect(colnames(CAS_DE_family$logFC_P), "_E_")]

# Subset of ASVs that are prevalant
temp.meta_RP <- temp.meta %>% filter(Compartment == "E")
id_RP <- as.character(temp.meta_RP$SampleID)
idx <- temp.meta_RP$Soil_Batch == "CAS11"
IDs_11 <- temp.meta_RP$SampleID[idx]


# Relative abundance table
asv.ra <- apply(temp.mat[,id_RP], 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv.ra <- apply(temp.mat[idx,id_RP], 2, function(x) x/sum(x))
idx <- colnames(asv.ra) %in% IDs_11
ra_RP_cut_11 <- asv.ra[, idx]
ra_RP_cut_13 <- asv.ra[, !idx]

# Compute prevelance
idx <- (rowSums(ra_RP_cut_11 > 0) > ncol(ra_RP_cut_11)/2) & (rowSums(ra_RP_cut_13 > 0) > ncol(ra_RP_cut_13)/2)
asv.prev <- asv.ra[idx,]
prev_taxa <- unique(all_exp_cut$taxa[match(row.names(asv.prev), all_exp_cut$taxa$asvb), "lineage_family"])


# Box-plot of the prev ASVs
df.prev <- as.data.frame(asv.prev) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$lineage_family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
group_by(tag, SampleID) %>%
summarise(RA = sum(RA)) %>%
ungroup(.) %>%
mutate(
  Soil_Batch = temp.meta$Soil_Batch[match(.$SampleID, temp.meta$SampleID)],
  Genotype = temp.meta$Genotype[match(.$SampleID, temp.meta$SampleID)],
) %>%
mutate(
  Soil_Batch = factor(Soil_Batch, 
    levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)]),
  Genotype = factor(Genotype, 
    levels = genotype$short[which(genotype$short %in% .$Genotype)]),
  tag = as.factor(tag)
  ) %>%
# group_by(tag, group) %>%
# summarise(cRA = mean(RA)) %>%
# spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

p <- df.prev %>%
ggplot(aes(x = Genotype, y = log10(RA+0.0001))) +
geom_point(aes(colour = Genotype), size = 0.5, shape = 2, position = position_jitterdodge(1)) +
geom_boxplot(aes(colour = Genotype), fill = c_white, alpha = 0.7, size = 1, outlier.alpha = 0) +
facet_grid(Soil_Batch ~ tag, scale = "fixed", space = "fixed", switch = "both") +
scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(df.prev$Genotype))]) +
theme_RTN +
theme(
  panel.spacing = unit(0.1, "lines"), 
  axis.text.x = element_text(size = 6, angle = 90, hjust = 0.9, vjust = 0.5),
  axis.text.y = element_text(size = 10), 
  strip.text.x = element_text(size = 2, angle = 90)) +
labs(x = "", y = " RA (%)")
ggsave(p, file = paste0(figs.out, "/RA_prevalant_E_family_level.png"), 
            dpi = 600, 
            units = img, 
            device = "png", 
            bg = "transparent", 
            width = 20, 
            height = 5, 
            limitsize = T)


idx <- (row.names(lfc) %in% prev_taxa)
lfc.prev <- lfc[idx, ]

ra.prev <- as.data.frame(asv.prev) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$lineage_family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = temp.meta$group[match(.$SampleID, temp.meta$SampleID)]) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(ra.prev)[-1] <- paste0("meanra_", colnames(ra.prev[,-1]), "_RA")
row.names(ra.prev) <- ra.prev$tag


ra <- as.data.frame(asv.ra) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$lineage_family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = temp.meta$group[match(.$SampleID, temp.meta$SampleID)]) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(ra)[-1] <- paste0("meanra_", colnames(ra[,-1]), "_RA")
row.names(ra) <- ra$tag

# Top 20 families
top_20 <- data.frame(mean_ra = rowSums(ra[,-1])) %>%
mutate(tag = row.names(ra)) %>%
arrange(desc(mean_ra)) %>%
.$tag %>% head(., n=20)

# box plot for top20s
df.t20 <- as.data.frame(asv.ra) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$lineage_family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
group_by(tag, SampleID) %>%
summarise(RA = sum(RA)) %>%
ungroup(.) %>%
mutate(
  Soil_Batch = temp.meta$Soil_Batch[match(.$SampleID, temp.meta$SampleID)],
  Genotype = temp.meta$Genotype[match(.$SampleID, temp.meta$SampleID)],
) %>%
mutate(
  Soil_Batch = factor(Soil_Batch, 
    levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)]),
  Genotype = factor(Genotype, 
    levels = genotype$short[which(genotype$short %in% .$Genotype)]),
  tag = as.factor(tag)
  ) %>%
filter(tag %in% top_20) %>%
# spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

p <- df.t20 %>%
ggplot(aes(x = Genotype, y = log10(RA+0.0001))) +
geom_point(aes(colour = Genotype), size = 0.5, shape = 2, position = position_jitterdodge(1)) +
geom_boxplot(aes(colour = Genotype), fill = c_white, alpha = 0.7, size = 1, outlier.alpha = 0) +
facet_grid(Soil_Batch ~ tag, scale = "free", space = "free_y", switch = "both") +
scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(df.prev$Genotype))]) +
theme_RTN +
theme(
  panel.spacing = unit(0.1, "lines"), 
  axis.text.x = element_text(size = 6, angle = 90, hjust = 0.9, vjust = 0.5),
  axis.text.y = element_text(size = 10), 
  strip.text.x = element_text(size = 2, angle = 90)) +
labs(x = "", y = " RA (%)")
ggsave(p, file = paste0(figs.out, "/RA_t20E_family_level.png"), 
            dpi = 600, 
            units = img, 
            device = "png", 
            bg = "transparent", 
            width = 20, 
            height = 5, 
            limitsize = T)

# Barplot for rank abundance
ra_table <- ra %>%
gather(key = "key", value = "vals", convert = FALSE, -tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "x2")) %>%
select(-x, -x1,-x2) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)])
	) %>%
filter(tag %in% top_20) %>%
mutate(tag = factor(tag, levels = top_20)) %>%
group_by(tag, Genotype) %>%
summarise(meanRA = mean(vals)) %>%
data.frame(., stringsAsFactors = FALSE)

p <- ra_table %>%
ggplot(aes(x = tag, y = 100*meanRA)) +
geom_bar(stat = "identity", position = "dodge", aes(fill = Genotype), colour = c_dark_grey) +
scale_fill_manual(values = genotype$colours[which(genotype$short %in% ra_table$Genotype)]) +
theme_RTN +
theme(
	panel.spacing = unit(0.5, "lines"), 
	axis.text.x = element_text(size = 6, angle = 90, hjust = 0.9, vjust = 0.5),
	axis.text.y = element_text(size = 10), 
	axis.title = element_text(size = 10)) +
labs(x = "Top 20 abundant families", y = "mean RA (%)")
ggsave(p, file = paste0(figs.out, "/rank_abundance_top20_family_E_family_level.png"), 
            dpi = 600, 
            units = img, 
            device = "png", 
            bg = "transparent", 
            width = 7, 
            height = 8, 
            limitsize = T)


# Heatmap to show the prevelant ASVs
ra.prev.hmap <- as.data.frame(asv.prev) %>% 
add_column(tag = row.names(.), .before = 1) %>%
mutate(tag = (all_exp_cut$taxa$lineage_family[match(row.names(.), all_exp_cut$taxa$asvb)])) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = paste("mean", temp.meta$group[match(.$SampleID, temp.meta$SampleID)], "ra", sep = "_")) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE, row.names = .$tag)

df_hmap <- cbind.data.frame(lfc.prev, ra.prev.hmap[row.names(lfc.prev), ]) %>%
gather(key = "key", value = "vals",-tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "value")) %>%
select(-x, -x1) %>%
spread(key = value, value = vals, fill = NA) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)]),
	sig = PValue <= 0.05,
	higher_tag = replace_na(all_exp_cut$taxa$Class[match(.$tag, all_exp_cut$taxa$lineage_family)], "Unclassified")
) %>%
# filter(Family != "Unclassified") %>%
data.frame(., stringsAsFactors = FALSE)

idt <- all_exp_cut$taxa$lineage_family[str_detect(all_exp_cut$taxa$lineage_family, "__Burkhol|__Xantho|__Rhizo|__Baci")]
taxa.mat_temp <- as.data.frame(asv.prev) %>% 
add_column(asv = row.names(.), .before = 1) %>%
add_column(tag = all_exp_cut$taxa[match(row.names(.), all_exp_cut$taxa$asvb), "lineage_family"], .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag, -asv) %>%
# filter(tag %in% ids) %>%
filter(tag %in% idt) %>%
data.frame(., stringsAsFactors = FALSE)

taxa.df <- taxa.mat_temp %>% 
cbind.data.frame(., temp.meta[match(.$SampleID, temp.meta$SampleID), c("Genotype", "Soil_Batch", "Experiment", "Replicate")])  %>%
cbind.data.frame(., all_exp_cut$taxa[match(.$tag, all_exp_cut$taxa$lineage_family), c("Phylum", "Class", "Order")]) %>%
mutate(new_Order = replace_na(Order, "Unclassified"))

taxa.df$new_Order <- ifelse(taxa.df$new_Order %in% names(col.pal), taxa.df$new_Order, "Rare_taxa")
taxa.df$new_Order <- factor(taxa.df$new_Order, levels = names(col.pal))

# Loop for the individual family barplots
mclapply(unique(as.character(taxa.df$new_Order)), function(x){

  # Subset the data
  temp.df <- taxa.df %>% filter(new_Order == x) %>%
  mutate(logRA = ifelse(!is.finite(log10(RA)), log10(RA+0.0001), log10(RA)))

  x.fig.out <- paste0(figs.out, "/", x)
  n_asvs <- length(unique(temp.df$asv))

  if(!dir.exists(paths = x.fig.out)){

        message(paste0("Directory created ", x.fig.out))
      dir.create(x.fig.out, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x.fig.out))

    }
  
  temp.df$Genotype <- factor(temp.df$Genotype, levels = genotype$short[which(genotype$short %in% temp.df$Genotype)])
  temp.df$Soil_Batch <- factor(temp.df$Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% temp.df$Soil_Batch)])


  # Plot Canvas
  p <- temp.df %>%
  group_by(Genotype, Soil_Batch, SampleID, tag) %>%
  summarise(RA = sum(RA)) %>%
  ggplot(aes(x = SampleID, y = RA*100)) +
  geom_bar(
    stat = "identity", 
    position = "stack",
    aes(fill = tag),
    colour = c_grey) +
  facet_grid(.~ Soil_Batch + Genotype, switch = "x", scales = "free_x", space = "free_x") +
  # scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(temp.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  scale_fill_brewer(palette = "Spectral") +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_blank(), 
          legend.position = "right",
          strip.text.x = element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="RA%", x="", colour="", shape = "")
    ggsave(p, filename = paste0(x.fig.out,"/", x, "_", "_barplot_E.png"),
           bg="transparent", 
           width=8, 
           height=4,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

  # Plot Canvas
  p <- temp.df %>%
  group_by(Genotype, Soil_Batch, Experiment, SampleID, tag) %>%
  summarise(RA = sum(RA)) %>%
  ungroup(.) %>%
  group_by(Genotype, Soil_Batch, tag) %>%
  summarise(RA = mean(RA)) %>%
  ggplot(aes(x = Genotype, y = RA*100)) +
  geom_bar(
    stat = "identity", 
    position = "stack",
    aes(fill = tag),
    colour = c_grey) +
  facet_grid(.~ Soil_Batch, switch = "x", scales = "free_x", space = "free_x") +
  # scale_colour_manual(values = genotype$colours[which(genotype$short %in% as.character(temp.df$Genotype))]) +
  scale_shape_manual(values = c(1:10)) +
  scale_fill_brewer(palette = "Spectral") +
  theme_RTN +
  theme(panel.spacing = unit(0.1, "lines"), 
          panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 4, hjust = 0.5, vjust = 0.5), 
          legend.position = "right",
          strip.text.x = element_text(angle = 90, size = 6, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 12),
          strip.text.y = element_text(size = 4, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 4, hjust = 0.5, vjust = 0.5)) +
    labs(y="RA%", x="", colour="", shape = "", fill = "")
    ggsave(p, filename = paste0(x.fig.out,"/", x, "_", "_mean_RA_barplot_E.png"),
           bg="transparent", 
           width=5, 
           height=3,
           units = img, 
           limitsize=F, 
           device = "png", 
           dpi = 600)

}, mc.cores = 4)

p <- df_hmap %>%
filter(Genotype != "Col") %>%
mutate(
  behaviour = as.factor(ifelse(!is.na(logFC), sign(logFC), -1)),
  ra = log10(ra),
  logFC = ifelse(abs(logFC) > 4, sign(logFC)*4, logFC)) %>%
ggplot(aes(x = Genotype, y = tag)) +
geom_point(alpha=0.8, aes(shape = behaviour, size = ra, colour = sig, fill = logFC), na.rm = FALSE) +
# geom_tile(alpha=1, aes(colour = sig), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(higher_tag ~ Soil_Batch, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_cudo_magenta, 
        mid = c_white,
        low = c_dark_green, 
        na.value = c_black,
        midpoint = 0, 
        breaks = c(-4, -2, 0, 2, 4), 
        limits = c(-4,4), 
        labels = paste0(c("<4", "-2", "0", "2", ">4"))) +
scale_shape_manual(values = c(`1` = 25, `-1` = 24)) +
scale_size_continuous(range = c(0.5,2)) +
scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_dark_grey)) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      axis.text.y=element_text(size = 2),
      strip.text.y = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/prevelant_heatmap_family_E.png"),
       bg="transparent", 
       width=3.5, 
       height=8,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

p <- df_hmap %>%
filter(tag %in% top_20) %>%
filter(Genotype != "Col") %>%
mutate(
  behaviour = as.factor(ifelse(!is.na(logFC), sign(logFC), -1)),
  ra = log10(ra+0.001),
  logFC = ifelse(abs(logFC) > 4, sign(logFC)*4, logFC)) %>%
ggplot(aes(y = Genotype, x = tag)) +
geom_point(alpha=0.8, aes(fill = logFC), size = 3, shape = 21, colour = c_black, na.rm = FALSE) +
# geom_tile(alpha=0.5, aes(fill = ra), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(Soil_Batch ~ higher_tag, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_dark_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = c_dark_grey,
        midpoint = 0, 
        breaks = c(-4, -2, 0, 2, 4), 
        limits = c(-4,4), 
        labels = paste0(c("≤-4", "-2", "0", "2", "≥4"))) +
# scale_shape_manual(values = c(`1` = 24, `-1` = 25)) +
# scale_size_continuous(range = c(1,4)) +
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
ggsave(p, filename = paste0(figs.out,"/prevelant_heatmap_t20family_E.png"),
       bg="transparent", 
       width=3.5, 
       height=7,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)


# BArplot
p <- df_hmap %>%
filter(tag %in% top_20) %>%
group_by(higher_tag, tag, Soil_Batch) %>%
summarise(RA = (sum(ra)), meanra = log10(mean(ra)+0.0001)) %>%
ggplot(aes(y = RA, x = tag)) +
geom_bar(position = "dodge", stat = "identity", aes(fill = Soil_Batch), colour = c_black) +
# geom_line(aes(y = meanra, x = tag, lty = Soil_Batch), colour = c_red) +
facet_grid(. ~ higher_tag, switch = "both", scales="free", space="free") +
scale_fill_manual(values = c(`CAS11` = c_white, `CAS13` = c_black)) +
scale_linetype_manual(values = c(`CAS11` = "dashed", `CAS13` = "solid")) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.y=element_text(angle = 0, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text.x=element_text(size = 4, angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/barplot_t20family_E.png"),
       bg="transparent", 
       width=3.5, 
       height=5,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

# LFC correlation
lfc_table <- lfc %>%
add_column(tag = row.names(.), .before = 1) %>%
gather(key = "key", value = "vals", convert = FALSE, -tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "value")) %>%
select(-x, -x1) %>%
spread(key = value, value = vals, fill = NA) %>%
select(-PValue) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)])
	) %>%
spread(key = Genotype, value = logFC, fill = NA) %>%
data.frame(., stringsAsFactors = FALSE)

# Correlagram
d <- (cor(lfc[,str_detect(colnames(lfc), "_logFC$") & str_detect(colnames(lfc), "^genotype")]))
d[lower.tri(d)] <- NA

cc_df <- as.data.frame(d) %>%
add_column(x = row.names(.)) %>%
gather(key = "y", value = "pcc", -x) %>%
separate(x, into = c("v1", "genotypeA", "v2", "Soil_BatchA", "v3"), sep = "_") %>%
separate(y, into = c("v4", "genotypeB", "v5", "Soil_BatchB", "v6"), sep = "_") %>%
select(genotypeA, Soil_BatchA, genotypeB, Soil_BatchB, pcc) %>%
mutate(
    genotypeA = factor(genotypeA, levels = genotype$short[which(genotypeA %in% genotype$short)]),
    genotypeB = factor(genotypeB, levels = genotype$short[which(genotypeB %in% genotype$short)]),
    Soil_BatchA = factor(Soil_BatchA, levels = c("CAS11", "CAS13")),
    Soil_BatchB = factor(Soil_BatchB, levels = c("CAS11", "CAS13"))
    ) %>%
data.frame(.)

# Plot the correlogram
p <- cc_df %>%
filter(Soil_BatchA == Soil_BatchB) %>%
ggplot(aes(x = genotypeA, y=genotypeB, fill = pcc)) +
geom_point(alpha=0.8, shape = 21, size = 10, colour = c_white, na.rm = TRUE) +
geom_text(aes(label = as.character(round(pcc, 2))), size = 3, colour = c_black, guide = FALSE) +
# geom_tile(alpha=1, aes(colour = sig), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(. ~ Soil_BatchA, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = NA,
        midpoint = 0, 
        breaks = c(-1, -0.5, 0, 0.5, 1), 
        limits = c(-1,1), 
        labels = paste0(c("<1", "-0.5", "0", "0.5", ">1"))
        ) +
theme_RTN +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      # strip.text.x = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      # axis.text.y=element_text(size = 2),
      # strip.text.y = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "", label = "")
ggsave(p, filename = paste0(figs.out,"/correlogram_heatmap_family_E.png"),
       bg="transparent", 
       width=3.5, 
       height= 3,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

# Pairwise Plot
for(i in c("nai1", "pyk10", "cyp", "myb")){

   for(j in c("nai1", "pyk10", "cyp", "myb")){

      temp <- lfc_table

      if(j != i && dim(temp)[1] > 2){

        # i <- "nai1"
        # j <- "pyk10"
      	message(paste0("combination: ", i, " and ", j))

      	temp$x <- temp[,i]
      	temp$y <- temp[,j]

      	stat <- cor.test(temp$x, temp$y) 
        til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", (stat$p.value))

        p <- temp %>%
        select(tag, Soil_Batch, x, y) %>%
        mutate(prevelant = tag %in% row.names(lfc.prev), note = tag %in% top_20) %>%
        ggplot(aes(x = x, y = y)) +
        ggtitle(paste0(paste(i,j, sep = "_"), til)) +
        geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
        geom_point(aes(
              fill = Soil_Batch, colour = prevelant), 
              shape = 21,  
              size = 2,
              na.rm = FALSE) +
        ggrepel::geom_text_repel(aes(label = note), 
            colour = c_black, 
            size = 2,  max.overlaps = 15,
            box.padding = unit(0.5, "lines"),
            point.padding = unit(0.3, "lines")
            ) +
        scale_fill_manual(values = as.character(soilbatch$colours[which(soilbatch$names %in% temp$Soil_Batch)])) +
        scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_grey)) +
        theme_RTN_MDS +
        theme(panel.spacing = unit(0.5, "lines"), plot.title = element_text(size=4)) +
        labs(x = i, y = j, fill = "")
        ggsave(p, file = paste0(figs.out, "/Parwise_E_family_level_", j,"_", i,"_logFC_tagged.png"), 
                    dpi = 600, 
                    units = img, 
                    device = "png", 
                    bg = "transparent", 
                    width = 4, 
                    height = 4, 
                    limitsize = T)

        p <- temp %>%
        select(tag, Soil_Batch, x, y) %>%
        mutate(prevelant = tag %in% row.names(lfc.prev), note = ifelse(tag %in% top_20, tag, "")) %>%
        ggplot(aes(x = x, y = y)) +
        ggtitle(paste0(paste(i,j, sep = "_"), til)) +
        geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
        geom_point(aes(
              fill = Soil_Batch, colour = prevelant), 
              shape = 21,  
              size = 2,
              na.rm = FALSE) +
        # ggrepel::geom_text_repel(aes(label = note), 
        #     colour = c_black, 
        #     size = 2,  max.overlaps = 15,
        #     box.padding = unit(0.5, "lines"),
        #     point.padding = unit(0.3, "lines")
        #     ) +
        scale_fill_manual(values = as.character(soilbatch$colours[which(soilbatch$names %in% temp$Soil_Batch)])) +
        scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_grey)) +
        theme_RTN_MDS +
        theme(panel.spacing = unit(0.5, "lines"), plot.title = element_text(size=4)) +
        labs(x = i, y = j, fill = "")
        ggsave(p, file = paste0(figs.out, "/Parwise_E_family_level_", j,"_", i,"_logFC.png"), 
                    dpi = 600, 
                    units = img, 
                    device = "png", 
                    bg = "transparent", 
                    width = 4, 
                    height = 4, 
                    limitsize = T)

		}

	}
}

# END OF SCRIPT
sessionInfo()

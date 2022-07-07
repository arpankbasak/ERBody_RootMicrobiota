# Summary of GLM, relative abundance and absolute abundance
# script by Arpan Kumar Basak

# Loading Dependencies and path
rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/parameters.R")
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/functions.R")
setwd(analysis.syncom)

# Loading required packages
pkgs <- c("tidyverse", "reshape2", "parallel", "broom", "RColorBrewer", "ggpubr", "vegan", "RColorBrewer", "ggtern", "patchwork")
lapply(pkgs, require, character.only = T)

# Load data
load("./data/GLM_analysis_spiked.RData")
# load("./data/sync_filtered_data.Rdata")
load("./data/sync_filtered_spike_data.Rdata")

# Absolute abundance
lfc_p <- de_analysis$lfc_p
db <- read.table(taxa.sync, header = TRUE, sep = "\t", as.is = TRUE, row.names = NULL)
db <- db[which(db[,1] %in% row.names(lfc_p)),] %>% group_by(row.names, Kingdom, Phylum, Class, Order, Family, Genus) %>%
summarise(n = n()) %>%
data.frame(., stringsAsFactors = FALSE, row.names = .$row.names)

db <- db[row.names(lfc_p),]


lfc_df <- lfc_p %>% 
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "key", value = "vals", -strain) %>%
separate(key, into = c("genotype", "dose", "time", "hypothesis", "x2"), convert = FALSE, sep = "_") %>%
spread(key = x2, value = vals, convert = FALSE) %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  significance = PValue < alpha & abs(logFC) > 1,
  family = sync.spiked.dat$taxa$family[match(.$strain, row.names(sync.spiked.dat$taxa))],
  hypothesis = factor(hypothesis, levels = c("genotype", "competency"))) %>%
data.frame(., stringsAsFactors = FALSE)

# Syncom spiked data
aa_mat <- sync.spiked.dat$cData[which(row.names(sync.spiked.dat$cData) %in% row.names(lfc_p)), which(colnames(sync.spiked.dat$cData) %in% sync.spiked.dat$metadata$SampleID)]
aa_mat <- aa_mat[row.names(lfc_p),]

aa_df <- as.data.frame(aa_mat) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(
  fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  AA = log2(AA+1)) %>%
group_by(strain, fixed) %>%
summarise(meanAA = mean(AA), sdAA = sd(AA)) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
mutate(
  genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  ) %>%
data.frame(., stringsAsFactors = FALSE)

aa_df$family <- db$Family[match(aa_df$strain, row.names(db))]

# Plot Canvas
# Heatmap -- AA
(fig_1a <- aa_df %>%
  filter(strain %in% unique(lfc_df$strain)) %>%
mutate(strain = factor(strain, levels = de_analysis$strain_clusters),
  # meanAA = ifelse(is.finite(log10(meanAA)), log10(meanAA), NA),
  family = as.factor(family),
  meanAA = ifelse(abs(meanAA) > 1, sign(meanAA) * 1, meanAA)) %>%
ggplot(aes(x = genotype, y = strain, fill = meanAA)) +
geom_raster(alpha=1) +
    scale_fill_gradient2(
      low=c_black,
      mid = c_yellow, 
      midpoint = 2, 
      high = c_cudo_skyblue, 
      na.value = c_black,
      breaks = c(0, 1, 2, 4, 6), 
      limits = c(0, 6), 
      labels = paste0(c("0","1", "2", "4", "6"))
      ) +
    facet_grid(family ~ dose + time, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 10),
          strip.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="mean(log2(reads/reads_spike))"))
    ggsave(fig_1a, filename = paste0("./figures/differential_analysis/spiked/hmap_AA_db.png"),
           bg="transparent", width=7, height=14, limitsize=F, device = "png", dpi = 600)

# HEatmap sample-wise aA
aa_df_all <- as.data.frame(aa_mat) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  AA = log2(AA+1)) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
separate(random, into = c("seq_batch", "experiment", "tech_replicate"), convert = FALSE, sep = "_") %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  seq_batch = as.factor(seq_batch),
  experiment = as.factor(experiment),
  tech_replicate = as.factor(tech_replicate)) %>%
data.frame(., stringsAsFactors = FALSE)

aa_df_all$family <- db$Family[match(aa_df_all$strain, row.names(db))]


# Selecting prevelant Strains
df_prev <- aa_df_all %>%
group_by(strain, dose, time) %>%
summarise(aa = (sum(AA))) %>%
data.frame(., stringsAsFactors = FALSE)

p <- df_prev %>%
ggplot(aes(x = aa)) +
geom_density() +
facet_wrap(dose ~ time, scales = "free") +
theme_RTN +
labs(x = "", y = "")

ggsave(p, 
  filename = paste0("./figures/differential_analysis/spiked/density_strains_db.png"),
           bg="transparent", units = img, 
           width=6, 
           height=3.5, limitsize=F, device = "png", dpi = 600)


df_prev <- aa_df_all %>%
group_by(strain) %>%
summarise(aa = (sum(AA))) %>%
data.frame(., stringsAsFactors = FALSE)

df_prev$cut_off <- (df_prev$aa/(sum(df_prev$aa))) > .01
strain_in <- df_prev$strain[df_prev$cut_off == TRUE]

p <- df_prev %>%
ggplot(aes(x = aa)) +
geom_density() +
# facet_wrap(dose ~ time, scales = "free") +
theme_RTN +
labs(x = "", y = "")

ggsave(p, 
  filename = paste0("./figures/differential_analysis/spiked/density_strains_all_db.png"),
           bg="transparent", units = img, 
           width=4, 
           height=2.5, limitsize=F, device = "png", dpi = 600)

(fig_1a1 <- aa_df_all %>%
filter(strain %in% unique(lfc_df$strain)) %>%
mutate(
  strain = factor(strain, levels = de_analysis$strain_clusters),
  family = as.factor(family),
  AA = ifelse(abs(AA) > 6, sign(AA) * 6, AA)) %>%
ggplot(aes(x = as.factor(SampleID), y = strain, fill = saturate(AA))) +
geom_raster(alpha=1) +
    scale_fill_gradient2(
      low=c_black,
      mid = c_yellow, 
      midpoint = 2, 
      high = c_cudo_skyblue, 
      na.value = c_black,
      breaks = c(0, 1, 2, 4, 6), 
      limits = c(0, 6), 
      labels = paste0(c("0","1", "2", "4", "6"))
      ) +
    facet_grid(family ~  seq_batch + experiment + dose + time + genotype + tech_replicate, 
      switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
          strip.text.y = element_text(size = 15, angle = 180, hjust = 0.5, vjust = 0.5),
          axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 10),
          legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="mean(log2(AA))"))
    ggsave(fig_1a1, filename = paste0("./figures/differential_analysis/spiked/hmap_ra_samplewise_db.png"),
           bg="transparent", units = img, width=25, height=16, limitsize=F, device = "png", dpi = 600)


# Heatmap -- LFC --genotype
(fig_1b <- lfc_df %>%
filter(hypothesis != "competency") %>%
mutate(strain = factor(strain, levels = de_analysis$strain_clusters),
  sig = as.factor(PValue < alpha),
  family = as.factor(family),
  logFC = ifelse(abs(logFC) > 2, sign(logFC) * 2, logFC)) %>%
ggplot(aes(x = genotype, y = strain, fill = saturate(logFC))) +
geom_raster(alpha=1) +
geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = sig)) +
    scale_fill_gradient2(
      low=c_cudo_magenta, 
      high = c_dark_green, 
      midpoint = 0.0, 
      mid = c_white,
      na.value = c_grey,
      breaks = c(-2, -1, 0, 1, 2), 
      limits = c(-3,3), 
      labels = paste0(c("-3","-1", "0", "1", "3"))
    ) +
    scale_colour_manual(values = c(`FALSE` = "transparent", `TRUE` = c_black), 
      labels = c("", paste0("FDR < ", alpha))) +
    facet_grid(family ~ dose + time, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 10),
          strip.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="log2(FC)vs Col-0 ", colour = ""))
    ggsave(fig_1b, filename = paste0("./figures/differential_analysis/spiked/hmap_lfc_db.png"), 
           bg="transparent", width=7, height=14, limitsize=F, device = "png", dpi = 600)

# Heatmap -- LFC --competency
(fig_1b1 <- lfc_df %>%
filter(hypothesis == "competency") %>%
mutate(strain = factor(strain, levels = de_analysis$strain_clusters),
  sig = as.factor(PValue < alpha),
  family = as.factor(family),
  logFC = ifelse(abs(logFC) > 2, sign(logFC) * 2, logFC)) %>%
ggplot(aes(x = genotype, y = strain, fill = saturate(logFC))) +
geom_raster(alpha=1) +
geom_tile(fill = NA, width = 0.95, height = 0.95, size = 0.6, aes(colour = sig)) +
    scale_fill_gradient2(
      low=c_cudo_magenta, 
      high = c_dark_green, 
      midpoint = 0.0, 
      mid = c_white,
      na.value = c_grey,
      breaks = c(-2, -1, 0, 1, 2), 
      limits = c(-3,3), 
      labels = paste0(c("-3","-1", "0", "1", "3"))
    ) +
    scale_colour_manual(values = c(`FALSE` = "transparent", `TRUE` = c_black), 
      labels = c("", paste0("FDR < ", alpha))) +
    facet_grid(family ~ dose + time, switch = "both", scales="free", space="free") +
    theme_RTN +
    theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA),
          axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
          strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 10),
          strip.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
          legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
    labs(y="", x="", fill="log2(FC)vs Col-0 ", colour = ""))
    ggsave(fig_1b1, filename = paste0("./figures/differential_analysis/spiked/hmap_lfc_competency_db.png"), 
           bg="transparent", width=7, height=14, limitsize=F, device = "png", dpi = 600)

# Strips for annotation
taxa_df <- sync.spiked.dat$taxa[match(row.names(sync.spiked.dat$taxa), row.names(lfc_p)),] %>%
mutate(class = factor(class, levels = taxonomy.default.colours$class[which(.$class %in% taxonomy.default.colours$class)]),
  family = as.factor(family),
  compartment = as.factor(compartment),
  strains = factor(row.names(.), levels = de_analysis$strain_clusters))

# Strip for root, soil and root and soil strains
(strip_comp <- taxa_df %>%
ggplot(aes(x = "", y = strains, fill = compartment)) +
geom_raster(alpha=1) +
scale_fill_manual(values = c(c_very_dark_green, c_dark_brown)) +
facet_grid(family ~ ., switch = "both", scales="free", space="free") +
theme_RTN +
theme(axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
      axis.text.y=element_text(size = 10),
      strip.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", colour = ""))
ggsave(strip_comp, filename = paste0("./figures/differential_analysis/spiked/strip_compartments_db.png"), 
           bg="transparent", width=3, height=14, limitsize=F, device = "png", dpi = 600)

# Strip for higher taxonomic rank -- class
(strip_taxa <- taxa_df %>%
ggplot(aes(x = "", y = strains, fill = class)) +
geom_raster(alpha=1) +
scale_fill_manual(values = taxonomy.default.colours$colour[which(as.character(taxonomy.default.colours$class) %in% as.character(taxa_df$class))]) +
facet_grid(family ~ ., switch = "both", scales="free", space="free") +
theme_RTN +
theme(axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
      axis.text.y=element_text(size = 10),
      strip.text.y = element_text(size = 2, angle = 180, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", colour = ""))
ggsave(strip_taxa, filename = paste0("./figures/differential_analysis/spiked/strip_taxa_db.png"), 
           bg="transparent", width=3, height=14, limitsize=F, device = "png", dpi = 600)


(strip_prev <- taxa_df %>%
mutate(mark = row.names(.) %in% strain_in) %>%
ggplot(aes(x = "", y = strains, fill = mark)) +
geom_raster(alpha=1) +
scale_fill_manual(values = c(`TRUE` = c_red, `FALSE` = c_white)) +
facet_grid(family ~ ., switch = "both", scales="free", space="free") +
theme_RTN +
theme(axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(angle = 90, size = 15, hjust = 0.5, vjust = 0.5),
      axis.text.y=element_text(size = 10),
      strip.text.y = element_text(size = 2, angle = 180, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 10, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", colour = ""))
ggsave(strip_prev, filename = paste0("./figures/differential_analysis/spiked/strip_prevelant_db.png"), 
           bg="transparent", width=3, height=14, limitsize=F, device = "png", dpi = 600)


# COmbined AA and LFC dataset for analysis
aa <- aa_df %>%
mutate(group = paste(genotype, dose, time, sep = "_")) %>%
select(strain, group, meanAA) %>%
spread(key = group, value = meanAA, convert = FALSE, fill = 0) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(aa) <- paste(colnames(aa), "aa", sep = "_")
lfc_p_g <- lfc_p[, str_detect(colnames(lfc_p), "_genotype_")]
colnames(lfc_p_g) <- str_replace_all(colnames(lfc_p_g), "_genotype_", "_")

df <- cbind.data.frame(lfc_p_g, aa[match(row.names(lfc_p_g), aa$strain_aa),]) %>%
gather(key = "x", value = "y", convert = FALSE, -strain_aa) %>%
separate(x, into = c("genotype", "dose", "time", "x"), sep = "_") %>%
spread(key = x, value = y, convert = FALSE, fill = 0) %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  # hypothesis = factor(hypothesis, levels = c("genotype", "competency")),
  family = sync.spiked.dat$taxa$family[match(aa_df$strain, row.names(sync.spiked.dat$taxa))],
  significance = PValue < alpha & abs(logFC) > 1,
  log_aa = aa
  ) %>%
data.frame(., stringsAsFactors = FALSE)

df_wt <- df %>%
filter(genotype == "Col") %>%
group_by(strain_aa, genotype, dose, time) %>%
summarise(mean_aa = mean(aa)) %>%
ungroup() %>%
mutate(group = paste(dose, time, sep = "_")) %>%
select(strain_aa, group, mean_aa) %>%
spread(key = group, value = mean_aa, convert = FALSE, fill = 0) %>%
data.frame(., stringsAsFactors = FALSE)

# Set X and Y range for Col-0
marker_x <- mean(df$log_aa[df$genotype == "Col"])
marker_y <- mean(df$logFC[df$genotype == "Col"])

(sct_plot <- df %>%
  mutate(text_target = strain_aa) %>%
  mutate(text_target = ifelse(significance == TRUE, as.character(text_target), "")) %>%
  ggplot(aes(x = log_aa, y = logFC)) +
  geom_hline(yintercept = c(-1,1), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_vline(xintercept = c(0), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_point(aes(fill = family, 
    # shape = genotype, 
    size = significance, 
    colour = significance, 
    alpha = significance), shape = 22) +
  # xlim(c(-21, 10)) +
  # ylim(c(-6, 6)) +
  ggrepel::geom_text_repel(aes(label = text_target)) +
  facet_grid(genotype ~ dose + time, 
    scale = "fixed", 
    switch = "both", 
    space = "free") +
  # scale_fill_brewer(palette = "Spectral") +
  # scale_shape_manual(values = c(22,23,24)) +
  scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
  theme_RTN +
   theme(legend.position = "top", 
    axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
   labs(x = "log2(read_strain/reads_spike)", y = "LFC", 
    colour = "", shape = "", fill = "", size = "")
  )
ggsave(sct_plot, "./figures/differential_analysis/spiked/pairwise_integrated_plot_db.png", 
    units = "in", 
    width = 12, 
    height = 8, 
    bg = "transparent")

# Correlation plots
(sct_plot_geno <- df %>%
  filter(!genotype %in% c("Col", "inoc")) %>%
  # group_by(strain_aa, dose, time, genotype) %>%
  select(-log_aa, -aa, -PValue) %>%
  spread(key = genotype, value = logFC, fill = 0) %>%
  mutate(text_target = strain_aa) %>%
  mutate(text_target = ifelse(significance == TRUE, as.character(text_target), "")) %>%
  ggplot(aes(x = cyp, y = pyk10)) +
  geom_hline(yintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
  # geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_smooth(method = glm, se = FALSE, colour = c_black, linewidth = 0.5) +
  geom_point(aes(fill = family, 
    # shape = genotype, 
    size = significance, 
    colour = significance, 
    alpha = significance), shape = 22) +
  # xlim(c(-6, 6)) +
  # ylim(c(-6, 6)) +
  ggrepel::geom_text_repel(aes(label = text_target)) +
  facet_grid(. ~ dose + time, 
    scale = "fixed", 
    switch = "both", 
    space = "free") +
  # scale_fill_brewer(palette = "Spectral") +
  # scale_shape_manual(values = c(22,23,24)) +
  scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
  theme_RTN +
   theme(
    legend.position = "top", 
    legend.text = element_text(size = 2),
    axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
   labs(x = "LFC_cyp", y = "LFC_pyk10", 
    colour = "", shape = "", fill = "", size = "")
  )
ggsave(sct_plot_geno "./figures/differential_analysis/spiked/pairwise_integrated_genotype_plot_db.png", 
    units = "in", 
    width = 7, 
    height = 5, 
    bg = "transparent")

(sct_plot_geno_aa <- df %>%
  filter(!genotype %in% c("Col", "inoc")) %>%
  select(-aa, -logFC, -PValue) %>%
  spread(key = genotype, value = log_aa, fill = 0) %>%
  mutate(text_target = strain_aa) %>%
  mutate(text_target = ifelse(significance == TRUE, as.character(text_target), "")) %>%
  ggplot(aes(x = cyp, y = pyk10)) +
  geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_vline(xintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_point(aes(fill = family, 
    # shape = genotype, 
    size = significance, 
    colour = significance, 
    alpha = significance), shape = 22) +
  # xlim(c(-21, 10)) +
  # ylim(c(-21, 10)) +
  ggrepel::geom_text_repel(aes(label = text_target)) +
  facet_grid(. ~ dose + time, 
    scale = "fixed", 
    switch = "both", 
    space = "free") +
  # scale_fill_brewer(palette = "Spectral") +
  # scale_shape_manual(values = c(22,23,24)) +
  scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
  theme_RTN +
   theme(
    legend.position = "top", 
    legend.text = element_text(size = 2),
    axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
   labs(x = "AA_cyp", y = "AA_pyk10", 
    colour = "", shape = "", fill = "", size = "")
  )
ggsave(sct_plot_geno_aa, "./figures/differential_analysis/spiked/pairwise_integrated_aa_genotype_plot_db.png", 
    units = "in", 
    width = 7, 
    height = 5, 
    bg = "transparent")

# Pairwise plot
df_temp <- lfc_df %>%
  filter(!genotype %in% c("Col", "inoc"), hypothesis == "genotype") %>%
  select(-PValue) %>%
  spread(key = genotype, value = logFC, fill = 0) %>%
  mutate(text_target = strain) %>%
  mutate(text_target = ifelse(significance == TRUE, as.character(text_target), ""),
    group = paste(dose, time, sep = "_"))

mclapply(as.character(unique(df_temp$group)),function(x){

  temp <- df_temp %>% filter(group == x, strain %in% strain_in) 
  # %>%
  # mutate(mean_wt = df_wt[match(.$strain_aa, df_wt$strain), 
  # 	which(str_replace_all(colnames(df_wt), "X", "") %in% x)]) %>%
  # # mutate(wt_size = cut(mean_wt, breaks = c(0, 0.01, 0.1, 1), labels = c(".01", ".1", "1")))

  stat <- cor.test((temp$cyp), (temp$pyk10), method = "pearson")
  til <- paste0("Spearman: ", round(stat$estimate, 4), "; p-value: ", c(stat$p.value))
  
  p <- temp %>%
  mutate(
    cyp = ifelse(abs(cyp) >= 1, sign(cyp)*1, cyp), 
    pyk10 = ifelse(abs(cyp) >= 1, sign(pyk10)*1, pyk10)) %>%
  ggplot(aes(x = (cyp), y = (pyk10))) +
  geom_hline(yintercept = c(0), lty = "solid", lwd = 1, colour = "darkgrey") +
  geom_vline(xintercept = c(0), lty = "solid", lwd = 1, colour = "darkgrey") +
  ggtitle(til) +
  geom_smooth(method = glm, se = FALSE, colour = c_black, lwd = 0.5) +
  geom_point(aes(fill = family, 
    # shape = genotype, 
    size = wt_size, 
    colour = significance, 
    alpha = significance), shape = 23) +
  # xlim(c(-6, 6)) +
  # ylim(c(-6, 6)) +
  # ggrepel::geom_text_repel(aes(label = text_target)) +
  # facet_grid(. ~ dose + time, 
  #   scale = "fixed", 
  #   switch = "both", 
  #   space = "free") +
  # scale_fill_brewer(palette = "Spectral") +
  # scale_shape_manual(values = c(22,23,24)) +
  scale_size_manual(values = c(`.01` = 1, 
  	`.1` = 2,
  	`1` = 3), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.8, `TRUE` = 1), guide = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
  theme_RTN_MDS +
   theme(
    legend.position = "top", 
    legend.text = element_text(size = 2),
    plot.title = element_text(size=4)
    ) +
   labs(x = "LFC_cyp", y = "LFC_pyk10", 
    colour = "", shape = "", fill = "", size = "")
ggsave(plot = p, paste0("./figures/differential_analysis/spiked/pairwise_integrated_genotype_", x,"_plot_db.png"), 
    units = "in", 
    device = "png",
    width = 4, 
    height = 5, 
    bg = "transparent")
  
  return(NULL)

}, mc.cores = 8)

# Spill Canvas in one aisle
mat_plot <- cowplot::plot_grid(
  
  # Features
  strip_taxa + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()), 

  strip_prev + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()),
  
  strip_comp + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()), 

  fig_1a + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()), 
  
  fig_1b + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()), 
  
  fig_1b1 + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()), 

  ncol=6,
  nrow = 1, 
  align='hv', 
  axis = 'tblr',
  labels=NULL, 
  # byrow = TRUE,
  rel_widths = c(.1, .1, .1, 1.5, 1, 1), 
  rel_heights = c(1, 1, 1, 1, 1, 1)
)
ggsave(mat_plot, filename = paste0(figs, "/differential_analysis/spiked//hmap_de_summary_db.png"), 
           bg="transparent", 
           units = "in",
           width=5, height=12, 
           limitsize=F, device = "png", dpi = 600)

# Spill Canvas in one aisle -- HUGE
mat_plot <- cowplot::plot_grid(
  
  # Features
  strip_taxa + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 

  strip_prev + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()),

  strip_comp + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 
  
  # Abundance
  fig_1a1 + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 
  fig_1a + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 
  
  # Hypothesis testing
  fig_1b + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 

  fig_1b1 + theme_void() + 
  theme(panel.spacing = unit(0.1, "lines"), panel.border = element_rect(fill="transparent", colour=NA), 
    legend.position = "none", 
    axis.text.x = element_blank(),
    strip.text.y = element_blank()), 

  ncol=6,
  nrow = 1, 
  align='hv', 
  axis = 'tblr',
  labels=NULL, 
  # byrow = TRUE,
  rel_widths = c(0.05, 0.05, 0.05, 6, .4, 0.4, .4), 
  rel_heights = c(1, 1, 1, 1, 1, 1, 1)
)
ggsave(mat_plot, filename = paste0(figs, "/differential_analysis/spiked/hmap_AIO_summary_db.png"), 
           bg="transparent", 
           units = "in",
           width=25, height=14, 
           limitsize=F, device = "png", dpi = 600)


# MAke ternary plot
sig_strain <- unique(df$strain_aa[df$significance == TRUE])
tern.df <- as.data.frame(aa_mat) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(
  fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)]) %>%
group_by(strain, fixed) %>%
summarise(AA = mean(AA)) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
filter(genotype != "inoc") %>%
spread(key = genotype, value = AA, fill = NA) %>%
mutate(
  dose = as.factor(dose),
  time = as.factor(time),
  sig = strain %in% sig_strain
  ) %>%
data.frame(., stringsAsFactors = FALSE)

# Plot
# require(ggtern)
# tern.df  %>% 
#         ggtern::ggtern(aes(x = cyp, y = Col, z = pyk10), na.rm = FALSE) +
#         # ggtern::tern_limits(T=0.1, L=10, R=4) +
#         geom_point(aes(size = Col, fill = sig, alpha = sig), shape = 21, colour = c_black,show.legend = TRUE, na.rm = FALSE) +
#         facet_wrap(dose~time, scale = "free") +
#         scale_size_continuous(range = c(0.3, 3)) +
#         scale_fill_manual(values = c(`TRUE` = c_red, `FALSE` = c_dark_grey), guide = FALSE) +
#         scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.6), guide = FALSE) +
#         labs(size = "") +
#         theme_bw() +
#         theme(legend.position      = c(0,1),
#               legend.justification = c(0, 1)) +
#         ggtern::ggsave(filename = paste0(figs,"/differential_analysis/spiked/ternary_diagram_db.png"),
#                width=7, height=7, units = "in", 
#                dpi = 300, 
#                device = "png", limitsize = FALSE)


# Within Order barplots


# END OF SCRIPT
sessionInfo()
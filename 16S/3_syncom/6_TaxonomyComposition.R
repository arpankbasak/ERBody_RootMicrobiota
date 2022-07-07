# Summary of GLM, relative abundance and absolute abundance
# screipt by Arpan Kumar Basak
# TODO::Make changes in y axis scale to log10

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
load("./data/sync_filtered_data.Rdata")
load("./data/sync_filtered_spike_data.Rdata")

# Absolute abundance
lfc_p <- de_analysis$lfc_p

lfc_df <- lfc_p %>% 
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "key", value = "vals", -strain) %>%
separate(key, into = c("genotype", "dose", "time", "hypothesis", "x2"), convert = FALSE, sep = "_") %>%
spread(key = x2, value = vals, convert = FALSE) %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  significance = PValue < alpha & abs(logFC) > 1,
  family = sync.spiked.dat$taxa$family[match(.$strain, row.names(sync.spiked.dat$taxa))]) %>%
data.frame(., stringsAsFactors = FALSE)

# Syncom spiked data
aa_mat <- sync.spiked.dat$cData[which(row.names(sync.spiked.dat$cData) %in% row.names(lfc_p)), which(colnames(sync.spiked.dat$cData) %in% sync.spiked.dat$metadata$SampleID)]
in_set <- sync.spiked.dat$metadata %>% filter(exudates_rep != 3) %>%
.$SampleID

aa_df <- as.data.frame(aa_mat[,in_set]) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)]) %>%
# group_by(strain, fixed) %>%
# summarise(AA = sum(AA)) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time)) %>%
data.frame(., stringsAsFactors = FALSE)

aa_df$genus <- sync.spiked.dat$taxa$genus[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$family <- sync.spiked.dat$taxa$family[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$class <- sync.spiked.dat$taxa$class[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$phylum <- sync.spiked.dat$taxa$phylum[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]

plot_list <- parallel::mclapply(c("phylum", "class", "family", "genus"), function(i){

    set.seed(seeder)
    # i <- "class"
    temp <- aa_df %>% 
    mutate(x = aa_df[,i]) %>%
    group_by(x, SampleID) %>%
    summarise(AA = sum(AA)) %>%
    spread(key = SampleID, value = AA, fill = 0, convert = FALSE) %>%
    data.frame(., stringsAsFactors = FALSE, row.names = .$x) %>%
    select(-x)
    
    n <- length(unique(aa_df[,i]))

    col.pal <- colorRampPalette(RColorBrewer::brewer.pal("Set2", n = 5))(nrow(temp))
    (plot_obj <- as.data.frame(temp) %>%
    add_column(x = row.names(.)) %>%
    gather(key = "SampleID", value = "AA", convert = FALSE,-x) %>%
    cbind.data.frame(sync.spiked.dat$metadata[match(.$SampleID, sync.spiked.dat$metadata$SampleID), -1]) %>%
    filter(biol_rep != 1, exudates_rep != 3) %>%
    group_by(x, dose, time, genotype, exudates_rep) %>%
    summarise(AA = mean(AA)) %>%
    # mutate(AA = ifelse(is.finite(log10(AA)), log10(AA+1), log10(AA+1))) %>%
    ggplot(aes(x = exudates_rep, y = log10(AA+1), fill = x)) +
    geom_bar(stat = "identity", size=0) +
    facet_grid(. ~ dose + time + genotype, 
      switch = "x", 
      scales = "free", 
      space = "free", 
      labeller = label_parsed) +
    scale_fill_manual(values = col.pal) +
    labs(y = "log10(Quantitative abundance+1)", x = "", fill = "") + 
    theme_RTN +
    theme(panel.spacing = unit(0.01, "lines"), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text(angle = 90, size = 10, vjust = .5, hjust = .5),
          legend.text = element_text(size = 10),
          legend.title = element_blank())) +
    ggsave(file = paste0(figs, "./taxonomy_composition/spiked/AA_", i, "_level.png"), 
      dpi = 600, 
      device = "png", 
      bg = "transparent", 
      width = 8, height = 6, units = img, limitsize = F)

    return(plot_obj)

}, mc.cores = 4)
    
names(plot_list) <- c("phylum", "class", "family", "genus")

# Spill Canvas in one aisle
cowplot::plot_grid(
  
  plot_list$phylum,
  # theme(
  #   panel.spacing = unit(0.01, "lines"),
  #   strip.text.x = element_blank(), 
  #   axis.text.x = element_blank()),
  plot_list$class,
  # theme(
  #   panel.spacing = unit(0.01, "lines"),
  #   strip.text.x = element_blank(), 
  #   axis.text.x = element_blank()),
  plot_list$family,
  # theme(
  #   panel.spacing = unit(0.01, "lines"),
  #   strip.text.x = element_blank(), 
  #   axis.text.x = element_blank()),
  plot_list$genus,
  # theme(
  #   panel.spacing = unit(0.01, "lines"),
  #   strip.text.x = element_blank(), 
  #   axis.text.x = element_blank()),
  ncol=1,
  nrow = 4, 
  align='hv', 
  axis = 'tblr',
  labels=NULL, 
  # byrow = TRUE,
  rel_widths = c(1, 1, 1, 1), 
  rel_heights = c(1, 1, 1, 1)
) +
ggsave(filename = paste0(figs, "/taxonomy_composition/spiked/taxonomy_summary.png"), 
           bg="transparent", 
           units = "in",
           width=8, height=25, 
           limitsize=F, device = "png", dpi = 600)


# END OF SCRIPT
sessionInfo()
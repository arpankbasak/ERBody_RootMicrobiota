# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
pkgs <- c("tidyverse", "reshape2", "vegan")

lapply(pkgs, require, character.only = T)

source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/parameters.R")
setwd(analysis.syncom)

load("./data/sync_filtered_spike_data.Rdata")

# ===== Beta Diversity =====
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/functions.R")

mat <- sync.spiked.dat$cData
meta <- sync.spiked.dat$metadata

# Path
figs <- paste0(figs, "/alpha_diversity/spiked/")
stats <- paste0(stats, "/alpha_diversity/spiked/")
out <- paste0(out, "//alpha_diversity/spiked/")

# ===== Alpha Diversity =====

# Calculate Shannon Indices
dat_si <- data.frame(index = diversity(t(mat), index = "shannon"))

# Calculate Simpson Indices
dat_sm <- data.frame(index = diversity(t(mat), index = "simpson"))

# Calculate Observed ASVs Indices
# dat_fi <- data.frame(index = fisher.alpha(t(mat)))

# Calculate Observed ASVs Indices
dat_oa <- data.frame(index = colSums(mat))

# Setting plotting parameters
title_oa <- "Observed Strains"
title_si <- "Shannon Diversity Index"
title_smp <- "Simpson Index"
# title_fi <- "Fisher Alpha"

df <- cbind.data.frame(dat_si,
dat_sm,
# dat_fi,
dat_oa
)

colnames(df) <- c("shannon", "simpson", 
  # "fishers", 
  "observed_strains")

plot_list <- parallel::mclapply(1:ncol(df), function(i){

    # i = 1
    # Compute statistics
    temp <- cbind.data.frame(meta, Index = df[match(meta$SampleID, row.names(df)), i])

    f <- as.formula(Index ~ 0 + fixed + random)

    mod <- glm(f = f, data = temp, family = Gamma)
    aov.mod <- aov(mod)
    tsd.mod <- as.data.frame(TukeyHSD(aov.mod)$fixed)

    #mod <- wilcox.test(Index ~ group, data = meta)
    #aov.mod <- aov(mod)
    #tsd.mod <- TukeyHSD(aov.mod)$group

    write.table(tsd.mod, paste(stats, "glm_", colnames(df)[i], ".txt", sep = ""), sep = "\t")

    # Plot Canvas
    (plot_obj <- temp %>% mutate(
      biol_rep = as.factor(biol_rep),
      exudates_rep = as.factor(as.character(exudates_rep))) %>% 
        ggplot(aes(x = time, y = Index)) +
        # ggtitle(title) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.2), 
            alpha = 0.4, 
            size = 1, 
            aes(colour = genotype, group = genotype, shape = exudates_rep)) +
        geom_boxplot(alpha = 0.8, outlier.colour = NA, size = 1, fill = NA, aes(colour = genotype)) +
        scale_y_continuous(labels = scales::scientific) +
        scale_color_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% as.character(temp$genotype))]) + 
        scale_shape_manual(values = c(17,18,19)) + 
        facet_grid(dose ~ genotype, scales = "free", space = "free", switch = "both") +
        theme_RTN +
        theme(
            axis.line.x = element_line(size = 2, colour = c_black),
            axis.line.y = element_line(size = 2, colour = c_black),
            axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5)) +
        labs(x = "", y = paste("", colnames(df)[i]), fill = "", colour = "", shape = "")) +
    ggsave(filename = paste0(figs, "observed_", colnames(df)[i],"_alpha_diversity.png"), 
      dpi = 600, 
      device = "png", 
      units = img,
      bg="transparent", 
      width = 4, height = 5)

    return(plot_obj)

}, mc.cores = 4)

names(plot_list) <- c("shannon", "simpson", 
  # "fishers", 
  "observed_strains")

# Plot composite
p1 <- cowplot::plot_grid(
  plot_list$shannon + theme(legend.position = "none"),
  plot_list$simpson + theme(legend.position = "none"),
  plot_list$fishers + theme(legend.position = "none"),
  plot_list$observed_strains + theme(legend.position = "none"),
  ncol=1,
  nrow = 4, 
  align='hv', 
  axis = 'tblr',
  labels= NULL, 
  byrow = TRUE,
  rel_widths = c(1, 1, 1, 1), 
  rel_heights = c(1, 1, 1, 1)
)

cowplot::plot_grid(
  cowplot::get_legend(plot_list$shannon),
  p1,
  ncol=1,
  nrow = 2, 
  align='hv', 
  axis = 'tblr',
  labels= NULL,
  byrow = TRUE,
  rel_widths = c(.3, 1), 
  rel_heights = c(.3, 1)
) + ggsave(filename = paste0(figs, "summary_alpha_diversity.png"), 
      dpi = 600, 
      device = "png", 
      units = img,
      bg="transparent", 
      width = 7, height = 16)

# END OF SCRIPT
sessionInfo()
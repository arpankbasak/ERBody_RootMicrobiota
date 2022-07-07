# Summary of abosolute abundance
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

# Syncom non-spiked data
aa_mat <- sync.spiked.dat$cData[which(row.names(sync.spiked.dat$cData) %in% row.names(lfc_p)), which(colnames(sync.spiked.dat$cData) %in% sync.spiked.dat$metadata$SampleID)]

aa_df <- as.data.frame(aa_mat) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(
  fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  tech_rep = sync.spiked.dat$metadata$tech_rep[match(.$SampleID, sync.spiked.dat$metadata$SampleID)]) %>%
# group_by(strain, fixed) %>%
# summarise(AA = sum(AA)) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
mutate(genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time)) %>%
data.frame(., stringsAsFactors = FALSE)

# Inoc
inoc_df <- aa_df %>% filter(genotype == "inoc") %>%
mutate(id = paste(strain, dose, random, tech_rep, sep = "_")) %>%
select(id, AA)

aa_df$genus <- sync.spiked.dat$taxa$genus[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$family <- sync.spiked.dat$taxa$family[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$class <- sync.spiked.dat$taxa$class[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]
aa_df$phylum <- sync.spiked.dat$taxa$phylum[match(aa_df$strain, row.names(sync.spiked.dat$taxa))]

# Plot Canvas
# HEatmap sample-wise RA
aa_df_all <- as.data.frame(aa_mat) %>%
add_column(strain = row.names(.), .before = 1) %>%
gather(
  key = "SampleID", 
  value = "AA", 
  convert = FALSE, -strain) %>%
mutate(fixed = sync.spiked.dat$metadata$fixed[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  random = sync.spiked.dat$metadata$random[match(.$SampleID, sync.spiked.dat$metadata$SampleID)],
  tech_rep = sync.spiked.dat$metadata$tech_rep[match(.$SampleID, sync.spiked.dat$metadata$SampleID)]) %>%
separate(fixed, into = c("genotype", "dose", "time"), convert = FALSE, sep = "_") %>%
separate(random, into = c("seq_batch", "experiment", "tech_replicate"), convert = FALSE, sep = "_") %>%
mutate(id = paste(strain, dose, random, tech_replicate, sep = "_"),
  time_s = inoc_df$AA[match(.$id, inoc_df$id)]) %>%
mutate(
  genotype = factor(genotype, levels = genotype.syncom$short), 
  dose = as.factor(dose),
  time = as.factor(time),
  seq_batch = as.factor(time),
  experiment = as.factor(experiment),
  tech_replicate = as.factor(paste(tech_replicate, tech_rep, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

aa_df_all$family <- sync.spiked.dat$taxa$family[match(aa_df_all$strain, row.names(sync.spiked.dat$taxa))]

sig_strains <- as.character(unique(lfc_df$strain[lfc_df$significance == TRUE]))

# for each family -- plot the bacterial abundance trend over time
aa_df_all %>%
mutate(sig = as.factor(strain %in% sig_strains)) %>%
ggplot(aes(x = as.numeric(as.character(time)), 
  y = log10(AA+1), 
  group = tech_replicate)) +
# geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .2), 
#   shape = 15, aes(colour = genotype, alpha = sig, size = sig)) +
# geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = NA, width = 0.6) +
# geom_hline(data = , aes()) +
# geom_line(data = aa_df, 
#   inherit.aes = FALSE, 
#   aes(x = as.numeric(as.character(time)), y = log10(AA), colour = genotype)) +
# geom_smooth(aes(colour = genotype), se = FALSE, method = "lm", lwd = 0.3) +
geom_line(aes(colour = genotype)) +
facet_wrap(family ~ dose, scale = "free", as.table = TRUE) +
scale_size_manual(values = c(`FALSE` = 0.2, `TRUE` = 1)) +
scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 0.8)) +
scale_colour_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% aa_df_all$genotype)]) +
# scale_shape_manual(c(21:24)) +
scale_x_continuous() +
theme_RTN +
labs(x = "", y = "log10(Reads/spike_reads)", colour = "", fill = "", shape = "") +
ggsave(filename = paste0(figs, "/individual_analysis/spiked/individual_growth_summary.png"),
           bg="transparent", units = "in", width=10, height=12, limitsize=F, device = "png", dpi = 600)

aa_df_all %>%
mutate(sig = as.factor(strain %in% sig_strains)) %>%
ggplot(aes(x = (as.factor(time)), y = log10(AA+1))) +
geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .2), 
  shape = 1, aes(colour = genotype, alpha = sig, size = sig)) +
geom_boxplot(outlier.alpha = 0, aes(colour = genotype), alpha = 0.6, fill = "white", width = 0.6) +
# geom_hline(data = , aes()) +
# geom_line(data = aa_df, 
#   inherit.aes = FALSE, 
#   aes(x = as.numeric(as.character(time)), y = log10(AA), colour = genotype)) +
# geom_smooth(aes(colour = genotype), se = FALSE, method = "lm", lwd = 0.3) +
# geom_line(aes(colour = genotype)) +
facet_wrap(family ~ dose, scale = "free", as.table = TRUE) +
scale_size_manual(values = c(`FALSE` = 0.2, `TRUE` = 1)) +
scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 0.8)) +
scale_colour_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% aa_df_all$genotype)]) +
# scale_shape_manual(c(21:24)) +
# scale_x_continuous() +
theme_RTN +
labs(x = "", y = "log10(Reads/spike_reads)", colour = "", fill = "", shape = "") +
ggsave(filename = paste0(figs, "/individual_analysis/spiked/individual_growth_summary_bxp.png"),
           bg="transparent", width=10, height=7, limitsize=F, device = "png", dpi = 600)

# Trend of competency of total microbial load
aa_df_all %>%
group_by(SampleID,
genotype,
dose,
time,
seq_batch,
experiment,
tech_replicate) %>%
summarise(load = sum(AA)) %>%
ggplot(aes(x = as.numeric(as.character(time)), y = log10(load), group = tech_replicate)) +
# geom_smooth(aes(colour = genotype), se = FALSE, method = "lm", lwd = 0.3) +
# geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .2), 
#   shape = 15, aes(colour = genotype)) +
# geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = NA, width = 0.6) +
# geom_hline(data = , aes()) +
# geom_line(data = aa_df, 
#   inherit.aes = FALSE, 
#   aes(x = as.numeric(as.character(time)), y = log10(AA), colour = genotype)) +
geom_line(aes(colour = genotype)) +
facet_wrap(dose + experiment ~., scale = "free", as.table = TRUE) +
scale_colour_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% aa_df_all$genotype)]) +
# scale_shape_manual(c(21:24)) +
scale_x_continuous() +
theme_RTN +
labs(x = "", y = "log10(Total Abundance)", colour = "", fill = "", shape = "") +
ggsave(filename = paste0(figs, "/individual_analysis/spiked/microbial_load_summary.png"),
           bg="transparent", width=10, height=7, limitsize=F, device = "png", dpi = 600)


# A ternary plot showing the strains that are commonly depleted or enriched


# END OF SCRIPT
sessionInfo()
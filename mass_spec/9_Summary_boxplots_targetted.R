#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for data pre processing
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(mass_spec)

# Loading required packages
pkgs <- c("tidyverse", "limma", "parallel")

lapply(pkgs, require, character.only = T)

# Loading feature data
# load(paste0("./data/lfc_DAMs_cologne.Rdata"))
load(paste0("./data/msms_data_msdial.Rdata"))

meta <- ms_data$meta
mat <- ms_data$mat[, meta$ID]
geno <- genotype[which(genotype$short %in% meta$genotype),]

target_bar <- read.delim2("./data/targetted_metabolome.txt", sep = "\t", header = TRUE)
# target_msd <- read.delim2("./data/msdial_targets.txt", sep = "\t", header = TRUE, row.names = 1)
# target_msd_mar <- read.delim2("./data/msdial_targets_subset_intensity.txt", sep = "\t", header = TRUE, row.names = 1)

# The previous set
target <- target_bar[, c("mzl", "rt", "metabolite_name", meta$ID)] %>%
gather(key = "sample", value = "intensity", convert = FALSE, -mzl,-rt, -metabolite_name) %>%
mutate(intensity = 2^as.numeric(intensity)) %>%
cbind.data.frame(., meta[match(.$sample, meta$ID), c("batch", "genotype")])

target_raw <- target_bar[, c("mzl", "rt", "metabolite_name", meta$ID)] %>%
gather(key = "sample", value = "intensity", convert = FALSE, -mzl,-rt, -metabolite_name) %>%
mutate(intensity = as.numeric(intensity)) %>%
cbind.data.frame(., meta[match(.$sample, meta$ID), c("batch", "genotype")])

target_col <- target %>% group_by(metabolite_name, batch, genotype) %>%
summarise(x = mean(intensity)) %>%
spread(key = genotype, value = x, fill = 0) %>%
mutate(id = paste(metabolite_name, batch, sep = "_")) %>%
data.frame(., stringsAsFactors = FALSE)

target_norm <- target %>%
mutate(id = paste(metabolite_name, batch, sep = "_")) %>%
mutate(
  intensity_col = intensity/target_col$Col[match(.$id, target_col$id)],
  Col = target_col$Col[match(.$id, target_col$id)])

df_target_col <- target_col %>% select(-pyk10,-cyp,-id) %>%
spread(key = batch, value = Col) %>%
mutate(batch_cor = Cologne_A/Krakow) %>%
data.frame(.)

target_norm$batch_cor <- df_target_col$batch_cor[match(target_norm$metabolite_name, df_target_col$metabolite_name)]

# Compute stats for the targetted metabolites
stat_list <- mclapply(unique(target_norm$metabolite_name), function(x){

  fit <- aov(lm(intensity_col ~ genotype + batch_cor, data = target_norm %>% filter(metabolite_name == x)))
  temp <- unlist(summary(fit))
  stat <- broom::tidy(TukeyHSD(fit, method = "BH")) %>%
  filter(term != "batch") %>%
  separate(comparison, into = c("A", "B"), sep = "-") %>%
  # separate(A, into = c("Genotype_A", "batch_A"), sep = ":") %>%
  # separate(B, into = c("Genotype_B", "batch_B"), sep = ":") %>%
  # filter(B == "Col") %>%
  mutate(anova.F = temp["F value1"], anova.p.genotype = temp["Pr(>F)1"], metabolite_name = x) %>%
  data.frame(., stringsAsFactors = FALSE)

  return(stat)
}, mc.cores = 2) %>% do.call(rbind.data.frame, .) %>% 
# mutate(adj.p.value = p.ajust(adj.p.value, method = "fdr")) %>%
arrange(desc(adj.p.value)) %>%
mutate(significant = ifelse(adj.p.value <= 0.05, "*", ""))

write.table(stat_list, "./statistics/targetted_metabolome_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Boxplots for the targetted metabolites
p <- target_norm %>%
# mutate(class = factor(target_bar$class[match(.$metabolite_name, target_bar$metabolite_name)], 
#   levels = c("aliphatic", "benzyl", "indole", "indolic compounds", "coumarins"))) %>%
ggplot(aes(x = genotype, y = intensity_col)) +
geom_point(
  size = 0.5, 
  # shape = 21, 
  position = position_jitterdodge(0.2), 
  alpha = 0.8, 
  aes(fill = genotype, shape = batch), 
  colour = c_black) +
geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.6) +
# geom_point(data = target_norm %>% filter(genotype == "Col"), aes(x = genotype, y = batch_cor), 
#   fill = c_white, 
#   shape = 22,
#   colour = c_black, alpha = 0.6,
#   size = 1, 
#   override.aes = TRUE) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    facet_wrap(.~ metabolite_name, scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    scale_shape_manual(values = c(21, 23)) +
    ylim(c(0,NA)) +
    theme_RTN +
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 2)) +
    labs(x = "", y = "Relative abundance normalised to Col-0", fill = "", colour = "")

ggsave(p, file = paste0(figs, "boxplot_targets_all_batch_norm.pdf"),
    units = "in",
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 600
  )

p <- target_norm %>%
# mutate(class = factor(target_bar$class[match(.$metabolite_name, target_bar$metabolite_name)], 
  # levels = c("aliphatic", "benzyl", "indole", "indolic compounds", "coumarins"))) %>%
ggplot(aes(x = genotype, y = intensity_col)) +
geom_point(
  size = 0.5, 
  # shape = 21, 
  position = position_jitterdodge(0.2), 
  alpha = 0.8, 
  aes(fill = genotype, shape = batch), 
  colour = c_black) +
geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.6) +
geom_point(data = target_norm %>% filter(genotype == "Col"), aes(x = genotype, y = (batch_cor), shape = batch), 
  fill = c_white, 
  colour = c_black, alpha = 0.6,
  size = 1, 
  override.aes = TRUE) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    facet_wrap(.~ metabolite_name, scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    scale_shape_manual(values = c(21, 23)) +
    ylim(c(0,NA)) +
    theme_RTN +
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 2)) +
    labs(x = "", y = "Relative abundance normalised to Col-0", fill = "", colour = "")

ggsave(p, file = paste0(figs, "boxplot_targets_all_col_marked.pdf"),
    units = "in",
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 600
  )


p <- target_raw %>%
# mutate(class = factor(target_bar$class[match(.$metabolite_name, target_bar$metabolite_name)], 
  # levels = c("aliphatic", "benzyl", "indole", "indolic compounds", "coumarins"))) %>%
ggplot(aes(x = genotype, y = 2^(intensity))) +
geom_point(
  size = 0.5, 
  # shape = 21, 
  position = position_jitterdodge(0.2), 
  alpha = 0.8, 
  aes(fill = genotype, shape = batch), 
  colour = c_black) +
    geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.6) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    facet_wrap(.~ metabolite_name, scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% target$genotype)]) +
    scale_shape_manual(values = c(21, 23)) +
    ylim(c(0,NA)) +
    theme_RTN +
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 2)) +
    labs(x = "", y = "Relative abundance", fill = "", colour = "")

ggsave(p, file = paste0(figs, "boxplot_targets_all_raw.pdf"),
    units = "in",
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 600
  )

# END OF SCRIPT
sessionInfo()
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

target_bar <- read.delim2("./data/targetted_metabolome_root_extract.txt", sep = "\t", header = TRUE)
# target_msd <- read.delim2("./data/msdial_targets.txt", sep = "\t", header = TRUE, row.names = 1)
# target_msd_mar <- read.delim2("./data/msdial_targets_subset_intensity.txt", sep = "\t", header = TRUE, row.names = 1)

# Create metadata

idx <- str_detect(colnames(target_bar), "^col|^pyk|^cyp")
meta <- data.frame(SampleID = colnames(target_bar)[idx], row.names = colnames(target_bar)[idx]) %>%
separate(SampleID, remove = FALSE, into = c("Genotype", "ID1", "ID2", "ID3"), sep = "_")

meta$Genotype[str_detect(meta$Genotype, "^col")] <- "Col"
meta$Genotype[str_detect(meta$Genotype, "^pyk")] <- "pyk10"
meta$Genotype[str_detect(meta$Genotype, "^cyp")] <- "cyp"

# The previous set
df <- cbind.data.frame(metabolite_name = target_bar[, "Metabolite.name"], target_bar[, idx]) %>%
gather(key = "SampleID", value = "RA", -metabolite_name, convert = FALSE) %>%
mutate(
  RA = as.numeric(RA), 
  RA_1k = RA*100, 
  log_RA = log2(RA+.001),
  Genotype = factor(meta$Genotype[match(as.character(.$SampleID), as.character(meta$SampleID))], levels = geno$short),
  group = paste(metabolite_name, Genotype, sep = "_")
  )

df_mean <- df %>%
filter(!is.na(RA)) %>%
group_by(group) %>%
summarise(RA = mean(RA)) %>%
mutate(
  RA = as.numeric(RA), 
  RA_1k = RA*100, 
  log_RA = log2(RA+.001)
) %>%
data.frame(., row.names = .$group)

df_mean_col <- df %>%
filter(!is.na(RA) & Genotype == "Col") %>%
group_by(metabolite_name) %>%
summarise(RA = mean(RA)) %>%
mutate(
  RA = as.numeric(RA), 
  RA_1k = RA*100, 
  log_RA = log2(RA+.001)
) %>%
data.frame(., row.names = .$metabolite_name)

df$mean_center_genotype <- df_mean$RA[match(df$group, row.names(df_mean))]
df$mean_center_col <- df_mean_col$RA[match(df$metabolite_name, row.names(df_mean_col))]
df$norm_col <- df$RA/df$mean_center_col



# Compute stats for the targetted metabolites
stat_list <- mclapply(unique(df$metabolite_name), function(x){

  # x = "IAA"
  fit <- aov(lm(RA ~ Genotype, data = df %>% filter(metabolite_name == x)))
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

write.table(stat_list, "./statistics/targetted_metabolome_all_root_extract_.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Boxplots for the targetted metabolites
p <- df %>%
mutate(mean_center_genotype = mean_center_genotype/mean_center_col) %>%
# mutate(class = factor(target_bar$class[match(.$metabolite_name, target_bar$metabolite_name)], 
#   levels = c("aliphatic", "benzyl", "indole", "indolic compounds", "coumarins"))) %>%
ggplot(aes(x = Genotype, y = norm_col)) +
geom_point(
  size = 0.5, 
  # shape = 21, 
  position = position_jitterdodge(0.2), 
  alpha = 0.8, 
  aes(fill = Genotype),
  shape = 21, 
  colour = c_black) +
geom_boxplot(outlier.alpha = 0, aes(colour = Genotype), fill = c_white, alpha = 0.6) +
# geom_point(
#   size = 1.2, 
#   # shape = 21, 
#   position = position_jitterdodge(0.2), 
#   alpha = 0.8, 
#   aes(
#     y = mean_center_genotype, 
#     fill = Genotype),
#   shape = 23, 
#   colour = c_black) +
# geom_point(data = target_norm %>% filter(genotype == "Col"), aes(x = genotype, y = batch_cor), 
#   fill = c_white, 
#   shape = 22,
#   colour = c_black, alpha = 0.6,
#   size = 1, 
#   override.aes = TRUE) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% df$Genotype)]) +
    facet_wrap(.~ metabolite_name, scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% df$Genotype)]) +
    # scale_shape_manual(values = c(21, 23)) +
    ylim(c(0,NA)) +
    theme_RTN +
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 2)) +
    labs(x = "", y = "Relative abundance normalised to Col-0", fill = "", colour = "")

ggsave(p, file = paste0(figs, "boxplot_targets_all_batch_norm_root_extract.pdf"),
    units = "in",
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 600
  )

# Boxplots for the targetted metabolites
p <- df %>%
mutate(mean_center_genotype = mean_center_genotype/mean_center_col) %>%
# mutate(class = factor(target_bar$class[match(.$metabolite_name, target_bar$metabolite_name)], 
#   levels = c("aliphatic", "benzyl", "indole", "indolic compounds", "coumarins"))) %>%
ggplot(aes(x = Genotype, y = RA)) +
geom_point(
  size = 0.5, 
  # shape = 21, 
  position = position_jitterdodge(0.2), 
  alpha = 0.8, 
  aes(fill = Genotype),
  shape = 21, 
  colour = c_black) +
geom_boxplot(outlier.alpha = 0, aes(colour = Genotype), fill = c_white, alpha = 0.6) +
# geom_point(
#   size = 1.2, 
#   # shape = 21, 
#   position = position_jitterdodge(0.2), 
#   alpha = 0.8, 
#   aes(
#     y = mean_center_genotype, 
#     fill = Genotype),
#   shape = 23, 
#   colour = c_black) +
# geom_point(data = target_norm %>% filter(genotype == "Col"), aes(x = genotype, y = batch_cor), 
#   fill = c_white, 
#   shape = 22,
#   colour = c_black, alpha = 0.6,
#   size = 1, 
#   override.aes = TRUE) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% df$Genotype)]) +
    facet_wrap(.~ metabolite_name, scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% df$Genotype)]) +
    # scale_shape_manual(values = c(21, 23)) +
    ylim(c(0,NA)) +
    theme_RTN +
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(size = 2)) +
    labs(x = "", y = "Relative abundance", fill = "", colour = "")

ggsave(p, file = paste0(figs, "boxplot_targets_all_batch_RA_root_extract.pdf"),
    units = "in",
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 600
  )

# END OF SCRIPT
sessionInfo()
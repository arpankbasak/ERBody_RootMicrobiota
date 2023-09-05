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
load(paste0("./data/lfc_DAMs_annotated.Rdata"))
load(paste0("./data/msms_data.Rdata"))
meta <- ms_data$meta

# Separate negative and positive ionisation
geno <- genotype[which(genotype$short %in% ms_data$meta$genotype),]

FC_threshold <- 0.05
alpha <- 0.05

# Carpent data
df <- as.data.frame(lfc_obj$logFC_P.wide)
id <- str_detect(colnames(df), "_logFC$")
d <- dist((df[,id]))
hc <- hclust(as.dist(d), "ward.D2")
lfc_sorted <- row.names(df)[hc$order]

lfc_df <- as.data.frame(lfc_obj$logFC_P.long) %>%
  mutate(
    Genotype = factor(contrast, levels = geno$short[which(geno$short %in% .$contrast)]),
    sig = PValue <= alpha,
    peak_id = factor(peak_id, levels = lfc_sorted)
    )

# Plot heatmap
lfc_df %>%
ggplot(aes(x = Genotype, y=peak_id, fill = saturate(logFC))) +
geom_raster() +
# facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
scale_fill_gradient2(midpoint = 0, high = c_green, mid = c_white, low = c_cudo_magenta) +
theme_RTN +
theme(axis.text.y = element_blank()) +
labs(x = "", y = "", fill = "logFC") +
ggsave(file = paste0(figs, "heatmap_lfc_annotated.png"),
units = "in",
device = "png",
width = 4.5,
height = 7.5
)

# stat <- cor.test(df$cyp_logFC, df$pyk10_logFC)
# til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", round(stat$p.value, 4))

df <- as.data.frame(lfc_obj$annotated_mat) %>% 
gather(key = "SampleID", value = "RA", convert = FALSE, -Compound) %>%
cbind.data.frame(., meta[match(.$SampleID, meta$given_id),]) %>%
mutate(
  Compound = factor(Compound, levels = lfc_sorted),
  log_ra = ifelse(!is.finite(log10(RA)), NA, log10(RA))
  ) %>%
data.frame(.)

df %>%
  ggplot(aes(x = SampleID, y=Compound, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  facet_grid(.~ batch + genotype, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    ) +
  labs(x = "", y = "", fill = "log10RA") +
  ggsave(file = paste0(figs, "heatmap_ra_annotated.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 7.5
  )

df %>%
group_by(Compound, genotype, batch) %>%
summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(mean_ra)), NA, log10(mean_ra))),
  scale_ra = scale(mean_ra)) %>%
  ggplot(aes(x = genotype, y=Compound, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  # scale_fill_gradient2(midpoint = 0, high = c_red, mid = c_white, low = c_blue, na.value = c_black) +
  facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    ) +
  labs(x = "", y = "", fill = "log10RA") +
  ggsave(file = paste0(figs, "heatmap_meanra_compound.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 7.5
  )

df %>%
mutate(
  # log_ra = (ifelse(!is.finite(log10(RA)), NA, log10(RA))),
  scale_ra = scale(RA)
  ) %>%
  ggplot(aes(x = genotype, y = log_ra)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8), 
    alpha = 0.9, shape = 21, colour = c_black, aes(fill = genotype)
    ) +
  geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.5, size = 0.2) +
  facet_grid(
    Compound ~ batch,  
    scale = "free_x", 
    switch = "both", 
    # nrow = length(unique(as.character(df$Compound))), 
    # ncol = length(unique(as.character(df$batch)))
    ) +
  scale_fill_manual(values = c(geno$colours, "darkgrey")) +
  scale_colour_manual(values = c(geno$colours, "darkgrey")) +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 3), 
    axis.text.x = element_text(angle = 90, size = 3),
    strip.text.x = element_text(angle = 90, size = 3),
    strip.text.y = element_text(angle = 90, size = 3)
    ) +
  labs(x = "", y = "", fill = "", colour = "") +
  ggsave(file = paste0(figs, "boxplot_ra_annotated.png"),
  units = "in",
  device = "png",
  width = 3,
  height = 45
  )
      

df %>%
  ggplot(aes(x = genotype, y = log_ra)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.8), 
    alpha = 0.9, colour = c_black, 
    aes(fill = genotype, shape = batch)
    ) +
  geom_boxplot(outlier.alpha = 0, aes(colour = genotype), 
    fill = c_white, 
    alpha = 0.5, size = 0.2) +
  facet_wrap(
    . ~ Compound,  
    scale = "free", 
    switch = "both") +
  scale_fill_manual(values = c(geno$colours, "darkgrey")) +
  scale_colour_manual(values = c(geno$colours, "darkgrey")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 3), 
    axis.text.x = element_text(angle = 90, size = 3),
    strip.text.x = element_text(angle = 0, size = 3),
    strip.text.y = element_text(angle = 0, size = 3)
    ) +
  labs(x = "", y = "", fill = "", colour = "") +
  ggsave(file = paste0(figs, "boxplot_ra_annotated_batch.png"),
  units = "in",
  device = "png",
  width = 16,
  height = 7
  )            

# END OF SCRIPT
sessionInfo()
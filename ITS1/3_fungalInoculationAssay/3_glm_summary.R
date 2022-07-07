# Set R path

# Script by
# @Arpan Kumar Basak
# Quantification of fungal inoculation on plants that are:
# - Treated with flg22 to trigger immune response 10nM
# - Mutants that lack in PYK10 myrosinases and Trp-pathway

# Script for summarising statistical analysis

rm(list = ls())

out <- "./output/"
stats <- "./statistics/"
figs <- "./figures/"

# Read data
source("./scripts/parameters.R")
setwd(analysis.fungali)
require(tidyverse)

# source("./scripts/analysis_statistics.R")
load("./data/data.RData")
load("./data/stats.RData")
load("./data/taxonomy.RData")
load("./data/stats_mock.RData")


# Load data
fit_df <- stats$stat_fit
fit_temp_strain <- stats$stat_strain
fit_temp <- stats_mock$stat_strain
strain_sorted_phy <- sort_list$strain_sorted_phy
strain_sorted <- sort_list$strain_sorted
df_prl <- dat$prl %>% filter(strain != "F131")
df <- dat$sfw %>% filter(strain != "F131")

out <- setdiff(unique(fit_df$strain), taxa_obj$strain_phylogeny)
labs <- c(as.character(taxa_obj$taxa$ID), out)
names(labs) <- c(taxa_obj$taxa$strain_ids, out)

labs <- labs[strain_sorted]

strain_sorted_phy <- c(as.character(taxa_obj$strain_phylogeny), out)
strain_labs_phy <- c(as.character(taxa_obj$taxa$ID), out)
strain_class <- c(as.character(taxa_obj$taxa$class), out)
names(strain_class) <- c(taxa_obj$taxa$strain_ids, out)
levs <- c(
	"Sordariomycetes",
	"Leotiomycetes",
	"Dothideomycetes",
	"Agaricomycetes",
	"Mortierellomycetes", 
	"mock")


# One heatmap to summarise everything
(hmap_phy_2 <- fit_temp_strain %>%
	filter(treatment != "flg22") %>%
	 ggplot(aes(
	 	x = genotype, 
	 	y = strain)) +
	 geom_raster(aes(fill = diff_sat)) +
	 geom_tile(aes(colour = FDR), 
	 	fill = NA, 
	 	width = 0.9, 
	 	height = 0.9, 
	 	size = 0.9) +
	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
	 scale_fill_gradient2(
	 	low = "darkmagenta", 
	 	high = "darkgreen", 
	 	mid = "white", 
	 	na.value = "lightgrey") +
	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
	 	guide = FALSE) +
	 theme_classic() +
	 theme(
	 	panel.spacing = unit(0.01, "lines"), 
	 	panel.border = element_rect(fill="transparent", colour=NA),
	 	legend.position = "top", 
	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
	ggsave("./figures/sfw_selection_phylogeny_vs_mockWT_2.png", 
		units = "in", 
		width = 3, 
		height = 7, 
		bg = "white")

(hmap_phy_mock_2 <- fit_temp %>%
	filter(treatment != "flg22") %>%
	 ggplot(aes(
	 	x = genotype, 
	 	y = strain)) +
	 geom_raster(aes(fill = diff_sat)) +
	 geom_tile(aes(colour = FDR), 
	 	fill = NA, 
	 	width = 0.9, 
	 	height = 0.9, 
	 	size = 0.9) +
	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
	 scale_fill_gradient2(
	 	low = "darkmagenta", 
	 	high = "darkgreen", 
	 	mid = "white", 
	 	na.value = "lightgrey") +
	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
	 	guide = FALSE) +
	 theme_classic() +
	 theme(
	 	panel.spacing = unit(0.01, "lines"), 
	 	panel.border = element_rect(fill="transparent", colour=NA),
	 	legend.position = "top", 
	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
	ggsave("./figures/hmap_phylogeny_vs_mockWT_2.png", 
		units = "in", 
		width = 2, 
		height = 7, 
		bg = "white")

# With the Mock control
# (hmap_sfw_phy_mock <- fit_group %>%
# 	 mutate(strain = factor(A_strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A_genotype),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$A_strain, df$strain)],
# 	 		tag = factor(replace_na(strain_class[match(.$A_strain, names(strain_class))], "Mock"), levels = levs),
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate) # saturate the negatives
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 facet_grid(tag~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
# 	ggsave("./figures/sfw_selection_phylogeny_vs_mockWT.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# # Selection basis heatmap -- SFW clustering
# (hmap_sfw_mock <- fit_group %>%
# 	 mutate(strain = factor(A_strain, levels = (strain_sorted_group), labels = (labs)),
# 	 		genotype = as.factor(A_genotype),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$A_strain, df$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
# 	ggsave("./figures/sfw_selection_vs_mockWT.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# # With phylogey
# (hmap_prl_phy_mock <- fit_group_prl %>%
# 	 mutate(strain = factor(A_strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A_genotype),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(adj.p.value < 0.05),
# 	 		cc = df_prl$cc[match(.$A_strain, df_prl$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(replace_na(estimate, -6)))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white",
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
# 	ggsave("./figures/prl_selection_phylogeny_vs_mockWT.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# # With clustering
# (hmap_prl_mock <- fit_group_prl %>%
# 	 mutate(strain = factor(A_strain, levels = (strain_sorted_group), labels = (labs)),
# 	 		genotype = as.factor(A_genotype),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(adj.p.value < 0.05),
# 	 		cc = df_prl$cc[match(.$A_strain, df_prl$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(replace_na(estimate, -6)))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white",
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
# 	ggsave("./figures/prl_selection_vs_mockWT.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")


# For the line marked for 845_mock
mock_levs <- 1000*c(mean(df$shoot_fresh_weight[df$strain == "mock" & df$genotype == "845"]))

# # Plot data
# ## Clustering
# (bxp_sfw <- df %>%
# 	 mutate(strain = factor(strain, levels = rev(strain_sorted_group), labels = rev(labs))) %>%
# 	 ggplot(aes(y = 1000*shoot_fresh_weight, x = group)) +
# 	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
# 	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
# 	 	aes(shape = random, 
# 	 		colour = group), 
# 	 	size = 0.6, alpha = 0.6) +
# 	 geom_boxplot(fill = NA, 
# 	 	aes(colour = group),
# 	 	outlier.alpha = 0, 
# 	 	alpha = 0.8) +
# 	 facet_wrap(strain ~ ., scale = "free_y") +
# 	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
# 	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "")) +
# 	ggsave("./figures/sfw_bxp.png", 
# 		units = "in", 
# 		width = 12, 
# 		height = 9, 
# 		bg = "white")

(bxp_sfw <- df %>%
	filter(group != "845_flg22") %>%		
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = 1000*shoot_fresh_weight, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "")) +
	ggsave("./figures/sfw_bxp_2.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")


# Phylogeny
(bxp_sfw_phy <- df %>%
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = 1000*shoot_fresh_weight, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "")) +
	ggsave("./figures/sfw_bxp_phylogeny.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

mock_levs <- c(mean(df_prl$prl[df_prl$strain == "mock" & df_prl$genotype == "845"]))

# (bxp_prl <- df_prl %>%
# 	 mutate(strain = factor(strain, levels = rev(strain_sorted_group), labels = rev(labs))) %>%
# 	 ggplot(aes(y = prl, x = group)) +
# 	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
# 	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
# 	 	aes(shape = random, 
# 	 		colour = group), 
# 	 	size = 0.6, alpha = 0.6) +
# 	 geom_boxplot(fill = NA, 
# 	 	aes(colour = group),
# 	 	outlier.alpha = 0, 
# 	 	alpha = 0.8) +
# 	 facet_wrap(strain ~ ., scale = "free_y") +
# 	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
# 	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
# 	ggsave("./figures/prl_bxp.png", 
# 		units = "in", 
# 		width = 12, 
# 		height = 9, 
# 		bg = "white")


# (bxp_prl <- df_prl %>%
# 	filter(group != "845_flg22") %>%	
# 	 mutate(strain = factor(strain, levels = rev(strain_sorted_group), labels = rev(labs))) %>%
# 	 ggplot(aes(y = prl, x = group)) +
# 	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
# 	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
# 	 	aes(shape = random, 
# 	 		colour = group), 
# 	 	size = 0.6, alpha = 0.6) +
# 	 geom_boxplot(fill = NA, 
# 	 	aes(colour = group),
# 	 	outlier.alpha = 0, 
# 	 	alpha = 0.8) +
# 	 facet_wrap(strain ~ ., scale = "free_y") +
# 	 scale_colour_manual(values = c("black","darkcyan", "indianred")) +
# 	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
# 	ggsave("./figures/prl_bxp_2.png", 
# 		units = "in", 
# 		width = 12, 
# 		height = 9, 
# 		bg = "white")

# Phylogeny
(bxp_prl_phy <- df_prl %>%
	# filter(group != "845_flg22") %>%	
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = prl, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
	ggsave("./figures/prl_bxp_phylogeny.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

(bxp_prl_phy <- df_prl %>%
	filter(group != "845_flg22") %>%	
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = prl, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
	ggsave("./figures/prl_bxp_phylogeny_2.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")


# A complete heatmap
# (hmap_phy_2 <- fit_temp_strain %>%
# 	filter(treatment != "flg22") %>%
# 	 ggplot(aes(
# 	 	x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = diff_sat)) +
# 	 geom_tile(aes(colour = FDR), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
# 	 scale_fill_gradient2(
# 	 	low = c_cudo_magenta, 
# 	 	high = c_dark_green, 
# 	 	mid = c_white, 
# 	 	na.value = c_black) +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(
# 	 	panel.spacing = unit(0.01, "lines"), 
# 	 	panel.border = element_rect(fill="transparent", colour=NA),
# 	 	legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
# 	ggsave("./figures/sfw-prl_selection_phylogeny_mock_2.png", 
# 		units = "in", 
# 		width = 2, 
# 		height = 7, 
# 		bg = "white")

(hmap_phy <- fit_temp_strain %>%
	# filter(treatment != "flg22") %>%
	 ggplot(aes(
	 	x = genotype, 
	 	y = strain)) +
	 geom_raster(aes(fill = diff_sat)) +
	 geom_tile(aes(colour = FDR), 
	 	fill = NA, 
	 	width = 0.9, 
	 	height = 0.9, 
	 	size = 0.9) +
	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
	 scale_fill_gradient2(
	 	low = c_cudo_magenta, 
	 	high = c_dark_green, 
	 	mid = c_white, 
	 	na.value = c_black) +
	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
	 	guide = FALSE) +
	 theme_classic() +
	 theme(
	 	panel.spacing = unit(0.01, "lines"), 
	 	panel.border = element_rect(fill="transparent", colour=NA),
	 	legend.position = "top", 
	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
	ggsave("./figures/sfw-prl_selection_phylogeny.png", 
		units = "in", 
		width = 2, 
		height = 7, 
		bg = "white")

(hmap_phy_2 <- fit_temp_strain %>%
	filter(treatment != "flg22") %>%
	 ggplot(aes(
	 	x = genotype, 
	 	y = strain)) +
	 geom_raster(aes(fill = diff_sat)) +
	 geom_tile(aes(colour = FDR), 
	 	fill = NA, 
	 	width = 0.9, 
	 	height = 0.9, 
	 	size = 0.9) +
	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
	 scale_fill_gradient2(
	 	low = c_cudo_magenta, 
	 	high = c_dark_green, 
	 	mid = c_white, 
	 	na.value = c_black) +
	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
	 	guide = FALSE) +
	 theme_classic() +
	 theme(
	 	panel.spacing = unit(0.01, "lines"), 
	 	panel.border = element_rect(fill="transparent", colour=NA),
	 	legend.position = "top", 
	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
	ggsave("./figures/sfw-prl_selection_phylogeny_2.png", 
		units = "in", 
		width = 2, 
		height = 7, 
		bg = "white")

(hmap_phy_mock_2 <- fit_temp %>%
	filter(treatment != "flg22") %>%
	 ggplot(aes(
	 	x = genotype, 
	 	y = strain)) +
	 geom_raster(aes(fill = diff_sat)) +
	 geom_tile(aes(colour = FDR), 
	 	fill = NA, 
	 	width = 0.9, 
	 	height = 0.9, 
	 	size = 0.9) +
	 facet_grid(class ~ val, space = "free", scale = "free", switch = "both") +
	 scale_fill_gradient2(
	 	low = c_cudo_magenta, 
	 	high = c_dark_green, 
	 	mid = c_white, 
	 	na.value = c_black) +
	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
	 	guide = FALSE) +
	 theme_classic() +
	 theme(
	 	panel.spacing = unit(0.01, "lines"), 
	 	panel.border = element_rect(fill="transparent", colour=NA),
	 	legend.position = "top", 
	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "", colour = "", fill = "log2(vs 845:mock)")) +
	ggsave("./figures/sfw-prl_selection_phylogeny_mock_2.png", 
		units = "in", 
		width = 2, 
		height = 7, 
		bg = "white")

# Within each strain
# Selection basis heatmap -- SFW Phylogeny
# (hmap_sfw_phy <- fit_temp_strain %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$strain, df$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate) # saturate the negatives
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/sfw_selection_phylogeny.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# (hmap_sfw_phy_2 <- fit %>%
# 	filter(A != "845_flg22") %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$strain, df$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate) # saturate the negatives
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/sfw_selection_phylogeny_2.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# Selection basis heatmap -- SFW clustering
# (hmap_sfw <- fit %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted), labels = (labs)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$strain, df$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/sfw_selection.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# (hmap_sfw_2 <- fit %>%
# 	filter(A != "845_flg22") %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted), labels = (labs)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(significance),
# 	 		cc = df$cc[match(.$strain, df$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white", 
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/sfw_selection_2.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# With phylogey
# (hmap_prl_phy <- fit_prl %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(adj.p.value < 0.05),
# 	 		cc = df_prl$cc[match(.$strain, df_prl$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white",
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/prl_selection_phylogeny.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# (hmap_prl_phy_2 <- fit_prl %>%
# 	filter(A != "845_flg22") %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted_phy), labels = (strain_labs_phy)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(adj.p.value < 0.05),
# 	 		cc = df_prl$cc[match(.$strain, df_prl$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white",
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/prl_selection_phylogeny_2.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")

# # With clustering
# (hmap_prl <- fit_prl %>%
# 	 mutate(strain = factor(strain, levels = (strain_sorted), labels = (labs)),
# 	 		genotype = as.factor(A),
# 	 		# A_treatment = as.factor(A_treatment),
# 	 		significance = as.factor(adj.p.value < 0.05),
# 	 		cc = df_prl$cc[match(.$strain, df_prl$strain)],
# 	 		estimate = ifelse(abs(estimate) > 2, sign(estimate) * 2, estimate)
# 	 	) %>%
# 	 ggplot(aes(x = genotype, 
# 	 	y = strain)) +
# 	 geom_raster(aes(fill = saturate(estimate))) +
# 	 geom_tile(aes(colour = significance), 
# 	 	fill = NA, 
# 	 	width = 0.9, 
# 	 	height = 0.9, 
# 	 	size = 0.9) +
# 	 # facet_grid(cc~., space = "free", scale = "free_y", switch = "y") +
# 	 scale_fill_gradient2(low = "darkmagenta", 
# 	 	high = "darkgreen", mid = "white",
# 	 	na.value = "darkmagenta") +
# 	 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA),
# 	 	guide = FALSE) +
# 	 theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
# 	 	axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "", y = "", colour = "", fill = "log2(vs 845_mock)")) +
# 	ggsave("./figures/prl_selection.png", 
# 		units = "in", 
# 		width = 4, 
# 		height = 9, 
# 		bg = "white")


# For the line marked for 845_mock
mock_levs <- 1000*c(mean(df$shoot_fresh_weight[df$strain == "mock" & df$genotype == "845"]))

# Plot data
## Clustering
(bxp_sfw <- df %>%
	 mutate(strain = factor(strain, levels = rev(strain_sorted), labels = rev(labs))) %>%
	 ggplot(aes(y = 1000*shoot_fresh_weight, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "")) +
	ggsave("./figures/sfw_bxp.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

# Phylogeny
(bxp_sfw_phy <- df %>%
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = 1000*shoot_fresh_weight, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "")) +
	ggsave("./figures/sfw_bxp_phylogeny.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

mock_levs <- c(mean(df_prl$prl[df_prl$strain == "mock" & df_prl$genotype == "845"]))

(bxp_prl <- df_prl %>%
	 mutate(strain = factor(strain, levels = rev(strain_sorted), labels = rev(labs))) %>%
	 ggplot(aes(y = prl, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
	ggsave("./figures/prl_bxp.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

# Phylogeny
(bxp_prl_phy <- df_prl %>%
	 mutate(strain = factor(strain, levels = rev(strain_sorted_phy), labels = rev(strain_labs_phy))) %>%
	 ggplot(aes(y = prl, x = group)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = random, 
	 		colour = group), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = group),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ ., scale = "free_y") +
	 scale_colour_manual(values = c(c_black, c_dark_grey, c_cudo_skyblue, c_red)) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14, 8, 13, 9, 12)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "")) +
	ggsave("./figures/prl_bxp_phylogeny.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")



# Compared to mock WT
# Regression analysis for mean estimates
# idx <- match(fit_mat_group$A_strain, fit_mat_group_prl$A_strain)
# fit_mat_group <- cbind.data.frame(fit_mat_group, fit_mat_group_prl[idx, -1])

# fit_df_group <- fit_mat_group %>%
# gather(key = "key", value = "vals", convert = FALSE, -A_strain) %>%
# separate(key, into = c("measure", "A"), sep = ":") %>%
# spread(key = measure, value = vals, convert = FALSE) %>%
# data.frame(., stringsAsFactors = FALSE)

# # 
# mock_lev <- fit_df_group %>% group_by(A_strain) %>%
# filter(A_strain == "mock") %>%
# ungroup() %>%
# select(-A_strain) %>%
# data.frame(., stringsAsFactors = FALSE)

# fit_df_group$behaviour_PRL <- ((fit_df_group$PRL - mock_lev$PRL[match(fit_df_group$A, mock_lev$A)]))
# # fit_df$behaviour_PRL <- ifelse(fit_df$PRL == 1, min(fit_df$behaviour_PRL), fit_df$behaviour_PRL)

# fit_df_group$behaviour_SFW <- ((fit_df_group$SFW - mock_lev$SFW[match(fit_df_group$A, mock_lev$A)]))
# fit_df$behaviour_SFW <- ifelse(fit_df$SFW == 1, min(fit_df$behaviour_SFW), fit_df$behaviour_SFW)

# Set X and Y range for mock
# marker_x <- mean(fit_df$PRL[fit_df$strain == "mock"])
# marker_y <- mean(fit_df$SFW[fit_df$strain == "mock"])
# marker_x <- 0
# marker_y <- 0
# cut_off <- 0.5

# # fit_df_group$lab_group <- paste(fit_df_group$strain, taxa_obj$taxa$ID[match(fit_df_group$A_strain, taxa_obj$taxa$strain_id)], sep = "_")

# # Qurdant analysis using MDA using estimates
# (scatter_plt_mock <- fit_df_group %>%
# 	separate(A, 
# 		remove = FALSE, into = c("genotype", "treatment"), 
# 		sep = "_", 
# 		convert = FALSE) %>%
# 	mutate(genotype = factor(as.character(genotype), levels = c("845", "992", "994")),
# 		treatment = factor(treatment, levels = c("mock", "flg22")),
# 		# significance = as.factor(p.value <= 0.05 | prl_pvalue <= 0.05)
# 		target = (abs(PRL) > cut_off & abs(SFW) > cut_off),
# 	) %>%
# 	mutate(text_target = ifelse(target == TRUE, as.character(A_strain), "")) %>%
# 	ggplot(aes(x = saturate(PRL), y = saturate(SFW))) +
# 	geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 	geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 	geom_point(aes(fill = A_strain, 
# 		# shape = genotype, 
# 		size = target, 
# 		colour = target, 
# 		alpha = target), shape = 22) +
# 	# ylim(c(-7, 7)) +
# 	# xlim(c(-4, 4)) +
# 	ggrepel::geom_text_repel(aes(label = text_target)) +
# 	facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 	# scale_fill_brewer(palette = "Spectral") +
# 	# scale_shape_manual(values = c(22,23,24)) +
# 	scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 	scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 	scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 	theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "log2(vs 845:mock) (Primary Root length)", y = "log2(vs 845:mock) (Shoot Fresh Weight)", 
# 	 	colour = "", shape = "", fill = "", size = "")
# 	) +
# ggsave("./figures/regression_mockWT.png", 
# 		units = "in", 
# 		width = 6, 
# 		height = 4, 
# 		bg = "white")

# (scatter_plt_behav_mock <- fit_df_group %>%
# 	separate(A, remove = FALSE, into = c("genotype", "treatment"), 
# 		sep = "_", 
# 		convert = FALSE) %>%
# 	mutate(genotype = factor(as.character(genotype), levels = c("845", "992", "994")),
# 		treatment = factor(treatment, levels = c("mock", "flg22")),
# 		# significance = as.factor(p.value <= 0.05 | prl_pvalue <= 0.05)
# 		target = (abs(PRL) > cut_off & abs(SFW) > cut_off),
# 	) %>%
# 	mutate(text_target = ifelse(target == TRUE, as.character(A_strain), "")) %>%
# 	ggplot(aes(x = saturate(behaviour_PRL), y = saturate(behaviour_PRL))) +
# 	geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 	geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 	geom_point(aes(fill = A_strain, 
# 		# shape = genotype, 
# 		size = target, 
# 		colour = target, 
# 		alpha = target), shape = 22) +
# 	# ylim(c(-7, 7)) +
# 	# xlim(c(-4, 4)) +
# 	ggrepel::geom_text_repel(aes(label = text_target)) +
# 	facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 	# scale_fill_brewer(palette = "Spectral") +
# 	# scale_shape_manual(values = c(22,23,24)) +
# 	scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 	scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 	scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 	theme_classic() +
# 	 theme(legend.position = "top", 
# 	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 	 labs(x = "log2(vs 845:mock) (Primary Root length)", y = "log2(vs 845:mock) (Shoot Fresh Weight)", 
# 	 	colour = "", shape = "", fill = "", size = "")
# 	) +
# ggsave("./figures/regression_behaviour_mockWT.png", 
# 		units = "in", 
# 		width = 6, 
# 		height = 4, 
# 		bg = "white")


# Within strain
marker_x <- 0
marker_y <- 0
cut_off <- 0.5

# Qurdant analysis using MDA using estimates
(scatter_plt <- fit_df %>%
	separate(A, remove = FALSE, into = c("genotype", "treatment"), 
		sep = "_", 
		convert = FALSE) %>%
	mutate(genotype = factor(as.character(genotype), levels = c("845", "992", "994")),
		treatment = factor(treatment, levels = c("mock", "flg22")),
		# significance = as.factor(p.value <= 0.05 | prl_pvalue <= 0.05)
		target = (abs(PRL) > cut_off & abs(SFW) > cut_off),
	) %>%
	mutate(text_target = ifelse(target == TRUE, as.character(strain), "")) %>%
	ggplot(aes(x = saturate(PRL), y = saturate(SFW))) +
	geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
	geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
	geom_point(aes(fill = strain, 
		# shape = genotype, 
		size = target, 
		colour = target, 
		alpha = target), shape = 22) +
	# ylim(c(-7, 7)) +
	# xlim(c(-4, 4)) +
	ggrepel::geom_text_repel(aes(label = text_target)) +
	facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
	# scale_fill_brewer(palette = "Spectral") +
	# scale_shape_manual(values = c(22,23,24)) +
	scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
	scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
	scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
	theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "log2(vs 845_mock) (Primary Root length)", y = "log2(vs 845_mock) (Shoot Fresh Weight)", 
	 	colour = "", shape = "", fill = "", size = "")
	) +
ggsave("./figures/regression.png", 
		units = "in", 
		width = 8, 
		height = 6, 
		bg = "white")


(scatter_plt_behav <- fit_df %>%
	separate(A, remove = FALSE, into = c("genotype", "treatment"), 
		sep = "_", 
		convert = FALSE) %>%
	mutate(genotype = factor(as.character(genotype), levels = c("845", "992", "994")),
		treatment = factor(treatment, levels = c("mock", "flg22")),
		# significance = as.factor(p.value <= 0.05 | prl_pvalue <= 0.05)
		target = (abs(PRL) > cut_off & abs(SFW) > cut_off),
	) %>%
	mutate(text_target = ifelse(target == TRUE, as.character(strain), "")) %>%
	ggplot(aes(x = saturate(behaviour_PRL), y = saturate(behaviour_PRL))) +
	geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
	geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
	geom_point(aes(fill = strain, 
		# shape = genotype, 
		size = target, 
		colour = target, 
		alpha = target), shape = 22) +
	# ylim(c(-7, 7)) +
	# xlim(c(-4, 4)) +
	ggrepel::geom_text_repel(aes(label = text_target)) +
	facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
	# scale_fill_brewer(palette = "Spectral") +
	# scale_shape_manual(values = c(22,23,24)) +
	scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
	scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
	scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
	theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "log2(vs 845_mock) (Primary Root length)", y = "log2(vs 845_mock) (Shoot Fresh Weight)", 
	 	colour = "", shape = "", fill = "", size = "")
	) +
ggsave("./figures/regression_behaviour.png", 
		units = "in", 
		width = 8, 
		height = 6, 
		bg = "white")

# ix <- abs(mean(fit_df$PRL[fit_df$strain == "mock"]))
# iy <- abs(mean(fit_df$SFW[fit_df$strain == "mock"]))

# Composite for the comparison against Mock treatment
# plot_obj_group <- parallel::mclapply(c(as.character(unique(fit_df_group$A))), function(i){

# 	res <- lapply(c(as.character(unique(fit_df$A))), function(j){

# 		# i = "845_flg22"
# 		# j = "992_mock"

# 		if(i != j){

# 			# PRL
# 			(plt_prl <- fit_df_group %>%
# 				select(-SFW, -behaviour_PRL, -behaviour_SFW) %>%
# 				filter(A %in% c(i,j)) %>%
# 				spread(key = A, value = PRL, fill = 0, convert = FALSE) %>%
# 				mutate(x = ifelse(abs(.[,2]) >= 1, sign(.[,2]) * 1, .[,2]), 
# 					y = ifelse(abs(.[,3]) >= 1, sign(.[,3]) * 1, .[,3])) %>%
# 				mutate(target = (abs(x) >= .2 & abs(y) >= .2),
# 					strain = factor(A_strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
# 				mutate(text_target = ifelse(abs(x) >= .2 & abs(y) >= .2, as.character(strain), "")) %>%
# 				ggplot(aes(x = saturate(x), y = saturate(y))) +
# 				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_point(aes(fill = strain, 
# 					# shape = genotype, 
# 					size = target, 
# 					colour = target, 
# 					alpha = target), shape = 22) +
# 				# ylim(c(-7, 7)) +
# 				# xlim(c(-4, 4)) +
# 				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
# 				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 				# scale_fill_brewer(palette = "Spectral") +
# 				# scale_shape_manual(values = c(22,23,24)) +
# 				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 				theme_classic() +
# 				 theme(legend.position = "top", 
# 				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 				 labs(x = paste0("PRL: LFC vs 845_mock ", i), 
# 				 	y = paste0("PRL: LFC vs 845_mock ", j), 
# 				 	colour = "", shape = "", fill = "", size = "")
# 				)

# 			# SFW
# 			(plt_sfw <- fit_df_group %>%
# 				select(-PRL, -behaviour_PRL, -behaviour_SFW) %>%
# 				filter(A %in% c(i,j)) %>%
# 				spread(key = A, value = SFW, fill = 0, convert = FALSE) %>%
# 				mutate(x = ifelse(abs(.[,2]) >= 2, sign(.[,2]) * 2, .[,2]), 
# 					y = ifelse(abs(.[,3]) >= 2, sign(.[,3]) * 2, .[,3]),
# 					strain = factor(A_strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
# 				mutate(target = ((abs(x) >= .5 & abs(y) >= .5))) %>%
# 				mutate(text_target = ifelse(abs(x) >= .5 & abs(y) >= .5, as.character(strain), "")) %>%
# 				ggplot(aes(x = saturate(x), y = saturate(y))) +
# 				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_point(aes(fill = strain, 
# 					# shape = genotype, 
# 					size = target, 
# 					colour = target, 
# 					alpha = target), shape = 22) +
# 				# ylim(c(-7, 7)) +
# 				# xlim(c(-4, 4)) +
# 				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
# 				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 				# scale_fill_brewer(palette = "Spectral") +
# 				# scale_shape_manual(values = c(22,23,24)) +
# 				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 				theme_classic() +
# 				 theme(legend.position = "top", 
# 				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 				 labs(x = paste0("SFW: LFC vs 845_mock ", i), 
# 				 	y = paste0("SFW: LFC vs 845_mock ", j), 
# 				 	colour = "", shape = "", fill = "", size = "")
# 				)

# 			plt <- list(plt_prl, plt_sfw)
# 			return(plt)
# 		}


# 	})
# 	return(res)


# }, mc.cores = 8)

# Composite figure for pairwise
plot_obj <- parallel::mclapply(c(as.character(unique(fit_df$A))), function(i){

	res <- lapply(c(as.character(unique(fit_df$A))), function(j){

		# i = "845_flg22"
		# j = "992_mock"

		if(i != j){

			# PRL
			(plt_prl <- fit_df %>%
				select(-SFW, -behaviour_PRL, -behaviour_SFW) %>%
				filter(A %in% c(i,j)) %>%
				spread(key = A, value = PRL, fill = 0, convert = FALSE) %>%
				mutate(x = ifelse(abs(.[,2]) >= 1, sign(.[,2]) * 1, .[,2]), 
					y = ifelse(abs(.[,3]) >= 1, sign(.[,3]) * 1, .[,3])) %>%
				mutate(target = (abs(x) >= .2 & abs(y) >= .2),
					strain = factor(strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
				mutate(text_target = ifelse(abs(x) >= .2 & abs(y) >= .2, as.character(strain), "")) %>%
				ggplot(aes(x = saturate(x), y = saturate(y))) +
				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_point(aes(fill = strain, 
					# shape = genotype, 
					size = target, 
					colour = target, 
					alpha = target), shape = 22) +
				# ylim(c(-7, 7)) +
				# xlim(c(-4, 4)) +
				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
				# scale_fill_brewer(palette = "Spectral") +
				# scale_shape_manual(values = c(22,23,24)) +
				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
				theme_classic() +
				 theme(legend.position = "top", 
				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
				 labs(x = paste0("PRL: LFC vs 845_mock ", i), 
				 	y = paste0("PRL: LFC vs 845_mock ", j), 
				 	colour = "", shape = "", fill = "", size = "")
				)

			# SFW
			(plt_sfw <- fit_df %>%
				select(-PRL, -behaviour_PRL, -behaviour_SFW) %>%
				filter(A %in% c(i,j)) %>%
				spread(key = A, value = SFW, fill = 0, convert = FALSE) %>%
				mutate(x = ifelse(abs(.[,2]) >= 2, sign(.[,2]) * 2, .[,2]), 
					y = ifelse(abs(.[,3]) >= 2, sign(.[,3]) * 2, .[,3])) %>%
				mutate(target = (abs(x) >= .5 & abs(y) >= .5),
					strain = factor(strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
				mutate(text_target = ifelse(abs(x) >= cut_off & abs(cut_off) >= cut_off, as.character(strain), "")) %>%
				ggplot(aes(x = saturate(x), y = saturate(y))) +
				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_point(aes(fill = strain, 
					# shape = genotype, 
					size = target, 
					colour = target, 
					alpha = target), shape = 22) +
				# ylim(c(-7, 7)) +
				# xlim(c(-4, 4)) +
				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
				# scale_fill_brewer(palette = "Spectral") +
				# scale_shape_manual(values = c(22,23,24)) +
				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
				theme_classic() +
				 theme(legend.position = "top", 
				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
				 labs(x = paste0("SFW: LFC vs 845_mock ", i), 
				 	y = paste0("SFW: LFC vs 845_mock ", j), 
				 	colour = "", shape = "", fill = "", size = "")
				)

			plt <- list(plt_prl, plt_sfw)
			return(plt)
		}


	})
	return(res)


}, mc.cores = 4)

names(plot_obj) <- c(as.character(unique(fit_df$A)))
# names(plot_obj_group) <- c(as.character(unique(fit_df$A)))


# Mock control
# plot_obj_behav_group <- parallel::mclapply(c(as.character(unique(fit_df_group$A))), function(i){

# 	res <- lapply(c(as.character(unique(fit_df_group$A))), function(j){

# 		# i = "845_flg22"
# 		# j = "992_mock"

# 		if(i != j){

# 			# PRL
# 			(plt_prl <- fit_df_group %>%
# 				select(-behaviour_SFW, -PRL, -SFW) %>%
# 				filter(A %in% c(i,j)) %>%
# 				spread(key = A, value = behaviour_PRL, fill = 0, convert = FALSE) %>%
# 				mutate(x = ifelse(abs(.[,2]) >= 1, sign(.[,2]) * 1, .[,2]), 
# 					y = ifelse(abs(.[,3]) >= 1, sign(.[,3]) * 1, .[,3])) %>%
# 				mutate(target = (abs(x) >= .2 & abs(y) >= .2),
# 					strain = factor(A_strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
# 				mutate(text_target = ifelse(abs(x) >= .2 & abs(y) >= .2, as.character(strain), "")) %>%
# 				ggplot(aes(x = saturate(x), y = saturate(y))) +
# 				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_point(aes(fill = strain, 
# 					# shape = genotype, 
# 					size = target, 
# 					colour = target, 
# 					alpha = target), shape = 22) +
# 				# ylim(c(-7, 7)) +
# 				# xlim(c(-4, 4)) +
# 				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
# 				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 				# scale_fill_brewer(palette = "Spectral") +
# 				# scale_shape_manual(values = c(22,23,24)) +
# 				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 				theme_classic() +
# 				 theme(legend.position = "top", 
# 				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 				 labs(x = paste0("PRL: LFC vs 845_mock ", i), 
# 				 	y = paste0("PRL: LFC vs 845_mock ", j), 
# 				 	colour = "", shape = "", fill = "", size = "")
# 				)

# 			# SFW
# 			(plt_sfw <- fit_df_group %>%
# 				select(-behaviour_PRL, -PRL, -SFW) %>%
# 				filter(A %in% c(i,j)) %>%
# 				spread(key = A, value = behaviour_SFW, fill = 0, convert = FALSE) %>%
# 				mutate(x = ifelse(abs(.[,2]) >= 2, sign(.[,2]) * 2, .[,2]), 
# 					y = ifelse(abs(.[,3]) >= 2, sign(.[,3]) * 2, .[,3])) %>%
# 				mutate(target = (abs(x) >= .5 & abs(y) >= .5),
# 					strain = factor(A_strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
# 				mutate(text_target = ifelse(abs(x) >= cut_off & abs(cut_off) >= cut_off, as.character(strain), "")) %>%
# 				ggplot(aes(x = saturate(x), y = saturate(y))) +
# 				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
# 				geom_point(aes(fill = strain, 
# 					# shape = genotype, 
# 					size = target, 
# 					colour = target, 
# 					alpha = target), shape = 22) +
# 				# ylim(c(-7, 7)) +
# 				# xlim(c(-4, 4)) +
# 				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
# 				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
# 				# scale_fill_brewer(palette = "Spectral") +
# 				# scale_shape_manual(values = c(22,23,24)) +
# 				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
# 				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
# 				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
# 				theme_classic() +
# 				 theme(legend.position = "top", 
# 				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
# 				 labs(x = paste0("SFW: LFC vs 845_mock ", i), 
# 				 	y = paste0("SFW: LFC vs 845_mock ", j), 
# 				 	colour = "", shape = "", fill = "", size = "")
# 				)

# 			plt <- list(plt_prl, plt_sfw)
# 			return(plt)
# 		}


# 	})
# 	return(res)


# }, mc.cores = 4)

# names(plot_obj_behav) <- c(as.character(unique(fit_df$A)))

# Growth promoting/detrimental
plot_obj_behav <- parallel::mclapply(c(as.character(unique(fit_df$A))), function(i){

	res <- lapply(c(as.character(unique(fit_df$A))), function(j){

		# i = "845_flg22"
		# j = "992_mock"

		if(i != j){

			# PRL
			(plt_prl <- fit_df %>%
				select(-behaviour_SFW, -PRL, -SFW) %>%
				filter(A %in% c(i,j)) %>%
				spread(key = A, value = behaviour_PRL, fill = 0, convert = FALSE) %>%
				mutate(x = ifelse(abs(.[,2]) >= 1, sign(.[,2]) * 1, .[,2]), 
					y = ifelse(abs(.[,3]) >= 1, sign(.[,3]) * 1, .[,3])) %>%
				mutate(target = (abs(x) >= .2 & abs(y) >= .2),
					strain = factor(strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
				mutate(text_target = ifelse(abs(x) >= .2 & abs(y) >= .2, as.character(strain), "")) %>%
				ggplot(aes(x = saturate(x), y = saturate(y))) +
				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_point(aes(fill = strain, 
					# shape = genotype, 
					size = target, 
					colour = target, 
					alpha = target), shape = 22) +
				# ylim(c(-7, 7)) +
				# xlim(c(-4, 4)) +
				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
				# scale_fill_brewer(palette = "Spectral") +
				# scale_shape_manual(values = c(22,23,24)) +
				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
				theme_classic() +
				 theme(legend.position = "top", 
				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
				 labs(x = paste0("PRL: LFC vs 845_mock ", i), 
				 	y = paste0("PRL: LFC vs 845_mock ", j), 
				 	colour = "", shape = "", fill = "", size = "")
				)

			# SFW
			(plt_sfw <- fit_df %>%
				select(-behaviour_PRL, -PRL, -SFW) %>%
				filter(A %in% c(i,j)) %>%
				spread(key = A, value = behaviour_SFW, fill = 0, convert = FALSE) %>%
				mutate(x = ifelse(abs(.[,2]) >= 2, sign(.[,2]) * 2, .[,2]), 
					y = ifelse(abs(.[,3]) >= 2, sign(.[,3]) * 2, .[,3])) %>%
				mutate(target = (abs(x) >= .5 & abs(y) >= .5),
					strain = factor(strain, levels = strain_sorted_phy, labels= strain_labs_phy)) %>%
				mutate(text_target = ifelse(abs(x) >= cut_off & abs(cut_off) >= cut_off, as.character(strain), "")) %>%
				ggplot(aes(x = saturate(x), y = saturate(y))) +
				geom_hline(yintercept = c(marker_x), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_vline(xintercept = c(marker_y), lty = "solid", lwd = 1, colour = "darkgrey") +
				geom_point(aes(fill = strain, 
					# shape = genotype, 
					size = target, 
					colour = target, 
					alpha = target), shape = 22) +
				# ylim(c(-7, 7)) +
				# xlim(c(-4, 4)) +
				ggrepel::geom_text_repel(aes(label = text_target), size = 2) +
				# facet_grid(.~ genotype, scale = "fixed", switch = "x", space = "free") +
				# scale_fill_brewer(palette = "Spectral") +
				# scale_shape_manual(values = c(22,23,24)) +
				scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2), guide = FALSE) +
				scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
				scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
				theme_classic() +
				 theme(legend.position = "top", 
				 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
				 labs(x = paste0("SFW: LFC vs 845_mock ", i), 
				 	y = paste0("SFW: LFC vs 845_mock ", j), 
				 	colour = "", shape = "", fill = "", size = "")
				)

			plt <- list(plt_prl, plt_sfw)
			return(plt)
		}


	})
	return(res)


}, mc.cores = 4)

names(plot_obj_behav) <- c(as.character(unique(fit_df$A)))
# names(plot_obj_behav_group) <- c(as.character(unique(fit_df_group$A)))



# Composite figure
# ## Get legend on seperate composite
# composite <- cowplot::plot_grid(

# 	bxp_sfw + theme(legend.position = "none"), 
# 	hmap_sfw + theme(legend.position = "none"), 
# 	hmap_prl + theme(legend.position = "none"), 
# 	bxp_prl + theme(legend.position = "none"),
# 	ncol = 4, nrow = 1,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(2, 0.6, 0.6, 2), 
# 	rel_heights = c(1, 1, 1, 1),
# 	scale = c(1, 1, 1, 1),
# 	greedy = FALSE

# )

# leg <- cowplot::plot_grid(
# 	cowplot::get_legend(bxp_sfw),
# 	cowplot::get_legend(hmap_sfw),
# 	cowplot::get_legend(hmap_prl),
# 	ncol = 3, 
# 	nrow = 1, 
# 	# byrow  = TRUE,
# 	# align = "hv",
# 	# axis = "tblr",
# 	rel_widths = c(1, 0.5, 0.5), 
# 	rel_heights = c(1, 1, 1),
# 	# scale = c(1, 1, 1),
# 	greedy = FALSE)

# cowplot::plot_grid(

# 	leg,
# 	composite,
# 	ncol = 1, nrow = 2, 
# 	byrow  = TRUE,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(0.8, 1), 
# 	rel_heights = c(.2, 1),
# 	scale = c(0.8, 1),
# 	hjust = 0.5, 
# 	vjust = 0.5,
# 	greedy = TRUE

# ) + ggsave( 
# 	filename = "./figures/summary_figure.png", 
# 	device = "png", 
# 	units = "in",
# 	width = 21, 
# 	height = 9, 
# 	bg = "white")

# Composite plot for phylogeny
# Composite figure
## Get legend on seperate composite
# composite <- cowplot::plot_grid(

# 	bxp_sfw_phy + theme(legend.position = "none"), 
# 	hmap_sfw_phy + theme(legend.position = "none"), 
# 	hmap_prl_phy + theme(legend.position = "none"), 
# 	bxp_prl_phy + theme(legend.position = "none"),
# 	ncol = 4, nrow = 1,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(2, 0.6, 0.6, 2), 
# 	rel_heights = c(1, 1, 1, 1),
# 	scale = c(1, 1, 1, 1),
# 	greedy = FALSE

# )

# leg <- cowplot::plot_grid(
# 	cowplot::get_legend(bxp_sfw_phy),
# 	cowplot::get_legend(hmap_sfw_phy),
# 	cowplot::get_legend(hmap_prl_phy),
# 	ncol = 3, 
# 	nrow = 1, 
# 	# byrow  = TRUE,
# 	# align = "hv",
# 	# axis = "tblr",
# 	rel_widths = c(1, 0.5, 0.5), 
# 	rel_heights = c(1, 1, 1),
# 	# scale = c(1, 1, 1),
# 	greedy = FALSE)

# cowplot::plot_grid(
	
# 	leg,
# 	composite,
# 	ncol = 1, nrow = 2, 
# 	byrow  = TRUE,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(0.8, 1), 
# 	rel_heights = c(.2, 1),
# 	scale = c(0.8, 1),
# 	hjust = 0.5, 
# 	vjust = 0.5,
# 	greedy = TRUE

# ) + ggsave( 
# 	filename = "./figures/summary_figure_phylogeny.png", 
# 	device = "png", 
# 	units = "in",
# 	width = 21, 
# 	height = 8, 
# 	bg = "white")

# Pairwise plot
pw_composite <- cowplot::plot_grid(

	plot_obj[["845_flg22"]][[2]][[1]] + theme(legend.position = "none"), 
	plot_obj[["845_flg22"]][[3]][[1]] + theme(legend.position = "none"),
	plot_obj[["992_mock"]][[3]][[1]] + theme(legend.position = "none"),
	plot_obj[["845_flg22"]][[2]][[2]] + theme(legend.position = "none"), 
	plot_obj[["845_flg22"]][[3]][[2]] + theme(legend.position = "none"),
	plot_obj[["992_mock"]][[3]][[2]] + theme(legend.position = "none"),
	ncol = 3, 
	nrow = 2, 
	byrow  = TRUE,
	align = "hv",
	axis = "tblr",
	rel_widths = c(1, 1, 1), 
	rel_heights = c(1, 1),
	# scale = c(1, 0.9, 0.9, 1),
	greedy = FALSE)

pw_composite <- cowplot::plot_grid(

	cowplot::get_legend(plot_obj[["845_flg22"]][[2]][[1]]),
	pw_composite,
	ncol = 1, nrow = 2, 
	byrow  = TRUE,
	align = "hv",
	axis = "tblr",
	rel_widths = c(1, 1), 
	rel_heights = c(.1, 1),
	scale = c(0.001, 1),
	greedy = FALSE

)

ggsave(pw_composite, 
	filename = "./figures/summary_regression_pwise.png", 
	device = "png", 
	units = "in",
	width = 8, 
	height = 6, 
	bg = "white")

pw_composite <- cowplot::plot_grid(

	plot_obj_behav[["845_flg22"]][[2]][[1]] + theme(legend.position = "none"), 
	plot_obj_behav[["845_flg22"]][[3]][[1]] + theme(legend.position = "none"),
	plot_obj_behav[["992_mock"]][[3]][[1]] + theme(legend.position = "none"),
	plot_obj_behav[["845_flg22"]][[2]][[2]] + theme(legend.position = "none"), 
	plot_obj_behav[["845_flg22"]][[3]][[2]] + theme(legend.position = "none"),
	plot_obj_behav[["992_mock"]][[3]][[2]] + theme(legend.position = "none"),
	ncol = 3, nrow = 2, 
	byrow  = TRUE,
	align = "hv",
	axis = "tblr",
	rel_widths = c(1, 1, 1), 
	rel_heights = c(1, 1),
	# scale = c(1, 0.9, 0.9, 1),
	greedy = FALSE

)

pw_composite <- cowplot::plot_grid(

	cowplot::get_legend(plot_obj_behav[["845_flg22"]][[2]][[1]]),
	pw_composite,
	ncol = 1, nrow = 2, 
	byrow  = TRUE,
	align = "hv",
	axis = "tblr",
	rel_widths = c(1, 1), 
	rel_heights = c(.1, 1),
	scale = c(0.2, 1),
	greedy = FALSE

)

ggsave(pw_composite, 
	filename = "./figures/summary_regression_behaviour_pwise.png", 
	device = "png", 
	units = "in",
	width = 8, 
	height = 6, 
	bg = "white")

# HEatmap with tree
# hmap_composite <- cowplot::plot_grid(

# 	# taxa_obj$tree_graph + theme(legend.position = "none"), 
# 	hmap_sfw_phy_2 + coord_fixed() + theme(legend.position = "none"), 
# 	hmap_prl_phy_2 + coord_fixed() + theme(legend.position = "none", axis.text.y = element_blank()), 
# 	ncol = 2, nrow = 1,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(0.5, 0.5), 
# 	rel_heights = c(1, 1),
# 	scale = c(1, 1),
# 	greedy = FALSE

# )

# leg <- 
# 	cowplot::plot_grid(
# 	cowplot::get_legend(hmap_sfw_phy_2),
# 	cowplot::get_legend(hmap_prl_phy_2),
# 	ncol = 2, 
# 	nrow = 1, 
# 	# byrow  = TRUE,
# 	# align = "hv",
# 	# axis = "tblr",
# 	rel_widths = c(0.5, 0.5), 
# 	rel_heights = c(1, 1),
# 	# scale = c(1, 1, 1),
# 	greedy = FALSE)

# cowplot::plot_grid(
	
# 	leg,
# 	hmap_composite,
# 	ncol = 1, 
# 	nrow = 2, 
# 	byrow  = TRUE,
# 	align = "hv",
# 	axis = "tblr",
# 	rel_widths = c(1, 1), 
# 	rel_heights = c(.2, 1),
# 	scale = c(0.8, 1),
# 	hjust = 0.5, 
# 	vjust = 0.5,
# 	greedy = TRUE

# ) + 
# ggsave( 
# 	filename = "./figures/hmap_phylogeny_2.png", 
# 	device = "png", 
# 	units = "in",
# 	width = 3, 
# 	height = 7, 
# 	bg = "white")

# HEatmap with tree --using aplot
require(aplot)
require(ggtree)
tree_phy <- taxa_obj$tree_graph
composite <- (hmap_phy_2 + theme(legend.position = "none", axis.text.x = element_blank())) %>% 
insert_left(plot = (tree_phy + 
	xlim_tree(1) + 
	theme(legend.position = "none")), width = 0.5) # %>%
# insert_right(plot = (hmap_phy_mock_2 + theme(legend.position = "none", axis.text.x = element_blank())))

ggsave(composite, 
	filename = "./figures/hmap_phylogeny.png", 
	device = "png", 
	units = "in",
	width = 2.5, 
	height = 7, 
	bg = "white")

composite <- (hmap_phy_mock_2 + theme(legend.position = "none", axis.text.x = element_blank())) %>% 
insert_left(plot = (tree_phy + 
	xlim_tree(1) + 
	theme(legend.position = "none")), width = 0.5)


# composite <- (tree_phy + xlim_tree(1) + theme(legend.position = "none")) + 
# ylim2(hmap_phy_2 + theme(legend.position = "none", axis.text.x = element_blank()))

ggsave(composite, 
	filename = "./figures/hmap_phylogeny_mock.png", 
	device = "png", 
	units = "in",
	width = 3.5, 
	height = 7, 
	bg = "white")



# END
sessionInfo()
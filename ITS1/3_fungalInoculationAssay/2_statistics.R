#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript

# Script by
# @Arpan Kumar Basak
# Quantification of fungal inoculation on plants that are:
# - Treated with flg22 to trigger immune response 10nM
# - Mutants that lack in PYK10 myrosinases and Trp-pathway

# Script for statistical analysis

rm(list = ls())

# Read data
source("./scripts/parameters.R")
setwd(analysis.fungali)
require(tidyverse)

load("./data/data.RData")
load("./data/taxonomy.RData")

# Load the data
df_prl <- dat$prl %>% filter(strain != "F131")
df <- dat$sfw %>% filter(strain != "F131")

set1_strains <- intersect(unique(df_prl$strain), unique(df$strain))

# Diagnose normality
require(qqplotr)
df_prl %>% 
	ggplot(aes(sample = prl)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free")
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/prl_diagnosis_qq.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

# Log2+1 transformed
df_prl %>% 
	mutate(prl = log2(prl+1)) %>%
	ggplot(aes(sample = prl)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/prl_diagnosis_qq_logmaxT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

df_prl %>% 
	mutate(prl = log2(prl+1e-12)) %>%
	ggplot(aes(sample = prl)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/prl_diagnosis_qq_logminT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

df_prl %>% 
	mutate(prl = sqrt(prl)) %>%
	ggplot(aes(sample = prl)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/prl_diagnosis_qq_srtT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

# For shoot freshweiight
df %>% 
	ggplot(aes(sample = shoot_fresh_weight)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/sfw_diagnosis_qq.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

# Log2+1 transformed
df %>% 
	mutate(shoot_fresh_weight = log2(shoot_fresh_weight+1)) %>%
	ggplot(aes(sample = shoot_fresh_weight)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/sfw_diagnosis_qq_logmaxT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

df %>% 
	mutate(shoot_fresh_weight = log2(shoot_fresh_weight+1e-12)) %>%
	ggplot(aes(sample = shoot_fresh_weight)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/diagnosis_qq_logminT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

df %>% 
	mutate(shoot_fresh_weight = sqrt(shoot_fresh_weight)) %>%
	ggplot(aes(sample = shoot_fresh_weight)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(.~ strain, as.table = TRUE, scale = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_classic() +
	theme(
		legend.position = "top", 
		legend.text= element_text(angle = 60, hjust = 0.5, vjust = 0.5),
		axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
	labs(x = "", y = "", colour = "") +
	ggsave("./figures/sfw_sfw_diagnosis_qq_srtT.png", 
		units = "in", 
		width = 7, 
		height = 7.5, 
		bg = "white")

# Conduct strain wise comparative satitstics using GLM
res <- parallel::mclapply(set1_strains, function(x){

	# x = "F131"
	temp <- df %>% filter(as.character(strain) == x) %>% 
	# filter(shoot_fresh_weight > 0) %>%
	mutate(shoot_fresh_weight = sqrt(1000*shoot_fresh_weight), 
		random = as.factor(as.character(rep)))
	

	mod <- lme4::lmer(shoot_fresh_weight ~ 0 + group + (1|random), data = temp)
	fit <- 
		broom::tidy(summary(multcomp::glht(mod,linfct=multcomp::mcp(group="Tukey")))) %>%
		separate(contrast, into = c("A", "B"), sep = " - ") %>%
		mutate(strain = x) %>%
		data.frame(.)

	return(fit)

}, mc.cores = 4)

# Diagnose normality

# Stats for primary root length
res_prl <- parallel::mclapply(set1_strains, function(x){

	# x = "F91"
	temp <- df_prl %>% filter(as.character(strain) == x) %>% 
	# filter(prl > 0) %>%
	mutate(prl = sqrt(prl), 
		random = as.factor(as.character(replicate)))
	

	mod <- lme4::lmer(prl ~ 0 + group + (1|random), data = temp)
	fit <- 
		broom::tidy(summary(multcomp::glht(mod,linfct=multcomp::mcp(group="Tukey")))) %>%
		separate(contrast, into = c("A", "B"), sep = " - ") %>%
		mutate(strain = x) %>%
		data.frame(.)

	return(fit)

}, mc.cores = 4)

# Concatenate statistical data
fit <- do.call(rbind.data.frame, res) %>%
	mutate(significance =  adj.p.value < 0.05, 
		FDR = p.adjust(adj.p.value, method = "BH")) %>%
	filter(B == "845_mock") %>%
	data.frame(., stringsAsFactors = FALSE)

fit_prl <- do.call(rbind.data.frame, res_prl) %>%
	mutate(significance =  adj.p.value < 0.05, FDR = p.adjust(adj.p.value, method = "BH")) %>%
	filter(B == "845_mock") %>%
	data.frame(., stringsAsFactors = FALSE)


# Wether there are strains that differ when compared to mock treatment
# mod_group <- aov(glm(sqrt(1000*shoot_fresh_weight) ~ 0 + strain:group + random, data = df))
# fit_group <- broom::tidy(TukeyHSD(mod_group, method = "BH")) %>%
# filter(term == "strain:group") %>%
# separate(contrast, into = c("A", "B"), sep = "-") %>%
# separate(A, into = c("A_strain", "A_group"), sep = ":") %>%
# separate(B, into = c("B_strain", "B_group"), sep = ":") %>%
# separate(A_group, into = c("A_genotype", "A_treatment"), sep = "_") %>%
# # separate(B_group, into = c("B_genotype", "B_treatment"), sep = "_") %>%
# filter(B_strain == "mock", B_group == "845_mock") %>%
# mutate(significance = adj.p.value < 0.05 & abs(estimate) > 1e-3) %>%
# data.frame(., stringsAsFactors = FALSE)

# mod_group_prl <- aov(glm(sqrt(prl) ~ 0 + strain:group + random, data = df_prl))
# fit_group_prl <- broom::tidy(TukeyHSD(mod_group_prl, method = "BH")) %>%
# filter(term == "strain:group") %>%
# separate(contrast, into = c("A", "B"), sep = "-") %>%
# separate(A, into = c("A_strain", "A_group"), sep = ":") %>%
# separate(B, into = c("B_strain", "B_group"), sep = ":") %>%
# separate(A_group, into = c("A_genotype", "A_treatment"), sep = "_") %>%
# # separate(B_group, into = c("B_genotype", "B_treatment"), sep = "_") %>%
# filter(B_strain == "mock", B_group == "845_mock") %>%
# mutate(significance = adj.p.value < 0.05 & abs(estimate) > .1) %>%
# data.frame(., stringsAsFactors = FALSE)

# Sort on the basis of mock control
# fit_mat_group <- fit_group %>% 
# 	select(estimate, A_genotype, A_treatment, A_strain) %>%
# 	mutate(A = paste(A_genotype, A_treatment, sep = "_")) %>%
# 	select(-A_treatment, -A_genotype) %>%
# 	spread(key = A, value = estimate, convert = FALSE, fill = 0) %>%
# 	data.frame(.)

# fit_mat_group_p <- fit_group %>% 
# 	select(adj.p.value, A_genotype, A_treatment, A_strain) %>%
# 	mutate(A = paste(A_genotype, A_treatment, sep = "_")) %>%
# 	select(-A_treatment, -A_genotype) %>%
# 	spread(key = A, value = adj.p.value, convert = FALSE, fill = 0) %>%
# 	data.frame(.)

# colnames(fit_mat_group) <- str_replace_all(colnames(fit_mat_group), "X", "SFW:")
# colnames(fit_mat_group_p) <- str_replace_all(colnames(fit_mat_group_p), "X", "pvalSFW:")

# fit_mat_group_prl <- fit_group_prl %>% 
# 	select(estimate, A_genotype, A_treatment, A_strain) %>%
# 	mutate(A = paste(A_genotype, A_treatment, sep = "_")) %>%
# 	select(-A_treatment, -A_genotype) %>%
# 	spread(key = A, value = estimate, convert = FALSE, fill = 0) %>%
# 	data.frame(.)

# fit_mat_group_prl_p <- fit_group %>% 
# 	select(adj.p.value, A_genotype, A_treatment, A_strain) %>%
# 	mutate(A = paste(A_genotype, A_treatment, sep = "_")) %>%
# 	select(-A_treatment, -A_genotype) %>%
# 	spread(key = A, value = adj.p.value, convert = FALSE, fill = 0) %>%
# 	data.frame(.)

# colnames(fit_mat_group_prl) <- str_replace_all(colnames(fit_mat_group_prl), "X", "PRL:")
# colnames(fit_mat_group_prl_p) <- str_replace_all(colnames(fit_mat_group_prl_p), "X", "pvalPRL:")

# idx <- match(fit_mat_group$A_strain, fit_mat_group_prl$A_strain)
# fit_mat <- cbind.data.frame(fit_mat_group, fit_mat_group_prl[idx, -1])

# # To make clusters of similar behaving strains considering both SFW and PRL
# row.names(fit_mat) <- fit_mat$A_strain
# d <- 1-cor(t(fit_mat[,-1]))

# hc <- hclust(as.dist(d), "average")
# strain_sorted_group <- row.names(fit_mat)[hc$order]
# strain_labs_group <- strain_sorted_group

# labs <- taxa_obj$taxa$ID
# labs <- c(as.character(taxa_obj$taxa$ID), out)
# names(labs) <- c(taxa_obj$taxa$strain_ids, out)

# labs <- labs[strain_sorted]

# # Sort by phylogeny
# out <- setdiff(fit_mat$A_strain, taxa_obj$strain_phylogeny)
# strain_sorted_phy <- c(as.character(taxa_obj$strain_phylogeny), out)
# strain_labs_phy <- c(as.character(taxa_obj$taxa$ID), out)
# strain_class <- c(as.character(taxa_obj$taxa$class), out)
# names(strain_class) <- c(taxa_obj$taxa$strain_ids, out)
# levs <- c(
# 	"Sordariomycetes",
# 	"Leotiomycetes",
# 	"Dothideomycetes",
# 	"Agaricomycetes",
# 	"Mortierellomycetes", 
# 	"mock")

# idx <- match(fit_mat_group$A_strain, fit_mat_group_prl$A_strain)
# fit_mat <- cbind.data.frame(fit_mat_group, fit_mat_group_prl[idx, -1])

# idx <- match(fit_mat_group_p$A_strain, fit_mat_group_prl_p$A_strain)
# fit_mat_p <- cbind.data.frame(fit_mat_group_p, fit_mat_group_prl_p[idx, -1])

# idx <- match(fit_mat$A_strain, fit_mat_p$A_strain)
# fit_mat <- cbind.data.frame(fit_mat, fit_mat_p[idx, -1])


# fit_temp <- cbind.data.frame(
# 	id = strain_sorted_phy, 
# 	strain = strain_labs_phy, 
# 	class = strain_class,
# 	fit_mat[match(strain_sorted_phy, fit_mat$A_strain),-1]
# 	) %>%
# gather(key = "key", value = "value", -id, -strain, -class) %>%
# separate(key, into = c("val", "group"), sep = ":", remove = FALSE, convert = FALSE) %>%
# mutate(temp = ifelse(grepl(.$val, pattern = "^pval"), "Pvalue", "Diff")) %>%
# mutate(val = str_replace_all(val, "pval", "")) %>%
# select(-key) %>%
# spread(key = temp, value = value, fill=NA) %>%
# mutate(FDR = Pvalue < 0.05 & abs(Diff) > 0.05,
# 	significance = Pvalue < 0.05,
# 	diff_sat = ifelse(abs(Diff) > 2, sign(Diff) * 2, Diff)) %>%
# separate(group, into = c("genotype", "treatment"), sep = "_", convert = FALSE, remove = FALSE) %>%
# mutate(
# 	genotype = factor(genotype, levels = c("845", "992", "994"), labels = c("Col", "pyk10", "cyp")),
# 	treatment = factor(treatment, levels = c("mock", "flg22")),
# 	strain = factor(strain, levels = (strain_labs_phy)),
# 	class = factor(class, levels = levs)) %>%
# data.frame(., stringsAsFactors = FALSE)



# Within strain comparison
fit_mat <- fit %>% 
	select(estimate, A, strain) %>%
	spread(key = A, value = estimate, convert = FALSE, fill = 0) %>%
	data.frame(.)
colnames(fit_mat) <- str_replace_all(colnames(fit_mat), "X", "SFW:")

fit_prl_mat <- fit_prl %>% 
	select(estimate, A, strain) %>%
	spread(key = A, value = estimate, convert = FALSE, fill = 0) %>%
	data.frame(.)
colnames(fit_prl_mat) <- str_replace_all(colnames(fit_prl_mat), "X", "PRL:")

fit_mat_p <- fit %>% 
	select(adj.p.value, A, strain) %>%
	spread(key = A, value = adj.p.value, convert = FALSE, fill = 0) %>%
	data.frame(.)
colnames(fit_mat_p) <- str_replace_all(colnames(fit_mat_p), "X", "pvalSFW:")

fit_prl_mat_p <- fit_prl %>% 
	select(adj.p.value, A, strain) %>%
	spread(key = A, value = adj.p.value, convert = FALSE, fill = 0) %>%
	data.frame(.)
colnames(fit_prl_mat_p) <- str_replace_all(colnames(fit_prl_mat_p), "X", "pvalPRL:")

idx <- match(fit_mat$strain, fit_prl_mat$strain)
fit_mat <- cbind.data.frame(fit_mat, fit_prl_mat[idx, -1])

# To make clusters of similar behaving strains considering both SFW and PRL
row.names(fit_mat) <- fit_mat$strain
matrix_lfc <- apply(fit_mat[,-1], 2, as.numeric)
row.names(matrix_lfc) <- fit_mat$strain
d <- 1-cor(t(matrix_lfc))

hc <- hclust(as.dist(d), "average")
strain_sorted <- row.names(matrix_lfc)[hc$order]
strain_labs <- strain_sorted

# labs <- taxa_obj$taxa$ID
out <- setdiff(fit_mat$strain, taxa_obj$strain_phylogeny)
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

row.names(matrix_lfc) <- str_replace_all(row.names(matrix_lfc), "\\.Wt", "")
save(matrix_lfc, file = "./data/lfc_matrix_prl.sfw.RData")


# idx <- match(fit_mat$strain, fit_mat_prl$strain)
# fit_mat <- cbind.data.frame(fit_mat, fit_mat_prl[idx, -1])

idx <- match(fit_mat_p$strain, fit_prl_mat_p$strain)
fit_mat_p <- cbind.data.frame(fit_mat_p, fit_prl_mat_p[idx, -1])

idx <- match(fit_mat$strain, fit_mat_p$strain)
fit_mat <- cbind.data.frame(fit_mat, fit_mat_p[idx, -1])


fit_temp_strain <- cbind.data.frame(
	id = strain_sorted_phy, 
	strain = strain_labs_phy, 
	class = strain_class,
	fit_mat[match(strain_sorted_phy, fit_mat$strain),-1]
	) %>%
gather(key = "key", value = "value", -id, -strain, -class) %>%
separate(key, into = c("val", "group"), sep = ":", remove = FALSE, convert = FALSE) %>%
mutate(temp = ifelse(grepl(.$val, pattern = "^pval"), "Pvalue", "Diff")) %>%
mutate(val = str_replace_all(val, "pval", "")) %>%
select(-key) %>%
spread(key = temp, value = value, fill=NA) %>%
mutate(FDR = Pvalue < 0.05 & abs(Diff) > 0.05,
	significance = Pvalue < 0.05,
	diff_sat = ifelse(abs(Diff) > 2, sign(Diff) * 2, Diff)) %>%
separate(group, into = c("genotype", "treatment"), sep = "_", convert = FALSE, remove = FALSE) %>%
mutate(
	genotype = factor(genotype, levels = c("845", "992", "994"), labels = c("Col", "pyk10", "cyp")),
	treatment = factor(treatment, levels = c("mock", "flg22")),
	strain = factor(strain, levels = (strain_labs_phy)),
	class = factor(class, levels = levs)) %>%
data.frame(., stringsAsFactors = FALSE)


# Within strain

# idx <- match(fit_mat$A_strain, fit_mat_group_prl$A_strain)
# fit_mat <- cbind.data.frame(fit_mat, fit_mat_group_prl[idx, -1])

fit_df <- fit_mat %>%
gather(key = "key", value = "vals", convert = FALSE, -strain) %>%
separate(key, into = c("measure", "A"), sep = ":") %>%
spread(key = measure, value = vals, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

# 
mock_lev <- fit_df %>% 
group_by(strain) %>%
filter(strain == "mock") %>%
ungroup() %>%
data.frame(., stringsAsFactors = FALSE)

fit_df$behaviour_PRL <- ((fit_df$PRL^2 - mock_lev$PRL[match(fit_df$A, mock_lev$A)]^2))
# fit_df$behaviour_PRL <- ifelse(fit_df$PRL == 1, min(fit_df$behaviour_PRL), fit_df$behaviour_PRL)

fit_df$behaviour_SFW <- ((fit_df$SFW^2 - mock_lev$SFW[match(fit_df$A, mock_lev$A)]^2))
# fit_df$behaviour_SFW <- ifelse(fit_df$SFW == 1, min(fit_df$behaviour_SFW), fit_df$behaviour_SFW)

sort_list <- list(
	strain_sorted_phy = strain_sorted_phy, 
	strain_sorted = strain_sorted
	)

stats <- list(
	stat_fit = fit_df, 
	stat_strain = fit_temp_strain
	)

save(list = c("stats", "sort_list"), file = "./data/stats.RData")

# END
sessionInfo()
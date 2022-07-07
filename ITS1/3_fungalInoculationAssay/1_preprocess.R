#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript

# Script by
# @Arpan Kumar Basak
# Quantification of fungal inoculation on plants that are:
# - Treated with flg22 to trigger immune response 10nM
# - Mutants that lack in PYK10 myrosinases and Trp-pathway

# Script for preprocessing the data

rm(list = ls())

# Update this into the parameter script
# Functions
saturate <- function(x){
    max <- quantile(x, .99)
    min <- quantile(x, .01)
    
    idx <- (x > max)
    x[idx] <- max
    
    idx <- x < min
    x[idx] <- min
    return(x)
}

col.pal <- c(
	`Agaricomycetes` = "indianred",
	`Sordariomycetes` = "steelblue",
	`Leotiomycetes` = "green",
	`Dothideomycetes` = "darkkhaki",
	`Mock` = "black",
	`Mortierellomycetes` = "darkgrey")

# Read data
setwd("/biodata/dep_psl/grp_psl/Arpan/fungal_inoculation/stage/fungal_inoculation/")
load("./data/taxonomy.RData")
require(tidyverse)

# REad data
df <- read_delim("./data/sfw_data.txt", delim = "\t") %>%
fill(strain, genotype, rep, .direction = "down") %>%
na.omit(.) %>%
filter(strain != "F131") %>%
data.frame(., stringsAsFactors = FALSE)

df_prl <- read_delim("./data/prl_data.txt", delim = "\t") %>%
fill(strain, Plate, genotype, replicate, .direction = "down") %>%
na.omit(.) %>%
filter(strain != "F131") %>%
mutate(Plate = str_replace_all(Plate, "\\.", "-")) %>%
data.frame(., stringsAsFactors = FALSE)


# Boxplot for the snijder chamber samples
target <- "Set3"
id <- str_detect(df$exp_id, target)
df_subset <- df[id,]
mock_levs <- 1000*c(mean(df_subset$shoot_fresh_weight[df_subset$strain == "mock" & df_subset$genotype == "845" & df_subset$rep == "1"]))
df_subset %>%
	mutate(
		strain = str_replace_all(str_replace_all(strain, "mock", "Mock"), "F46a", "F46CA"),
		genotype = str_replace_all(genotype, "1039", "845")
		) %>%
	mutate(
		strain = as.factor(strain),
		genotype = factor(as.character(genotype), 
			levels = c("845", "845f", "992", "994"), 
				labels = c("Col-0", "Col-0 flg22", "pyk10bglu21", "cyp79b2b3")),
		rep = factor(rep, levels = c("1", "2"), labels = c("Cabinet", "Chamber")),
		exp_id = as.factor(exp_id)
			) %>%
	 ggplot(aes(y = 1000*shoot_fresh_weight, x = genotype)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = rep, 
	 		colour = genotype), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = genotype),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ rep, scale = "free_y") +
	 scale_colour_manual(values = c("black", "darkgrey", "darkcyan", "indianred")) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Shoot fresh weight (mg)", colour = "", shape = "") +
	ggsave("./figures/sfw_bxp_snijder.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")

id <- str_detect(df_prl$exp_id, target)
df_subset <- df_prl[id,]
mock_levs <- c(mean(df_subset$prl[df_subset$strain == "mock" & df_subset$genotype == "845" & df_subset$rep == "1"]))
df_subset %>%
	mutate(
		strain = str_replace_all(str_replace_all(strain, "mock", "Mock"), "F46a", "F46CA"),
		genotype = str_replace_all(genotype, "1039", "845")
		) %>%
	mutate(
		strain = as.factor(strain),
		genotype = factor(as.character(genotype), 
			levels = c("845", "845f", "992", "994"), 
				labels = c("Col-0", "Col-0 flg22", "pyk10bglu21", "cyp79b2b3")),
		replicate = factor(replicate, levels = c("1", "2"), labels = c("Cabinet", "Chamber")),
		exp_id = as.factor(exp_id)
			) %>%
	 ggplot(aes(y = prl, x = genotype)) +
	 geom_hline(yintercept = mock_levs, lwd = 1, lty = "solid", alpha = 0.8) +
	 geom_point(position = position_jitterdodge(jitter.width = .02, dodge.width = .2), 
	 	aes(shape = replicate, 
	 		colour = genotype), 
	 	size = 0.6, alpha = 0.6) +
	 geom_boxplot(fill = NA, 
	 	aes(colour = genotype),
	 	outlier.alpha = 0, 
	 	alpha = 0.8) +
	 facet_wrap(strain ~ replicate, scale = "free_y") +
	 scale_colour_manual(values = c("black", "darkgrey", "darkcyan", "indianred")) +
	 scale_shape_manual(values = c(0, 15, 1, 16, 2, 17, 3, 18, 4, 19, 6, 21, 7, 14)) +
	 theme_classic() +
	 theme(legend.position = "top", 
	 	axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5),
	 	strip.text.y = element_text(size = 2, hjust = 0.5, vjust = 0.5)) +
	 labs(x = "", y = "Primary root length (mm)", colour = "", shape = "") +
	ggsave("./figures/prl_bxp_snijder.png", 
		units = "in", 
		width = 12, 
		height = 9, 
		bg = "white")


# plant count data by Jana
plant_count <- read.table("./data/plant_count.tsv", sep="\t", header=T, row.names=NULL, fill=T, stringsAsFactors=F)

# for shoot fresh weight
df <- parallel::mclapply(1:nrow(plant_count), function(x){
	# extract corresponding data
	target_strain   <- plant_count$Strain[x]
	target_rep      <- plant_count$Replicate[x]
	target_genotype <- plant_count$Genotype[x]

	idx <- df$strain == target_strain & df$genotype == target_genotype & df$rep == target_rep

	if(sum(idx) > 0) {
		# when there is at least one value
		temp <- df[idx, c("strain", "genotype", "rep", "shoot_fresh_weight")]

		# count the number of missing values
		missing <- nrow(temp) - plant_count$total[x]

		if(missing > 0){
			# make additional rows with zero values
			add <- data.frame(
				strain   = target_strain,
				genotype = target_genotype,
				rep      = target_rep,
				shoot_fresh_weight = rep(0, missing)
			)
			temp <- rbind(temp, add)
		}

	} else {
		# when there are no data recorded, make data with zeros
		temp <- data.frame(
			strain   = target_strain,
			genotype = target_genotype,
			rep      = target_rep,
			shoot_fresh_weight = rep(0, plant_count$total[x])
		)
	}

	# remove suspicous plates
	if(plant_count$doubts[x] == "X") temp <- NULL

	# return corrected data
	return(temp)

}, mc.cores = 4) %>% do.call(rbind.data.frame, .)


# for primay root length
df_prl <- parallel::mclapply(1:nrow(plant_count), function(x){
	# extract corresponding data
	target_strain   <- plant_count$Strain[x]
	target_rep      <- plant_count$Replicate[x]
	target_genotype <- plant_count$Genotype[x]

	idx <- df_prl$strain == target_strain & df_prl$genotype == target_genotype & df_prl$replicate == target_rep

	if(sum(idx) > 0) {
		# when there is at least one value
		temp <- df_prl[idx, c("strain", "genotype", "replicate", "prl")]

		# count the number of missing values
		missing <- nrow(temp) - plant_count$total[x]

		if(missing > 0){
			# make additional rows with zero values
			add <- data.frame(
				strain    = target_strain,
				genotype  = target_genotype,
				replicate = target_rep,
				prl       = rep(0, missing)
			)
			temp <- rbind(temp, add)
		}

	} else {
		# when there are no data recorded, make data with zeros
		temp <- data.frame(
			strain    = target_strain,
			genotype  = target_genotype,
			replicate = target_rep,
			prl       = rep(0, plant_count$total[x])
		)
	}

	# remove suspicous plates
	if(plant_count$doubts[x] == "X") temp <- NULL

	# return corrected data
	return(temp)

}, mc.cores = 4) %>% do.call(rbind.data.frame, .)


# Data carpentary
df$treatment <- "mock"
df$treatment[str_detect(df$genotype, "f$")] <- "flg22"
df_prl$treatment <- "mock"
df_prl$treatment[str_detect(df_prl$genotype, "f$")] <- "flg22"

# ---> Note that in Set4 the strain 68 will be 66 and it is changed here <--
idx <- str_detect(df$exp_id, "Set4") & str_detect(df$strain, "F68")
df$strain[idx] <- "F66" 

idx <- str_detect(df_prl$exp_id, "Set4") & str_detect(df_prl$strain, "F68")
df_prl$strain[idx] <- "F66" 

df$genotype <- str_replace_all(df$genotype, "f", "")
df$genotype <- str_replace_all(df$genotype, "1039", "845")
df$strain   <- str_replace_all(df$strain, "^M", "m")
df$strain   <- str_replace_all(df$strain, "F46CA", "F46c")
df$strain   <- str_replace_all(df$strain, "F46[A|a]", "F46c")
# df$strain   <- str_replace_all(df$strain, "-", "\\.")
df$strain   <- str_replace_all(df$strain, "Ct-Wt", "Ct")
df$strain   <- str_replace_all(df$strain, "Ci-Wt", "Ci")
df$strain   <- str_replace_all(df$strain, "\\s+", "")

df_prl$genotype <- str_replace_all(df_prl$genotype, "f", "")
df_prl$genotype <- str_replace_all(df_prl$genotype, "1039", "845")
df_prl$genotype <- str_replace_all(df_prl$genotype, "854", "845")
df_prl$genotype <- str_replace_all(df_prl$genotype, "991", "992")
df_prl$strain   <- str_replace_all(df_prl$strain, "^M", "m")
df_prl$strain   <- str_replace_all(df_prl$strain, "F46CA", "F46c")
df_prl$strain   <- str_replace_all(df_prl$strain, "F46[A|a]", "F46c")
# df_prl$strain   <- str_replace_all(df_prl$strain, "-", "\\.")
df_prl$strain   <- str_replace_all(df_prl$strain, "Ct-Wt", "Ct")
df_prl$strain   <- str_replace_all(df_prl$strain, "Ci-Wt", "Ci")
df_prl$strain   <- str_replace_all(df_prl$strain, "\\s+", "")

# Check the strain names in CC
out <- rev(unique(df_prl$strain)[which(!unique(df_prl$strain) %in% taxa_obj$taxa$strain_ids)])
df_prl$cc <- "MPI_hq"
df_prl$cc[df_prl$strain %in% out] <- "other"

df$cc <- "MPI_hq"
df$cc[df$strain %in% out] <- "other"

df$plateid <- paste(df$exp_id, df$strain, df$treatment, df$genotype, df$rep, sep = "_")
# df_prl$plateid <- paste(df_prl$exp_id, df_prl$strain, df_prl$treatment, df_prl$genotype, df_prl$replicate, sep = "_")
df_prl$plateid <- paste(df_prl$strain, df_prl$treatment, df_prl$genotype, df_prl$replicate, sep = "_")

df$group <- factor(paste(df$genotype, df$treatment, sep = "_"), 
	levels = c("845_mock", "845_flg22",  "992_mock",  "994_mock") )
df_prl$group <- factor(paste(df_prl$genotype, df_prl$treatment, sep = "_"), 
	levels = c("845_mock", "845_flg22",  "992_mock",  "994_mock") )

df$random <- as.factor(paste(df$exp_id, df$rep, sep = "_"))
df_prl$random <- as.factor(paste(df_prl$exp_id, df_prl$replicate, sep = "_"))

df$random <- as.factor(df$rep)
df_prl$random <- as.factor(df_prl$replicate)

# set1_strains <- intersect(unique(df$strain[df$exp_id == "Exp2Set1"]),
# 	unique(df$strain[df$exp_id == "Exp1Set12"]))

# df <- df %>% filter(strain %in% set1_strains)

set1_strains <- intersect(unique(df_prl$strain), unique(df$strain))

# Multiariate analysis
df$strain <- factor(df$strain, levels = c("mock", set1_strains[set1_strains != "mock"]))
df_prl$strain <- factor(df_prl$strain, levels = c("mock", set1_strains[set1_strains != "mock"]))
df_prl <- df_prl %>% filter(strain %in% set1_strains)
df <- df %>% filter(strain %in% set1_strains)

# Save in one object
dat <- list(prl = df_prl, sfw = df)

save(list = "dat", file = "./data/data.RData")

# END
sessionInfo()
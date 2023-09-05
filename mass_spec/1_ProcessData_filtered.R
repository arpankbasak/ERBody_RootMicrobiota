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
pkgs <- c("tidyverse", "reshape2", "sva")

lapply(pkgs, require, character.only = T)

# Loading feature data
# ms_tab <- read.csv2(paste0("./data/raw_msms.csv"), header = T, stringsAsFactors = F, row.names = 1)
ms_tab <- read.csv2(paste0("./data/root_exudates_120721_combined_FOR STATISTICS.csv"), header = T, stringsAsFactors = F, row.names = 1)
# colnames(ms_tab)[1] <- "MZ"

sample_table <- data.frame(given_id = colnames(ms_tab), given_geno = as.character(ms_tab[1,]))
sample_table$batch <- "Krakow"
sample_table$batch[str_detect(sample_table$given_id, "A_")] <- "Cologne_A"
sample_table$batch[str_detect(sample_table$given_id, "B_")] <- "Cologne_B"

sample_table$genotype[str_detect(sample_table$given_geno, "col0")] <- "Col"
sample_table$genotype[str_detect(sample_table$given_geno, "cyp|CYP")] <- "cyp"
sample_table$genotype[str_detect(sample_table$given_geno, "pyk10")] <- "pyk10"
sample_table$genotype[str_detect(sample_table$given_geno, "nai1")] <- "nai1"
sample_table$genotype[str_detect(sample_table$given_geno, "myb34")] <- "myb"
sample_table$genotype[str_detect(sample_table$given_geno, "ms")] <- "MS"
sample_table$genotype[str_detect(sample_table$given_geno, "qc")] <- "qc"

sample_table$group <- paste(sample_table$genotype, sample_table$batch, sep = ":")
levs <- lapply(c("Cologne_A", "Cologne_B", "Krakow"), function(x) paste(unique(as.character(sample_table$genotype)), x, sep = ":")) %>% unlist(.)

# Create the peak dataframe
# mz_tab <- data.frame(mz = ms_tab$MZ) %>%
# add_column(pid = row.names(.), .before = 1) %>%
# mutate(mz_id = str_replace_all(mz, "/", ":")) %>%
# separate(mz_id, into = c("Mass", "RT", "mode"), sep = ":", remove = FALSE) %>%
# mutate(Mass = as.numeric(Mass),
#   RT = as.numeric(RT),
#   pid = paste(mode, pid, sep=""))

# The data here is not transformed in any dimension so the ROWIDs are identical
# mat <- apply(ms_tab[,-1], 2, as.numeric)
# row.names(mat) <- mz_tab$pid
# colnames(mat) <- str_replace_all(colnames(mat), "\\.raw\\.Peak\\.area", "")

# Create sample dataframe
meta <- sample_table %>%
mutate(
	genotype = factor(genotype, levels = c("Col", "nai1", "pyk10", "cyp", "myb", "MS", "qc")), 
	group = factor(group, levels = levs[levs %in% as.character(group)]),
  replicate = as.factor(given_id),
  batch = factor(batch, levels = c("Cologne_A", "Cologne_B", "Krakow"))
  ) %>%
filter(genotype %in% c("Col", "cyp", "pyk10", "nai1", "myb") & batch != "Cologne_B")

# Normalize matrices
# mat_norm <- t(apply(mat, 1, function(x) (x-mean(x))))
mat <- apply(ms_tab[-1,], 2, as.numeric)
row.names(mat) <- row.names(ms_tab[-1,])
colSums(mat)

mat <- mat[,meta$given_id]
# idx <- mat < 10e-17
# mat[idx] <- 10e-17 # Check if needed or not
mat_lognorm <- log2(mat+1e-9)
mat_norm <- apply(mat_lognorm, 2, function(x) (x - mean(x))/sd(x))
summary(mat_norm)


# Use glog transformation
# mat_lognorm <- log2(mat + sqrt(mat^2 + lambda[4]))
# mat_lognorm <- t(apply(mat_lognorm, 1, function(x) saturate(x)))



# Make R objects
ms_data <- list(

  metadata = meta,
  mat = mat,
  mat_norm = mat_norm,
  mat_lognorm = mat_lognorm
  # mz_tab = mz_tab

  )

# Save objects of clusters
save(list = "ms_data", 
  file = paste0("./data/msms_data.Rdata"))

# !END OF SCRIPT
sessionInfo()
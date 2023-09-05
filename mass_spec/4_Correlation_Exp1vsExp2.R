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
load(paste0("./data/msms_data.Rdata"))
adat <- read.delim2("./data/strehmel_rootexudate_metabolome_adata.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# Separate negative and positive ionisation
meta <- ms_data$meta 
# %>% filter(genotype %in% c("Col", "pyk10", "cyp"))
mat <- ms_data$mat[, meta$given_id]

geno <- genotype[which(genotype$short %in% meta$genotype),]
levs <- lapply(unique(meta$batch), function(x) {
  paste(geno$short, x, sep = "_")
}) %>% unlist(.)

design <- meta
row.names(design) <- design$given_id
design <- design %>% 
filter(genotype %in% c("Col", "pyk10", "cyp")) %>%
mutate(group = factor(paste(genotype, batch, sep = "_"), levels = levs))

design$genotype <- droplevels(design$genotype)
design$group <- droplevels(design$group)
design$batch <- droplevels(design$batch)

## Correlation among the experiments overall

mat_1 <- mat[,design$given_id[design$batch == "Cologne_A"]]
mat_2 <- mat[,design$given_id[design$batch == "Krakow"]]
df <- data.frame(id = row.names(mat)) %>% mutate(
  Cologne_A = log10(mat_1[match(.$id, row.names(mat_1))] + 1e-5), 
  Krakow = log10(mat_2[match(.$id, row.names(mat_2))] + 1e-5)
  )
ct <- cor.test(df$Cologne_A, df$Krakow, method = "pearson")

# Subsetting Col
mat_1 <- mat[,design$given_id[design$group == "Col_Cologne_A"]]
mat_2 <- mat[,design$given_id[design$group == "Col_Krakow"]]
df_col <- data.frame(id = row.names(mat)) %>% mutate(
  Cologne_A = log10(mat_1[match(.$id, row.names(mat_1))] + 1e-5), 
  Krakow = log10(mat_2[match(.$id, row.names(mat_2))] + 1e-5)
  )
# ct <- cor.test(df_col$Cologne_A, df_col$Krakow, method = "pearson") # .485

# Subsetting pyk10
mat_1 <- mat[,design$given_id[design$group == "pyk10_Cologne_A"]]
mat_2 <- mat[,design$given_id[design$group == "pyk10_Krakow"]]
df_pyk10 <- data.frame(id = row.names(mat)) %>% mutate(
  Cologne_A = log10(mat_1[match(.$id, row.names(mat_1))] + 1e-5), 
  Krakow = log10(mat_2[match(.$id, row.names(mat_2))] + 1e-5)
  )
# ct <- cor.test(df_pyk10$Cologne_A, df_pyk10$Krakow, method = "pearson") # .0.49

# Subsetting cyp
mat_1 <- mat[,design$given_id[design$group == "cyp_Cologne_A"]]
mat_2 <- mat[,design$given_id[design$group == "cyp_Krakow"]]
df_cyp<- data.frame(id = row.names(mat)) %>% mutate(
  Cologne_A = log10(mat_1[match(.$id, row.names(mat_1))] + 1e-5), 
  Krakow = log10(mat_2[match(.$id, row.names(mat_2))] + 1e-5)
  )
ct <- cor.test(df_cyp$Cologne_A, df_cyp$Krakow, method = "pearson") # 0.79

# Projecting each genotype one by one in one canvas
p <- df_col %>%
ggplot(aes(x = Cologne_A, y = Krakow)) +
ggtitle(paste0("r.sq = ", round(ct$estimate, 4), "; p.value = ", ifelse(ct$p.value <= 0.001, "<=.0001", ct$p.value))) +
geom_point(fill = genotype$colours[genotype$short == "Col"], colour = "black", size = 1, shape = 21, alpha = 0.4) +
geom_point(data = df_pyk10, fill = genotype$colours[genotype$short == "pyk10"], colour = "black", size = 1, shape = 21, alpha = 0.4) +
geom_point(data = df_cyp, fill = genotype$colours[genotype$short == "cyp"], colour = "black", size = 1, shape = 21, alpha = 0.4) +
scale_fill_manual(values = genotype$colours[which(genotype$short %in% df$genotype)]) +
geom_abline(colour = "darkgrey") +
theme_RTN_MDS +
theme(text = element_text(size = 8)) +
labs(x = "log10(Exp1_Ra)", y = "log10(Exp2_Ra)")

ggsave(file = paste0(figs, "cor_experiment_geno_coloured.png"),
units = "in",
device = "png",
width = 5,
height = 4.5
)

# Projecting all in one
p <- df %>%
ggplot(aes(x = Cologne_A, y = Krakow)) +
ggtitle(paste0("r.sq = ", round(ct$estimate, 4), "; p.value = ", ifelse(ct$p.value <= 0.001, "<=.0001", ct$p.value))) +
geom_point(data = df[,], fill = genotype$colours[genotype$short == "Col"], colour = "black", size = 1, shape = 21, alpha = 0.8) +
geom_abline(colour = "darkred") +
theme_RTN_MDS +
theme(text = element_text(size = 8)) +
labs(x = "Exp1", y = "Exp2")

ggsave(file = paste0(figs, "cor_experiment.png"),
units = "in",
device = "png",
width = 5,
height = 4.5
)

# END OF SCRIPT
sessionInfo()
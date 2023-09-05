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
pkgs <- c("tidyverse", "reshape2", "caret", "vegan", "mda")

lapply(pkgs, require, character.only = T)
dat <- read.delim2("./data/from_AP/splsda_score.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Loading feature data
load(paste0("./data/msms_data.Rdata"))

# Separate negative and positive ionisation
# mat_p <- ms_data$mat_norm[str_detect(row.names(ms_data$mat_norm), "^p"),]
# mat_n <- ms_data$mat_norm[str_detect(row.names(ms_data$mat_norm), "^n"),]
meta <- ms_data$meta %>% select(given_id, genotype, batch, group)
row.names(meta) <- meta$given_id

splsda <- meta %>% cbind.data.frame(., dat[match(str_replace_all(.$given_id, "X", ""), dat$X), c(2:4)])

temp <- splsda %>% group_by(group) %>%
summarise(m1 = mean(as.numeric(comp.1)), m2 = mean(as.numeric(comp.2)), m3 = mean(as.numeric(comp.3))) %>%
data.frame(., stringsAsFactors = FALSE)

geno <- genotype[which(genotype$short %in% meta$genotype),]

sp12 <- splsda %>%
cbind.data.frame(., temp[match(as.character(.$group), as.character(temp$group)), c("m1", "m2")]) %>%
ggplot(aes(x = as.numeric(comp.1), y= as.numeric(comp.2))) +
# ggtitle(ti) +
geom_point(aes(fill = genotype, shape = batch), colour = c_black, size = 3) +
# geom_hline(yintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
# geom_vline(xintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
geom_segment(aes(x = as.numeric(m1), y = as.numeric(m2),
          xend = as.numeric(comp.1), yend = as.numeric(comp.2), 
          colour = genotype), 
        lty = "solid",
        lwd = 0.3) +
scale_fill_manual(values = c(geno$colours)) +
scale_shape_manual(values = c(21, 23, 24)) +
scale_colour_manual(values = c(geno$colours)) +
theme_RTN_MDS +
theme(plot.title = element_text(size = 6)) +
labs(x = paste0("Component1"), 
  y = paste0("Component2"), colour = "", fill = "", shape = "")

ggsave(sp12, file = paste0(figs, "splsda_draw12.png"),
  width = cpcoa.box[1],
  height = cpcoa.box[2]
  )


sp23 <- splsda %>%
cbind.data.frame(., temp[match(as.character(.$group), as.character(temp$group)), c("m3", "m2")]) %>%
ggplot(aes(x = as.numeric(comp.3), y= as.numeric(comp.2))) +
# ggtitle(ti) +
geom_point(aes(fill = genotype, shape = batch), colour = c_black, size = 3) +
# geom_hline(yintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
# geom_vline(xintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
geom_segment(aes(x = as.numeric(m3), y = as.numeric(m2),
          xend = as.numeric(comp.3), yend = as.numeric(comp.2), 
          colour = genotype), 
        lty = "solid",
        lwd = 0.3) +
scale_fill_manual(values = c(geno$colours)) +
scale_shape_manual(values = c(21, 23, 24)) +
scale_colour_manual(values = c(geno$colours)) +
theme_RTN_MDS +
theme(plot.title = element_text(size = 6)) +
labs(x = paste0("Component3"), 
  y = paste0("Component2"), colour = "", fill = "", shape = "")

ggsave(sp23, file = paste0(figs, "splsda_draw23.png"),
  width = cpcoa.box[1],
  height = cpcoa.box[2]
  )




# geno <- genotype[which(genotype$short %in% meta$genotype),]

# preproc.param <- mat.cut %>% 
#   preProcess(method = c("center", "scale"))


# mat.norm <- preproc.param %>% predict(mat)
# mat.norm <- as.data.frame(t(mat.norm)) %>% add_column(genotype = meta$genotype[match(row.names(.), row.names(meta))])

# # MDA analysis
# mda_mod <- mda(genotype ~., data = mat.norm)


# # RDA analysis
# rda_mod <- rda(genotype ~., data = mat.norm)















save(list = "ranks", file = paste("./data/top_loadings.Rdata", sep = ""))


# !END OF SCRIPT
sessionInfo()
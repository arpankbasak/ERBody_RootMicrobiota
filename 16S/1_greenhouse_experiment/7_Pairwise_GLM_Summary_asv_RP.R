#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for GLM based clustering of ASVs
# Removing ASVs of importance and obtaining the most contributing taxa
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "RColorBrewer", "parallel", "cowplot")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))
load(paste0("./data/CAS_DE_asv.RData"))

# REad the taxonomy table from the ATSPhereCC
# taxa_at_cc <- read.table("./output/taxonomy_syncom.txt", sep = "\t", header = TRUE, as.is = TRUE)

# What are we anlysing
pro <- c("gh_GG")

# pro <- names(meta$datasets)
figs.out <- paste0(figs, paste0("/gh/differential_analysis/uspiked/asv/", pro))
stats.out <- paste0(stats, paste0("/gh/differential_analysis/uspiked/asv/", pro))
out.out <- paste0(out, paste0("/gh/differential_analysis/uspiked/asv/", pro))

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

# Plot Canvas
set.seed(seeder)

temp.meta <- rbind.data.frame(meta$datasets[["CAS"]]) %>%
filter(Compartment != "rhizosphere", !SampleID %in% outlier)
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$Replicate))
temp.meta$Experiment <- factor(as.character(temp.meta$Experiment))
temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# Use short terms for factoring
temp.meta$Genotype <- genotype$short[match(temp.meta$Genotype, genotype$names)]
temp.meta$Compartment <- compartment$short[match(temp.meta$Compartment, compartment$names)]
temp.meta$group <- paste(temp.meta$Genotype, temp.meta$Compartment, temp.meta$Soil_Batch, sep = "_")

# LFC matrix
lfc <- CAS_DE_asv$logFC_P[,!str_detect(colnames(CAS_DE_asv$logFC_P), "_E_")]

# Subset of ASVs that are prevalant
temp.meta_E <- temp.meta %>% filter(Compartment == "RP")
id_E <- as.character(temp.meta_E$SampleID)
idx <- temp.meta_E$Soil_Batch == "CAS11"
IDs_11 <- temp.meta_E$SampleID[idx]


# Relative abundance table
asv.ra <- apply(temp.mat[,id_E], 2, function(x) x/sum(x))
idx <- rowSums(asv.ra * 100 > threshold) >= 1
asv.ra <- apply(temp.mat[idx,id_E], 2, function(x) x/sum(x))
idx <- colnames(asv.ra) %in% IDs_11
ra_E_cut_11 <- asv.ra[, idx]
ra_E_cut_13 <- asv.ra[, !idx]
​
# Compute prevelance
idx <- (rowSums(ra_E_cut_11 > 0) > ncol(ra_E_cut_11)/2) & (rowSums(ra_E_cut_13 > 0) > ncol(ra_E_cut_13)/2)
asv.prev <- asv.ra[idx,]
​idx <- (row.names(lfc) %in% row.names(asv.prev))
lfc.prev <- lfc[idx, ]

ra.prev <- as.data.frame(asv.prev) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$Family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = temp.meta$group[match(.$SampleID, temp.meta$SampleID)]) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(ra.prev)[-1] <- paste0("meanra_", colnames(ra.prev[,-1]), "_RA")
row.names(ra.prev) <- ra.prev$tag


ra <- as.data.frame(asv.ra) %>% 
mutate(tag = replace_na(all_exp_cut$taxa$Family[match(row.names(.), all_exp_cut$taxa$asvb)], "Unclassified")) %>%
# add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = temp.meta$group[match(.$SampleID, temp.meta$SampleID)]) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

colnames(ra)[-1] <- paste0("meanra_", colnames(ra[,-1]), "_RA")
row.names(ra) <- ra$tag

# Top 20 families
top_20 <- data.frame(mean_ra = rowSums(ra[,-1])) %>%
mutate(tag = row.names(ra)) %>%
filter(tag != "Unclassified") %>%
arrange(desc(mean_ra)) %>%
.$tag %>% head(., n=20)

# Barplot for rank abundance
ra_table <- ra %>%
gather(key = "key", value = "vals", convert = FALSE, -tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "x2")) %>%
select(-x, -x1,-x2) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)])
	) %>%
filter(tag %in% top_20) %>%
mutate(tag = factor(tag, levels = top_20)) %>%
group_by(tag, Genotype) %>%
summarise(meanRA = mean(vals)) %>%
data.frame(., stringsAsFactors = FALSE)

p <- ra_table %>%
ggplot(aes(x = tag, y = 100*meanRA)) +
geom_bar(stat = "identity", position = "dodge", aes(fill = Genotype), colour = c_dark_grey) +
scale_fill_manual(values = genotype$colours[which(genotype$short %in% ra_table$Genotype)]) +
theme_RTN +
theme(
	panel.spacing = unit(0.5, "lines"), 
	axis.text.x = element_text(size = 6, angle = 90, hjust = 0.9, vjust = 0.5),
	axis.text.y = element_text(size = 10), 
	axis.title = element_text(size = 10)) +
labs(x = "Top 20 abundant families", y = "mean RA (%)")
ggsave(p, file = paste0(figs.out, "/rank_abundance_top20_family_E.png"), 
            dpi = 600, 
            units = img, 
            device = "png", 
            bg = "transparent", 
            width = 7, 
            height = 4, 
            limitsize = T)


# Heatmap to show the prevelant ASVs
ra.prev.hmap <- as.data.frame(asv.prev) %>% 
add_column(tag = row.names(.), .before = 1) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -tag) %>%
mutate(group = paste("mean", temp.meta$group[match(.$SampleID, temp.meta$SampleID)], "ra", sep = "_")) %>%
group_by(tag, group) %>%
summarise(cRA = mean(RA)) %>%
spread(key = group, value = cRA, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE, row.names = .$tag)

df_hmap <- lfc.prev %>% 
cbind.data.frame(., ra.prev.hmap[match(row.names(.), ra.prev.hmap$tag), ]) %>%
gather(key = "key", value = "vals",-tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "value")) %>%
select(-x, -x1) %>%
# group_by(tag, Genotype, Soil_Batch, value) %>%
spread(key = value, value = vals, fill = NA) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)]),
	sig = PValue <= 0.05,
	higher_tag = replace_na(all_exp_cut$taxa$Family[match(.$tag, all_exp_cut$taxa$asvb)], "Unclassified")
) %>%
filter(higher_tag != "Unclassified") %>%
data.frame(., stringsAsFactors = FALSE)

p <- df_hmap %>%
filter(Genotype != "Col") %>%
mutate(
  behaviour = as.factor(ifelse(!is.na(logFC), sign(logFC), -1)),
  ra = log10(ra),
  logFC = ifelse(abs(logFC) > 4, sign(logFC)*4, logFC)) %>%
ggplot(aes(x = Genotype, y = tag)) +
geom_point(alpha=0.8, aes(shape = behaviour, size = ra, colour = sig, fill = logFC), na.rm = FALSE) +
# geom_tile(alpha=1, aes(colour = sig), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(higher_tag ~ Soil_Batch, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_cudo_magenta, 
        mid = c_white,
        low = c_dark_green, 
        na.value = c_black,
        midpoint = 0, 
        breaks = c(-4, -2, 0, 2, 4), 
        limits = c(-4,4), 
        labels = paste0(c("<4", "-2", "0", "2", ">4"))) +
scale_shape_manual(values = c(`1` = 25, `-1` = 24)) +
scale_size_continuous(range = c(0.5,2)) +
scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_dark_grey)) +
theme_RTN_MDS +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      axis.text.x=element_text(angle = 90, size = 6, face = "bold", hjust = 1, vjust = 0.5), 
      strip.text.x = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      axis.text.y=element_text(size = 2),
      strip.text.y = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "")
ggsave(p, filename = paste0(figs.out,"/prevelant_heatmap_asvs.png"),
       bg="transparent", 
       width=4, 
       height= 20,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)


# LFC correlation
lfc_table <- lfc %>%
add_column(tag = row.names(.), .before = 1) %>%
gather(key = "key", value = "vals", convert = FALSE, -tag) %>%
separate(key, into = c("x", "Genotype", "x1", "Soil_Batch", "value")) %>%
select(-x, -x1) %>%
spread(key = value, value = vals, fill = NA) %>%
select(-PValue) %>%
mutate(
	Genotype = factor(Genotype, levels = genotype$short[which(genotype$short %in% .$Genotype)]),
	Soil_Batch = factor(Soil_Batch, levels = soilbatch$names[which(soilbatch$names %in% .$Soil_Batch)])
	) %>%
spread(key = Genotype, value = logFC, fill = NA) %>%
data.frame(., stringsAsFactors = FALSE)

# Correlagram
d <- (cor(lfc[,str_detect(colnames(lfc), "_logFC$") & str_detect(colnames(lfc), "^genotype")]))
d[lower.tri(d)] <- NA

cc_df <- as.data.frame(d) %>%
add_column(x = row.names(.)) %>%
gather(key = "y", value = "pcc", -x) %>%
separate(x, into = c("v1", "genotypeA", "v2", "Soil_BatchA", "v3"), sep = "_") %>%
separate(y, into = c("v4", "genotypeB", "v5", "Soil_BatchB", "v6"), sep = "_") %>%
select(genotypeA, Soil_BatchA, genotypeB, Soil_BatchB, pcc) %>%
mutate(
    genotypeA = factor(genotypeA, levels = genotype$short[which(genotypeA %in% genotype$short)]),
    genotypeB = factor(genotypeB, levels = genotype$short[which(genotypeB %in% genotype$short)]),
    Soil_BatchA = factor(Soil_BatchA, levels = c("CAS11", "CAS13")),
    Soil_BatchB = factor(Soil_BatchB, levels = c("CAS11", "CAS13"))
    ) %>%
data.frame(.)

# Plot the correlogram
p <- cc_df %>%
filter(Soil_BatchA == Soil_BatchB) %>%
ggplot(aes(x = genotypeA, y=genotypeB, fill = pcc)) +
geom_point(alpha=0.8, shape = 21, size = 10, colour = c_white, na.rm = TRUE) +
geom_text(aes(label = as.character(round(pcc, 2))), size = 3, colour = c_black, guide = FALSE) +
# geom_tile(alpha=1, aes(colour = sig), fill = NA, width = 1, height = 1, size = 0.5) +
facet_grid(. ~ Soil_BatchA, switch = "both", scales="free", space="free") +
scale_fill_gradient2(
        high=c_green, 
        mid = c_white,
        low = c_cudo_magenta, 
        na.value = NA,
        midpoint = 0, 
        breaks = c(-1, -0.5, 0, 0.5, 1), 
        limits = c(-1,1), 
        labels = paste0(c("<1", "-0.5", "0", "0.5", ">1"))
        ) +
theme_RTN +
theme(panel.spacing = unit(0.1, "lines"), 
      panel.border = element_rect(fill="transparent", colour=NA),
      # strip.text.x = element_text(size = 6, hjust = 0.5, vjust = 0.5),
      # axis.text.y=element_text(size = 2),
      # strip.text.y = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = 6, hjust = 0.5, vjust = 0.5)) +
labs(y="", x="", fill="", shape = "", size = "", colour = "", label = "")
ggsave(p, filename = paste0(figs.out,"/correlogram_heatmap_asvs_RP.png"),
       bg="transparent", 
       width=3.5, 
       height= 3,
       units = img, 
       limitsize=F, 
       device = "png", 
       dpi = 600)

# Pairwise Plot
for(i in c("nai1", "pyk10", "cyp", "myb")){

   for(j in c("nai1", "pyk10", "cyp", "myb")){

      temp <- lfc_table
      if(j != i && dim(temp)[1] > 2){

        # i <- "nai1"
        # j <- "pyk10"
      	message(paste0("combination: ", i, " and ", j))


      	temp$x <- temp[,i]
      	temp$y <- temp[,j]

      	stat <- cor.test(temp$x, temp$y) 
        til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", (stat$p.value))

        p <- temp %>%
        select(tag, Soil_Batch, x, y) %>%
        mutate(prevelant = tag %in% row.names(lfc.prev)) %>%
        ggplot(aes(x = x, y = y)) +
        ggtitle(paste0(paste(i,j, sep = "_"), til)) +
        geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
        geom_point(aes(
              fill = Soil_Batch, colour = prevelant), 
              shape = 21, 
              size = 2,
              na.rm = FALSE) +
        scale_colour_manual(values = c(`TRUE` = c_black, `FALSE` = c_grey)) +
        scale_fill_manual(values = as.character(soilbatch$colours[which(soilbatch$names %in% temp$Soil_Batch)])) +
        theme_RTN_MDS +
        theme(panel.spacing = unit(0.5, "lines"), plot.title = element_text(size=4)) +
        labs(x = i, y = j, fill = "")
        ggsave(p, file = paste0(figs.out, "/Parwise_RP_", j,"_", i,"_logFC.png"), 
                    dpi = 600, 
                    units = img, 
                    device = "png", 
                    bg = "transparent", 
                    width = 4, 
                    height = 4, 
                    limitsize = T)

		}

	}
}

# END OF SCRIPT
sessionInfo()

#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript
# Script for Exploratory data analysis
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
path <- "/netscratch/dep_psl/grp_psl/Arpan/analysis/"
source(paste0(path, "/manifest/parameters.R"))
source(paste0(path, "/manifest/functions.R"))
setwd(analysis_combat.16s)

# Loading required packages
pkgs <- c("tidyverse", "vegan", "RColorBrewer", "parallel")
lapply(pkgs, require, character.only = T)
load(paste0("./data/metadata_asvs.RData"))
load(paste0("./data/all_experiments_filtered_asvs.RData"))

# What are we anlysing
p <- "st"

figs.out <- paste0(figs, "/gh/taxonomy_composition/uspiked/")
stats.out <- paste0(stats, "/gh/taxonomy_composition/uspiked/")
out.out <- paste0(out, "/gh/taxonomy_composition/uspiked/")

mclapply(c(figs.out,stats.out,out.out), function(x){

    if(!dir.exists(paths = x)){

        message(paste0("Directory created ", x))
      dir.create(x, recursive = TRUE)

    } else{
      
      message(paste0("Directory exists ", x))

    }

}, mc.cores = 4)

# Sub path
plot_path <- paste0(figs, "/taxonomy_composition/uspiked/")
stat_path <- paste0(stats, "/taxonomy_composition/uspiked/")
out_path <- paste0(out, "/taxonomy_composition/uspiked/")

lapply(c(plot_path,stat_path,out_path), function(x){

    exp_path <- paste0(x, p)
    if(!dir.exists(paths = exp_path)){

        message(paste0("Directory created ", exp_path))
        dir.create(exp_path, recursive = FALSE)

    } else{
      
      message(paste0("Directory exists ", exp_path))

    }

})

plot_path <- paste0(plot_path, p)
stat_path <- paste0(stat_path, p)
out_path <- paste0(out_path, p)


# Prepare the ASV and metadata
temp.meta <- (meta$datasets[["STREX"]])
temp.meta <- data.frame(meta$design_strex[match(temp.meta$SampleID, meta$design_strex$SampleID),])
idx <- temp.meta$SampleID
temp.mat <- all_exp_cut$asv
idx <- which(colnames(temp.mat) %in% idx)
temp.mat <- temp.mat[,idx]
temp.meta <- temp.meta[which(temp.meta$SampleID %in% colnames(temp.mat)),]
# temp.meta$Soil_Batch[temp.meta$Soil_Batch == ""] <- temp.meta$soil[temp.meta$Soil_Batch == ""]
temp.meta$Replicate <- as.factor(as.character(temp.meta$replicate))
temp.meta$tech_rep <- as.factor(as.character(temp.meta$tech_rep))
temp.meta$Experiment <- factor(as.character(temp.meta$experiment))
temp.meta$Treatment <- factor(as.character(temp.meta$Treatment), levels = c("Bulk soil", "Exudate", "Extract"))
# temp.meta$Run <- as.factor(as.character(temp.meta$Run))

# temp.mat <- temp.mat[common_in_CAS,]
idg <- which(genotype$short %in% as.character(temp.meta$genotype))
# idsb <- which(soilbatch$name %in% as.character(temp.meta$Soil_Batch))
# idc <- which(compartment$name %in% as.character(temp.meta$Compartment))

temp.meta$Genotype <- factor(as.character(temp.meta$genotype), levels = genotype$short[idg])

mclapply(c("Bulk soil", "Exudate", "Extract"), function(x){


		temp.metadata <- temp.meta %>%
        filter(Treatment == x)
        
        idg <- which(genotype$short %in% as.character(temp.metadata$genotype))
        # idsb <- which(soilbatch$name %in% as.character(temp.metadata$Soil_Batch))
        temp.metadata$Genotype <- factor(as.character(temp.metadata$genotype), levels = genotype$short[idg])
        

	    idx <- which(colnames(temp.mat) %in% temp.metadata$SampleID)
	    temp.mat.ra <- temp.mat[,idx]
	    temp.mat.ra <- apply(temp.mat.ra, 2, function(x) x/sum(x))
	    taxa.mat <- all_exp_cut$taxa
	    taxa_level <- "Family"
	    n <- 10

	    taxa.mat <- taxa.mat[which(taxa.mat$asvb %in% row.names(temp.mat.ra)),]
        taxa <- taxa.mat
        temp <- temp.metadata
        feature <- temp.mat.ra
        
        taxa$tag <- taxa[, taxa_level]
        row.names(taxa) <- taxa[,grepl(colnames(taxa), pattern = "^asv", ignore.case = T)]
        
        # Build RA table
        ra_table <- cbind.data.frame(feature, 
        tag = taxa[row.names(feature), "tag"]) %>% 
            gather(key = "SampleID", value = "ra", convert = F, -tag) %>% 
            group_by(tag, SampleID) %>% 
            summarise(RA = sum(ra)) %>% 
            data.frame(.,stringsAsFactors = F)
        
        # Extract all the tags and sort
        all_tags <- cbind.data.frame(feature, 
        tag = taxa[row.names(feature), "tag"]) %>% 
            gather(key = "SampleID", value = "ra", convert = F, -tag) %>% 
            group_by(tag) %>% 
            summarise(RA = sum(ra)) %>% 
            arrange(desc(RA)) %>% 
            select(tag) %>% 
            unlist(.,use.names = F)

        all_tags <- str_replace_na(all_tags, "Undetermined")
        
        # Set the top aggragated tags
        top_tags <- head(all_tags, n = n+1)
        
        ra_table$tag <- str_replace_na(ra_table$tag, "Undetermined")
        ra_table$tag[!ra_table$tag %in% top_tags] <- "rare_taxa"
        idx <- match(ra_table$SampleID, temp$SampleID)
        
        # Join metadata
        ra_table <- cbind.data.frame(ra_table, temp[idx, 
                                         c("Experiment", "tech_rep", "Genotype", "Replicate")])
        
        
        # Set plotting parametes
        # Experiment <- Experiment[which(Experiment$names %in% temp.metadata$experiment),]
        #soiltype <- soiltype[which(soiltype$names %in% meta$soil),]
        # genotype <- genotype[which(genotype$name %in% temp.metadata$genotype),]
        
        # Set colour pallete --- make a contrasting colour palette for n 
        col.pal <- paint_taxa(taxa_label = c("Undetermined", "rare_taxa", top_tags[!top_tags == "Undetermined"]),
            others = c("Undetermined", "rare_taxa"), 
            palette = "Set1", 
            colour_n = 5)
        
        ra_table$tag <- factor(ra_table$tag, levels = c("Undetermined", "rare_taxa", top_tags[!top_tags == "Undetermined"]))

        # levels(ra_table$Compartment) <- compartment$short[match(as.character(levels(ra_table$Compartment)), compartment$name)]

        # Grouping factors
        
        # Plot Canvas
        plot <- ra_table %>% mutate(SampleID = as.factor(SampleID)) %>%
            ggplot2::ggplot(aes(x= SampleID, y = RA, fill = tag)) +
            ggtitle(paste("Top ",n," - ", taxa_level, sep = "")) +
            geom_bar(stat = "identity", size=0) +
            facet_grid(. ~ Experiment + Genotype + tech_rep, 
                switch = "x", scales = "free", space = "free") +
            scale_fill_manual(values = col.pal) +
            labs(y = "Aggregated Relative abundance", x = "", fill = "") + 
            theme_RTN +
            theme(panel.spacing = unit(0.01, "lines"), 
                axis.text.x = element_blank(),
                  axis.text.y = element_text(size = 10),
                  axis.title = element_text(size = 8, hjust = 0.5, vjust = 0.5),
                  strip.text.x = element_text(angle = 90, size = 8, vjust = .5, hjust = .5),
                  legend.text = element_text(size = 8),
                  legend.title = element_blank()) +
            ggsave(file = paste0(plot_path, "/RA_", taxa_level,"_ST", x,".png"), 
                dpi = 600, 
                device = "png", 
                bg = "transparent", 
                width = 8, 
                units = img,
                height = 4, 
                limitsize = F)


}, mc.cores = 4)

# END
sessionInfo()
#!/biodata/dep_psl/grp_psl/Arpan/anaconda3/envs/gsl_analysis/bin/Rscript

rm(list = ls())

# Script by
# @Arpan Kumar Basak

# Read data
setwd("/biodata/dep_psl/grp_psl/Arpan/fungal_inoculation/stage/fungal_inoculation/")
require(tidyverse)
require(ggtree)

# REad data
tree <- ape::read.tree("./data/SpeciesTree_rooted_node_labels.newick")
taxa <- read.table("./data/taxonomy.txt", sep = "\t", header = TRUE, as.is = TRUE) %>%
mutate(taxa_id = paste0(str_replace_all(species, " ", "_"), "_MPI")) %>%
# select(phylum, class, order, family, genus, species, taxa_id) %>%
# group_by(phylum, class, order, family, genus, species, taxa_id) %>%
# summarise(n = n()) %>%
data.frame(., stringsAsFactors = FALSE)

# Prepare the metadata
metadata <- cbind.data.frame(
	ID = tree$tip.label, taxa[match(as.character(tree$tip.label), taxa$jgi_id), ]) %>%
separate(id, into = c("strain", "source", "host", "strain_id"), remove = FALSE, sep = "-") %>%
mutate(strain_ids = str_replace_all(strain_id, "^0+", "F")) %>%
data.frame(., stringsAsFactors = FALSE)

# Group the Class corresponding to the strain
strain_sorted_phy <- metadata$strain_ids
grp <- list(
	Agaricomycetes = as.character(metadata$ID[str_detect(metadata$class, "Agaricomycetes")]),
	Dothideomycetes = as.character(metadata$ID[str_detect(metadata$class, "Dothideomycetes")]),
	Leotiomycetes = as.character(metadata$ID[str_detect(metadata$class, "Leotiomycetes")]),
	Sordariomycetes = as.character(metadata$ID[str_detect(metadata$class, "Sordariomycetes")]))
	# Mortierellomycetes = as.character(metadata$ID[str_detect(metadata$class, "Mortierellomycetes")]))

tree_grouped <- groupOTU(tree, .node = grp)
# tree_grouped$node.label <- as.character(metadata$class)

# Plot the phylogenetic tree
(taxa_tree <- 
ggtree(tree_grouped, size = 0.7, aes(colour = group)) +
geom_tiplab(align = TRUE, linetype = "dotted", colour = "black", size = 4) +
xlim_tree(1.5) +
# geom_nodelab(aes(label=group)) +
scale_colour_manual(values = c(
	"black",
	`Agaricomycetes` = "indianred",
	`Dothideomycetes` = "darkkhaki",
	`Leotiomycetes` = "green",
	`Sordariomycetes` = "steelblue"
	# `Mortierellomycetes` = "deepskyblue",
	)) +
theme(legend.position = "top", legend.text = element_text(size = 4)) +
labs(colour = "")) +
ggsave("./figures/tree41fungi_CtCi.png", 
	units = "in", 
	# scale = 0.3,
	width = 3, 
	height = 7, 
	bg = "white",
	device = "png", limitsize = FALSE)

taxa_obj <- list(
	taxa = metadata,
	strain_phylogeny = strain_sorted_phy,
	strain_tree = tree,
	tree_graph = taxa_tree)

# Save objects
save(list = "taxa_obj", file = "./data/taxonomy.RData")


# END
sessionInfo()
#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscript
# Ryohei Thomas Nakano

#
rm(list=ls())

#
library(dplyr,    quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(ggplot2,  quietly=T, warn.conflicts=F)


#
data    <- "/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/data/"
scripts <- "/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/"
fig     <- "/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/figures/RTN_growth/"
source(paste(scripts, "parameters.R", sep=""))

# load
count  <- read.table(paste(data, "merged_strain_table.txt", sep=""), header=T, row.names=1,    sep="\t", comment.char="", stringsAsFactors=F) %>% as.matrix
design <- read.table(paste(data, "metadata_table.txt", sep=""),      header=T, row.names=NULL, sep="\t", comment.char="", stringsAsFactors=F)


genotype.syncom <- data.frame(
    names=c("Col-0", "pyk10bglu21", "cyp79b2b3", "inoculum"),
    short=c("Col", "pyk", "cyp", "inoculum"),
    shapes=c(8, 15, 16, 5),
    colours=c(c_black, c_cudo_skyblue, c_red, c_dark_grey),
    stringsAsFactors=F)

# group
design$group <- paste(design$tech_rep, design$biol_rep, design$dose, sep="_")
design$geno_exu_group <- paste(design$genotype, design$exudates_rep, design$group, sep="_")


# remove samples whose spike reads are less than 10
spike_pos <- rownames(count) == "DH5alpha"
idx <- count[spike_pos,] < 10
count <- count[, !idx]

aa <- apply(count, 2, function(x) x[!spike_pos]/x[spike_pos])

# plot overall growth
# total bacterial load
aaSum <- colSums(aa)

# merge meta data
idx <- match(names(aaSum), paste("Sample_", design$SampleID, sep=""))
aaSum <- data.frame(aa=aaSum, design[idx,], stringsAsFactors=F)

aaSum$genotype <- factor(aaSum$genotype, levels=genotype.syncom$short)

# direct plot
p <- ggplot(aaSum, aes(x=factor(time), y=log2(aa), colour=genotype)) +
	geom_point(position=position_jitterdodge(jitter.width=.25), alpha=.5, size=.5) +
	geom_boxplot(fill=c_white, outlier.shape=NA) +
	scale_colour_manual(values=genotype.syncom$colours) +
	facet_grid( ~ dose + biol_rep, labeller=label_both) +
	theme_RTN
ggsave(p, file=paste(fig, "overall_growth.box.pdf", sep=""), width=10, height=4, bg="transparent")





# give inoculum values to each genotype
rel_aa <- lapply(genotype.syncom$short[1:3], function(x){
	lapply(dosage.syncom$names, function(y){
		lapply(1:3, function(i){
			lapply(1:4, function(j){ # for each biological replicate
				lapply(1:8, function(k){ # for each tech_replicate

					idx <- which(aaSum$genotype == x          & aaSum$dose == y & aaSum$exudates_rep == i & aaSum$biol_rep ==j  & aaSum$tech_rep == k)
					
					if(length(idx)){
						temp <- aaSum[idx,]

						idx <- which(aaSum$genotype == "inoculum" & aaSum$dose == y &                           aaSum$biol_rep ==j  & aaSum$tech_rep == k)
						if(length(idx)==0){
							aa_0 <- NA
						} else {
							aa_0 <- aaSum$aa[idx]
						}

						temp_0 <- temp[1,]
						temp_0$time <- 0
						temp_0$aa <- aa_0

						temp <- rbind(temp, temp_0)
						temp$inoc_aa <- aa_0

						return(temp)

					} else {

						return(NULL)

					}

				}) %>% do.call(rbind, .)
			}) %>% do.call(rbind, .)
		}) %>% do.call(rbind, .)
	}) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


# relative to time 0
rel_aa$rel_aa <- log2(rel_aa$aa / rel_aa$inoc_aa)
rel_aa$geno_exu_group <- paste(rel_aa$genotype, rel_aa$exudates_rep, rel_aa$group, sep="_")

# remove biol_rep 4
# idx <- rel_aa$biol_rep != 4
# rel_aa <- rel_aa[idx,]

# sort
rel_aa$genotype <- factor(rel_aa$genotype, levels=genotype.syncom$short)


# plot
p <- ggplot(rel_aa, aes(x=time, y=rel_aa, colour=genotype)) +
	geom_hline(yintercept=0, size=.2, colour=c_black, linetype="dashed") +
	
	geom_line(aes(group=geno_exu_group), alpha=.3, size=.2) +
	scale_colour_manual(values=genotype.syncom$colours) +
	facet_grid(biol_rep ~ dose) +
	theme_RTN
ggsave(p, file=paste(fig, "overall_growth.relatve_to_0hpi.pdf", sep=""), width=5.5, height=4, bg="transparent")


# Grouped Plot
group_df <- rel_aa %>%
filter(!is.na(inoc_aa) & biol_rep %in% c(2,3)) %>%
group_by(dose, genotype, time, biol_rep) %>%
summarise(mean_aa = mean(aa)) %>%
mutate(group = as.factor(paste(genotype, biol_rep, sep = "_"))) %>%
data.frame(.)


p <- rel_aa %>%
filter(!is.na(inoc_aa) & biol_rep %in% c(2,3)) %>%
ggplot(aes(x = time, y = log2(aa+.01), colour = genotype, fill = genotype)) +
geom_point(shape = 21, colour = c_black, position = position_jitterdodge(6), size = 0.4, alpha = 0.3) +
geom_line(data = group_df, aes(group=group, x = time, y = log2(mean_aa+.01)), alpha=.6, size=.4) +
ylim(c(-8,8)) +
	scale_colour_manual(values=genotype.syncom$colours) +
	scale_fill_manual(values=genotype.syncom$colours) +
	facet_grid(biol_rep ~ dose, scale = "free", space = "free") +
	theme_RTN
ggsave(p, file=paste(fig, "grouped_overall_growth.png", sep=""), 
		device = "png",
		units = "in",
		dpi = 600,
		width=5.5, height=4, bg="transparent")




# ANOVA for each dose and time
stat_df <- lapply(c("0.005", "0.05"), FUN = function(d){

	stat_list <- lapply(c("0", "24", "72"), FUN = function(t){

		# d = "0.005"
		# t = "24"

		temp <- rel_aa %>%
				filter(time == t, dose == d, !is.na(inoc_aa) & biol_rep %in% c(2,3))

		aov_obj <- broom::tidy((aov(rel_aa ~ genotype, temp))) %>%
		filter(term == "genotype") %>%
		mutate(dose = d, time = t) %>%
		data.frame(.)
	
	}) %>% do.call(rbind.data.frame,.)

	return(stat_list)

}) %>% do.call(rbind.data.frame,.) %>% 
na.omit(.) %>% 
mutate(significance = p.value <= 0.05) %>% 
select(-term)


# Save output
write.table(stat_df, 
	file = "./statistics/growth_curve_statistics_anova.txt", 
	quote = F, 
	row.names = F, 
	sep = "\t")

# stat_obj_anova <- broom::tidy(aov_obj) %>%
# mutate(adj.p = p.adjust(p.value, "BH")) %>%
# data.frame(.)

# write.table(stat_obj_anova, 
# 	file = "./statistics/differential_analysis/growth_curve_statistics_anova_summary.txt", 
# 	quote = F, 
# 	row.names = F, 
# 	sep = "\t")

# plot
idx <- rel_aa$time != 0
p <- ggplot(rel_aa[idx,], aes(x=factor(time), y=rel_aa, colour=genotype)) +
	geom_hline(yintercept=0, size=.2, colour=c_black, linetype="dashed") +
	geom_point(position=position_jitterdodge(jitter.width=.25), alpha=.5, size=.5) +
	geom_boxplot(fill=c_white, outlier.shape=NA) +
	scale_colour_manual(values=genotype.syncom$colours) +
	facet_grid( ~ dose+biol_rep, labeller=label_both) +
	theme_RTN
ggsave(p, file=paste(fig, "overall_growth.relatve_to_0hpi.box.pdf", sep=""), width=6.5, height=4, bg="transparent")

# plot
p <- ggplot(rel_aa %>% filter(biol_rep != 4), aes(x=time, y=log2(aa), colour=genotype)) +
	geom_line(aes(group=geno_exu_group), alpha=.3, size=.2) +
	scale_colour_manual(values=genotype.syncom$colours) +
	facet_grid(biol_rep ~ dose) +
	theme_RTN
ggsave(p, file=paste(fig, "overall_growth.png", sep=""), 
	device = "png", dpi = 600, units = "in",
	width=5.5, height=4, bg="transparent")



# END
sessionInfo()



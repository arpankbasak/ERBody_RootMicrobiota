# Script for estimating beta diversity using absolute abundance
# @Arpan Kumar Basak

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1, mc.cores = 12)

# Loading Dependencies and path
pkgs <- c("tidyverse", "reshape2", "vegan")

lapply(pkgs, require, character.only = T)

source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/parameters.R")
setwd(analysis.syncom)

load("./data/sync_filtered_spike_data.Rdata")

# ===== Beta Diversity =====
source("/biodata/dep_psl/grp_psl/Arpan/syncom/analysis/scripts/functions.R")

mat <- sync.spiked.dat$cData
meta <- sync.spiked.dat$metadata

figs <- paste0(figs, "/beta_diversity/spiked/pcoa/")
stats <- paste0(stats, "/beta_diversity/spiked/pcoa/")
out <- paste0(out, "/beta_diversity/spiked/pcoa/")

# Overall beta diversity
mat.ra <- apply(mat, 2, function(x) x/sum(x))
d <- vegan::vegdist(t(mat), method = "bray")
mds_obj <- cmdscale(as.dist(d), k = 3, eig = TRUE)
variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)

# Plot overall beta diversity
# Prepare MDS table
mds_df <- meta %>%
mutate(MDS1 = mds_obj$points[.$SampleID,1],
MDS2 = mds_obj$points[.$SampleID,2],
MDS3 = mds_obj$points[.$SampleID,3]
)
id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

# Plot canvas
plot_obj <-  mds_df %>%
ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
    # ggtitle(paste0("Dose: ", i)) +
    geom_point(alpha = 0.6, size = 3, 
        aes(colour = genotype, 
            group = genotype,
            shape = time),
        show.legend = T) +
    scale_color_manual(values = genotype.syncom$colours[id]) + 
    scale_shape_manual(values = c(2, 15, 17)) +
    theme_RTN_MDS +
    labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
         y = paste0("PCoA - 2: ",variance[2],"%")) +
ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered.png"), 
       dpi = 600, 
       device = "png", 
       bg = "transparent", 
       units = img, 
       width = cpcoa.box[1], height = cpcoa.box[2])

# Grouped
plot_obj <-  mds_df %>%
ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
    # ggtitle(paste0("Dose: ", i)) +
    geom_point(alpha = 0.6, size = 3, 
        aes(colour = genotype, 
            group = genotype,
            shape = as.factor(tech_rep)),
        show.legend = T) +
    facet_grid(biol_rep + dose ~ time, space = "free", scale = "free", switch = "both") +
    scale_color_manual(values = genotype.syncom$colours[id]) + 
    scale_shape_manual(values = c(1:8)) +
    theme_RTN_MDS +
    theme(strip.text.x = element_text(angle = 90),
      strip.text.y = element_text(angle = 90)) +
    labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
         y = paste0("PCoA - 2: ",variance[2],"%")) +
ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered_grouped.png"), 
       dpi = 600, 
       device = "png", 
       bg = "transparent", 
       units = img, 
       width = cpcoa.box[1]*2, height = cpcoa.box[2]*3)

# Subsetted beta diversity
res <- parallel::mclapply(as.character(unique(meta$dose)), function(i){

    set.seed(seeder)
    # i <- "0.05"
    temp <- meta %>%
    filter(dose == i) %>%
    data.frame(., stringsAsFactors = FALSE)

    temp.mat <- mat[, temp$SampleID]
    temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

    # d <- 1-cor(temp.mat)
    d <- vegan::vegdist(t(temp.mat), method = "bray")
    mds_obj <- cmdscale(as.dist(d), k = 3, eig = TRUE)
    variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)

    message(paste0("Computing for dose ", i, ":"))
    # Plot beta diversity for individual dosage
    # Prepare MDS table
    mds_df <- temp %>%
    mutate(MDS1 = mds_obj$points[.$SampleID,1],
        MDS2 = mds_obj$points[.$SampleID,2],
        MDS3 = mds_obj$points[.$SampleID,3]
        )
    id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

    # Plot canvas
    plot_obj <-  mds_df %>%
    ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
    geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
    geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            ggtitle(paste0("Dose: ", i)) +
            geom_point(alpha = 0.6, size = 3, 
                aes(colour = genotype, 
                    group = genotype,
                    shape = time),
                show.legend = T) +
            scale_color_manual(values = genotype.syncom$colours[id]) + 
            scale_shape_manual(values = c(2, 15, 17)) +
            theme_RTN_MDS +
            labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
                 y = paste0("PCoA - 2: ",variance[2],"%")) +
    ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered_(", i,")_dose.png"), 
               dpi = 600, 
               device = "png", 
               bg = "transparent", 
               units = img, 
               width = cpcoa.box[1], height = cpcoa.box[2])


    # Without inoculum
    temp <- meta %>%
    filter(dose == i, genotype != "inoc") %>%
    data.frame(., stringsAsFactors = FALSE)

    temp.mat <- mat[, temp$SampleID]
    temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

    d <- vegan::vegdist(t(temp.mat), method = "bray")
    mds_obj <- cmdscale(as.dist(d), k = 3, eig = TRUE)
    variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)

    message(paste0("Computing for dose without inoculum ", i, ":"))
    # Plot beta diversity for individual dosage
    # Prepare MDS table
    mds_df <- temp %>%
    mutate(MDS1 = mds_obj$points[.$SampleID,1],
        MDS2 = mds_obj$points[.$SampleID,2],
        MDS3 = mds_obj$points[.$SampleID,3]
        )
    id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

    # Plot canvas
    plot_obj <-  mds_df %>%
    ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
    geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
    geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            ggtitle(paste0("Dose: ", i)) +
            geom_point(alpha = 0.6, size = 3, 
                aes(colour = genotype, 
                    group = genotype,
                    shape = time),
                show.legend = T) +
            scale_color_manual(values = genotype.syncom$colours[id]) + 
            scale_shape_manual(values = c(2, 15, 17)) +
            theme_RTN_MDS +
            labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
                 y = paste0("PCoA - 2: ",variance[2],"%")) +
    ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered_(", i,")_dose_no_inoc.png"), 
               dpi = 600, 
               device = "png", 
               bg = "transparent", 
               units = img, 
               width = cpcoa.box[1], height = cpcoa.box[2])

    set.seed(seeder)
    plot_obj <- list()

    obj <- list()
    length(obj) <- 2
    names(obj) <- c("24", "72")

    # For individual time points and inoculum
    for(j in c(as.character(unique(meta$time)))) {

        if(j != "0") {

            message(paste0("Computing for time ", j, ":"))
            # j = "24"
            temp <- meta %>%
            filter(dose == i & time == j) %>%
            data.frame(., stringsAsFactors = FALSE)

            temp.mat <- mat[, temp$SampleID]
            temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

            # d <- 1-cor(temp.mat)
            d <- vegan::vegdist(t(temp.mat), method = "bray")

            mds_obj <- cmdscale(as.dist(d), k = 3, eig = TRUE)
            variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)

            id <- which(genotype.syncom$names %in% as.character(mds_df$genotype))

            # Plot beta diversity for individual time points
            # Prepare MDS table
            mds_df <- temp %>%
            mutate(MDS1 = mds_obj$points[.$SampleID,1],
                MDS2 = mds_obj$points[.$SampleID,2],
                MDS3 = mds_obj$points[.$SampleID,3]
                )
            id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

            # Plot canvas
            plot_obj <-  mds_df %>%
            ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
            geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
                    ggtitle(paste0("Dose: ", i, " ; Time: ", j)) +
                    geom_point(alpha = 0.6, size = 3, 
                        aes(fill = genotype, 
                            group = genotype,
                            shape = biol_rep),
                        colour = c_black,
                        show.legend = T) +
                    scale_fill_manual(values = genotype.syncom$colours[id]) + 
                    scale_shape_manual(values = c(21:24)) +
                    theme_RTN_MDS +
                    labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
                         y = paste0("PCoA - 2: ",variance[2],"%")) +
            ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered_(", i,")_dose_(", j, ")_time.png"), 
                       dpi = 600, 
                       device = "png", 
                       bg = "transparent", 
                       units = img, 
                       width = cpcoa.box[1], height = cpcoa.box[2])


            # Compute pairwise BC stats relative to Col-0
            stat_df <- as.data.frame(as.matrix(d)) %>%
            add_column(X = row.names(.), .before = 1) %>%
            gather(key = "Y", value = "vals", convert = FALSE, -X) %>%
            mutate(Y_genotype = as.character(temp$genotype)[match(.$Y, temp$SampleID)],
                dose = i, 
                time = j) %>%
            filter(Y_genotype == "Col") %>%
            data.frame(., stringsAsFactors = FALSE)

            obj[[j]] <- list(plot = plot_obj, stats = stat_df)


        } else {

            message(paste0("Computing for time ", j, ":"))
            temp <- meta %>%
            filter(dose == i, time == j) %>%
            data.frame(., stringsAsFactors = FALSE)

            temp.mat <- mat[, temp$SampleID]
            temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))
            # d <- 1-cor(temp.mat)
            d <- vegan::vegdist(t(temp.mat), method = "bray")
            mds_obj <- cmdscale(as.dist(d), k = 3, eig = TRUE)
            variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)
            
            # Plot beta diversity for inocula
            id <- which(genotype.syncom$names %in% as.character(mds_df$genotype))

            # Plot beta diversity for individual time points
            # Prepare MDS table
            mds_df <- temp %>%
            mutate(MDS1 = mds_obj$points[.$SampleID,1],
                MDS2 = mds_obj$points[.$SampleID,2],
                MDS3 = mds_obj$points[.$SampleID,3]
                )
            id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

            # Plot canvas
            plot_obj <-  mds_df %>%
            ggplot(aes(x = saturate(MDS1), y = saturate(MDS2))) +
            geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
                    ggtitle(paste0("Dose: ", i, " ; Time: ", j)) +
                    geom_point(alpha = 0.6, size = 3, 
                        aes(fill = genotype, 
                            group = genotype,
                            shape = biol_rep),
                        colour = c_black,
                        show.legend = T) +
                    scale_fill_manual(values = genotype.syncom$colours[id]) + 
                    scale_shape_manual(values = c(21:24)) +
                    theme_RTN_MDS +
                    labs(x = paste0("PCoA - 1: ",variance[1],"%"), 
                         y = paste0("PCoA - 2: ",variance[2],"%")) +
            ggsave(filename = paste0(figs, "overall_beta_diversity_pcoafiltered_(", i,")_dose_(", j, ")_time.png"), 
                       dpi = 600, 
                       device = "png", 
                       bg = "transparent", 
                       units = img, 
                       width = cpcoa.box[1], height = cpcoa.box[2])

            # Compute stats for inocula
            stat_df <- as.data.frame(as.matrix(d)) %>%
            add_column(X = row.names(.), .before = 1) %>%
            gather(key = "Y", value = "vals", convert = FALSE, -X) %>%
            mutate(Y_genotype = as.character(temp$genotype)[match(.$Y, temp$SampleID)],
                dose = i, 
                time = j) %>%
            # filter(Y_genotype == "Col") %>%
            data.frame(., stringsAsFactors = FALSE)
            
            obj[[j]] <- list(plot = plot_obj, stats = stat_df)
        }
        
    }
    
    return(obj)

}, mc.cores = 4)

# For every dose 

# Compute paiwise PCC
stat_df <- lapply(c(1:length(res)),
    function(x){
        
        df <- rbind.data.frame(res[[x]][[1]]$stats, res[[x]][[2]]$stats)
        return(df)
})

stat_df <- do.call(rbind.data.frame, stat_df) %>%
filter(Y_genotype == "Col" & X != Y) %>%
cbind.data.frame(meta[match(.$X, meta$SampleID), c("genotype", "random")])

mark1 <- mean(stat_df$vals[stat_df$genotype == "Col" & stat_df$dose == "0.05"])
mark2 <- mean(stat_df$vals[stat_df$genotype == "Col" & stat_df$dose == "0.005"])

# Plot canvas
(bxp <- stat_df %>% ggplot(aes(x = genotype, y = vals)) +
geom_hline(yintercept = c(mark1, mark2), 
    colour = "darkgrey", 
    lty = "solid", 
    alpha = 0.6, lwd = 0.8) +
geom_jitter(aes(colour = genotype), alpha = 0.6, size = 1, shape = 2, position = position_jitterdodge(jitter.width = 0.2)) +
geom_boxplot(outlier.alpha = 0, size = 0.8, fill = c_white, aes(colour = genotype)) +
facet_grid(dose + time ~., switch = "y", scales = "free_x") +
scale_fill_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% stat_df$genotype)], 
    label = geno.label[which(genotype.syncom$short %in% stat_df$genotype)]
    ) +
scale_colour_manual(values = genotype.syncom$colours[which(genotype.syncom$short %in% stat_df$genotype)], 
    label = geno.label[which(genotype.syncom$short %in% stat_df$genotype)], guide=FALSE
    ) +
theme_RTN +
theme(
    panel.spacing = unit(0.01, "lines"),
    strip.text.x = element_text(hjust = 0.5, vjust = 0.5), 
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(angle = 60)
    ) +
labs(x = "", y = "Euc(log2(reads/reads_spike))", colour = "",  fill = "") +
coord_flip()) +
ggsave(filename = paste0(figs, "Aitchison_boxplot.png"), 
           bg = "transparent", units = img, 
           width=5, height=10, limitsize=F, device = "png", dpi = 600)


# Make a plot composite
p1 <- cowplot::plot_grid(
  res[[1]][[1]]$plot,
  res[[1]][[2]]$plot,
  res[[2]][[1]]$plot,
  res[[2]][[2]]$plot,
  ncol=2,
  nrow = 2, 
  align='hv', 
  axis = 'tblr',
  labels=NULL, 
  byrow = TRUE,
  rel_widths = c(1,1), 
  rel_heights = c(1, 1)
) 

cowplot::plot_grid(
  p1, bxp,
  ncol=2,
  nrow = 1, 
  align='hv', 
  axis = 'tblr',
  labels=NULL, 
  # byrow = TRUE,
  rel_widths = c(1.5,1), 
  rel_heights = c(1.5, 1),
  scale = c(1, 0.5)
) +
ggsave(filename = paste0(figs, "PCOA_spiked_summary.png"), 
           bg="transparent", 
           units = "in",
           width=16, height=9, 
           limitsize=F, device = "png", dpi = 600)


# END OF SCRIPT
sessionInfo()
# Script for conducting constrained ordination
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
meta <- sync.spiked.dat$metadata %>% filter(biol_rep != 1, exudates_rep != 3)

figs <- paste0(figs, "/beta_diversity/spiked/cpcoa/")
stats <- paste0(stats, "/beta_diversity/spiked/cpcoa/")
out <- paste0(out, "/beta_diversity/spiked/cpcoa/")

# Compute individual cPCOA factors
# comb <- list(
#   fixed = formula(t(mat) ~ fixed + Condition(random)),
#   random = formula(t(mat) ~ random + Condition(fixed)),
#   genotype = formula(t(mat) ~ genotype + Condition(random + time)),
#   genotype_time = formula(t(mat) ~ genotype + time + Condition(random)),
#   time = formula(t(mat) ~ time + Condition(random + genotype)),
#   dose = formula(t(mat) ~ dose + Condition(random + genotype + time)),
#   dose_time = formula(t(mat) ~ dose + time + Condition(random + genotype))
# )


# # cPCoA including inoculum
# cpcoa_overall <- parallel::mclapply(as.character(names(comb)), function(i){

#   f <- comb[[i]]
#   cca_obj <- capscale(formula = f, 
#                         meta, 
#                         distance = "bray", 
#                         sqrt.dist = T, 
#                         add = F, 
#                         na.action = na.exclude)
  
#   # cca_obj <- cca(formula = f, 
#   #                       meta, 
#   #                       method = "pearson")

#   # Conduct PERMANOVA
#   set.seed(seeder)
#   p.nova <- anova.cca(cca_obj)
  
#   # Extract variables from stat objects
#   eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
#   pval <- p.nova$`Pr(>F)`[1]
#   chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
#   variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
#                          row.names = c("total", "constrianed", "unconstrained"))
  
#   write.table(variable, paste0(stats, "Overall_CPCoA_(", i,")_stats.txt"), 
#       sep = "\t", quote = F)


#   # Format P-Value for representation
#   variance <- format(100 * variable$proportion[2], digits = 4)
#   tf <- paste(as.character(f)[2], as.character(f)[1], as.character(f)[3])
#   ti <- paste0("Constrained by: ", i, "; [", tf ,"] :  Constrained = ", variance,"% p = ", pval, 
#                ifelse(pval < 0.05, "*", " (NS)"))
  

#   # Plot beta diversity for individual time points
#   # Prepare MDS table
#   mds_df <- meta %>%
#   mutate(CCA1 = as.numeric(cca_obj$CCA$wa[match(.$SampleID, row.names(cca_obj$CCA$wa)),1]),
#     CCA2 = as.numeric(cca_obj$CCA$wa[match(.$SampleID, row.names(cca_obj$CCA$wa)),2]))

#   id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

#   if(i!="random"){

#       # Plot canvas
#     (plot_obj <-  mds_df %>%
#     ggplot2::ggplot(aes(x= CCA1, y = CCA2)) +
#             ggtitle(ti) +
#             geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
#             geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
#             geom_point(alpha = 0.6, size = 3 , 
#               aes(shape = biol_rep, 
#                 fill = genotype), colour = c_black) +
#             # stat_ellipse(type = "norm", level = 0.95, linetype = 2, alpha = 0.5, size = 0.8, 
#             #     aes(group = genotype)) +
#             scale_fill_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
#             scale_colour_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
#             #scale_fill_manual(values = col.pal, labels = geno.label) + 
#             scale_shape_manual(values = c(21,23,24,25)) +
#             theme_RTN_MDS +
#             labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
#                  y = paste0("cPCoA - 2: ",eig[2],"%"))) +
#     ggsave(filename = paste0(figs, "overall_beta_diversity_cpcoa_(factor_", i,").png"), 
#                dpi = 600, 
#                device = "png", 
#                bg = "transparent", 
#                units = img, 
#                width = cpcoa.box[1], 
#                height = cpcoa.box[2])
#   } else{

#     # Plot canvas --random only
#     (plot_obj <-  mds_df %>%
#       mutate(primer = (Fwd_barcode)) %>%
#     ggplot2::ggplot(aes(x= CCA1, y = CCA2)) +
#             ggtitle(ti) +
#             geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
#             geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
#             geom_point(alpha = 0.6, size = 3 , 
#               aes(
#                 shape = biol_rep, 
#                 # size = primer, 
#               fill = primer), colour = c_black) +
#             # stat_ellipse(type = "norm", level = 0.95, linetype = 2, alpha = 0.5, size = 0.8, 
#             #     aes(group = genotype)) +
#             # scale_fill_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
#             # scale_colour_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
#             #scale_fill_manual(values = col.pal, labels = geno.label) + 
#             scale_shape_manual(values = c(21,23,24,25)) +
#             theme_RTN_MDS +
#             labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
#                  y = paste0("cPCoA - 2: ",eig[2],"%"))) +
#     ggsave(filename = paste0(figs, "overall_beta_diversity_cpcoa_(factor_", i,").png"), 
#                dpi = 600, 
#                device = "png", 
#                bg = "transparent", 
#                units = img, 
#                width = cpcoa.box[1]*3, 
#                height = cpcoa.box[2]*3)
#   }
    

#   return(plot_obj)
# }, mc.cores = 4)

# Subsetted constrained beta diversity
# Subsetted by dose and time --for each subset make with and without inocula
res <- parallel::mclapply(as.character(unique(meta$dose)), function(i){

    set.seed(seeder)
    plot_obj <- list()

    obj <- list()
    length(obj) <- 2
    names(obj) <- c("24", "72")

    # For individual dose with inoculum
    temp <- meta %>%
    filter(dose == i) %>%
    data.frame(., stringsAsFactors = FALSE)

    row.names(temp) <- temp$SampleID
    temp.mat <- mat[, temp$SampleID]
    temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

    f <- formula(t(temp.mat) ~ genotype*time + Condition(random))
    cca_obj <- capscale(formula = f, 
                data = temp, 
                distance = "bray", 
                sqrt.dist = T, 
                add = F, 
                na.action = na.exclude)

    # Conduct PERMANOVA
    set.seed(seeder)
    p.nova <- anova.cca(cca_obj)
    
    # Extract variables from stat objects
    eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
    pval <- p.nova$`Pr(>F)`[1]
    chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
    variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                           row.names = c("total", "constrianed", "unconstrained"))
    
    write.table(variable, paste0(stats, "CPCoA_(", i,")_dose_genotype.time_stats.txt", sep = ""), 
        sep = "\t", quote = F)
    
    # Format P-Value for representation
    variance <- format(100 * variable$proportion[2], digits = 4)
    tf <- paste(as.character(f)[2], as.character(f)[1], as.character(f)[3])
    ti <- paste0("Dosage: ", i, "; [", tf ,"] :  Genotype = ", variance,"% p = ", pval, 
                 ifelse(pval < 0.05, "*", " (NS)"))
    

    # Plot beta diversity for individual time points
    # Prepare MDS table
    mds_df <- temp %>%
    mutate(CCA1 = as.numeric(cca_obj$CCA$wa[.$SampleID,1]),
        CCA2 = as.numeric(cca_obj$CCA$wa[.$SampleID,2])
        )

    id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

    # Plot canvas
    plot_obj <-  mds_df %>%
    ggplot2::ggplot(aes(x= CCA1, y = CCA2, colour = genotype)) +
            ggtitle(ti) +
            geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_point(alpha = 0.6, size = 3 , 
              aes(shape = biol_rep, fill = genotype), colour = c_black) +
            # stat_ellipse(type = "norm", level = 0.95, linetype = 2, alpha = 0.5, size = 0.8, 
            #     aes(colour = genotype)) +
            scale_fill_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
            scale_colour_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
            #scale_fill_manual(values = col.pal, labels = geno.label) + 
            scale_shape_manual(values = c(21,23,24,25)) +
            theme_RTN_MDS +
            labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
                 y = paste0("cPCoA - 2: ",eig[2],"%")) +
    ggsave(filename = paste0(figs, "overall_beta_diversity_cpcoafiltered_(", i,")_dose_genotype.time.png"), 
               dpi = 600, 
               device = "png", 
               bg = "transparent", 
               units = img, 
               width = cpcoa.box[1], height = cpcoa.box[2])


    # For individual dose without inoculum
    temp <- meta %>%
    filter(dose == i, genotype != "inoc") %>%
    data.frame(., stringsAsFactors = FALSE)

    row.names(temp) <- temp$SampleID
    temp.mat <- mat[, temp$SampleID]
    temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))
    row.names(temp) <- temp$SampleID
    f <- formula(t(temp.mat) ~ genotype*time + Condition(random))
    cca_obj <- capscale(formula = f, 
                data = temp, 
                distance = "bray", 
                sqrt.dist = T, 
                add = F, 
                na.action = na.exclude)

    # Conduct PERMANOVA
    set.seed(seeder)
    p.nova <- anova.cca(cca_obj)
    
    # Extract variables from stat objects
    eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
    pval <- p.nova$`Pr(>F)`[1]
    chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
    variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                           row.names = c("total", "constrianed", "unconstrained"))
    
    write.table(variable, paste0(stats, "CPCoA_(", i,")_dose_genotype.time_stats_no_inoc.txt", sep = ""), 
        sep = "\t", quote = F)
    
    # Format P-Value for representation
    variance <- format(100 * variable$proportion[2], digits = 4)
    tf <- paste(as.character(f)[2], as.character(f)[1], as.character(f)[3])
    ti <- paste0("Dosage: ", i, "; [", tf ,"] :  Genotype = ", variance,"% p = ", pval, 
                 ifelse(pval < 0.05, "*", " (NS)"))
    

    # Plot beta diversity for individual time points
    # Prepare MDS table
    mds_df <- temp %>%
    mutate(CCA1 = as.numeric(cca_obj$CCA$wa[.$SampleID,1]),
      CCA2 = as.numeric(cca_obj$CCA$wa[.$SampleID,2])
        )

    id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

    # Plot canvas
    plot_obj <-  mds_df %>%
    ggplot2::ggplot(aes(x= CCA1, y = CCA2, colour = genotype)) +
            ggtitle(ti) +
            geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
            geom_point(alpha = 0.6, size = 3 , 
              aes(shape = exudates_rep, fill = genotype), colour = c_black) +
            # stat_ellipse(type = "norm", level = 0.95, linetype = 2, alpha = 0.5, size = 0.8, 
            #     aes(group = genotype)) +
            scale_color_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
            scale_fill_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
            #scale_fill_manual(values = col.pal, labels = geno.label) + 
            scale_shape_manual(values = c(21,23,24,25)) +
            theme_RTN_MDS +
            labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
                 y = paste0("cPCoA - 2: ",eig[2],"%")) +
    ggsave(filename = paste0(figs, "overall_beta_diversity_cpcoafiltered_(", i,")_dose_genotype.time_no_inoc.png"), 
               dpi = 600, 
               device = "png", 
               bg = "transparent", 
               units = img, 
               width = cpcoa.box[1], height = cpcoa.box[2])


    # For individual time points and inoculum
    for(j in c(as.character(unique(meta$time)))) {

        if(j != "0") {

            # i = "0.05";j = "24";k="cyp"
            temp <- meta %>%
            filter(dose == i & time == j) %>%
            # mutate(
            #   group = as.factor(paste(genotype, exudates_rep, biol_rep, sep = "_")),
            #   random1 = as.factor(paste(lib, sep = "_"))) %>%
            data.frame(., stringsAsFactors = FALSE)
            row.names(temp) <- temp$SampleID
            temp.mat <- mat[, temp$SampleID]
            temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

            # f <- formula(t(temp.mat) ~ group + Condition(random1))
            f <- formula(t(temp.mat) ~ genotype + Condition(random))
            cca_obj <- capscale(formula = f, 
                        data = temp, 
                        distance = "bray", 
                        sqrt.dist = T, 
                        add = F, 
                        na.action = na.exclude)
            
            # cca_obj <- cca(formula = f, 
            #             temp, 
            #             method = "pearson")

            # Conduct PERMANOVA
            set.seed(seeder)
            p.nova <- anova.cca(cca_obj)
            
            # Extract variables from stat objects
            eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
            pval <- p.nova$`Pr(>F)`[1]
            chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
            variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                                   row.names = c("total", "constrianed", "unconstrained"))
            
            write.table(variable, paste0(stats, "CPCoA_(", i,")_dose_(", j, ")_time_stats.txt", sep = ""), 
                sep = "\t", quote = F)
            
            # Format P-Value for representation
            variance <- format(100 * variable$proportion[2], digits = 4)
            tf <- paste(as.character(f)[2], as.character(f)[1], as.character(f)[3])
            ti <- paste0("Dosage: ", i, "; Time: ", j, "; [", tf ,"] :  Genotype = ", variance,"% p = ", pval, 
                         ifelse(pval < 0.05, "*", " (NS)"))
            

            markers <- as.data.frame(cca_obj$CCA$centroids) %>%
            mutate(group = str_replace_all(row.names(.), "genotype", "")) %>%
            # separate(group, into = c("A", "B"), sep = ":", convert = FALSE, remove = FALSE) %>%
            # mutate(A = factor(A)) %>%
            select(group, CAP1, CAP2)

            # Plot beta diversity for individual time points
            # Prepare MDS table
            mds_df <- temp %>%
            mutate(CCA1 = as.numeric(cca_obj$CCA$wa[.$SampleID,1]),
                CCA2 = as.numeric(cca_obj$CCA$wa[.$SampleID,2]),
                cen1 = markers$CAP1[match(as.character(genotype), markers$group)],
                cen2 = markers$CAP2[match(as.character(genotype), markers$group)]
                )

            id <- which(genotype.syncom$short %in% as.character(mds_df$genotype))

            # Plot canvas
            plot_obj <-  mds_df %>%
            ggplot2::ggplot(aes(x= CCA1, y = CCA2, colour = genotype)) +
                    ggtitle(ti) +
                    # geom_hline(yintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
                    # geom_vline(xintercept = 0.0, linetype = "solid", colour = c_grey, alpha = 0.3, size = 0.2) +
                    geom_segment(aes(xend = (CCA1), yend = (CCA2),
                        x = (cen1), y = (cen2), 
                        colour = genotype), lty = "solid", 
                    lwd = 0.3) +
                    geom_point(alpha = 0.6, size = 3 , 
                      aes(shape = exudates_rep, fill = genotype), colour = c_black) +
                    # stat_ellipse(type = "norm", level = 0.95, linetype = 2, alpha = 0.5, size = 0.8, 
                    #     aes(group = genotype)) +
                    scale_color_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
                    scale_fill_manual(values = genotype.syncom$colours[id], labels = geno.label[id]) + 
                    #scale_fill_manual(values = col.pal, labels = geno.label) + 
                    scale_shape_manual(values = c(21,23,24,25)) +
                    theme_RTN_MDS +
                    labs(x = paste0("cPCoA - 1: ",eig[1],"%"), 
                         y = paste0("cPCoA - 2: ",eig[2],"%")) +
            ggsave(filename = paste0(figs, "overall_beta_diversity_cpcoafiltered_(", i,")_dose_(", j, ")_time.png"), 
                       dpi = 600, 
                       device = "png", 
                       bg = "transparent", 
                       units = img, 
                       width = cpcoa.box[1], height = cpcoa.box[2])


            # Compute pairwise PERMANOVA with Col-0 as intercept
            pwnova <- lapply(c("pyk10", "cyp"), function(k){

                # Filter dataset for pairwise multi-variate constrained ordination
                temp <- meta %>%
                filter(dose == i, time == j, genotype %in% c("Col", k)) %>%
                data.frame(., stringsAsFactors = FALSE)
                row.names(temp) <- temp$SampleID
                temp.mat <- mat[, temp$SampleID]
                temp.mat <- apply(temp.mat, 2, function(x) x/sum(x))

                # d <- vegdist(t(temp.mat.ra), "bray")
                # d <- as.dist(1-cor(temp.mat))
                ad_obj <- adonis2(t(temp.mat) ~ genotype + random,
                    sqrt.dist = TRUE, 
                    method = "bray",
                    data = temp,
                    strata = temp$genotype,
                    permutations = 1000)

                stat_sum <- broom::tidy(ad_obj)
                
                temp <- stat_sum %>% 
                dplyr::filter(term == "genotype") %>%
                  mutate(
                    con.var.explained = .$SumOfSqs/stat_sum$SumOfSqs[which(stat_sum$term == "Total")],
                    genotype = k,
                    time = j,
                    dose = i) %>%
                  data.frame(., stringsAsFactors = FALSE)
                return(temp)
            })

            df_pwnova <- do.call(rbind.data.frame, pwnova)
            obj[[j]] <- list(plot = plot_obj, stats = df_pwnova)



        } else {

            message(paste0("Not applicable for inoculum ", j))
        }
        
    }

    return(obj)

}, mc.cores = 12)

# Make variance explained barplots
stat_df <- lapply(c(1:length(res)),
    function(x){
        
        df <- rbind.data.frame(res[[x]][[1]]$stats, res[[x]][[2]]$stats)
        return(df)
})

stat_df <- do.call(rbind.data.frame, stat_df) %>%
mutate(significance = p.value <= 0.05, 
    genotype = factor(genotype, levels = c("pyk10", "cyp")),
    dose = as.factor(dose),
    time = as.factor(time),
    variance_explained = 100*con.var.explained) %>%
select(-term)

write.table(stat_df, 
    paste0(stats,"sync_pw_CPCOA.txt"), sep = "\t", quote = FALSE)


# Plot canvas
ln <- median(as.numeric(stat_df$variance_explained))

# Barplot
(bp <- stat_df %>%
mutate(tag = paste0(as.character(round(variance_explained, 2)), "%,\n", as.character(round(p.value, 4)))) %>%
ggplot(aes(x = genotype, y = variance_explained)) +
# geom_hline(yintercept = ln, 
#     colour = "darkgrey", 
#     lty = "solid", 
#     alpha = 0.6, lwd = 0.8) +
geom_bar(stat = "identity", aes(colour = significance, 
    fill = genotype, alpha = significance)) +
geom_text(aes(label = tag, alpha = significance), nudge_y = 1, size = 1) +
facet_grid(dose + time ~., space = "free", scale = "free_y") +
scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 1), guide = FALSE) +
scale_fill_manual(values = genotype.syncom$colour[which(genotype.syncom$short %in% stat_df$genotype)], 
    guide = FALSE) +
coord_flip() +
theme_RTN +
theme() +
labs(y = "Variance explained %", x = "", colour = "", fill = "pvalue", alpha = "")) +
ggsave(filename = paste0(figs, "pw_PNOVA_genotype.png"), 
           bg="transparent", 
           units = "in",
           width=3, height=6, 
           limitsize=F, device = "png", dpi = 600)

# Make CPCOA plot composite
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
  p1, bp,
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
ggsave(filename = paste0(figs, "cPCOA_spiked_summary.png"), 
           bg="transparent", 
           units = "in",
           width=16, height=9, 
           limitsize=F, device = "png", dpi = 600)

# Overall cPCOA for diagnosis
# p1 <- cowplot::plot_grid(
#   cpcoa_overall$fixed + theme(legend.position = "none"),
#   cpcoa_overall$genotype + theme(legend.position = "none"),
#   cpcoa_overall$genotype_time + theme(legend.position = "none"),
#   cpcoa_overall$time + theme(legend.position = "none"),
#   cpcoa_overall$dose + theme(legend.position = "none"),
#   cpcoa_overall$dose_time + theme(legend.position = "none"),
#   ncol=3,
#   nrow = 3, 
#   align='hv', 
#   axis = 'tblr',
#   labels= names(comb), 
#   byrow = TRUE,
#   rel_widths = c(1, 1, 1), 
#   rel_heights = c(1, 1, 1)
# )

# cowplot::plot_grid(
#   cowplot::get_legend(cpcoa_overall$fixed), 
#   p1,
#   ncol=1,
#   nrow = 2, 
#   align='hv', 
#   axis = 'tblr',
#   labels=NULL, 
#   # byrow = TRUE,
#   rel_widths = c(0.5, 2), 
#   rel_heights = c(0.5, 2),
#   scale = c(1, 1)
# ) +
# ggsave(filename = paste0(figs, "overall_spiked_cPCOA_summary.png"), 
#       bg = "transparent", 
#       units = "in",
#       width = 16, 
#       height = 9, 
#       limitsize = F, 
#       device = "png", 
#       dpi = 600)

# END OF SCRIPT
sessionInfo()
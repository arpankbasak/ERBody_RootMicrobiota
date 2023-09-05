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
load(paste0("./data/lfc_DAMs_cologne.Rdata"))
load(paste0("./data/msms_data.Rdata"))

meta <- ms_data$meta
mat <- ms_data$mat[, meta$given_id]
geno <- genotype[which(genotype$short %in% meta$genotype),]


FC_threshold <- 2
alpha <- 0.05

# mclapply(c(names(DAM_glm)), function(x){

  # x = "mat_p"
  df <- as.data.frame(lfc_obj$logFC_P.wide)
  id <- str_detect(colnames(df), "_logFC$")
  d <- dist((df[,id]))
  hc <- hclust(as.dist(d), "ward.D2")
  lfc_sorted <- row.names(df)[hc$order]
  # mat <- as.data.frame(DAM_glm[[x]]$voom_mat)

  # # PCA 
  # d <- 1-cor((mat))
  # mds_obj <- cmdscale(d, eig = TRUE)
  # # hc <- hclust(d, "ward.D2")
  # variance <- round(100*mds_obj$eig/sum(mds_obj$eig), 2)

  # # Make data frame
  # mds_df <- ms_data$meta %>%
  # cbind.data.frame(., MDS1 = mds_obj$points[match(.$id, row.names(mds_obj$points)),1],
  #  MDS2 = mds_obj$points[match(.$id, row.names(mds_obj$points)),2]) %>%
  # mutate(genotype = factor(genotype, levels = geno$short))

  # # Plot Canvas
  # (mdsn <- mds_df %>%
  # ggplot(aes(x = saturate(MDS1), y= saturate(MDS2))) +
  # geom_point(aes(fill = genotype), colour = c_black, shape = 21, size = 5) +
  # geom_hline(yintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
  # geom_vline(xintercept = 0, lwd = 0.5, lty = "solid", colour = c_grey) +
  # scale_fill_manual(values = geno$colours) +
  # theme_RTN_MDS +
  # theme() +
  # labs(x = paste0("PCA1:", variance[1], "%"), 
  #   y = paste0("PCA2:", variance[2], "%"), colour = "")) +
  # ggsave(file = paste0(figs, "PCA_", x,"_mode.png"),
  #   width = cpcoa.box[1],
  #   height = cpcoa.box[2],
  #   units = img, 
  #   bg = "transparent"
  #   )



  # Carpent data
  lfc_df <- as.data.frame(lfc_obj$logFC_P.long) %>%
  mutate(
    Genotype = factor(contrast, levels = geno$short[which(geno$short %in% .$contrast)]),
    sig = PValue <= alpha,
    peak_id = factor(peak_id, levels = lfc_sorted))

  # Plot heatmap
  lfc_df %>%
  ggplot(aes(x = Genotype, y=peak_id, fill = saturate(logFC))) +
  geom_raster() +
  # facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  scale_fill_gradient2(midpoint = 0, high = c_green, mid = c_white, low = c_cudo_magenta) +
  theme_RTN +
  theme(axis.text.y = element_blank()) +
  labs(x = "", y = "", fill = "logFC") +
  ggsave(file = paste0(figs, "heatmap_lfc.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 7.5
  )

  lfc_df_sig <- lfc_df %>% filter(sig == TRUE)
  lfc_df_sig %>%
  ggplot(aes(x = Genotype, y=peak_id, fill = saturate(logFC))) +
  geom_raster() +
  # facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  scale_fill_gradient2(midpoint = 0, high = c_green, mid = c_white, low = c_cudo_magenta) +
  theme_RTN +
  theme(axis.text.y = element_blank()) +
  labs(x = "", y = "", fill = "logFC") +
  ggsave(file = paste0(figs, "heatmap_lfc_sig.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 4.5
  )



  stat <- cor.test(df$cyp_logFC, df$pyk10_logFC)
  til <- paste0("Pearson: ", round(stat$estimate, 4),"; t-value: ", round(stat$statistic, 4), "; p-value: ", round(stat$p.value, 4))

  # Scatterplot
  p <- df %>%
  mutate(sig = cyp_PValue <= alpha | pyk10_PValue <= alpha) %>%
  ggplot(aes(x = cyp_logFC, y = pyk10_logFC)) +
  ggtitle(paste0(til)) +
  geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
  geom_point(aes(fill = sig), 
              shape = 21, 
              colour = c_black, 
              alpha = 0.75,
              na.rm = FALSE) +
  scale_fill_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
  theme_RTN +
  theme(panel.spacing = unit(0.5, "lines"), 
                                  axis.text.x = element_text(size = 12),
                                  axis.text.y = element_text(size = 12),
                                  # axis.title = element_text(size = 6, hjust = 0.5, vjust = 0.5),
                                  strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
                                  legend.text = element_text(size = 2),
                                  legend.title = element_blank(),
                                  axis.line.y = element_line(size = 1),
                                  axis.line.x = element_line(size = 1),
                                  axis.ticks.y = element_line(size = 1),
                                  axis.ticks.x = element_line(size = 1),,
                                  plot.title = element_text(size=4)
                                  )
  ggsave(p, file = paste0("./figures/Parwise_metabolome_cyp_pyk_logFC.png"), 
                    dpi = 600, 
                    units = img, 
                    device = "png", 
                    width = 4, 
                    height = 5, 
                    limitsize = T)

  # fdr_line <- -log10(alpha)

  # # Volcanoplot
  # lfc_df %>%
  # ggplot(aes(x = logFC, y = -log10(PValue))) +
  # geom_hline(yintercept = c(fdr_line), 
  #               lty = "solid", 
  #               colour = c_red, 
  #               lwd = 1) +
  # geom_point(aes(
  #             fill = sig, 
  #             size = AveExpr), 
  #             shape = 21, 
  #             colour = c_black, 
  #             alpha = 0.75,
  #             na.rm = FALSE) +
  # facet_grid(.~ Genotype, scale = "free_x", space = "free", switch = "x") +
  # scale_fill_manual(values = c(`FALSE` = c_dark_grey, `TRUE` = c_black), guide = FALSE) +
  # scale_size_continuous(range = c(0.2,2)) +
  # theme_RTN +
  # theme(panel.spacing = unit(0.5, "lines"), 
  #                                 axis.text.x = element_text(size = 12),
  #                                 axis.text.y = element_text(size = 12),
  #                                 # axis.title = element_text(size = 6, hjust = 0.5, vjust = 0.5),
  #                                 strip.text.y = element_text(angle = 180, size = 3, vjust = .5, hjust = .5),
  #                                 legend.text = element_text(size = 2),
  #                                 legend.title = element_blank(),
  #                                 axis.line.y = element_line(size = 1),
  #                                 axis.line.x = element_line(size = 1),
  #                                 axis.ticks.y = element_line(size = 1),
  #                                 axis.ticks.x = element_line(size = 1),,
  #                                 plot.title = element_text(size=4)
  #                                 ) +
  # ggsave(file = paste0("./figures/volcano_metabolome_",x,"_logFC.png"), 
  #                   dpi = 600, 
  #                   units = img, 
  #                   device = "png", 
  #                   width = 6, 
  #                   height = 4, 
  #                   limitsize = T)

  # Number of differentially abundant metabolites
  df <- as.data.frame(lfc_obj$logFC_P.wide)[unique(as.character(lfc_df_sig$peak_id)),]
  enriched <- list(
      # nai1 = row.names(df)[
      #     which(sign(df[,str_detect(colnames(df), "nai1_logFC")]) == 1 )],
      pyk10 = row.names(df)[
          which(sign(df[,str_detect(colnames(df), "pyk10_logFC")]) == 1)],
      cyp = row.names(df)[
          which(sign(df[,str_detect(colnames(df), "cyp_logFC")]) == 1)]
      # myb = row.names(df)[
      # which(sign(df[,str_detect(colnames(df), "myb_logFC")]) == 1)]

  )

  # Enriched
  VennDiagram::venn.diagram(x = enriched,
                            category.names = names(enriched),
                            filename = paste(figs, "/DAMs-enrichment_VenD.png", sep = ""),
                            imagetype = "png",
                            sigdigs = 2,
                            hyper.test = TRUE, 
                            lower.tail = TRUE,
                            output = T,
                            height = 500,
                            width = 500,
                            resolution = 600,
                            lwd = 3,
                            lty = "blank",
                            fill = geno$colours[which(geno$short %in% names(enriched))],
                            cex = 0.2,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 0.2,
                            cat.fontface = "bold",
                            cat.default.pos = "outer", 
                            main =  paste0("Differential Metabolites: (|LFC| >", FC_threshold,
                                          " & FDR < ", alpha,") - vs Col-0"), 
  #                           print.mode = "percent",
                            main.cex = 0.2

  )
  depleted <- list(
     # nai1 = row.names(df)[
     #      which(sign(df[,str_detect(colnames(df), "nai1_logFC")]) == -1)],
      pyk10 = row.names(df)[
          which(sign(df[,str_detect(colnames(df), "pyk10_logFC")]) == -1)],
      cyp = row.names(df)[
          which(sign(df[,str_detect(colnames(df), "cyp_logFC")]) == -1)]
      # myb = row.names(df)[
      # which(sign(df[,str_detect(colnames(df), "myb_logFC")]) == -1)]
  )
  # Depleted
  VennDiagram::venn.diagram(x = depleted,
                            category.names = names(depleted),
                            filename = paste(figs, "/DAMs-depleted_VenD.png", sep = ""),
                            imagetype = "png",
                            sigdigs = 2,
                            hyper.test = TRUE, 
                            lower.tail = TRUE,
                            output = T,
                            height = 500,
                            width = 500,
                            resolution = 600,
                            lwd = 3,
                            lty = "blank",
                            fill = geno$colours[which(geno$short %in% names(depleted))],
                            cex = 0.2,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 0.2,
                            cat.fontface = "bold",
                            cat.default.pos = "outer", 
                            main =  paste0("Depleted Metabolites (|LFC| >", FC_threshold,
                                          " & FDR < ", alpha,") -vs Col-0"), 
  #                           print.mode = "percent",
                            main.cex = 0.2

)

# d <- 1-cor(t(mat))
# hc <- hclust(as.dist(d), "ward.D2")
# compound_clustered <- row.names(d)[hc$order]

# TArgetted peaks
# target <- c("322.09335", "247.05484", "318.30164", "237.11282", "160.03337", "306.06242")
target <- c(
  RA = "162.9323", 
  ICA = "160.0334", 
  Camalexin = "199.0248", 
  IAA = "174.0496", 
  Glucobrassicin = "477.0439", 
  Neoglucobrassicin = "477.0445",
  ICN = "169.0432",
  MI3AN = "185.0738",
  I3AA = "173.0748",
  Glucoberteroin = "434.0738"
  )

df <- as.data.frame(mat) %>% 
add_column(mz = row.names(.)) %>%
gather(key = "SampleID", value = "RA", convert = FALSE, -mz) %>%
cbind.data.frame(., meta[match(.$SampleID, meta$given_id),]) %>%
mutate(
  mz = str_replace_all(mz, levels = lfc_sorted),
  log_ra = ifelse(!is.finite(log10(RA)), NA, log10(RA))
  ) %>%
data.frame(.)

# plot(, pch=20 , cex=1.5 , col="#69b3a2")
p <- ggpairs(as.data.frame(log10(mat)+1e-9))
ggsave(p, filename = paste(figs, "/correlation_between_samples.png", sep = ""), 
  device = "png",
  width = 30, 
  height = 30, 
  bg = c_white, 
  units = "in", 
  dpi=72)


parallel::mclapply(target, function(x){

  temp <- df %>%
  separate(mz, into = c("MZ", "RT", "mode"), sep = ".\\/", remove = FALSE) %>%
  mutate(MZ = round(as.numeric(MZ), 4)) %>%
  filter(between(MZ,as.numeric(x)+0.001,as.numeric(x)-0.001))


  if(nrow(temp) != 0){

  
    p <- temp %>%
    # group_by(mz, genotype, batch) %>%
    # summarise(mean_ra = mean(RA)) %>%
    # mutate(log_ra = (ifelse(!is.finite(log10(mean_ra)), NA, log10(mean_ra)))) %>%
    data.frame(.) %>%
    ggplot(aes(x = genotype, y=log_ra)) +
    ggtitle(paste(names(target)[target == x], "Signal", x, sep = "_")) +
    geom_point(size = 1, shape = 21, position = position_jitterdodge(2), alpha = 0.8, aes(fill = genotype), colour = c_black) +
    geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.6) +
    scale_fill_manual(values= genotype$colours[which(genotype$short %in% temp$genotype)]) +
    facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
    scale_colour_manual(values= genotype$colours[which(genotype$short %in% temp$genotype)]) +
    theme_RTN +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, size = 6)) +
    labs(x = "", y = "", fill = "logFC")
    ggsave(p, file = paste0(figs, "boxplot_ra_", names(target)[target == x],".png"),
    units = "in",
    device = "png",
    width = 4,
    height = 3
  )

  } else{
    message("None ", names(target)[target == x]," signal detected")
  }

}, mc.cores = 4)


p <- df %>%
  ggplot(aes(x = SampleID, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  facet_grid(.~ batch + genotype, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(
    axis.text.y = element_blank(), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    ) +
  labs(x = "", y = "", fill = "log10RA")
  ggsave(p, file = paste0(figs, "heatmap_ra.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 7.5
  )

p <- df %>%
  filter(as.character(mz) %in% target) %>%
  ggplot(aes(x = SampleID, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  facet_grid(.~ genotype + batch, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, size = 6)) +
  labs(x = "", y = "", fill = "logFC") +
  ggsave(p, file = paste0(figs, "heatmap_ra_diagnosis.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 3.5
  )

p <- df %>%
group_by(mz, genotype, batch) %>%
summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(mean_ra)), NA, log10(mean_ra))),
  scale_ra = scale(mean_ra)) %>%
  ggplot(aes(x = genotype, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  # scale_fill_gradient2(midpoint = 0, high = c_red, mid = c_white, low = c_blue, na.value = c_black) +
  facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(
    axis.text.y = element_blank(), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    ) +
  labs(x = "", y = "", fill = "log10RA") +
  ggsave(file = paste0(figs, "heatmap_meanra.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 7.5
  )

p <- df %>%
filter(as.character(mz) %in% target) %>%
group_by(mz, genotype, batch) %>%
summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(mean_ra)), NA, log10(mean_ra))),
  scale_ra = scale(mean_ra)) %>%
  ggplot(aes(x = genotype, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  # scale_fill_gradient2(midpoint = 0, high = c_red, mid = c_white, low = c_blue, na.value = c_black) +
  facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, size = 6)) +
  labs(x = "", y = "", fill = "log10RA")
  ggsave(p, file = paste0(figs, "heatmap_meanra_daignosis.png"),
  units = "in",
  device = "png",
  width = 4.5,
  height = 3.5
  )

df %>%
filter(as.character(mz) %in% target) %>%
# group_by(mz, genotype, batch) %>%
# summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(RA)), NA, log10(RA))),
  scale_ra = scale(RA)) %>%
  ggplot(aes(x = SampleID, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  # scale_fill_gradient2(midpoint = 0, high = c_red, mid = c_white, low = c_blue, na.value = c_black) +
  facet_grid(.~ batch + genotype, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    ) +
  labs(x = "", y = "", fill = "log10RA") +
  ggsave(file = paste0(figs, "heatmap_ra_daignosis.png"),
  units = "in",
  device = "png",
  width = 7.5,
  height = 4.5
  )

df %>%
filter(mz %in% unique(as.character(lfc_df_sig$peak_id))) %>%
group_by(mz, genotype, batch) %>%
summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(mean_ra)), NA, log10(mean_ra))),
  scale_ra = scale(mean_ra)) %>%
  ggplot(aes(x = genotype, y=mz, fill = log_ra)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = -5, high = c_red, mid = c_yellow, low = c_black, na.value = c_black) +
  # scale_fill_gradient2(midpoint = 0, high = c_red, mid = c_white, low = c_blue, na.value = c_black) +
  facet_grid(.~ batch, space = "free", scale = "free", switch = "x") +
  theme_RTN +
  theme(axis.text.y = element_blank()) +
  labs(x = "", y = "", fill = "log10RA") +
  ggsave(file = paste0(figs, "heatmap_meanra_sig.png"),
  units = "in",
  device = "png",
  width = 2.5,
  height = 4.5
  )

df %>%
filter(as.character(mz) %in% target) %>%
# group_by(mz, genotype, batch) %>%
# summarise(mean_ra = mean(RA)) %>%
mutate(log_ra = (ifelse(!is.finite(log10(RA)), NA, log10(RA))),
  scale_ra = scale(RA)) %>%
  ggplot(aes(x = genotype, y = log_ra)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8), 
    alpha = 0.9, shape = 21, colour = c_black, aes(fill = genotype)
    ) +
  geom_boxplot(outlier.alpha = 0, aes(colour = genotype), fill = c_white, alpha = 0.5, size = 0.2) +
  facet_wrap(
    mz ~ batch,  
    scale = "free", 
    switch = "both", nrow = length(target), ncol = length(unique(as.character(df$batch)))) +
  scale_fill_manual(values = c(geno$colours, "darkgrey")) +
  scale_colour_manual(values = c(geno$colours, "darkgrey")) +
  theme_RTN +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text.x = element_text(angle = 90, size = 10),
    strip.text.y = element_text(angle = 0, size = 10)
    ) +
  labs(x = "", y = "", fill = "", colour = "") +
  ggsave(file = paste0(figs, "boxplot_ra_daignosis.png"),
  units = "in",
  device = "png",
  width = 6,
  height = 7
  )

# }, mc.cores = 8)                     

# END OF SCRIPT
sessionInfo()
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
meta <- ms_data$meta %>% 
filter(genotype %in% c("Col", "pyk10", "cyp") & batch == "Krakow")
mat <- ms_data$mat_lognorm[, meta$given_id]

geno <- genotype[which(genotype$short %in% meta$genotype),]

design <- meta
row.names(design) <- design$given_id
design <- design %>% select(genotype, batch)

model <- model.matrix(~ 0 + genotype + batch, data = design)
colnames(model) <- str_replace(colnames(model), "genotype", "")
model <- model[,which(colSums(model) != 0)]

# Build contrast
# contrast.mat <- limma::makeContrasts(
# # Effect of genotype
#   # Col = (Col - MS),
#   # nai1 = (nai1 - Col),
#   pyk10 = (pyk10 - Col),
#   cyp = (cyp - Col),
#   # myb = (myb - Col),
# levels = colnames(model)
# )

# contrast.names <- attr(contrast.mat, "dimnames")$Contrasts
# n <- length(contrast.names)

# Stabilising variance - max intensity
# mat.voom <- limma::voom(mat, lib.size = colMeans(mat), design = model)
# fit_obj <- limma::lmFit(mat, design = model)
# fit_list <- lapply(contrast.names, 
# function(x) {
# fit <- contrasts.fit(fit_obj, contrast=contrast.mat[, which(contrast.names == x)])
# fit <- eBayes(fit)
# res <- limma::topTable(fit, adjust.method = p.adj.method, number = nrow(mat))
# return(res)
# })

# Alternate linear model
# fit_obj <- limma::lmFit(mat, design = model)
fit_list <- mclapply(row.names(mat), 
function(x) {
	
	temp_dat <- design %>%
	mutate(y = t(mat)[match(row.names(.), names(t(mat)[,x])),x])

	fit <- lm(y ~ 0 + genotype, data = temp_dat)
	anova <- aov(fit)
	
	stat <- broom::tidy(TukeyHSD(anova, p.adjust = "BH")) %>%
	separate(term, into = c("A", "B")) %>%
	filter(B == "Col") %>%
	select(A, estimate, adj.p.value) %>%
	data.frame(.)

	return(stat)
}, mc.cores = 8)
names(fit_list) <- row.names(mat)

fit_mat_diff <- do.call(rbind.data.frame, fit_list) %>%
mutate(given_id = str_replace(row.names(.), "\\.[1-2]$", "")) %>%
select(A, estimate, given_id) %>%
spread(key = "A", value = c("estimate")) %>%
data.frame(., row.names = 1)
colnames(fit_mat_diff) <- paste(colnames(fit_mat_diff), "logFC", sep = "_")

fit_mat_p <- do.call(rbind.data.frame, fit_list) %>%
mutate(given_id = str_replace(row.names(.), "\\.[1-2]$", "")) %>%
select(A, adj.p.value, given_id) %>%
spread(key = "A", value = c("adj.p.value")) %>%
data.frame(., row.names = 1)
colnames(fit_mat_p) <- paste(colnames(fit_mat_p), "PValue", sep = "_")

fit_mat <- cbind.data.frame(fit_mat_diff, fit_mat_p[row.names(fit_mat_diff),])
fit_mat_sig <- fit_mat[which(fit_mat$cyp_PValue <= 0.05 | fit_mat$pyk10_PValue <= 0.05),]

# Long format
# LogFC and PValue tables

#names(logFC_P.list) <- contrast.names
logFC_P.long <- fit_mat_sig %>%
add_column(given_id = row.names(.), .before = 1) %>%
gather(key = "key", value = "value", convert = FALSE, -given_id) %>%
separate(key, into = c("genotype", "vals"), sep = "_") %>%
spread(key = "vals", value = "value") %>%
mutate(FDR = p.adjust(PValue, "fdr")) %>%
mutate(significance = FDR <= 0.05 & abs(logFC) >= 1e-3) %>%
data.frame(.)

message(paste0("Number of DAMs at FDR ≤ ", sum(logFC_P.long$FDR <= 0.05)))
dim(logFC_P.long[(logFC_P.long$FDR <= 0.05),])

# Wide format
# LogFC and PValue tables
logFC_P.wide <- fit_mat %>%
mutate(
	significance = 
	(p.adjust(cyp_PValue, "fdr") <= 0.05 & p.adjust(pyk10_PValue, "fdr") <= 0.05)
)

stat_cc <- cor.test(fit_mat$cyp_logFC, fit_mat$pyk10_logFC)
til <- paste0("Pearson: ", round(stat_cc$estimate, 4),"; t-value: ", round(stat_cc$statistic, 4), "; p-value: ", round(stat_cc$p.value, 4), "; no of significant metabolites = ", sum(logFC_P.long$FDR <= 0.05))
fit_mat %>%
  mutate(sig = cyp_PValue <= alpha | pyk10_PValue <= alpha) %>%
  ggplot(aes(x = cyp_logFC, y = pyk10_logFC)) +
  ggtitle(paste0(til)) +
  geom_smooth(method = glm, se = FALSE, lwd = 0.75, colour = c_black) +
  geom_point(aes(fill = sig), 
              shape = 21, 
              size = 1,
              colour = c_black, 
              alpha = 0.75,
              na.rm = FALSE) +
  scale_fill_manual(values = c(`FALSE` = c_white, `TRUE` = c_dark_red), guide = FALSE) +
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
                                  ) +
  ggsave(file = paste0("./figures/Parwise_metabolome_cyp_pyk_logFC_krakow.png"), 
                    dpi = 600, 
                    units = img, 
                    device = "png", 
                    width = 5, 
                    height = 4, 
                    limitsize = T)


write.table(logFC_P.long, 
paste(stats, "/logFC_P_long_annotated.txt", sep = ""), 
sep = "\t")

write.table(logFC_P.wide, 
paste(stats, "/logFC_P_wide_annotated.txt", sep = ""), 
sep = "\t")

lfc_obj <- list(
logFC_P.long = logFC_P.long, 
logFC_P.wide = logFC_P.wide
)
save(list = "lfc_obj", file = paste("./data/lfc_DAMs_cologne.Rdata", sep = ""))


# write.table(subset, 
# paste(stats, "/abundance_annotated.txt", sep = ""), 
# sep = "\t")

# annotate_obj <- list(mat = subset, lfc = logFC_P.wide, lfc_df = logFC_P.long)
# save(list = "annotate_obj", file = paste("./data/annotated_objects.Rdata", sep = ""))


# Non parametric version --DUNTEST


# DAM_glm <- mclapply(names(ms_assay), function(x) {

#   # x = "mat_p"
#   mat <- ms_assay[[x]]

#   # Make factors
#   design <- ms_data$metadata
#   row.names(design) <- design$id
#   design$genotype <- factor(design$genotype, levels = geno$short)
#   design$replicate <- as.factor(design$replicate)

#   # Build design matrix
#   model <- model.matrix(~ 0 + genotype, data = design)
#   colnames(model) <- str_replace(colnames(model), "genotype", "")
#   model <- model[,which(colSums(model) != 0)]

#   # Build contrast
#   contrast.mat <- limma::makeContrasts(
#       # Effect of genotype
#       pyk10 = (pyk10 - Col),
#       cyp = (cyp - Col),
#   levels = colnames(model))

#   contrast.names <- attr(contrast.mat, "dimnames")$Contrasts
#   n <- length(contrast.names)

#   # Stabilising variance - max intensity
#   # mat.voom <- limma::voom(mat, lib.size = colMeans(mat), design = model)
#   fit_obj <- limma::lmFit(mat, design = model)
#   fit_list <- lapply(contrast.names, 
#                      function(x) {
#                          fit <- contrasts.fit(fit_obj, contrast=contrast.mat[, which(contrast.names == x)])
#                      fit <- eBayes(fit)
#                      res <- limma::topTable(fit, adjust.method = p.adj.method, number = nrow(mat))
#                      return(res)
#                          })
                          
#   names(fit_list) <- contrast.names

#   # Long format
#   # LogFC and PValue tables
#   logFC_P.list <- lapply(1:n, function(x) {
#       table <- as.data.frame(fit_list[[x]])
#       table <- table[,c("logFC", "AveExpr", "adj.P.Val")]
#       table$peak_id <- row.names(table)
#       table$contrast <- contrast.names[x]
#       #colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
#       return(table)
#   })
#   #names(logFC_P.list) <- contrast.names
#   logFC_P.long <- do.call(rbind.data.frame, logFC_P.list)
#   colnames(logFC_P.long) <- str_replace_all(colnames(logFC_P.long), "adj.P.Val", "PValue")

#   # Wide format
#   # LogFC and PValue tables
#   logFC_P.list <- lapply(1:n, function(x) {
#       table <- as.data.frame(fit_list[[x]])
#       table <- table[,c("logFC", "AveExpr", "adj.P.Val")]
#       #table$peak_id <- row.names(table)
#       #table$contrast <- contrast.names[x]
#       colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
#       return(table)
#   })

#   logFC_P.wide <- do.call(data.frame, logFC_P.list)
#   colnames(logFC_P.wide) <- str_replace_all(colnames(logFC_P.wide), "adj.P.Val", "PValue")

#   write.table(logFC_P.long, 
#               paste(stats, "/logFC_P_long_", x,".txt", sep = ""), 
#               sep = "\t")

#   write.table(logFC_P.wide, 
#               paste(stats, "/logFC_P_wide_", x,".txt", sep = ""), 
#               sep = "\t")                        

#   # Make everyting in one list
#   obj <- list(
#     logFC_P = logFC_P.wide,
#     LFC_df = logFC_P.long
#     )

#   return(obj)


# }, mc.cores = 8)

# names(DAM_glm) <- names(ms_assay)

# save(list = "DAM_glm",
#      file = paste("./data/DAMs_objects.Rdata", sep = ""))                     

# END OF SCRIPT
sessionInfo()
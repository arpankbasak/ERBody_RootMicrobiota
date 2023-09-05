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

design <- meta
row.names(design) <- design$given_id
dsign <- design %>% select(genotype, batch)

model <- model.matrix(~ 0 + genotype + batch, data = design)
colnames(model) <- str_replace(colnames(model), "genotype", "")
model <- model[,which(colSums(model) != 0)]

# Build contrast
contrast.mat <- limma::makeContrasts(
# Effect of genotype
  # Col = (Col - MS),
  nai1 = (nai1 - Col),
  pyk10 = (pyk10 - Col),
  cyp = (cyp - Col),
  myb = (myb - Col),
levels = colnames(model)
)

contrast.names <- attr(contrast.mat, "dimnames")$Contrasts
n <- length(contrast.names)

# Stabilising variance - max intensity
# mat.voom <- limma::voom(mat, lib.size = colMeans(mat), design = model)
fit_obj <- limma::lmFit(mat, design = model)
fit_list <- lapply(contrast.names, 
function(x) {
fit <- contrasts.fit(fit_obj, contrast=contrast.mat[, which(contrast.names == x)])
fit <- eBayes(fit)
res <- limma::topTable(fit, adjust.method = p.adj.method, number = nrow(mat))
return(res)
})

names(fit_list) <- contrast.names

# Long format
# LogFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
table <- as.data.frame(fit_list[[x]])
table <- table[,c("logFC", "AveExpr", "adj.P.Val")]
table$peak_id <- row.names(table)
table$contrast <- contrast.names[x]
#colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
return(table)
})
#names(logFC_P.list) <- contrast.names
logFC_P.long <- do.call(rbind.data.frame, logFC_P.list)
colnames(logFC_P.long) <- str_replace_all(colnames(logFC_P.long), "adj.P.Val", "PValue")
message(paste0("Number of DAMs at FDR ≤ ", sum(logFC_P.long$PValue <= alpha)))

# Wide format
# LogFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
table <- as.data.frame(fit_list[[x]])
table <- table[,c("logFC", "AveExpr", "adj.P.Val")]
#table$peak_id <- row.names(table)
#table$contrast <- contrast.names[x]
colnames(table) <- paste(contrast.names[x], colnames(table), sep="_")
return(table)
})

logFC_P.wide <- do.call(data.frame, logFC_P.list)
colnames(logFC_P.wide) <- str_replace_all(colnames(logFC_P.wide), "adj.P.Val", "PValue")

# Annotate
# logFC_P.wide$annotation <- temp$annotate[match(row.names(logFC_P.wide), temp$pid)]
# logFC_P.long$annotation <- temp$annotate[match(logFC_P.long$peak_id, temp$pid)]
# subset$annotation <- temp$annotate[match(row.names(subset), temp$pid)]

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
save(list = "lfc_obj", file = paste("./data/lfc_DAMs.Rdata", sep = ""))


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
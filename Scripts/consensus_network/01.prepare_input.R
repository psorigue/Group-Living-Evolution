library(DESeq2)
library(matrixStats)

# -----------------------------
# 0. Paths and parameters
# -----------------------------
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/04.counts/metadata_datasetR.txt"

region <- "TL"
#species_list <- c("Ornatipinnis", "Meeli", "Multifasciatus")
species_list <- c("Meeli", "Multifasciatus")

median_counts <- 10
# Do only with 3 populations
var_quantile <- 0.2 # Percentile to reduce low-variability genes  

path_cons <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/06.WGCNA/consensus/", region, "/RData/")

samples_to_remove <- c("Or49DE")

# -----------------------------
# 1. Load data
# -----------------------------
cts <- read.csv(file_cts, sep = "\t", header = TRUE, row.names = 1)

meta <- read.csv(file_metadata, sep = "\t", header = TRUE)
rownames(meta) <- meta$id

# -----------------------------
# 2. Split data by species
# -----------------------------
cts_list <- list()
meta_list <- list()

for (spp in species_list) {
  
  meta_fil <- meta[meta$region == region & meta$species == spp, ]
  meta_fil <- meta_fil[!rownames(meta_fil) %in% samples_to_remove, ]
  
  cts_fil <- cts[, rownames(meta_fil)]
  
  cts_list[[spp]] <- cts_fil
  meta_list[[spp]] <- meta_fil
}

# Save meta_list for later use
save(meta_list, file = paste0(path_cons, "2spp_meta_list.RData"))


# -----------------------------
# 3. Keep only shared genes
# -----------------------------
common_genes <- Reduce(intersect, lapply(cts_list, rownames))
cts_list <- lapply(cts_list, function(x) x[common_genes, ])

# -----------------------------
# 4. Expression filtering (per species, then intersect)
# -----------------------------
keep_list <- lapply(cts_list, function(x) {
  rowMedians(as.matrix(x)) >= median_counts
})

keep_genes <- Reduce(intersect, lapply(names(cts_list), function(spp) {
  rownames(cts_list[[spp]])[keep_list[[spp]]]
}))

cts_list <- lapply(cts_list, function(x) x[keep_genes, ])

# -----------------------------
# 5. Normalize each species separately
# -----------------------------
expr_list <- list()

for (spp in species_list) {
  
  cat("Processing:", spp, "\n")
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts_list[[spp]],
    colData = meta_list[[spp]],
    design = ~ 1
  )
  
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)
  
  expr <- t(assay(vsd))  # samples x genes
  
  expr_list[[spp]] <- expr
}

# -----------------------------
# 6. CONSENSUS VARIANCE FILTERING 
# -----------------------------

# 1. Compute variance per species
var_list <- lapply(expr_list, function(expr) {
  colVars(expr)
})

# 2. Combine variances across species
var_min <- Reduce(pmin, var_list)

# 3. Define threshold
var_threshold <- quantile(var_min, probs = var_quantile)

# 4. Keep genes variable in ALL species
keep_genes_var <- names(var_min)[var_min > var_threshold]

cat("Genes before variance filter:", ncol(expr_list[[1]]), "\n")
cat("Genes after variance filter:", length(keep_genes_var), "\n")

# 5. Apply to all species (preserves identical gene set)
expr_list <- lapply(expr_list, function(expr) {
  expr[, keep_genes_var]
})

# -----------------------------
# 7. Save individual matrices
# -----------------------------
for (spp in species_list) {
  expr <- expr_list[[spp]]
  save(expr, file = paste0(path_cons, "2spp_norm_cts_", spp, ".RData"))
}

# -----------------------------
# 8. Build multiExpr 
# -----------------------------
multiExpr <- lapply(expr_list, function(x) list(data = x))

# Save for WGCNA
save(multiExpr, file = paste0(path_cons, "/2spp_multiExpr_input.RData"))



meta <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/10.phydget/metadata_datasetR.txt",
                 sep = "\t")
meta$group_phenotype <- paste(meta$species, meta$phenotype, sep = "_")
cts <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/10.phydget/counts_all_samples_datasetR.txt",
                sep = "\t")

# Remove outlier
meta_fil <- meta[meta$id != "Or49DE",]
cts$Or49DE <- NULL

# Remove "gene-" string
cts$X <- gsub(pattern = "gene-", replacement = "", cts$X)

# Split into brain regions
# 1. Find all column positions that contain the pattern "TL"
tl_indices <- grep("TL", colnames(cts))
de_indices <- grep("DE", colnames(cts))

# 2. Combine the 1st column index with the "TL" indices
# using unique() prevents duplication if the 1st column also happens to contain "TL"
selected_columns_tl <- unique(c(1, tl_indices))
selected_columns_de <- unique(c(1, de_indices))

# 3. Subset your count table
tl_filtered_by_col_counts <- cts[, selected_columns_tl]
de_filtered_by_col_counts <- cts[, selected_columns_de]


# For TL region
rownames(tl_filtered_by_col_counts) <- tl_filtered_by_col_counts[, 1]
tl_counts_numeric <- tl_filtered_by_col_counts[, -1]

# For DE region
rownames(de_filtered_by_col_counts) <- de_filtered_by_col_counts[, 1]
de_counts_numeric <- de_filtered_by_col_counts[, -1]

# ==========================================
# 2. METADATA SPLITTING & ALIGNMENT
# ==========================================
# Synchronize and extract metadata for TL samples
meta_tl <- meta_fil[match(colnames(tl_counts_numeric), meta_fil$id), ]
groups_tl <- factor(meta_tl$group_phenotype) # Extracts your real Solitary/Group factors

# Synchronize and extract metadata for DE samples
meta_de <- meta_fil[match(colnames(de_counts_numeric), meta_fil$id), ]
groups_de <- factor(meta_de$group_phenotype)


# ==========================================
# 3. LOW-COUNT FILTERING VIA EDGER
# ==========================================
library(edgeR)

# --- Process TL Region ---
dge_tl <- DGEList(counts = tl_counts_numeric)

# The filter automatically detects smallest group size (4) using 'groups_tl'
keep_genes_tl <- filterByExpr(
  dge_tl, 
  group = groups_tl,
  min.count = 10,       # Requires ~10 raw reads per sample (adjusted for depth)
  min.total.count = 15  # Minimum total reads across all TL samples to keep a gene
)
tl_counts_filtered <- tl_counts_numeric[keep_genes_tl, ]


# --- Process DE Region ---
dge_de <- DGEList(counts = de_counts_numeric)

keep_genes_de <- filterByExpr(
  dge_de, 
  group = groups_de,
  min.count = 10,       
  min.total.count = 15  
)
de_counts_filtered <- de_counts_numeric[keep_genes_de, ]


# ==========================================
# 4. FINAL CONSOLE SANITY CHECK
# ==========================================
# Check how many genes survived the filtering step in each region
print(paste("Remaining genes in TL matrix:", nrow(tl_counts_filtered)))
print(paste("Remaining genes in DE matrix:", nrow(de_counts_filtered)))


# Save
write.table(tl_counts_filtered, "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/10.phydget/TL_cts.tsv",
            sep = "\t", col.names = NA, row.names = T, quote = F)
write.table(de_counts_filtered, "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/10.phydget/DE_cts.tsv",
            sep = "\t", col.names = NA, row.names = T, quote = F)

# This script performs differential gene expression analysis using DESeq2. 
# It reads in count data and metadata, filters out low-count genes, 
# and fits a DESeq model to identify differentially expressed genes 
# between specified groups. The results are then merged with gene 
# descriptions and saved to a text file for downstream analysis.

{
  library(DESeq2) # version 1.50.2
  library(dplyr) # version 1.1.4
  library(clusterProfiler) # version 4.18.4
  library(AnnotationHub) # version 4.0.0
}

# Set region
region <- "TL"

# Set paths and read files.
home <- path.expand("~")
file_cts <- file.path(home, "03.Mapping and Counts", "counts_all_samples.txt") 
file_metadata <- file.path(home, "03.Mapping and Counts", "metadata.txt")
path_deg <- file.path(home, "04.Differential Gene Expression", region)
file_gene_info <- file.path(home, "Ref_genome", "gene_info_full.txt")


# INDEX
# 1. Prepare data and fit DESeq model
# 2. Differential gene expression analysis
# 3. GO and KEGG enrichment analysis
  # 3.1. Obtain gene ids for enrichment
  # 3.2. GO enrichment
  # 3.3. KEGG enrichment


# 1. Prepare data and fit DESeq model
#------------------------------------

# Read counts and metadata
cts  <- read.csv(file_cts, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file_metadata, sep = "\t", header = T, row.names = 1)

# Remove outlier
samples_to_remove <- "Or49DE"
meta <- meta[!rownames(meta) %in% samples_to_remove, ]
cts  <- cts[, rownames(meta)]

# Ensure metadata matches count matrix
stopifnot(all(rownames(meta) == colnames(cts)))

# Filter metadata and counts
meta_reg <- meta[meta$region == region, ]
cts_reg  <- cts[, rownames(meta_reg)]

# Create combined group
meta_reg$group <- factor(paste(meta_reg$species, meta_reg$phenotype, sep = "_"))
meta_reg$group <- droplevels(meta_reg$group)

# Check what groups exist
print(table(meta_reg$group))

# Build DDS
dds <- DESeqDataSetFromMatrix(
  countData = cts_reg,
  colData   = meta_reg,
  design    = ~ group
)

# Remove low-count genes
smallestGroupSize <- min(table(dds$species)) # Recommended by DESeq2 developers
keep_genes <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep_genes, ]

# Fit DESeq model
dds_obj <- DESeq(dds_filtered)
# This model creates five groups according to spp and phenotype: 
# Ornatipinnis_S, Meeli_S, Meeli_GL, Multifasciatus_S, Multifasciatus_GL.


# 2. DEG analysis
#-----------------

# Define species and groups to compare
ref_spp <- "Ornatipinnis_S" # Species_Phenotype "starting point" for fold change comparison
dir_spp <- "Meeli_S" # Species_Phenotype direction foldChange
comp_name <- "OrS_vs_MeS" # Name for output files

# DEG analysis
res <- DESeq2::results(
  dds_obj,
  contrast = c("group", dir_spp, ref_spp) 
)

# Order results by adjusted p-value
res_ord <- res[order(res$padj),]

# Join gene description in the results table
gene_info <- read.csv(file_gene_info, sep = "\t")

merged_df <- as.data.frame(res_ord) %>%
  mutate(gene = sub("^gene-", "", rownames(res_ord))) %>%
  left_join(
    gene_info[, c("symbol_gtf", "name_datasets")],
    by = c("gene" = "symbol_gtf")
  ) %>%
  rename(gene_description = name_datasets) %>%
  select(gene, gene_description, everything())

# Write results
write.table(merged_df, 
            file = file.path(path_deg, paste0("DEG_", comp_name, ".txt")),
            sep = "\t",
            quote = F,
            row.names = F)


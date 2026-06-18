{
  library(DESeq2)
  library(dplyr)
  #library(apeglm)
  library(ashr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
}



# 1. Read and create DDS object
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/metadata_datasetR.txt"
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/05.DEG/"

file_cts <- "../../to_work/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "../../to_work/04.counts/metadata_datasetR.txt"
path_deg <- "../../to_work/05.DEG/"

# 2. Read counts and metadata
cts  <- read.csv(file_cts, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file_metadata, sep = "\t", header = T, row.names = 1)

# 3. Remove outlier BEFORE creating DDS
samples_to_remove <- "Or49DE"
meta <- meta[!rownames(meta) %in% samples_to_remove, ]
cts  <- cts[, rownames(meta)]

# Ensure metadata matches count matrix
stopifnot(all(rownames(meta) == colnames(cts)))


# 3. Create initial DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = meta,
  design    = ~ phenotype 
)


# 4. Filter and subset
# Define comparison parameters
species_sl  <- "Ornatipinnis"  # solitary species
species_gl  <- "Meeli"         # group species
comp_name   <- "OrS_vs_MeGL"
region <- "DE"

# Remove low-count genes
smallestGroupSize <- 4
keep_genes <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep_genes, ]


# Subset samples by species, phenotype, and region
meta_subset <- meta[
  ((meta$species == species_gl & meta$phenotype == "GL") |
     (meta$species == species_sl & meta$phenotype == "S")) &
    meta$region == region, 
]

dds_subset <- dds_filtered[, rownames(meta_subset)]


# 5. Relevel phenotype (control = "S")
dds_subset$phenotype <- relevel(dds_subset$phenotype, ref = "S")

# 6. Save pre-fitted DDS
save(dds_subset, file = paste0(path_deg, region, "/DDS_files/DDSdata_", comp_name, ".RData"))

# 7. Fit DESeq model
dds_obj <- DESeq(dds_subset)

# 8. Save fitted DDS object
save(dds_obj, file = paste0(path_deg, region, "/DDS_files/DDSobj_", comp_name, ".RData"))
  
  
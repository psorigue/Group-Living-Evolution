# This script performs PCA and outlier detection on gene counts prior to Differential Gene Expression analysis.

# Load libraries
{
  library(DESeq2) # version 1.50.2
  library(dplyr) # version 1.1.4
  library(ggplot2) # version 4.0.1
}


# Set region for analysis
region <- "DE"

# INDEX
# 1. PCA all samples by region
  # 1.1. Prepare data
  # 1.2. Make and plot PCA
# 2. Outlier detection

# 1. PCA all samples by region
#-----------------------------

  # 1.1. Prepare data
  #-------------------

# Set paths and read files.
home <- path.expand("~")
file_cts <- file.path(home, "03.Mapping_and_Counts", "counts_all_samples.txt") 
file_metadata <- file.path(home, "03.Mapping_and_Counts", "metadata.txt")
path_deg <- file.path(home, "04.Differential_Gene_Expression")

# Read files. The datasets are already curated 
cts <- read.csv(file = file_cts, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file = file_metadata, sep = "\t", header = T, row.names = 1)

# Add group variable for PCA
meta$group <- paste(meta$phenotype, meta$species, sep = "_")

# Filter region
meta_reg <- meta[meta$region == region, ]

# Remove outliers from metadata -> This step only after first outlier detection,
# which is done on the unfiltered dataset. Skip this step for the first PCA, 
# and then remove outliers for the second PCA. The outliers are identified
# in the "2. Outlier detection" section below.
outliers <- "Or49DE"
meta_reg <- meta_reg[!rownames(meta_reg) %in% outliers, ]

# Subset counts to match filtered metadata
cts_reg <- cts[, rownames(meta_reg)]

# Arrange order of meta and cts to match, and check
cts_reg <- cts_reg[, rownames(meta_reg)]
all(rownames(meta_reg) == colnames(cts_reg))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts_reg,
                              colData = meta_reg,
                              design = ~ phenotype + species) 

# Prefilter low counts following developer’s recommendation
smallestGroupSize <- min(table(dds$species)) # Recommended by developer
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Create dds object
dds_obj <- DESeq(dds)

# Variance-stabilized transformation
vsd <- vst(dds_obj)

  # 1.2. Make and plot PCA
  #-----------------------
  
file_plot <- file.path(path_deg, "PCA_outliers", paste0(region, "_unfiltered.pdf"))

# Compute PCA scores
pc_scores <- prcomp(t(assay(vsd)))

# Add group info to PCA scores
sample_info <- as.data.frame(colData(vsd))  # metadata
pca_data <- cbind(pc_scores$x, sample_info)

# Variance explained by each PC
percentVar <- (prcomp(t(assay(vsd)))$sdev^2) / sum(prcomp(t(assay(vsd)))$sdev^2)

# Define which phenotype gets filled ellipses
filled_pheno <- "GL"

# Plot 2 PCAs with ggplot2, adding filled ellipses for the specified phenotype
# and outline ellipses for all groups
p <- ggplot(pca_data, aes(PC1, PC2,
                          color = species,
                          shape = phenotype)) +
  
  # Points
  geom_point(size = 3) +
  
  # Filled ellipses
  stat_ellipse(
    data = subset(pca_data, phenotype == filled_pheno),
    aes(group = interaction(species, phenotype),
        fill = species),
    geom = "polygon",
    type = "t",
    level = 0.95,
    alpha = 0.3,
    color = NA           
  ) +
  
  # Outline ellipses
  stat_ellipse(
    data = subset(pca_data),
    aes(group = interaction(species, phenotype)),
    geom = "path",
    type = "t",
    level = 0.95,
    linewidth = 1
  ) +
  
  # Manual colors (ONLY ONCE)
  scale_color_manual(values = c(
    "Meeli" = "#104E8B",
    "Ornatipinnis" = "#CD1076",
    "Multifasciatus" = "#8B4500"
  )) +
  
  scale_fill_manual(values = c(
    "Meeli" = "#104E8B",
    "Ornatipinnis" = "#CD1076",
    "Multifasciatus" = "#8B4500"
  )) +
  
  xlab(paste0("PC1: ", round(100 * percentVar[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * percentVar[2], 1), "% variance")) +
  
  ggtitle("Diencephalon") +
  
  # Legends
  labs(
    color = "Species",
    fill = NULL,
    shape = "Phenotype"
  ) +
  
  theme_minimal()

p

ggsave(plot = p, filename = file_plot, device = "pdf", width = 5, height = 5)



# 2. Outlier detection
#---------------------
file_out <- file.path(path_deg, "PCA_outliers", paste0(region, "_outlier_zscore.txt"))

# Run PCA on variance-stabilized data (samples in rows, genes in columns)
pc <- prcomp(t(assay(vsd)))

# Calculate proportion of variance explained by each PC
percentVar <- pc$sdev^2 / sum(pc$sdev^2)

# Select number of PCs to keep (e.g., first 4)
k <- 4
pc_scores <- pc$x[, 1:k]

# Add group info
pca_data <- as.data.frame(pc_scores)
pca_data$group <- meta_reg$group
pca_data$Sample <- rownames(pca_data)

# Group-wise z-score based outlier detection
results <- pca_data %>%
  group_by(group) %>%
  do({
    df <- .
    
    # Extract PCs as matrix
    pcs <- as.matrix(df[, paste0("PC", 1:k)])
    
    # Standardize PCs (z-scores)
    pcs_z <- scale(pcs)
    
    # Euclidean distance from group centroid (origin in z-space)
    dist_z <- sqrt(rowSums(pcs_z^2))
    
    # Flag outliers: distance > 3 SD
    outlier_flag <- dist_z > 3
    
    data.frame(
      Sample = df$Sample,
      Zdist = dist_z,
      BeyondThreshold = outlier_flag
    )
  }) %>%
  ungroup()

# Save results
write.table(results, file = file_out, sep = "\t", quote = FALSE, 
            col.names = TRUE, row.names = FALSE)

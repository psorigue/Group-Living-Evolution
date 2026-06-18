############################################################
# WGCNA Consensus Network Analysis
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

region <- "TL"

path_cons <- paste0(
  "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/06.WGCNA/consensus/",  region, "/RData/"
)

# Load metadata
load(paste0(path_cons, "2spp_meta_list.RData"))


# Load libraries
library(WGCNA)
library(ggplot2)
library(gridExtra)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# -----------------------------
# 2. Load data
# -----------------------------

load(paste0(path_cons, "/2spp_multiExpr_input.RData")) # loads variable called 'multiExpr'
# multiExpr = list of expression datasets
# Load metadata
load(paste0(path_cons, "2spp_meta_list.RData"))

str(multiExpr)

# -----------------------------
# 3. Quality Control
# -----------------------------

# Multi-set QC (correct function for multiExpr)
gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)

if (!gsg$allOK) {
  message("Removing bad genes/samples...")
  
  for (set in seq_along(multiExpr)) {
    multiExpr[[set]]$data <- multiExpr[[set]]$data[
      gsg$goodSamples[[set]],
      gsg$goodGenes
    ]
  }
}

# Check that gene names match across datasets
stopifnot(
  all.equal(
    colnames(multiExpr[[1]]$data),
    colnames(multiExpr[[2]]$data),
    colnames(multiExpr[[3]]$data)
  )
)

# -----------------------------
# 4. Sample Clustering (Outlier detection)
# -----------------------------

par(mfrow = c(1, length(multiExpr)))

for (i in seq_along(multiExpr)) {
  sampleTree <- hclust(
    dist(multiExpr[[i]]$data),
    method = "average"
  )
  
  plot(
    sampleTree,
    main = names(multiExpr)[i],
    xlab = "",
    sub = ""
  )
}

# -----------------------------
# 5. Soft-threshold selection
# -----------------------------

powers <- 1:20

# Check for all 3 datasets and decide a consensus power
sft = pickSoftThreshold(
  multiExpr[[1]]$data,
  powerVector = powers,
  corFnc = "bicor", networkType = "signed hybrid",
  corOptions = list(use = 'p', maxPOutliers = 0.1) # Default robust settings
)
wgcna_powers_plot(sft$fitIndices)

# Choose power manually
softPower <- 7
mch <- 0.4

# -----------------------------
# 6. Consensus Module Detection
# -----------------------------

consensus_MuMe <- blockwiseConsensusModules(
  multiExpr,
  corType="bicor",
  maxPOutliers = 0.1,
  power = softPower,
  networkType = "signed hybrid",
  
  minModuleSize = 30,
  deepSplit = 1,
  mergeCutHeight = mch,
  
  pamRespectsDendro = FALSE,
  maxBlockSize = 20000,
  
  TOMType = "signed", 
  TOMDenom = "min",
  numericLabels=FALSE,
  saveTOMs=F,
  
  nThreads=23,
  verbose = 3
)

# Save network
save(consensus_MuMe, file = paste0(path_cons, "consensus_MuMe_network.RData"))
load(file = paste0(path_cons, "consensus_MuMe_network.RData"))

net <- consensus_MuMe

net$multiMEs
length(table(net$colors))



# -----------------------------
# 7. Module Visualization
# -----------------------------
# How many modules?
length(unique(net$colors))
# How many genes per module?
moduleSizes <- as.data.frame(table(net$colors))
colnames(moduleSizes) <- c("Module", "GeneCount")
moduleSizes[order(moduleSizes$GeneCount),] 

moduleColors <- net$colors
moduleColors <- as.character(moduleColors)
table(moduleColors)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors,
  "Consensus modules",
  dendroLabels = FALSE,
  hang = 0.03
)


head(moduleColors)
# -----------------------------
# 8. Module Eigengenes
# -----------------------------

MEs <- net$multiMEs

# Inspect eigengenes per dataset
names(MEs)
# Example:
# MEs[[1]], MEs[[2]], etc.

# -----------------------------
# 9. Module Preservation
# -----------------------------
multiColor <- lapply(multiExpr, function(x) moduleColors)

# Run module preservation analysis. Tests whether modules defined in one dataset (reference) are preserved in others.
mp <- modulePreservation(
  multiData = multiExpr,        # list of expression datasets
  multiColor = multiColor,      # module assignments for each dataset
  referenceNetworks = c(1,2),        
  networkType = "signed hybrid",
  nPermutations = 500,     # IF IT WORKS, RE-RUN WITH 500     # permutation test (increase for publication-quality)
  verbose = 3
)

# Save results to avoid recomputation (this step is computationally expensive)
save(mp, file = paste0(path_cons, "/RData/mp_MuMe.RData"))
load(file = paste0(path_cons, "/RData/mp_MuMe.RData")) # loads variable 'mp'

# Access preservation statistics (Z scores)
mp$preservation$Z

# Example: Reference 2, Test 3
z_facultative_comp <- mp$preservation$Z[[1]][[2]]

head(z_facultative_comp)

z_facultative_comp$Zsummary
# Interpretation:
# Z-summary > 10  → highly preserved
# Z-summary 2–10 → moderately preserved
# Z-summary < 2  → not preserved



# CREATE TABLE WITH ALL Z-SUMMARY VALUES
# 1. Initialize a list to hold the individual data frames
all_results <- list()

# 2. Define your reference and test indices (1, 2, 3)
refs <- c(1, 2)
tests <- c(1, 2)

# 3. Loop through each reference-test pair
for (r in refs) {
  for (t in tests) {
    if (r != t) {
      # 1. Extract data frames
      z_pres_table   <- as.data.frame(mp$preservation$Z[[r]][[t]])
      obs_pres_table <- as.data.frame(mp$preservation$observed[[r]][[t]])
      z_qual_table   <- as.data.frame(mp$quality$Z[[r]][[t]])
      
      # 2. Identify the number of modules (excluding 'gold')
      # This makes the percentile accurate even if module counts change
      n_mods <- sum(rownames(z_pres_table) != "gold")
      
      # 3. Create the data frame
      z_table <- data.frame(
        Ref_Species     = names(multiExpr)[r],
        Test_Species    = names(multiExpr)[t],
        Module_Name     = rownames(z_pres_table),
        Module_Size     = z_pres_table[, "moduleSize"],
        Zsummary_Qual   = z_qual_table[, "Zsummary.qual"],
        Zsummary_Pres   = z_pres_table[, "Zsummary.pres"],
        MedianRank_Pres = obs_pres_table[, "medianRank.pres"]
      )
      
      # 4. Remove 'gold' module BEFORE calculating percentile
      z_table <- z_table[z_table$Module_Name != "gold", ]
      
      # 5. Add percentile correctly using the current table's rank
      z_table$Rank_Percentile <- (z_table$MedianRank_Pres / n_mods) * 100
      
      all_results[[paste0(r, "_vs_", t)]] <- z_table
    }
  }
}

# 6. Bind and Clean up
master_table <- do.call(rbind, all_results)
rownames(master_table) <- NULL

# 7. Final Column Selection
master_table <- master_table[, c("Ref_Species", "Test_Species", "Module_Name", "Module_Size",
                                 "Zsummary_Pres", "Zsummary_Qual", "MedianRank_Pres", "Rank_Percentile")]

# 8. Save
write.table(master_table, paste0(path_cons, "consensus_MuMe_scores.txt"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# VISUALIZATION
# Define the reference and test set indices
# For example: Ref 1 (Solitary) vs Test 2 (Social A)
ref = 1
test = 2

# Extract the data for the plot
plot_data = mp$preservation$Z[[ref]][[test]]
mod_sizes = plot_data$moduleSize
z_scores = plot_data$Zsummary.pres

# Create the plot
plot(mod_sizes, z_scores, 
     pch = 19, col = rownames(plot_data), 
     main = paste("Preservation:", names(multiExpr)[ref], "vs", names(multiExpr)[test]),
     xlab = "Module Size", ylab = "Z-summary",
     ylim = c(0, max(z_scores) + 2))

# Add significance thresholds
abline(h = 2, col = "blue", lty = 2)   # Weak to moderate preservation threshold
abline(h = 10, col = "red", lty = 2)   # Strong preservation threshold

# Add module labels to the points
text(mod_sizes, z_scores, labels = rownames(plot_data), cex = 0.7, pos = 3)



# -----------------------------
# 10. Module–Trait Relationships (Optional)
# -----------------------------

# Make metadata
meta <- read.csv(file = file_metadata, sep = "\t", header = T)
rownames(meta) <- meta$id
#Filter out samples
samples_to_remove <- "Or49DE"
meta_fil <- meta[!row.names(meta) %in% samples_to_remove, ]
# Filter region
meta_reg <- meta_fil[meta_fil$region == region, ]

# 1. Prepare Trait Data correctly
datTraits <- list()

for (s in seq_along(multiExpr)) {
  
  spp <- names(multiExpr)[s]
  
  expr_samples <- rownames(multiExpr[[s]]$data)
  meta_spp <- meta_list[[spp]]
  
  # Ensure perfect alignment
  meta_spp <- meta_spp[expr_samples, ]
  
  stopifnot(all(rownames(meta_spp) == expr_samples))
  
  trait <- as.numeric(as.factor(meta_spp$phenotype))
  
  datTraits[[s]] <- list(
    data = data.frame(Target_Trait = trait)
  )
}

# Sanity checks
all(rownames(multiExpr[[2]]$data) == rownames(meta_list[[names(multiExpr)[2]]]))

table(is.na(datTraits[[2]]$data))


# Correlations
# Initialize lists to store results
all_cor_results <- list()
MEs <- net$multiMEs

for (s in refs) {

  set_name <- names(multiExpr)[s]
  
  # 1. Calculate raw Correlation and P-values
  # Ensure your trait is numeric: as.numeric(datTraits[[s]]$data[,1])
  cor_vector <- as.vector(cor(MEs[[s]]$data, datTraits[[s]]$data, use = "p"))
  raw_p_values <- as.vector(corPvalueStudent(cor_vector, nrow(multiExpr[[s]]$data)))
  
  # 2. Apply Multiple Testing Correction (FDR / Benjamini-Hochberg)
  # This adjusts for the 163 tests performed in this species
  adj_p_values <- p.adjust(raw_p_values, method = "BH")
  
  # 3. Create the data frame for this species
  cor_df <- data.frame(
    Species = set_name,
    Module_Name = substring(colnames(MEs[[s]]$data), 3), # Removes "ME"
    Correlation = cor_vector,
    P_value_Raw = raw_p_values,
    P_value_Adj = adj_p_values
  )
  
  # 4. Add a significance label based on the adjusted p-value
  cor_df$Is_Significant_Adj <- ifelse(cor_df$P_value_Adj < 0.05, "YES", "no")
  
  all_cor_results[[s]] <- cor_df
}

# 5. Combine and Merge with your Preservation Table
master_cor_table <- do.call(rbind, all_cor_results)

final_summary <- merge(master_table, master_cor_table, 
                       by.x = c("Ref_Species", "Module_Name"), 
                       by.y = c("Species", "Module_Name"))

# 6. Save the final robust table
write.table(final_summary, "Final_Analysis_FDR_Corrected.txt", 
            sep="\t", row.names=F, quote=F)


# 1. Merge all MEs into one big matrix
# We bind the 'data' part of each species in the multiMEs list
all_MEs_merged <- do.call(rbind, lapply(MEs, function(x) x$data))

# 2. Merge all Traits into one big vector
# We do the same for the traits list we created earlier
all_traits_merged <- do.call(rbind, lapply(datTraits, function(x) x$data))

# 3. Calculate the Global Correlation
# This calculates one R and one P-value for each of the 163 modules
global_cor <- cor(all_MEs_merged, all_traits_merged, use = "p")
global_p   <- corPvalueStudent(global_cor, nrow(all_MEs_merged))

# 4. Create the Global Summary Table
global_results <- data.frame(
  Module_Name = substring(rownames(global_cor), 3),
  Global_Correlation = as.vector(global_cor),
  Global_P_value_Raw = as.vector(global_p),
  Global_P_value_Adj = p.adjust(as.vector(global_p), method = "BH")
)

# 5. Filter for Significant Global Modules
global_winners <- global_results[global_results$Global_P_value_Adj < 0.05, ]







# Cross-species trait

names(multiExpr)
sapply(multiExpr, function(x) nrow(x$data))
# Extract eigengenes for each dataset
MEs_Obli <- net$multiMEs[[1]]$data
MEs_Fac1 <- net$multiMEs[[1]]$data
MEs_Fac2 <- net$multiMEs[[2]]$data

spp1 <- "Meeli"
spp2 <- "Multifasciatus"

meta <- read.csv(file = file_metadata, sep = "\t", header = T)
rownames(meta) <- meta$id
samples_to_remove <- "Or49DE"
meta_fil <- meta[!row.names(meta) %in% samples_to_remove, ]
meta_fil$phenotype <- ifelse(meta_fil$phenotype == "GL", 1, 0)
meta_reg1 <- meta_fil[meta_fil$region == region & meta_fil$species == spp1, ]
meta_reg2 <- meta_fil[meta_fil$region == region & meta_fil$species == spp2, ]

traits_Fac1 <- meta_reg1$phenotype
traits_Fac2 <- meta_reg2$phenotype

# 1. Correlate Fac1 MEs with Fac1 Phenotype (0=Solitary, 1=Group)
cor_Fac1 = cor(MEs_Fac1, traits_Fac1, use = "p")
cor_Fac2 = cor(MEs_Fac2, traits_Fac2, use = "p")

# Compute significance
moduleTraitPvalue_Fac1 = as.data.frame(
  corPvalueStudent(cor_Fac1, nSamples = nrow(MEs_Fac1))
)

moduleTraitPvalue_Fac2 = as.data.frame(
  corPvalueStudent(cor_Fac2, nSamples = nrow(MEs_Fac2))
)

# Find modules significant in BOTH (e.g., p < 0.05) -> TRY CORRECTING FDR
plasticity_modules = rownames(moduleTraitPvalue_Fac1)[
  moduleTraitPvalue_Fac1$V1 < 0.05 & 
    moduleTraitPvalue_Fac2$V1 < 0.05 & 
    sign(corPvalueStudent(cor_Fac1, nSamples = nrow(MEs_Fac1))) == sign(corPvalueStudent(cor_Fac2, nSamples = nrow(MEs_Fac2)))
]

print(plasticity_modules)



# Check a module
moduleColors <- net$colors

# Gene names from the chosen dataset
geneNames <- colnames(multiExpr[[1]]$data)
# Genes in mistyrose module
module <- "firebrick4"
moduleGenes <- geneNames[moduleColors == module]
moduleGenes <- gsub(moduleGenes, pattern = "gene-", replacement = "")



# Link trait relevance with preservation
# Extract preservation stats:
stats_obli <- mp$preservation$Z$ref.Ornatipinnis$spec.set1

# Get Z-summary only for your modules of interest
preservation_in_obligate <- stats_obli[
  plasticity_modules,
  "Zsummary.pres",
  drop = FALSE
]

# Sort modules from least to most preserved
preservation_in_obligate <- preservation_in_obligate[
  order(preservacio_en_obligat$Zsummary.pres),
  , drop = FALSE
]

print(preservation_in_obligate)


# Visualization
combined_cors = cbind(cor_Fac1, cor_Fac2)

labeledHeatmap(
  Matrix = combined_cors, 
  xLabels = c("Fac1_Sociality", "Fac2_Sociality"),
  yLabels = names(MEs_Fac1),
  colors = blueWhiteRed(50)
)





# Get the genes of a module
cyan, darkorange, darkred lightgreen

# Module assignment vector (one color per gene)
moduleColors <- net$colors

# Gene names (same order as moduleColors)
geneNames <- colnames(multiExpr[[1]]$data)

# Extract genes belonging to the cyan module
genes <- geneNames[moduleColors == "red"]
genes <- gsub(genes, pattern = "gene-", replacement = "")

# Inspect
head(cyan_genes)
length(cyan_genes)





# -----------------------------
# 11. Helper Function
# -----------------------------

wgcna_powers_plot <- function(sft_indices){
  
  # Scale-free topology fit
  p1 <- ggplot(sft_indices,
               aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.05) +
    geom_hline(yintercept = 0.8, color = "red") +
    labs(
      x = "Power",
      y = "Scale-free topology fit (R²)"
    ) +
    theme_classic()
  
  # Mean connectivity
  p2 <- ggplot(sft_indices,
               aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.7) +
    labs(
      x = "Power",
      y = "Mean connectivity"
    ) +
    theme_classic()
  
  grid.arrange(p1, p2, nrow = 2)
  
  return(sft_indices)
}

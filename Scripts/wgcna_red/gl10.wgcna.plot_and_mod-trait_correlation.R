{
  library(stringr)
  library(WGCNA)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
  library(dplyr)
  library(tidyverse)
  library(tidyr)
}

# Path
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/metadata_datasetR.txt"

#Region:
region <- "TL"
path_wgcna_data <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/RData/")
path_wgcna <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/")


# Load data
net_name <- "MeSMeGL_TL"

comp <- "MeSMeGL"
load(paste0(path_wgcna_data, paste0("norm_cts_", comp, ".RData"))) # Variable named norm_cts
load(paste0(path_wgcna_data, net_name, ".RData"))


####### Make datTraits object
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/metadata_datasetR.txt"
meta <- read.csv(file = file_metadata, sep = "\t", header = T)
rownames(meta) <- meta$id
#Filter out samples
samples_to_remove <- "Or49DE"
meta_fil <- meta[!row.names(meta) %in% samples_to_remove, ]
# Filter region
meta_reg <- meta_fil[meta_fil$region == region & meta$species == "Multifasciatus", ]
# Make spp columns
meta_reg$phenotype <- ifelse(meta_reg$phenotype == "GL", 1, 0)
# Make datTraits
datTraits <- meta_reg$phenotype
# Save for later access
file_name <- paste0("datTraits_", comp)
save(datTraits, file = paste0(path_wgcna_data, file_name, ".RData"))

# Change network variable name
net <- MuSMuGL_TL
mod <- net$MEs
length(mod)

# Set output folder
folder_out <- paste0(path_wgcna, "/", net_name, "/")
dir.create(folder_out, showWarnings = F)


# How many modules?
length(unique(net$colors))
# How many genes per module?
moduleSizes <- as.data.frame(table(net$colors))
colnames(moduleSizes) <- c("Module", "GeneCount")
moduleSizes[order(moduleSizes$GeneCount),] 

# 1. Plot dendrogram
# Plot the dendrogram and the module colors before and after merging underneath 
plotDendroAndColors(net$dendrograms[[1]],  #order
                    cbind(net$unmergedColors, net$colors), 
                    c("unmerged", "merged"),  dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.03, guideHang = 0.05, main = "Weighted Gene Correlation Network Analysis")

# 2. Plot heatmap of sample-to-sample correlations
mod <- net$MEs
mod_ord <- orderMEs(mod) # Reorder modules so similar modules are next to each other
module_order = names(mod_ord) %>% gsub("ME","", .)
mod_ord$genotype = row.names(mod_ord) # Add treatment names
# tidy data
mME = mod_ord %>%
  pivot_longer(-genotype) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
# Plot
p <- mME %>% ggplot(., aes(x=genotype, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
p
# Save plot
ggsave(plot = p, file = paste0(folder_out, "sample_sample_corr_", net_name, ".pdf"), device = "pdf", width = 10, height = 20)



# 3. Quantify module-trait correlations
moduleTraitCor <- cor(mod, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(norm_cts))
merged_val <- merge(moduleTraitCor,moduleTraitPvalue, by = "row.names", all = TRUE )

# Convert p-values matrix to vector, adjust, then reshape back
pval_vec <- as.vector(moduleTraitPvalue)

# FDR correction (Benjamini-Hochberg)
pval_adj_vec <- p.adjust(pval_vec, method = "fdr")

# Put back into matrix with same dimensions
moduleTraitFDR <- matrix(pval_adj_vec,
                         nrow = nrow(moduleTraitPvalue),
                         ncol = ncol(moduleTraitPvalue))

rownames(moduleTraitFDR) <- rownames(moduleTraitPvalue)
colnames(moduleTraitFDR) <- colnames(moduleTraitPvalue)

# Merge everything
merged_val <- data.frame(
  module = rownames(moduleTraitCor),
  moduleTraitCor,
  moduleTraitFDR
)

colnames(merged_val) <- c("module", "GL_corr", "Ornatipinnis_corr", "Multifasciatus_corr", "Meeli_corr", "GL_padj", "Ornatipinnis_padj", "Multifasciatus_padj", "Meeli_padj")

# Join number of genes
num <- as.data.frame(table(net$colors))
colnames(num) <- c("module", "num_genes")
num$module <- paste0("ME", num$module)
merge_num <- merge(merged_val, num, by = "module")
merge_num_ord <- merge_num[order(abs(merge_num$GL_corr), decreasing = T),]

# Save
write.table(merge_num_ord,
            paste0(folder_out, "corr_module_trait_", net_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)




# 4. Plot heatmap module-trait correlations
textMatrix <- paste(base::signif(moduleTraitCor, 2), "\n(",
                    base::signif(moduleTraitPvalue, 1), ")", sep="")
textMatrix <- base::signif(moduleTraitCor, 2)
dim(textMatrix) <- dim(moduleTraitCor)
sizeGrWindow(10,6)
par(mar=c(6,8.5,3,3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = c("GL", "Orna", "Multi", "Meeli"), 
               yLabels = names(mod), xLabelsAdj = 0.52, ySymbols = names(mod),
               colorLabels = F, colors = colorRampPalette(c("#4472C4", "white", "#BF9000"))(50), zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = F, cex.text = 0.5,xLabelsPosition = "bottom", xLabelsAngle = 0,cex.lab.x = 1,
               cex.legendLabel = 0.5,
               main = paste("Module-trait correlation"))
# Save manually: 7 x 20. Name: "module_trait_hmap_" + net name




# 5. Lollipop chart modules
# Load data
data <- read.csv(paste0(folder_out, "corr_module_trait_", net_name, ".txt"),
                 sep = "\t",
                 header = TRUE)

# Reorder modules ONLY by GL correlation
data <- data %>%
  mutate(module = fct_reorder(module, abs(GL_corr), .desc = TRUE))

# Plot
p <- data %>%
  # Filter significant GL modules
  filter(GL_padj < 0.05) %>%
  
  # Long format
  pivot_longer(
    cols = c(GL_corr, Ornatipinnis_corr, Multifasciatus_corr, Meeli_corr),
    names_to = "variable",
    values_to = "value"
  ) %>%
  
  ggplot(aes(
    x = module,
    y = abs(value),
    color = variable,
    alpha = variable == "GL_corr" 
  )) +
  
  # Lollipop sticks
  geom_segment(
    aes(xend = module, y = 0, yend = abs(value)),
    color = "grey70"
  ) +
  
  # Points
  geom_point(aes(size = num_genes)) +
  
  # Transparency control
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.35),  # GL = solid, others faded
    guide = "none"
  ) +
  
  # Colors
  scale_color_manual(
    name = "Variable",
    values = c(
      "GL_corr" = "#104E8B",
      "Ornatipinnis_corr" = "#8B4500",
      "Multifasciatus_corr" = "#CD1076",
      "Meeli_corr" = "gold"
    ),
    labels = c(
      "GL_corr" = "Group-Living",
      "Ornatipinnis_corr" = "Ornatipinnis",
      "Multifasciatus_corr" = "Multifasciatus",
      "Meeli_corr" = "Meeli"
    )
  ) +
  
  # Labels
  labs(
    x = "Module",
    y = "Correlation (absolute value)",
    title = "Module–Trait Correlations"
  ) +
  
  coord_cartesian(ylim = c(0.3, NA)) +
  
  # Size scaling
  scale_size_continuous(
    name = "Number of genes",
    range = c(4, 12)
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Show
p

# Save
ggsave(plot = p,
       filename = paste0(folder_out, "module_lollipop.pdf"),
       device = "pdf",
       height = 4,
       width = 10)






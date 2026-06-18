{
  library(stringr)
  library(WGCNA)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
  library(dplyr)
}

#Region:
region <- "TL"
path_wgcna_data <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/RData/")
path_wgcna <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/")


# Load data
net_name <- "OrMuMe_TL"
comp <- "OrMuMe"
load(paste0(path_wgcna_data, paste0("datTraits_", comp, ".RData"))) # Variable named datTraits
load(paste0(path_wgcna_data, paste0("norm_cts_", comp, ".RData"))) # Variable named norm_cts
load(paste0(path_wgcna_data, net_name, ".RData"))

# Change network variable name
net <- OrMuMe_TL
mod <- net$MEs
length(mod)

# Set output folder
folder_out <- paste0(path_wgcna, "/", net_name, "/")
dir.create(folder_out, showWarnings = F)

# 1. Prepare data
#Define the trait of interest
trait_of_int <- as.data.frame(datTraits)
#colnames(trait_of_int) <- "group"
#Change module names (remove 'ME' from the begining)
modNames <- substring(names(mod), 3)



# 2. Gene - Module membership - Correlation between gene expression and ME 
# Correlations and p-values
# Gene–module membership correlations
geneModuleMembership <- as.data.frame(cor(norm_cts, mod, use = "p"))
MMPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneModuleMembership), nrow(norm_cts))
)

# Rename columns
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

# Add gene IDs
geneModuleMembership$gene <- rownames(geneModuleMembership)
MMPvalue$gene <- rownames(MMPvalue)

# Module assignment
dt <- data.frame(
  gene = names(net$colors),
  module = paste0("ME", net$colors)
)

# Merge everything
MM_dt <- dt %>%
  left_join(geneModuleMembership, by = "gene") %>%
  left_join(MMPvalue, by = "gene") %>%
  rowwise() %>%
  mutate(
    MM_corr = get(paste0("MM", substring(module, 3))),
    MM_pval = get(paste0("p.MM", substring(module, 3)))
  ) %>%
  ungroup() %>%
  dplyr::select(gene, module, MM_corr, MM_pval)

# FDR correction (global across all genes)
MM_dt <- MM_dt %>%
  mutate(MM_padj = p.adjust(MM_pval, method = "fdr")) %>%
  select(-MM_pval)

# Write output
write.table(MM_dt,
            file = paste0(folder_out, "module_membership_", net_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)


# 3. Gene Significance - Correlation value of each gene with trait of interest
# Gene–trait correlations
trait_of_int <- trait_of_int$phenotype # Select only GL as trait of interest
geneTraitSignificance <- as.data.frame(cor(norm_cts, trait_of_int, use = "p"))
colnames(geneTraitSignificance) <- "corrGS_GL"

# P-values
GSPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneTraitSignificance), nrow(norm_cts))
)
colnames(GSPvalue) <- "GS_GL_pval"

# Merge correlation + p-values
geneTraitSignificance_mer <- merge(
  geneTraitSignificance,
  GSPvalue,
  by = "row.names"
)
colnames(geneTraitSignificance_mer)[1] <- "gene"

# Module assignment
gene_module_key <- tibble::enframe(net$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))

# Combine all
GeneSignificance <- left_join(
  gene_module_key,
  geneTraitSignificance_mer,
  by = "gene"
)

# 🔹 FDR correction (global across all genes, per trait column)
GeneSignificance <- GeneSignificance %>%
  mutate(
    GS_GL_padj = p.adjust(GS_GL_pval, method = "fdr")
  ) %>%
  select(-GS_GL_pval)

# Order by gene significance
GeneSignificance_ord <- GeneSignificance %>%
  arrange(desc(abs(corrGS_GL)))

# Write output
write.table(GeneSignificance_ord, file = paste0(folder_out, "gene-module_trait_significance_", net_name, ".txt"),
            sep = "\t",
            col.names = T,
            quote = F,
            row.names = F)

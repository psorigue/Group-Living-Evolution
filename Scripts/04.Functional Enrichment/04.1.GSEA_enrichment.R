# This script performs functional enrichment analysis (GO and KEGG) 
# on the results of differential gene expression analysis. 
# It uses the clusterProfiler package to conduct 
# Gene Set Enrichment Analysis (GSEA) based on the 'stat' values from 
# DESeq2 results. The script reads in the DEG results, prepares a ranked 
# gene list, and then performs GSEA for both GO terms and KEGG pathways.

{
  library(clusterProfiler) # version 4.18.4
  library(AnnotationHub) # version 4.0.0
}

# Set region
region <- "TL"

# Set paths and read files.
home <- path.expand("~")
path_deg <- file.path(home, "04.Differential Gene Expression", region)
path_enrich <- file.path(home, "05.Functional Enrichment", region, "GSEA_datasets")
file_gene_info <- file.path(home, "Ref_genome", "gene_info_full.txt")

# Download tilapia functional annotation data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]

# INDEX
# 1. Read DEG results
# 2. Gene Set Enrichment Analysis (GSEA) -> Gene Ontology (GO)
# 3. Gene Set Enrichment Analysis (GSEA) -> KEGG Pathways


# 1. Read DEG results
#--------------------
# Set comparison name. In this case, Ornatipinnis Solitary vs Multi Group-Living
comp_name <- "OrS_vs_MuGL"
# Read DEG results and gene info file
ds_deg <- read.csv(file = file.path(path_deg, paste0("DEG_", comp_name, ".txt")),
                   sep = "\t", header = T, row.names = 1)
ds_info_genes <- read.csv(file = file.path(home, "Ref_genome", "gene_info_full.txt"),
                          sep = "\t", header = T)
ds_info_genes_sub <- ds_info_genes[,c("id", "symbol_gtf")] # Subset only 2 cols of interest


# 2. Gene Set Enrichment Analysis (GSEA) -> Gene Ontology (GO)
#------------------------------------------------------
# Read gene list for GSEA keeping the 'stat' column from DESeq2 results
all_genes <- data.frame(symbol_gtf = rownames(ds_deg),
                        stat = ds_deg$stat)

# Include Entrez gene ID and remove NAs
all_genes_mer <- merge(all_genes, ds_info_genes_sub, by = "symbol_gtf", all.x = TRUE)
all_genes_mer <- na.omit(all_genes_mer)

# Create named vector
geneList <- all_genes_mer$stat
names(geneList) <- all_genes_mer$id

# Sort decreasing (IMPORTANT for later interpretation of results)
geneList <- sort(geneList, decreasing = TRUE)

# Perform GSEA for GO terms
gsea_go <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Oni.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# Simplify GO terms to remove redundancy and order by adjusted p-value
gsea_go_s_df <- as.data.frame(simplify(gsea_go))
gsea_go_s_df_ord <- gsea_go_s_df[order(gsea_go_s_df$p.adjust),]

# Write GSEA GO results
write.table(gsea_go_s_df_ord,
            file = file.path(path_enrich, paste0("GSEA_GO_", comp_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)


# 3. Gene Set Enrichment Analysis (GSEA) -> KEGG Pathways
#--------------------------------------------------------
# Perform GSEA for KEGG pathways. Uses same geneList as GO
gsea_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "onl", # KEGG code for Nile tilapia
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# Convert to data frame and order by adjusted p-value
gsea_kegg_df <- as.data.frame(gsea_kegg)
gsea_kegg_df <- gsea_kegg_df[order(gsea_kegg_df$p.adjust), ]

# Write GSEA KEGG results
write.table(gsea_kegg_df,
            file = file.path(path_enrich, paste0("GSEA_KEGG_", comp_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)


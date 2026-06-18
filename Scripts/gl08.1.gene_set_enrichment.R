library(clusterProfiler)
library(AnnotationHub)


# Download tilapia data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb


# Path
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"

# Load data
region <- "DE"
comp_name <- "OrS_vs_MuGL"

ds_deg <- read.csv(file = paste0(path_deg, region, "/DEG/DEG_", comp_name, ".txt"),
                   sep = "\t",
                   header = T, row.names = 1)
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)
ds_info_genes_sub <- ds_info_genes[,c("id", "symbol_gtf")] # Subset only 2 cols of interest


# GSEA

# Merge ALL genes (not filtered)
all_genes <- data.frame(symbol_gtf = rownames(ds_deg),
                        stat = ds_deg$stat)



# Merge with IDs
all_genes_mer <- merge(all_genes, ds_info_genes_sub, by = "symbol_gtf", all.x = TRUE)

# Remove NA
all_genes_mer <- na.omit(all_genes_mer)

# Create named vector
geneList <- all_genes_mer$stat
names(geneList) <- all_genes_mer$id

# Sort decreasing (IMPORTANT)
geneList <- sort(geneList, decreasing = TRUE)

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
gsea_go_s <- simplify(gsea_go)
gsea_go_s_df <- as.data.frame(gsea_go_s)
gsea_go_s_df_ord <- gsea_go_s_df[order(gsea_go_s_df$p.adjust),]

gsea_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "onl",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
gsea_kegg_df <- as.data.frame(gsea_kegg)
gsea_kegg_df <- gsea_kegg_df[order(gsea_kegg_df$p.adjust), ]

write.table(as.data.frame(gsea_go_s_df_ord),
            file = paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name, ".txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(as.data.frame(gsea_kegg),
            file = paste0(path_deg, region, "/enrichment/GSEA/GSEA_KEGG_", comp_name, ".txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

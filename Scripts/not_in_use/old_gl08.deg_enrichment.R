{
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  library(gprofiler2)
  library(AnnotationHub)
  library(ggplot2)
}


# Path
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"

# Load data
region <- "TL"
comp_name <- "MuS_vs_MuGL"

ds_deg <- read.csv(file = paste0(path_deg, region, "/DEG/DEG_", comp_name, ".txt"),
                   sep = "\t",
                   header = T, row.names = 1)
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)

# From the DEG, obtain gene ids
# Filter genes by p-adj value
padj_threshold <- 0.05
foldch_threshold <- 1

{
  ds_deg_fil <- na.omit(ds_deg[ds_deg$padj < padj_threshold,])
  ds_deg_fil2 <- ds_deg_fil[abs(ds_deg_fil$log2FoldChange) > foldch_threshold,]
  genes <- rownames(ds_deg_fil2)
}

# Download tilapia data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb


# Convert data
new_df <- data.frame(genes) # Create df
colnames(new_df) <- "symbol_gtf"
#new_df$symbol_gtf <- gsub(new_df$symbol_gtf, pattern = "gene-", replacement = "")
ds_info_genes_sub <- ds_info_genes[,c("id", "symbol_gtf")] # Subset only 2 cols of interest
new_df_mer <- merge(new_df, ds_info_genes_sub, by = "symbol_gtf", all.x = T) # Merge
ids <- new_df_mer$id
sum(is.na(ids))
ids <- na.omit(ids)


# Go analysis
ego <- enrichGO(
  gene          = ids,        # vector of Entrez IDs
  OrgDb         = org.Oni.eg.db,
  keyType       = "ENTREZID",      # other options: "SYMBOL", etc.
  ont           = "all",            # "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego

dotplot(ego, showCategory = 30, font.size = 9)
barplot(ego, showCategory = 30, font.size = 9)
cnetplot(ego, showCategory = 30, font.size = 9)

# Save plot
folder_enr <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/05.DEG/", region, "/enrichment/plots/")
ggsave(dotplot(ego, showCategory = 30, font.size = 9), file = paste0(folder_enr, "/GOdot_", comp_name, ".pdf"),
       device = "pdf", width = 7, height = 7)
ggsave(cnetplot(ego, showCategory = 30, font.size = 9), file = paste0(folder_enr, "/GOnet_", comp_name, ".pdf"),
       device = "pdf", width = 7, height = 7)

# Save data frame
ego_df <- as.data.frame(ego)
ego_df_ord <- ego_df[order(ego_df$p.adjust),]

# Write enrichment data
write.table(ego_df_ord, file = paste0(path_deg, region, "/enrichment/GO_", comp_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)


# Simplify
ego_s <- simplify(ego)
# Save simplified plot
folder_enr <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/05.DEG/", region, "/enrichment/plots/")
ggsave(dotplot(ego_s, showCategory = 30, font.size = 9), file = paste0(folder_enr, "/GOdot_simplif_", comp_name, ".pdf"),
       device = "pdf", width = 7, height = 7)
ggsave(cnetplot(ego_s, showCategory = 30, font.size = 9), file = paste0(folder_enr, "/GOnet_simplif_", comp_name, ".pdf"),
       device = "pdf", width = 7, height = 7)
# Save data frame
ego_s_df <- as.data.frame(ego_s)
ego_s_df_ord <- ego_s_df[order(ego_s_df$p.adjust),]
# Write simplified enrichment data
write.table(ego_s_df_ord, file = paste0(path_deg, region, "/enrichment/GOsimplif_", comp_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)


# KEGG
ekegg <- enrichKEGG(
  gene         = ids,
  organism     = "onl",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# Convert + clean
ekegg_df <- as.data.frame(ekegg)
ekegg_df <- ekegg_df[order(ekegg_df$p.adjust), ]

# Save
write.table(
  ekegg_df,
  file = paste0(path_deg, region, "/enrichment/KEGG_", comp_name, ".txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

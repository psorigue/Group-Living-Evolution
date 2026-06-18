{
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  #library(gprofiler2)
  #BiocManager::install("AnnotationHub")
  library(AnnotationHub)
}




# Download tilapia data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb


# Variables
net_name <- "OrMuMe_TL"
comp <- "OrMuMe"
region <- "TL"

# Load files
path_wgcna <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/"
df_module_corr <- read.csv(file = paste0(path_wgcna,
                                         region, 
                                         "/", net_name, "/corr_module_trait_", net_name, ".txt"),
                           sep = "\t",
                           header = T)
df_gene_module <- read.csv(file = paste0(path_wgcna,
                                         region, 
                                         "/", net_name, "/gene-module_trait_significance_", net_name, ".txt"),
                          sep = "\t",
                          header = T)
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)



module <- "MEdeeppink1"
genes <- df_gene_module[df_gene_module$module %in% module,]$gene

# Set out folder
folder_enrich <- paste0(path_wgcna, region, "/", net_name, "/", module, "/enrichment/") ; dir.create(folder_enrich, showWarnings = F)

# Convert data
new_df <- data.frame(genes) # Create df
colnames(new_df) <- "symbol_gtf"
new_df$symbol_gtf <- gsub(new_df$symbol_gtf, pattern = "gene-", replacement = "")
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
  ont           = "ALL",            # "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego

dotplot(ego, showCategory = 30)
cnetplot(ego, showCategory = 30)

ego_df <- as.data.frame(ego)
ego_df_fil <- ego_df[, c("ONTOLOGY", "ID", "Description", "GeneRatio", "RichFactor", "p.adjust")]
ego_df_fil_ord <- ego_df_fil[order(ego_df_fil$p.adjust),]

# Write enrichment data
write.table(ego_df_fil_ord, file = paste0(folder_enrich, "/GO_", module, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Save plots
folder_plot <- paste0(folder_enrich, "plots/") ; dir.create(folder_plot, showWarnings = F)
ggsave(dotplot(ego, showCategory = 30), file = paste0(folder_plot, "/GO_dot_", net_name, ".pdf"), width = 5, height = 7, device = "pdf")
ggsave(cnetplot(ego, showCategory = 30), file = paste0(folder_plot, "/GO_net_", net_name, ".pdf"), width = 10, height = 10, device = "pdf")


# If I need to simplify
ego_s <- clusterProfiler::simplify(ego)
dotplot(ego_s, showCategory = 30)
cnetplot(ego_s, showCategory = 30)

ego_s_df <- as.data.frame(ego_s)
ego_s_df_fil <- ego_s_df[, c("ONTOLOGY", "ID", "Description", "GeneRatio", "RichFactor", "p.adjust")]
ego_s_df_fil_ord <- ego_s_df_fil[order(ego_s_df_fil$p.adjust),]

# Write enrichment data
write.table(ego_s_df_fil_ord, file = paste0(folder_enrich, "/GOsimplif_", module, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Save plots
folder_plot <- paste0(folder_enrich, "plots/") ; dir.create(folder_plot, showWarnings = F)
ggsave(dotplot(ego_s, showCategory = 30, font.size = 9), file = paste0(folder_plot, "/GOsimplif_dot_", net_name, ".pdf"), width = 5, height = 7, device = "pdf")
ggsave(cnetplot(ego_s, showCategory = 30), file = paste0(folder_plot, "/GOsimplif_net_", net_name, ".pdf"), width = 10, height = 10, device = "pdf")










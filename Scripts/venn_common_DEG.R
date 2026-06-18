library(ggVennDiagram)
library(UpSetR)
library(AnnotationHub)


region <- "TL"

# Read DEG results for each comparison

comparisons <- c("OrS_vs_MuGL", "OrS_vs_MeGL", "MuS_vs_MeGL", "MeS_vs_MuGL")



# Store DEG tables
deg_tables <- list()

# Store significant genes
gene_lists <- list()

for (comp_name in comparisons) {
  
  file <- paste0(
    "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/",
    region,
    "/DEG/DEG_",
    comp_name,
    ".txt"
  )
  
  dt <- read.csv(file, sep = "\t", header = TRUE)
  
  # Keep significant genes
  gene_lists[[comp_name]] <- dt[dt$padj < 0.05 & abs(dt$log2FoldChange) > 1, "gene"]
  
  # Store full table
  deg_tables[[comp_name]] <- dt
}



ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "white", high = "steelblue")


# Another way to visualize it
upset(fromList(gene_lists))



# Common genes across all comparisons
common_genes_all <- Reduce(intersect, gene_lists)
common_genes <- common_genes_all[!is.na(common_genes_all)]

# Annotation
gene_info <- read.csv(
  "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
  sep = "\t",
  header = TRUE
)

gene_info_red <- gene_info[, c("symbol_gtf", "name_datasets")]

# Base dataframe
common_names <- gene_info_red[
  match(common_genes, gene_info_red$symbol_gtf),
]

# Replace symbol column with original vector
common_names$symbol_gtf <- common_genes

# Manually add names (extracted from NCBI annotation)
if (region == "TL") {
  common_names[common_names$symbol_gtf == "LOC102082837", "name_datasets"] <- 
    "uncharacterized LOC102082837"
  
  common_names[common_names$symbol_gtf == "LOC109201417", "name_datasets"] <- 
    "uncharacterized LOC109201417"
} 
if (region == "DE") {
    common_names[common_names$symbol_gtf == "nbisL1-trna-16", "name_datasets"] <- 
      "tRNA"
}



# Add fold-change direction for each comparison
for (comp_name in comparisons) {
  
  dt <- deg_tables[[comp_name]]
  
  # Keep only common genes
  dt_sub <- dt[dt$gene %in% common_genes, ]
  
  # Match order
  dt_sub <- dt_sub[
    match(common_names$symbol_gtf, dt_sub$gene),
  ]
  
  # Add log2FC
  common_names[[paste0(comp_name, "_log2FC")]] <- dt_sub$log2FoldChange
  
  # Add direction
  common_names[[paste0(comp_name, "_direction")]] <-
    ifelse(
      dt_sub$log2FoldChange > 0,
      "UP",
      "DOWN"
    )
}

direction_cols <- grep("_direction$", names(common_names), value = TRUE)

common_names$overall_direction <- apply(
  common_names[, direction_cols],
  1,
  function(x) {
    x <- na.omit(x)
    
    if (length(x) == 0) return(NA)
    
    if (all(x == "UP")) {
      "UP"
    } else if (all(x == "DOWN")) {
      "DOWN"
    } else {
      "MIXED"
    }
  }
)

# Save table
write.table(common_names, paste0(
  "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/",
  region,
  "/Venn_interspp_common/common_genes.txt"),
  sep = "\t", quote = F, col.names = T, row.names = F)

# What genes are not MIXED?
genes_not_mixed <- common_names[common_names$overall_direction != "MIXED",]$symbol_gtf

common_names[common_names$overall_direction != "MIXED",]




# GO
genes <- common_names$symbol_gtf

# Download tilapia data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb

# Convert data
new_df <- data.frame(genes) # Create df
colnames(new_df) <- "symbol_gtf"
#new_df$symbol_gtf <- gsub(new_df$symbol_gtf, pattern = "gene-", replacement = "")
ds_info_genes_sub <- gene_info[,c("id", "symbol_gtf")] # Subset only 2 cols of interest
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



## ---- CHANGE PATHS
dotplot(ego, showCategory = 30, font.size = 9)
barplot(ego, showCategory = 30, font.size = 9)
cnetplot(ego, showCategory = 30, font.size = 9)

# Save plot
folder_enr <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/", region, "/enrichment/plots/")
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
dotplot(ego_s, showCategory = 30, font.size = 9)
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

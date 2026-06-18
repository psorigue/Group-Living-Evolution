


region <- "TL"

genes_AS_dt <- read.csv(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/08.isoform_analysis/AS-DEG_common/AS_ids_",
                            region, ".txt"), sep = "\t", header = TRUE)


# Read DEG results for each comparison
comparisons <- c("OrS_vs_MuS", "OrS_vs_MeS", "OrS_vs_MuGL", "OrS_vs_MeGL", "MuS_vs_MeGL", "MeS_vs_MuGL")

# Store significant genes
deg_tables <- list()

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


# Comupte common genes Ornatipinnis and facultative removing phylogenetic signal
# OrS_MuGL minus OrS_MuS
genes_Mu <- setdiff(gene_lists$OrS_vs_MuGL, gene_lists$OrS_vs_MuS)
genes_Me <- setdiff(gene_lists$OrS_vs_MeGL, gene_lists$OrS_vs_MeS)



# With gene symbol
genes <- gene_lists$MuS_vs_MeGL
intersect(genes_AS_dt$gene_symbol,genes_Me)
intersect(genes_Me,genes_AS_dt$gene_symbol)


# With gene_id
genes <- gene_lists$MeS_vs_MuGL
genes <- genes_Me


new_df <- data.frame(genes[!is.na(genes)]) # Create df
colnames(new_df) <- "symbol_gtf"
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)
ds_info_genes_sub <- ds_info_genes[,c("id", "symbol_gtf")] # Subset only 2 cols of interest
ds_info_genes_sub <- na.omit(ds_info_genes_sub)
new_df_mer <- merge(new_df, ds_info_genes_sub, by = "symbol_gtf", all.x = T) # Merge
# fill in "LOCs"
new_df_mer$id <- ifelse(
  is.na(new_df_mer$id) & grepl("^LOC", new_df_mer$symbol_gtf),
  sub("^LOC", "", new_df_mer$symbol_gtf),
  new_df_mer$id
)
new_df_mer[is.na(new_df_mer$id),]
ids <- new_df_mer$id
sum(is.na(ids))
ids <- na.omit(ids)

intersect(genes_AS_dt$gene_id,ids)
intersect(ids,genes_AS_dt$gene_id)

# This script computes semantic similarity between GO terms enriched 
# in the two comparisons, clusters them based on similarity, 
# and identifies shared clusters between datasets. It also annotates 
# clusters with representative GO terms and summarizes the number of 
# GO terms per cluster and dataset.

{
    library(GOSemSim) # version 2.36.0
    library(GO.db) # version 3.22.0
    library(AnnotationHub) # version 4.0.0
    library(AnnotationDbi) # version 1.72.0
    library(pheatmap) # version 1.0.13
    library(tidyr) # version 1.3.2
    library(dplyr) # version 1.1.4
}

# Set region
region <- "TL"

# Set paths and read files.
home <- path.expand("~")
path_datasets <- file.path(home, "05.Functional_Enrichment", region, "GSEA_datasets")
out_dir <- file.path(home, "05.Functional_Enrichment", region, "intraspp_common", "GO_similarity")

# Download tilapia functional annotation data
ah <- AnnotationHub()
org.Oni.eg.db <- ah[["AH119811"]]


# INDEX
# 1. Helper function to get GO term names from GO IDs.
# 2. Load GSEA GO results for both comparisons.
# 3. Compute semantic similarity matrix for all GO terms together.
# 4. Cluster GO terms based on similarity.
# 5. Annotate dataset origin for each GO term.
# 6. Identify clusters shared between datasets and add adjusted p-values and Normalized Enrichment Score.
# 7. Get representative GO term for each shared cluster.
# 8. Summarize number of GO terms per cluster and dataset, and add representative.



# 1. Helper function to get GO term names from GO IDs
# ---------------------------------------------------
get_go_terms <- function(go_ids) {
  sapply(go_ids, function(x) {
    term <- GOTERM[[x]]
    if (is.null(term)) return(NA)
    Term(term)
  })
}

# 2. Load GSEA GO results for both comparisons
# --------------------------------------------
# Set comparison names
comp_name1 <- "MuS_vs_MuGL"
comp_name2 <- "MeS_vs_MeGL"
set1 <- "Multifasciatus"
set2 <- "Meeli"

# Read GSEA GO results
ds_deg1 <- read.csv(file.path(path_datasets, paste0("GSEA_GO_", comp_name1, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1)
ds_deg2 <- read.csv(file.path(path_datasets, paste0("GSEA_GO_", comp_name2, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1)
# Extract GO IDs from both datasets
go_list1 <- rownames(ds_deg1)
go_list2 <- rownames(ds_deg2)


# 3. Compute semantic similarity matrix for all GO terms together
# ---------------------------------------------------------------
# Get all unique GO terms
all_terms <- unique(c(go_list1, go_list2))

# Compute semantic similarity matrix using Wang method
sim_all <- mgoSim(all_terms, all_terms,
                  semData = semData,
                  measure = "Wang",
                  combine = NULL)

# Save similarity matrix
write.table(sim_all, file = file.path(out_dir, "semantic_similarity_matrix.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

# Save heatmap -> High dimensions needed
pdf(file.path(out_dir, "semantic_similarity_heatmap.pdf"), width = 30, height = 30)
pheatmap(sim_all)
dev.off()


# 4. Cluster GO terms based on similarity
# ---------------------------------------
hc <- hclust(as.dist(1 - sim_all), method = "average")

# Cut tree into clusters
h <- 0.8
cluster_assign <- cutree(hc, h = h)

# Inspect number of categories per module
table(cluster_assign)

# Save cluster assignments
cluster_df <- data.frame(
  GO = names(cluster_assign),
  cluster = cluster_assign
)
cluster_df$Description <- get_go_description(cluster_df$GO)


# 5. Annotate dataset origin for each GO term
# -------------------------------------------
cluster_df$set <- ifelse(cluster_df$GO %in% go_list1 & cluster_df$GO %in% go_list2, "both",
                         ifelse(cluster_df$GO %in% go_list1, set1, set2))
# Write output
write.table(cluster_df, file.path(out_dir, "GO_clusters_all.txt"), sep = "\t", row.names = FALSE, quote = F)


# 6. Identify clusters shared between datasets and add 
# adjusted p-values and Normalized Enrichment Score.
# ----------------------------------------------------
# Identify clusters that contain GO terms from both datasets
shared_clusters_flag <- with(cluster_df, tapply(set, cluster, function(x) {
  length(unique(x)) > 1
}))

# Get GO terms from shared clusters and their descriptions
shared_cluster_ids <- names(shared_clusters_flag[shared_clusters_flag])
shared_df <- cluster_df[cluster_df$cluster %in% shared_cluster_ids, ]
shared_df$Description <- get_go_description(shared_df$GO)

# Add adjusted p-values and Normalized Enrichment Score from initial datasets
padj1 <- ds_deg1$p.adjust
names(padj1) <- rownames(ds_deg1)
nes1 <- ds_deg1$NES
names(nes1) <- rownames(ds_deg1)

padj2 <- ds_deg2$p.adjust
names(padj2) <- rownames(ds_deg2)
nes2 <- ds_deg2$NES
names(nes2) <- rownames(ds_deg2)

# Create a combined column for adjusted p-values and NES for shared clusters
shared_df_all <- shared_df %>%
  rowwise() %>%
  mutate(
    padj = if (set == "both") {
      paste0(set1, ":", signif(padj1[GO], 3), ";",
             set2, ":", signif(padj2[GO], 3))
    } else if (set == set1) {
      as.character(signif(padj1[GO], 3))
    } else {
      as.character(signif(padj2[GO], 3))
    },
    
    NES = if (set == "both") {
      paste0(set1, ":", signif(nes1[GO], 3), ";",
             set2, ":", signif(nes2[GO], 3))
    } else if (set == set1) {
      as.character(signif(nes1[GO], 3))
    } else {
      as.character(signif(nes2[GO], 3))
    }
  ) %>%
  ungroup()

# Order by cluster and select relevant columns
shared_df_ord <- shared_df_all[order(shared_df_all$cluster),
                           c("cluster", "GO", "Description", "set", "padj", "NES")]
# Save shared clusters
write.table(shared_df_ord,
            file.path(out_dir, "GO_clusters_shared.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# 7. Get representative GO term for each shared cluster
# -----------------------------------------------------
# Function to get representative GO term for a cluster based on average semantic similarity
get_representative <- function(cluster_id, sim_matrix, clusters) {
  terms <- names(clusters[clusters == cluster_id])
  sub_sim <- sim_matrix[terms, terms, drop = FALSE]
  mean_sim <- rowMeans(sub_sim, na.rm = TRUE)
  rep_term <- names(which.max(mean_sim))
  return(rep_term)
}

# Get representative GO term for each shared cluster
cluster_ids <- unique(shared_df$cluster)
representatives <- sapply(cluster_ids, function(cl) {
  get_representative(cl, sim_all, cluster_assign)
})

# Create a data frame with cluster IDs and their representative GO terms
cluster_labels <- data.frame(
  cluster = cluster_ids,
  representative_GO = representatives
)

# Add GO term descriptions for representative GO terms
cluster_labels$Description <- get_go_description(cluster_labels$representative_GO)


# 8. Summarize number of GO terms per cluster and dataset, and add representative
# --------------------------------------------------------------------------------
cluster_summary <- aggregate(GO ~ cluster + set, data = shared_df, length)
cluster_wide <- cluster_summary %>%
  pivot_wider(
    names_from = set,
    values_from = GO,
    values_fill = 0   # fills missing combinations with 0
  )

final_table <- merge(cluster_wide, cluster_labels, by = "cluster", all.x = TRUE)

final_table <- final_table %>%
  select(cluster, representative_GO, Description, both, Meeli, Multifasciatus)

write.table(final_table, file.path(out_dir, "cluster_summary_GOcounts.txt"), sep = "\t", quote = F, row.names = FALSE)

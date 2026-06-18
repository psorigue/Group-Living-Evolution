# This script computes semantic similarity between GO terms enriched 
# in all interspecies comparisons, clusters them based on similarity, 
# and identifies shared clusters between datasets. 


# Load libraries
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
path_datasets <- file.path(home, "05.Functional Enrichment", region, "GSEA_datasets")
out_dir <- file.path(home, "05.Functional Enrichment", region, "interspp")

# Download tilapia functional annotation data
ah <- AnnotationHub()
org.Oni.eg.db <- ah[["AH119811"]]

# Build semantic data (REQUIRED for mgoSim)
semData <- godata(annoDb = org.Oni.eg.db, ont = "BP")


# INDEX
# 1.1. Helper function to get GO term names from GO IDs.
# 1.2. Helper function to compare two datasets and categorize GO terms.
# 2. Load GO lists
# 3. Apply comparison function and save results
# 4. Compute semantic similarity matrix for all GO terms together
# 5. Clustering GO terms based on semantic similarity
# 6. Annotate dataset origin for each GO term
# 7. Identify clusters shared between datasets 
# and add adjusted p-values and Normalized Enrichment Score




# 1.1. Helper function to get GO term names from GO IDs
# ---------------------------------------------------
get_go_terms <- function(go_ids) {
  sapply(go_ids, function(x) {
    term <- GOTERM[[x]]
    if (is.null(term)) return(NA)
    Term(term)
  })
}

# 1.2. Helper function to compare two datasets and categorize GO terms
# --------------------------------------------------------------------
compare_states <- function(gsea_sol, gsea_grp) {
  
  merged <- full_join(
    gsea_sol %>% dplyr::select(ID, NES_sol = NES),
    gsea_grp %>% dplyr::select(ID, NES_grp = NES),
    by = "ID"
  ) %>%
    mutate(
      in_sol = !is.na(NES_sol),
      in_grp = !is.na(NES_grp),
      
      category = case_when(
        in_grp & !in_sol ~ "group_specific",
        in_sol & !in_grp ~ "solitary_specific",
        in_grp & in_sol & sign(NES_grp) != sign(NES_sol) ~ "direction_switch",
        in_grp & in_sol ~ "shared_same_direction"
      ),
      
      delta_NES = NES_grp - NES_sol
    )
  
  return(merged)
}



# 2. Load GO lists
# ----------------
comp_name1 <- "OrS_vs_MuS"
comp_name2 <- "OrS_vs_MuGL"
comp_name3 <- "OrS_vs_MeS"
comp_name4 <- "OrS_vs_MeGL"
set1 <- "Multi_solitary"
set2 <- "Multi_group"
set3 <- "Meeli_solitary"
set4 <- "Meeli_group"

ds_deg1 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name1, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)
ds_deg2 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name2, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)
ds_deg3 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name3, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)
ds_deg4 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name4, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)

ds_deg1 <- ds_deg1 %>% mutate(ID = rownames(ds_deg1))
ds_deg2 <- ds_deg2 %>% mutate(ID = rownames(ds_deg2))
ds_deg3 <- ds_deg3 %>% mutate(ID = rownames(ds_deg3))
ds_deg4 <- ds_deg4 %>% mutate(ID = rownames(ds_deg4))



# 3. Apply comparison function and save results
# ---------------------------------------------
merged_Mu <- compare_states(ds_deg1, ds_deg2)
merged_Me <- compare_states(ds_deg3, ds_deg4)

write.table(merged_Mu, file = file.path(out_dir, "Multifasciatus_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(merged_Me, file = file.path(out_dir, "Meeli_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# 4. Compute semantic similarity matrix for all GO terms together
# ---------------------------------------------------------------
all_terms <- unique(c(merged_Mu$ID, merged_Me$ID))

sim_all <- mgoSim(all_terms, all_terms,
                  semData = semData,
                  measure = "Wang",
                  combine = NULL)

# Save similarity matrix
write.table(sim_all, file = file.path(out_dir, "semantic_similarity_matrix.txt"), 
            sep = "\t", quote = F, row.names = T, col.names = NA)

# Plot heatmap
pdf(file.path(out_dir, "semantic_similarity_heatmap.pdf"), width = 30, height = 30)
pheatmap(sim_all)
dev.off()



# 5. Clustering GO terms based on semantic similarity
# ---------------------------------------------------
hc <- hclust(as.dist(1 - sim_all), method = "average")

# Adjust height depending on resolution
cluster_assign <- cutree(hc, h = 0.75)

# Inspect number of categories per module
table(cluster_assign)

# Save cluster assignments
cluster_df <- data.frame(
  GO = names(cluster_assign),
  cluster = cluster_assign
)
cluster_df$Description <- get_go_terms(cluster_df$GO)



# 6. Annotate dataset origin for each GO term
# -------------------------------------------
cluster_df <- cluster_df %>%
  # Join Multifasciatus (Mu)
  left_join(
    merged_Mu %>%
      select(ID, NES_sol, NES_grp, category, delta_NES) %>%
      rename(
        NES_sol_Mu = NES_sol,
        NES_grp_Mu = NES_grp,
        category_Mu = category,
        delta_NES_Mu = delta_NES
      ),
    by = c("GO" = "ID")
  ) %>%
  
  # Join Meeli (Me)
  left_join(
    merged_Me %>%
      select(ID, NES_sol, NES_grp, category, delta_NES) %>%
      rename(
        NES_sol_Me = NES_sol,
        NES_grp_Me = NES_grp,
        category_Me = category,
        delta_NES_Me = delta_NES
      ),
    by = c("GO" = "ID")
  )

write.table(cluster_df, file.path(out_dir, "GO_clusters_all.txt"), sep = "\t", row.names = FALSE, quote = F)


# 7. Identify clusters shared between datasets 
# and add adjusted p-values and Normalized Enrichment Score
# ---------------------------------------------------------
shared_clusters_flag <- tapply(cluster_df$GO, cluster_df$cluster, function(go_terms) {
  
  subset <- cluster_df[cluster_df$GO %in% go_terms, ]
  
  any(subset$category_Mu == "group_specific", na.rm = TRUE) &
    any(subset$category_Me == "group_specific", na.rm = TRUE)
  
})

# Get GO terms from shared clusters and their descriptions
shared_cluster_ids <- names(shared_clusters_flag[shared_clusters_flag])
shared_df <- cluster_df[cluster_df$cluster %in% shared_cluster_ids, ]

# Add adjusted p-values for each dataset
padj1 <- setNames(ds_deg1$p.adjust, rownames(ds_deg1))
padj2 <- setNames(ds_deg2$p.adjust, rownames(ds_deg2))
padj3 <- setNames(ds_deg3$p.adjust, rownames(ds_deg3))
padj4 <- setNames(ds_deg4$p.adjust, rownames(ds_deg4))

# Add adjusted p-values for each dataset to shared_df
shared_df <- shared_df %>%
  mutate(
    padj_Mu_sol = padj1[GO],
    padj_Mu_grp = padj2[GO],
    
    padj_Me_sol = padj3[GO],
    padj_Me_grp = padj4[GO]
  )

# Reorder columns for clarity and save results
shared_df_ord <- shared_df %>%
  arrange(cluster) %>%
  select(cluster, GO, Description,
         starts_with("category"),
         starts_with("NES"),
         starts_with("delta"),
         starts_with("padj"))

write.table(shared_df_ord,
            file.path(out_dir, "GO_clusters_shared.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


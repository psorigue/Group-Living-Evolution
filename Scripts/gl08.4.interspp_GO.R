############################
### 1. LOAD LIBRARIES
############################
{
  library(clusterProfiler)
  library(GOSemSim)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(AnnotationHub)
  library(pheatmap)

}

ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb

# Build semantic data (REQUIRED for mgoSim)
semData <- godata(annoDb = org.Oni.eg.db, ont = "BP")

#### HELPER FUNCTION to extract GO descriptions from terms
get_go_description <- function(go_ids) {
  sapply(go_ids, function(x) {
    term <- GOTERM[[x]]
    if (is.null(term)) return(NA)
    Term(term)
  })
}


############################
### 3. INPUT PARAMETERS
############################
region <- "TL"

path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"

out_dir <- file.path(path_deg, region, "/enrichment/GSEA/Interspp/")
dir.create(out_dir, showWarnings = F, recursive = T)

comp_name1 <- "OrS_vs_MuS"
comp_name2 <- "OrS_vs_MuGL"
comp_name3 <- "OrS_vs_MeS"
comp_name4 <- "OrS_vs_MeGL"
set1 <- "Multi_solitary"
set2 <- "Multi_group"
set3 <- "Meeli_solitary"
set4 <- "Meeli_group"

############################
### 4. LOAD GO LISTS
############################
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

############################
### 2. FUNCTION: COMPARE STATES
############################
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

############################
### 3. APPLY TO SPECIES AND SAVE
############################
merged_Mu <- compare_states(ds_deg1, ds_deg2)
merged_Me <- compare_states(ds_deg3, ds_deg4)

write.table(merged_Mu, file = file.path(out_dir, "Multifasciatus_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(merged_Me, file = file.path(out_dir, "Meeli_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


############################
### 6. OPTIONAL: SUMMARIES
############################
# Count categories
cat("B category counts:\n")
print(table(merged_Mu$category))

cat("\nC category counts:\n")
print(table(merged_Me$category))



############################
### 6. COMPUTE SEMANTIC SIMILARITY
############################
all_terms <- unique(c(merged_Mu$ID, merged_Me$ID))

sim_all <- mgoSim(all_terms, all_terms,
                  semData = semData,
                  measure = "Wang",
                  combine = NULL)

# Save similarity matrix
write.table(sim_all, file = file.path(out_dir, "semantic_similarity_matrix.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

# Plot heatmap
pdf(file.path(out_dir, "/semantic_similarity_heatmap.pdf"), width = 30, height = 30)
pheatmap(sim_all)
dev.off()




############################
### 7. CLUSTERING
############################
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
cluster_df$Description <- get_go_description(cluster_df$GO)





############################
### 8. ANNOTATE DATASET ORIGIN
############################
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


############################
### 9. IDENTIFY SHARED CLUSTERS
############################
shared_clusters_flag <- tapply(cluster_df$GO, cluster_df$cluster, function(go_terms) {
  
  subset <- cluster_df[cluster_df$GO %in% go_terms, ]
  
  any(subset$category_Mu == "group_specific", na.rm = TRUE) &
    any(subset$category_Me == "group_specific", na.rm = TRUE)
  
})

shared_cluster_ids <- names(shared_clusters_flag[shared_clusters_flag])


shared_df <- cluster_df[cluster_df$cluster %in% shared_cluster_ids, ]


padj1 <- setNames(ds_deg1$p.adjust, rownames(ds_deg1))
padj2 <- setNames(ds_deg2$p.adjust, rownames(ds_deg2))
padj3 <- setNames(ds_deg3$p.adjust, rownames(ds_deg3))
padj4 <- setNames(ds_deg4$p.adjust, rownames(ds_deg4))

shared_df <- shared_df %>%
  mutate(
    padj_Mu_sol = padj1[GO],
    padj_Mu_grp = padj2[GO],
    
    padj_Me_sol = padj3[GO],
    padj_Me_grp = padj4[GO]
  )

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

############################
### 10. REPRESENTATIVE GO TERM (MEDOID)
############################
get_representative <- function(cluster_id, sim_matrix, clusters) {
  
  terms <- names(clusters[clusters == cluster_id])
  
  sub_sim <- sim_matrix[terms, terms, drop = FALSE]
  
  mean_sim <- rowMeans(sub_sim, na.rm = TRUE)
  
  rep_term <- names(which.max(mean_sim))
  
  return(rep_term)
}

cluster_ids <- unique(shared_df$cluster)

representatives <- sapply(cluster_ids, function(cl) {
  get_representative(cl, sim_all, cluster_assign)
})

cluster_labels <- data.frame(
  cluster = cluster_ids,
  representative_GO = representatives
)

cluster_labels$Description <- get_go_description(cluster_labels$representative_GO)

write.table(cluster_labels,
            file.path(out_dir, "GO_representatives.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


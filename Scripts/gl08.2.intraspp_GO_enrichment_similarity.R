############################
### 1. LOAD LIBRARIES
############################
library(GOSemSim)
library(GO.db)
library(AnnotationHub)
library(AnnotationDbi)
library(pheatmap)
library(tidyr)
library(dplyr)

############################
### 2. LOAD ANNOTATION DB
############################
ah <- AnnotationHub()
org.Oni.eg.db <- ah[["AH119811"]]

# Build semantic data (REQUIRED for mgoSim)
semData <- godata(annoDb = org.Oni.eg.db, ont = "BP")


############################
### 3. INPUT PARAMETERS
############################
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"
path_deg <- "../../to_work/"
region <- "DE"

comp_name1 <- "MuS_vs_MuGL"
comp_name2 <- "MeS_vs_MeGL"
set1 <- "Multifasciatus"
set2 <- "Meeli"

# Output folder
out_dir <- file.path(path_deg, region, "/enrichment/GSEA/Intraspp/GO_similarity/")
dir.create(out_dir, showWarnings = F, recursive = T)


############################
### 4. LOAD GO LISTS
############################
ds_deg1 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name1, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)

ds_deg2 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_GO_", comp_name2, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)

go_list1 <- rownames(ds_deg1)
go_list2 <- rownames(ds_deg2)


#### HELPER FUNCTION to extract GO descriptions from terms
get_go_description <- function(go_ids) {
  sapply(go_ids, function(x) {
    term <- GOTERM[[x]]
    if (is.null(term)) return(NA)
    Term(term)
  })
}


############################
### 5. CLEAN GO TERMS
############################
# Replace obsolete GO terms
#go_list1 <- gsub("GO:0006082", "GO:0008152", go_list1)
#go_list2 <- gsub("GO:0006082", "GO:0008152", go_list2)


############################
### 6. COMPUTE SEMANTIC SIMILARITY
############################
all_terms <- unique(c(go_list1, go_list2))

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
cluster_df$set <- ifelse(cluster_df$GO %in% go_list1 & cluster_df$GO %in% go_list2, "both",
                         ifelse(cluster_df$GO %in% go_list1, set1, set2))


write.table(cluster_df, file.path(out_dir, "GO_clusters_all.txt"), sep = "\t", row.names = FALSE, quote = F)


############################
### 9. IDENTIFY SHARED CLUSTERS
############################
shared_clusters_flag <- with(cluster_df, tapply(set, cluster, function(x) {
  length(unique(x)) > 1
}))

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

shared_df_ord <- shared_df_all[order(shared_df_all$cluster),
                           c("cluster", "GO", "Description", "set", "padj", "NES")]
# Save shared clusters
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


############################
### 11. CLUSTER SUMMARY
############################
cluster_summary <- aggregate(GO ~ cluster + set, data = shared_df, length)
cluster_wide <- cluster_summary %>%
  pivot_wider(
    names_from = set,
    values_from = GO,
    values_fill = 0   # fills missing combinations with 0
  )

final_table <- merge(cluster_wide, cluster_labels, by = "cluster", all.x = TRUE)

final_table <- final_table %>%
  dplyr::select(cluster, representative_GO, Description, both, Meeli, Multifasciatus)

write.table(final_table, file.path(out_dir, "cluster_summary_GOcounts.txt"), sep = "\t", quote = F, row.names = FALSE)

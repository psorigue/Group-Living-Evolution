############################
### 1. LOAD LIBRARIES
############################
library(clusterProfiler)

############################
### 2. LOAD KEGG ANNOTATION
############################
kegg_info <- clusterProfiler::download_KEGG("onl")
kegg_map <- kegg_info$KEGGPATHID2NAME

############################
### 3. INPUT PARAMETERS
############################
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"
region <- "DE"

comp_name1 <- "MuS_vs_MuGL"
comp_name2 <- "MeS_vs_MeGL"
set1 <- "Multifasciatus"
set2 <- "Meeli"

out_dir <- file.path(path_deg, region, "/enrichment/GSEA/Intraspp/")
dir.create(out_dir, showWarnings = FALSE, recursive = T)

############################
### 4. LOAD KEGG LISTS
############################
ds_deg1 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_KEGG_", comp_name1, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)

ds_deg2 <- read.csv(paste0(path_deg, region, "/enrichment/GSEA/GSEA_KEGG_", comp_name2, ".txt"),
                    sep = "\t", header = TRUE, row.names = 1)

kegg_list1 <- rownames(ds_deg1)
kegg_list2 <- rownames(ds_deg2)

############################
### 5. FIND COMMON PATHWAYS
############################
common_kegg <- intersect(kegg_list1, kegg_list2)

############################
### 6. ADD DESCRIPTIONS
############################
get_kegg_description <- function(kegg_ids, kegg_map) {
  
  df <- data.frame(from = kegg_ids)
  
  merged <- merge(df, kegg_map,
                  by = "from",
                  all.x = TRUE,
                  sort = FALSE)
  
  merged$to
}

common_df <- data.frame(
  KEGG = common_kegg,
  Description = get_kegg_description(common_kegg, kegg_map)
)

############################
### 7. ADD STATISTICS
############################
common_df$p_adj_set1 <- ds_deg1[common_kegg, "p.adjust"]
common_df$p_adj_set2 <- ds_deg2[common_kegg, "p.adjust"]

common_df$NES_set1 <- ds_deg1[common_kegg, "NES"]
common_df$NES_set2 <- ds_deg2[common_kegg, "NES"]

colnames(common_df) <- c("KEGG", "Description",
                         paste0(set1, "_padj"),
                         paste0(set2, "_padj"),
                         paste0(set1, "_NES"),
                         paste0(set2, "_NES"))

############################
### 8. SAVE OUTPUT
############################
write.table(common_df,
            file.path(out_dir, "KEGG_common.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

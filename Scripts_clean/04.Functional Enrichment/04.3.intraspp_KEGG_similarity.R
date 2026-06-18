# This script identifies common KEGG pathways between two GSEA comparisons.
# It reads GSEA KEGG results for both comparisons, finds common KEGG pathways,
# retrieves their descriptions, and compiles a summary table with adjusted 
# p-values and Normalized Enrichment Scores (NES) for both datasets. 

library(clusterProfiler)

# Set region
region <- "TL"

# Set paths and read files.
home <- path.expand("~")
path_datasets <- file.path(home, "05.Functional Enrichment", region, "GSEA_datasets")
out_dir <- file.path(home, "05.Functional Enrichment", region, "intraspp_common")

# Download tilapia functional annotation data
kegg_info <- clusterProfiler::download_KEGG("onl")
kegg_map <- kegg_info$KEGGPATHID2NAME


# INDEX
# 1. Load GSEA KEGG results for both comparisons.
# 2. Find common KEGG pathways between datasets.
# 3. Get descriptions for common KEGG pathways.
# 4. Add adjusted p-values and NES for both datasets.
# 5. Save results.


# 1. Load GSEA KEGG results for both comparisons
# --------------------------------------------
# Set comparison names
comp_name1 <- "MuS_vs_MuGL"
comp_name2 <- "MeS_vs_MeGL"
set1 <- "Multifasciatus"
set2 <- "Meeli"

# Read GSEA KEGG results
ds_deg1 <- read.csv(file.path(path_datasets, paste0("GSEA_KEGG_", comp_name1, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1)
ds_deg2 <- read.csv(file.path(path_datasets, paste0("GSEA_KEGG_", comp_name2, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1)
# Extract KEGG IDs from both datasets
kegg_list1 <- rownames(ds_deg1)
kegg_list2 <- rownames(ds_deg2)


# 2. Find common KEGG pathways between datasets
# ---------------------------------------------
common_kegg <- intersect(kegg_list1, kegg_list2)


# 3. Get descriptions for common KEGG pathways
# --------------------------------------------
# Function to get KEGG descriptions
get_kegg_description <- function(kegg_ids, kegg_map) {
  df <- data.frame(from = kegg_ids)
  merged <- merge(df, kegg_map,
                  by = "from",
                  all.x = TRUE,
                  sort = FALSE)
  merged$to
}

# Create a data frame for common KEGG pathways with descriptions
common_df <- data.frame(
  KEGG = common_kegg,
  Description = get_kegg_description(common_kegg, kegg_map)
)

# 4. Add adjusted p-values and NES for both datasets
# --------------------------------------------------
common_df$p_adj_set1 <- ds_deg1[common_kegg, "p.adjust"]
common_df$p_adj_set2 <- ds_deg2[common_kegg, "p.adjust"]

common_df$NES_set1 <- ds_deg1[common_kegg, "NES"]
common_df$NES_set2 <- ds_deg2[common_kegg, "NES"]

colnames(common_df) <- c("KEGG", "Description",
                         paste0(set1, "_padj"),
                         paste0(set2, "_padj"),
                         paste0(set1, "_NES"),
                         paste0(set2, "_NES"))

# 5. Save results
# ---------------
write.table(common_df,
            file.path(out_dir, "KEGG_common.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
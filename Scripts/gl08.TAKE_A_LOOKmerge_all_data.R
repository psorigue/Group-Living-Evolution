#install.packages("readxl")
library(readxl)
library(dplyr)
library(purrr)

path_gl <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/Expression_genes_databases/"

# Check sheets file
excel_sheets("path/to/your/file.xlsx")

# Example: Read a sheet by name
data <- read_excel("path/to/your/file.xlsx", sheet = "SheetName", skip = 1) # Skip -> skip rows


# Read gene info
file_gene_info <- paste0(path_gl, "All_genes_db.txt")
ds_gene_info <- read.csv(file = file_gene_info, sep = "\t", header = T)
ds_gene_info_red <- ds_gene_info[,c("id", "symbol", "name")]

# 1. Read PhyDGET
file_phydget <- paste0(path_gl, "PhyDGET.xlsx")
excel_sheets(file_phydget)
# all_spp model
ds_phydget_tips <- read_excel(file_phydget, sheet = "all_spp_corrected (tips)")
ds_phydget_tips_red <- ds_phydget_tips[,c("id", "BF")]
colnames(ds_phydget_tips_red) <- c("id", "PhyDGET_TipsComb")
# all_pb_npb
ds_phydget_nodes_pb_no_pb <- read_excel(file_phydget, sheet = "all_pb_npb (nodes)")
ds_phydget_nodes_pb_red <- ds_phydget_nodes_pb_no_pb[,c("id", "BF.PB")]
colnames(ds_phydget_nodes_pb_red) <- c("id", "PhyDGET_NodesPB")
ds_phydget_nodes_nopb_red <- ds_phydget_nodes_pb_no_pb[,c("id", "BF.nPB")]
colnames(ds_phydget_nodes_nopb_red) <- c("id", "PhyDGET_NodesNoPB")
ds_phydget_nodes_comb_red <- ds_phydget_nodes_pb_no_pb[,c("id", "BF.combined")]
colnames(ds_phydget_nodes_comb_red) <- c("id", "PhyDGET_NodesComb")
# single

# 2. Read DEG
file_deg <- paste0(path_gl, "DEG.xlsx")
excel_sheets(file_deg)
ds_OrMe <- read_excel(file_deg, sheet = "OrS_vs_MeGL", skip = 1) ; ds_OrMe_red <- ds_OrMe[,c("id", "logFC", "padj")] ; colnames(ds_OrMe_red) <- c("id", "DEG_OrMe_LFC", "DEG_OrMe_padj")
ds_OrMu <- read_excel(file_deg, sheet = "OrS_vs_MuGL3", skip = 1) ; ds_OrMu_red <- ds_OrMu[,c("id", "logFC", "padj")] ; colnames(ds_OrMu_red) <- c("id", "DEG_OrMu_LFC", "DEG_OrMu_padj")
ds_MeMe <- read_excel(file_deg, sheet = "MeS_vs_MeGL", skip = 1) ; ds_MeMe_red <- ds_MeMe[,c("id", "logFC", "padj")] ; colnames(ds_MeMe_red) <- c("id", "DEG_MeMe_LFC", "DEG_MeMe_padj")
ds_MuMu <- read_excel(file_deg, sheet = "MuS_vs_MuGL3", skip = 1) ; ds_MuMu_red <- ds_MuMu[,c("id", "logFC", "padj")] ; colnames(ds_MuMu_red) <- c("id", "DEG_MuMu_LFC", "DEG_MuMu_padj")

# 3. Read WGCNA
file_wgcna <- paste0(path_gl, "WGCNA_groupliving_gene-module.txt")
ds_wgcna <- read.csv(file = file_wgcna, sep = "\t", header = T)


# 4. Merge all
# Combine datasets using reduce and left_join
datasets_list <- list(ds_phydget_tips_red, ds_phydget_nodes_pb_red, ds_phydget_nodes_nopb_red, ds_phydget_nodes_comb_red, ds_OrMe_red, ds_OrMu_red, ds_MeMe_red, ds_MuMu_red)
merged_data <- reduce(datasets_list, full_join, by = "id")
merged_data_plus <- right_join(ds_gene_info_red, merged_data, by = "id")
merged_final <- left_join(merged_data_plus, ds_wgcna, by = "symbol")

write.table(merged_final, file = paste0(path_gl, "overview_expression_genes.txt"), 
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)



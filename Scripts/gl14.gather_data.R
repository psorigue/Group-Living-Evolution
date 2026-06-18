# Take from Epi last script in WGCNA

net_name <- "OrMuMe_s_010"
path_wgcna <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06n.WGCNA/TL/", net_name, "/")

deg_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06n.WGCNA/TL/DEG_list_join.txt"


# Process DEG file
deg <- read.csv(deg_file, sep = "\t", header = T)
deg_coll <- deg %>%
  group_by(gene) %>%
  summarise(
    across(c(comp, padj), ~ paste(unique(.), collapse = ","))
  )


# Load modules correlation file
file_cm <- paste0(path_wgcna, "/corr_module_trait_", net_name, ".txt")
cor_mod <- read.csv(file_cm, sep = "\t", header = T)
# Take significan modules
min_r <- 0.6
modules_dt <- cor_mod[abs(cor_mod$GL_corr) > min_r,c("module", "GL_corr", "num_genes")]
module_names <- modules_dt$module

# Load gene significance file
file_gs <- paste0(path_wgcna, "/gene-module_trait_significance_", net_name, ".txt")
gs_gen <- read.csv(file_gs, sep = "\t", header = T)

# Load module membership file
file_mm <- paste0(path_wgcna, "/module_membership_", net_name, ".txt")
mm_gen <- read.csv(file_mm, sep = "\t", header = T)

# Build dataset
dt <- gs_gen[gs_gen$module %in% module_names,c("gene", "module", "corrGS_GL")]
dt_mer <- merge(dt, modules_dt, by = "module")
colnames(dt_mer) <- c("module", "gene", "GScorr", "corr_mod_trait", "num_gen_mod")

# Merge module membership
dt_mer1 <- merge(dt_mer, MM_dt[,c("gene", "MM_corr")], by = "gene")
dt_mer1$gene <- gsub(dt_mer1$gene, pattern = "gene-", replacement = "")

# Is differentially expressed?
dt_mer2 <- merge(dt_mer1, deg_coll[,c("gene", "comp", "padj")], by = "gene", all.x = T)
colnames(dt_mer2) <- c("gene", "module", "GScorr", "corr_mod_trait", "num_gen_mod", "MMcorr", "comp_deg", "padj_deg")

# Is a TF?
tf_df <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/Oreochromis_niloticus_TF.txt",
                  sep = "\t",
                  header = T)
tf_df_fil <- tf_df[,c("Symbol", "Family")]
tf_df_fil2 <- tf_df_fil %>% distinct(Symbol, .keep_all = TRUE) # Only keeps first occurrence of each gene
colnames(tf_df_fil2) <- c("gene", "family_tf")
dt_mer3 <- left_join(dt_mer2, tf_df_fil2, by = "gene")
colnames(dt_mer3)
dt_mer3 <- dt_mer3[,c("gene", "module", "GScorr", "MMcorr", "corr_mod_trait", "num_gen_mod", "comp_deg", "padj_deg", "family_tf")]



# Check genes that are TFs
dt_mer3[!is.na(dt_mer3$family_tf),]



# Write dataset
# Order by padj value
write.table(dt_mer3, file = paste0(path_wgcna, "/../overview_", net_name, "_sign_modules.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

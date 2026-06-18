{
  library(stringr)
  library(WGCNA)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
  library(dplyr)
  library(tidyverse)
}


# Path
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/metadata_datasetR.txt"
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)

#Region:
region <- "TL"
path_wgcna_data <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/RData/")
path_wgcna <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/")


# Load data
net_name <- "OrMuMe_TL"
comp <- "OrMuMe"
load(paste0(path_wgcna_data, paste0("datTraits_", comp, ".RData"))) # Variable named datTraits
load(paste0(path_wgcna_data, paste0("norm_cts_", comp, ".RData"))) # Variable named norm_cts
load(paste0(path_wgcna_data, net_name, ".RData"))

# Change network variable name
net <- OrMuMe_TL
mod <- net$MEs
length(mod)

# Set output folder
folder_out <- paste0(path_wgcna, "/", net_name, "/")
dir.create(folder_out, showWarnings = F)


# 1. Prepare data
# Key table
gene_module_key <- tibble::enframe(net$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))
#Define the trait of interest
trait_of_int <- trait_of_int$phenotype
#Change module names (remove 'ME' from the begining)
modNames <- substring(names(mod), 3)


# 2. Prepare data individual module
module_name <- "MEdeeppink1"
folder_ind <- paste0(folder_out, module_name, "/")
dir.create(folder_ind, showWarnings = F)
ME_of_interest <- gene_module_key %>% dplyr::filter(module == module_name) 
dim(ME_of_interest)


# 3. Write in a file: gene names, gene significance, module membership, and gene description
## 1: gene names
genes <- ME_of_interest$gene
gene_df <- data.frame(genes) ; colnames(gene_df) <- "symbol_gtf"
gene_df$symbol_gtf <- gsub(gene_df$symbol_gtf, pattern = "gene-", replacement = "")
gene_df <- merge(gene_df, ds_info_genes[,c("symbol_gtf", "name_datasets")], by = "symbol_gtf", all.x = T)
colnames(gene_df) <- c("gene", "description")
## 2: Extract Module membership of genes of the module of interest
###Relating the Gene-tait value with the module key data frame
geneTraitSignificance <- read.csv(paste0(folder_out, "gene-module_trait_significance_", net_name, ".txt"), 
                                  sep = "\t")
geneTraitSignificance_fil <- geneTraitSignificance[geneTraitSignificance$module == module_name, ]

gene_module_key_mer1 <- merge(gene_module_key, geneTraitSignificance_fil[,c("gene", "corrGS_GL")], by = "gene")
geneModuleMembership <- read.csv(paste0(folder_out, "module_membership_", net_name, ".txt"),
                                 sep = "\t")
gene_module_key_fin <- left_join(gene_module_key_mer1, geneModuleMembership[,c("gene", "MM_corr")], by = "gene")
colnames(gene_module_key_fin) <- c("gene", "module", "GScorr", "MMcorr")
## Merge gene info and module info
gene_module_key_fin$gene <- gsub(gene_module_key_fin$gene, pattern = "gene-", replacement = "")
merged_info <- merge(gene_module_key_fin, gene_df, by = "gene", all.x = T)
## Order by MM
merged_info_ord <- merged_info[order(merged_info$MMcorr, decreasing = T),]
# Write output
write.table(merged_info_ord, file = paste0(folder_ind, "GL_genes_GS_MM_description.txt"), sep = "\t", row.names = F, quote = F, col.names = T)
# If I want to filter by Gene Significance value
#module_of_interest_fil <- module_of_interest[abs(module_of_interest$GScorr) > 0.75,]


# 4. Plot GS vs MM -> Save it in file
pdf(paste0(folder_ind, "MM_vs_GS_", module_name, ".pdf"), width = 7, height = 7)
verboseScatterplot(merged_info$MMcorr,
                   merged_info$GScorr,
                   xlab = paste0("Module Membership in ", module_name, " module"),
                   ylab = "gene significance for trait of interest",
                   main = paste("MM vs GS\n"), abline = TRUE, 
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = substring(module_name, 3))
dev.off()



# 5. Plot expression profile
{
  module_color <- substring(module_name, 3)
  meta <- read.csv(file_metadata,
                   sep = "\t", header = TRUE)
  
  # Filter by region
  meta <- meta[meta$region == region,]
  
  # set custom order for the group variable
  meta$phenotype <- factor(meta$phenotype,
                       levels = c("S", "GL"))
  
  # 1. Pull out list of genes in that module
  module_df <- geneTraitSignificance_fil[, c("gene", "module")]
  module_df <- as.data.frame(module_df)
  row.names(module_df) <- module_df$gene
  submod <- module_df %>% subset(module %in% module_name)
  
  # 2. Get normalized expression for those genes
  subexpr <- t(norm_cts)[submod$gene, ]
  
  # 3. Turn into long format and center per gene
  submod_df <- data.frame(subexpr) %>%
    mutate(gene = row.names(.)) %>%
    pivot_longer(-gene, names_to = "name", values_to = "value") %>%
    group_by(gene) %>%
    mutate(value_centered = value - mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(module = module_df[gene, ]$module)
  
  # 4. Filter for single module
  submod_single <- submod_df %>% filter(module == module_name)
  
  # 5. Merge sample metadata and reorder x-axis by phenotype
  submod_single <- submod_single %>%
    left_join(dplyr::select(meta, id, phenotype, species),
              by = c("name" = "id"))
  
  # define sample order based on phenotype (levels set above)
  sample_order <- submod_single %>%
    distinct(name, phenotype) %>%
    arrange(phenotype) %>%
    pull(name)
  
  submod_single <- submod_single %>%
    mutate(name = factor(name, levels = sample_order))
  
  # 6. Compute mean trajectory per sample
  summary_single <- submod_single %>%
    group_by(name) %>%
    summarise(median_value = median(value_centered, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(name = factor(name, levels = sample_order))
  
  # 7. Define background shading by phenotype
  phenotype_df <- submod_single %>%
    distinct(name, phenotype) %>%
    mutate(x = as.numeric(name)) %>%
    group_by(phenotype) %>%
    summarise(xmin = min(x) - 0.5,
              xmax = max(x) + 0.5,
              ymin = -Inf,
              ymax = Inf,
              .groups = "drop")
  
  # 8. Plot
  p <- ggplot(submod_single, aes(x = name, y = value_centered, group = gene)) +
    # background shading
    geom_rect(data = phenotype_df,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax, fill = phenotype),
              inherit.aes = FALSE, alpha = 0.6) +
    scale_fill_manual(values = c("S" = "lightblue", "GL" = "lightpink")) +
    # gene-level lines
    geom_line(color = module_color, alpha = 0.2) +
    
    # mean trajectory
    geom_line(data = summary_single,
              aes(x = name, y = median_value, group = 1),
              color = "black", linewidth = 1.2, inherit.aes = FALSE) +
    # horizontal zero line
    geom_hline(yintercept = 0, color = "black",
               linetype = "dashed", linewidth = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    labs(x = "Samples (ordered by phenotype)",
         y = "Normalized expression (centered at 0)",
         title = module_name) +
    coord_cartesian(ylim = c(-1, 1))
  
  p
  
  # 9. Save plot
  ggsave(plot = p,
         filename = paste0(folder_ind, "expr_profile_", module_name, ".pdf"),
         height = 7, width = 10, device = "pdf")
  
}

p

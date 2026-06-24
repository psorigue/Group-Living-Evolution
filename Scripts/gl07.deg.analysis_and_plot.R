
{
  library(DESeq2)
  library(EnhancedVolcano)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(ggplot2)
}


region <- "TL"

# Path
path_deg <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_solitary/05.DEG/"

# Load DDSobject
load(paste0(path_deg, region, "/DDS_files/DDSobj_all.RData")) # loads variable called dds_obj

ref_spp <- "Ornatipinnis_S" # Species reference
dir_spp <- "Multifasciatus_GL" # Species direction foldChange
comp_name <- "OrS_vs_MuGL"


# DEG analysis
res <- DESeq2::results(
  dds_obj,
  contrast = c("group", dir_spp, ref_spp) 
)
res_ord <- res[order(res$padj),]


# Summary results
summary(res)
sum(res$padj < 0.05, na.rm=TRUE) # How many adjusted p-values were less than 0.05?

# Information about which variables and tests were used
mcols(res)$description

# Join gene description
gene_info <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                      sep = "\t", header = T)
gene_info_red <- gene_info[,c("name_datasets", "symbol_gtf")]
res_ord$gene <- gsub(rownames(res_ord), pattern = "gene-", replacement = "")
merged_df <- as.data.frame(res_ord) %>% left_join(gene_info_red, by = c("gene" = "symbol_gtf"))
merged_df <- merged_df[,c(7,8,1:6)] # Reorder columns
colnames(merged_df)[2] <- "gene_description"

# Write results
write.table(merged_df, 
            file = paste0(path_deg, region, "/DEG/DEG_", comp_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F)



# Plots

# Volcano plot
# Convert to data frame and keep gene names
df <- as.data.frame(res)
df$gene <- rownames(df)

# Remove NA adjusted p-values
df <- df %>% filter(!is.na(padj))

# Add -log10(padj)
df$negLog10Padj <- -log10(df$padj)

# Define significance categories
df <- df %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "Up-regulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

p <- ggplot(df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = sig), size = 2) +
  
  # Threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # Colors
  scale_color_manual(
    values = c(
      "Up-regulated" = "#8B4500",
      "Down-regulated" = "gold2",
      "Not Significant" = "darkgray"
    ),
    labels = c(
      "Up-regulated" = expression(italic("N.multifasciatus")),
      "Down-regulated" = expression(italic("L.ornatipinnis")),
      "Not Significant" = "None"
    ),
    breaks = c(
      "Up-regulated",
      "Down-regulated",
      "Not Significant"
    )
  ) +
  
  labs(
    title = "Telencephalon",
    x = "Log2 Fold Change",
    y = "-log10(adj. p-value)",
    color = "Up-Regulation"
  ) +
  
  theme_minimal() +
  
  # ✅ Add axis lines
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  #theme(legend.position = "none") +
  
  coord_flip()

p
path_plot <- paste0(path_deg, region, "/volcanos/",
                    comp_name, "legend_small.pdf")
ggsave(plot = p, filename = path_plot, device = "pdf", width = 3, height = 2)




# 2. Counts of a particular gene
gene <- "gene-LOC100703747"
d <- plotCounts(dds_obj, gene= gene, intgroup="group", returnData=TRUE)
ggplot(d, aes(x=group, y=count, fill = group)) + 
  geom_violin() + 
  scale_y_log10()


# Make plots of counts 4 different genes:
library(ggplot2)
library(patchwork)

gene1 <- "LOC102083013" ; descr1 <- merged_df[merged_df$gene == gene1,"name_datasets"]
gene2 <- "LOC100696153" ; descr2 <- merged_df[merged_df$gene == gene2,"name_datasets"]
gene3 <- "LOC112847478" ; descr3 <- merged_df[merged_df$gene == gene3,"name_datasets"]
gene4 <- "LOC100710458" ; descr4 <- merged_df[merged_df$gene == gene4,"name_datasets"]

d1 <- plotCounts(dds_obj, gene=paste0("gene-",gene1) , intgroup="phenotype", returnData=TRUE)
d2 <- plotCounts(dds_obj, gene=paste0("gene-",gene2) , intgroup="phenotype", returnData=TRUE)
d3 <- plotCounts(dds_obj, gene=paste0("gene-",gene3), intgroup="phenotype", returnData=TRUE)
d4 <- plotCounts(dds_obj, gene=paste0("gene-",gene4), intgroup="phenotype", returnData=TRUE)

p1 <- ggplot(d1, aes(x = phenotype, y = count, fill = phenotype)) +
  geom_violin() + scale_y_log10() + ggtitle(descr1)

p2 <- ggplot(d2, aes(x = phenotype, y = count, fill = phenotype)) +
  geom_violin() + scale_y_log10() + ggtitle(descr2)

p3 <- ggplot(d3, aes(x = phenotype, y = count, fill = phenotype)) +
  geom_violin() + scale_y_log10() + ggtitle(descr3)

p4 <- ggplot(d4, aes(x = phenotype, y = count, fill = phenotype)) +
  geom_violin() + scale_y_log10() + ggtitle(descr4)

# Combine into 2x2 grid with patchwork
p_all <- (p1 | p2) /
  (p3 | p4)

p_all

ggsave(p_all, file = paste0(path_deg, region, "/count_plots/", comp_name, ".pdf"), width = 8, height = 7, device = "pdf")





#### Count plot for 3 conditions


#### Double volcano
p1 <- EnhancedVolcano(res,
                      #lab = rownames(res),
                      lab = NA,
                      #title = paste(comp_name, region, sep = " "),
                      title = "Multifasciatus",
                      FCcutoff = 1.5,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-10, 10)) #+ coord_flip()
p1
res2 <- results(dds_obj2)

p2 <- EnhancedVolcano(res2,
                      #lab = rownames(res2),
                      lab = NA,
                      #title = paste(comp_name, region, sep = " "),
                      title = "Meeli",
                      FCcutoff = 1.5,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-10, 10)) #+ coord_flip()
p2
grid.arrange(p1, p2,
             ncol=2,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))

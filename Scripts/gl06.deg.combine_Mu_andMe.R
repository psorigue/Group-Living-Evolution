
# CHECK GENES

# Filter metadata and counts
meta_r <- meta[meta$region == region, ]
meta_r <- meta_r[meta_r$species != "Ornatipinnis",]
cts_r  <- cts[, rownames(meta_r)]


# Check what groups exist
print(table(meta_r$group))

# Build DDS
dds <- DESeqDataSetFromMatrix(
  countData = cts_r,
  colData   = meta_r,
  design    = ~ species + phenotype
)

# Remove low-count genes
smallestGroupSize <- min(table(dds$species)) # Recommended by developer
keep_genes <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep_genes, ]


# 7. Fit DESeq model
dds_obj <- DESeq(dds_filtered)


res <- results(
  dds_obj,
  contrast = c("phenotype", "S", "GL")
)

res_ord <- res[order(res$padj),]
as.data.frame(res_ord)

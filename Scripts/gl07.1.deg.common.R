{
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  #library(gprofiler2)
  #BiocManager::install("AnnotationHub")
  library(AnnotationHub)
}



# Download tilapia data
ah <- AnnotationHub()
query(ah, "Oreochromis niloticus")
org.Oni.eg.db <- ah[["AH119811"]]  # Oreochromis niloticus OrgDb

# Read list
de_list <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/05.DEG/DE_interspp_common.txt"
tl_list <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/05.DEG/TL_interspp_common.txt"
de <- read.csv(de_list, header = F)
tl <- read.csv(tl_list, header = F)

ids <- de$V1
ids <- na.omit(ids)

# Go analysis
ego <- enrichGO(
  gene          = ids,        # vector of Entrez IDs
  OrgDb         = org.Oni.eg.db,
  keyType       = "ENTREZID",      # other options: "SYMBOL", etc.
  ont           = "ALL",            # "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego

dotplot(ego, showCategory = 30)

dotplot(simplify(ego), showCategory = 30)

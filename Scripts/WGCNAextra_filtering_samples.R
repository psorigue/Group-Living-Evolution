# This script concatenates the three WGCNA scripts with functions to filter out samples.



##### WGCNA1 #####

library(WGCNA)
library(SummarizedExperiment)
library(DESeq2)
library(gridExtra)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


# Path
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/metadata_datasetR.txt"
path_wgcna <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/"

# Variables to define
region <- "TL" # "TL" or "DE"
phenotypes <- c("S", "GL")
species <- c("Meeli", "Multifasciatus")
median_counts <- 5

# Read counts and colData dataset 
cts <- read.csv(file = file_cts, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file = file_metadata, sep = "\t", header = T)
rownames(meta) <- meta$id

# Filter samples in meta
meta_fil <- meta[meta$region == region & meta$phenotype %in% phenotypes & meta$species %in% species, ]
meta_fil <- meta_fil[,c("phenotype", "species")]
meta_fil$phenotype <- factor(meta_fil$phenotype, levels = phenotypes) # Categorize phenotype column for later use
meta_fil$species <- factor(meta_fil$species, levels = species) # Categorize species column for later use
nrow(meta_fil)

# Apply filter to counts
cts_fil <- cts[, rownames(meta_fil)]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts_fil,
                              colData = meta_fil,
                              design = ~ phenotype) # phenotype or species

# Remove samples if needed
#samples_to_remove <- c("Me33TL", "Or49TL", "Mu7DE", "Mu40TL", "Or48TL")
#samples_to_remove <- c("Me33TL")
dds_rm <- dds[, !colnames(dds) %in% samples_to_remove]
meta_fil_rm <- meta_fil[!row.names(meta_fil) %in% samples_to_remove, ]

# Run DESeq
dds_obj <- DESeq(dds_rm)
dds_obj <- DESeq(dds)

# Filter by low counts
keep <- rowMedians(counts(dds_obj)) >= median_counts
dds_obj_lc <- dds_obj[keep,]


# If needed, filter by varance genes. N is the number to keep
#n <- 2000 # Number of most variable genes to keep
{
  row_var <- rowVars(counts(dds_obj_  ))
  top_n_indices <- order(row_var, decreasing = TRUE)[1:n]
  keep <- rep(FALSE, nrow(counts(dds_obj_  )))
  keep[top_n_indices] <- TRUE
  dds_obj_var <- dds_obj_  [keep,]
}


# Quality control WGCNA
# 1 - Outlier detection and removal funtion
wgcna_good_samples <- function(dds) {
  dataExpr <- t(counts(dds, normalized=TRUE))
  # the data already goes transposed as necessary for wgcna
  gsg <- goodSamplesGenes(dataExpr, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", 
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  gsg$allOK
}
# 2 - Run the detection function. Proceed with the analysis only if the answer is TRUE
wgcna_good_samples(dds_obj_lc) 

# VST - data needs to be in sample to gene format
norm_cts <- vst(dds_obj_lc) %>% assay() %>% t()

# Create design matrix
## Binarise categorical variables
datTraits <- binarizeCategoricalColumns(meta_fil, convertColumns = c("phenotype", "species"), includePairwise = FALSE,
                                        includeLevelVsAll = TRUE, minCount = 1)
rownames(datTraits) <- rownames(meta_fil)
colnames(datTraits) <- c("GL", "Multi")


# Pick Threshold
# Function necessary to visualize soft power thresholds
wgcna_powers_plot <- function(sft_indices){
  # R2 values =SFT.R.sq (maximum) + main connectivity (minimum)
  a1 <- ggplot(sft_indices, aes(Power, SFT.R.sq, label=Power)) +
    geom_point()+
    geom_text(nudge_y = 0.05) +
    geom_hline(yintercept = 0.8, color = "red")+
    labs(x="Power", y="Scale free topology model fit, signed R^2")+
    theme_classic()
  a2 <- ggplot(sft_indices, aes(Power, mean.k., label=Power)) +
    geom_point()+
    geom_text(nudge_y = 0.7) +
    labs(x="Power", y="Mean Connectivity")+
    theme_classic()
  grid.arrange(a1, a2, nrow=2)
  sft_indices
}

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from=12, to=50, by=2))

# Call the notwork topology function
sft <- pickSoftThreshold(norm_cts, dataIsExpr = TRUE, networkType = "signed", 
                         RsquaredCut = 0.85, powerVector = power)
wgcna_powers_plot(sft$fitIndices)



# Construct network
picked_power <- 8
network_type <- "signed hybrid"
dsplit <- 3

# convert matrix to numeric
norm_cts[] <- sapply(norm_cts, as.numeric)

temp_cor <- cor #assign correlation function to function, so WGCNA uses its own function and not the default cor() from R
cor <- WGCNA::cor


### MAKE UNIQUE NAME ON NETWORK
Sys.time()
TL_MuMe_1 = blockwiseModules(norm_cts,  # <= input here
                                           # Network construction arguments: correlation options
                                           corType="pearson",
                                           
                                           # == Adjacency Function ==
                                           power = picked_power,     # power determined by soft threshold 
                                           networkType = network_type,   # Unsigned = negative and positively correlated genes group together, if network is signed, the neg. and pos. correlated genes are stores in separated modules 
                                           
                                           # == Tree and Block Options ==
                                           deepSplit = dsplit, # this goes between 1-4, the biger the value, more refined the clustering more clusters in total
                                           pamRespectsDendro = F,
                                           detectCutHeight = 0.995, 
                                           minModuleSize = 20,       # min number of genes in a module
                                           maxBlockSize = 30000,     # number of genes by block, to save computational resource, one can separate the genes is block of 4k-5k. I prefer that they are analyzed together so my block size > my # genes. If decide to proceed with blocks, there is a penalty threshold blockSizePenaltyPower = 5
                                           
                                           # == Module Adjustments ==
                                           reassignThreshold = 1e-6,
                                           mergeCutHeight = 0.15, #correlation of 0.80 between the clusters to be merged | Default = 0.15
                                           
                                           # == TOM == Archive the run results in TOM file (saves time)
                                           TOMType = "signed", 
                                           TOMDenom = "min",
                                           numericLabels=FALSE,
                                           saveTOMs=F,
                                           #saveTOMFileBase = "ER", # Input the name of your TOM file
                                           
                                           # == Output Options
                                           verbose=5, 
                                           nThreads=30, 
)

cor <- temp_cor


# Save network
net_name <- "TL_all_sign_18_4_15_010"
save(TL_all_sign_18_4_15_010, file = paste0(path_wgcna, net_name, ".RData"))
Sys.time()

length(TL_all_sign_18_4_15_010$MEs)



##### WGCNA2 #####
{
  library(stringr)
  library(WGCNA)
  #library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
}

# Path
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/metadata_datasetR.txt"
path_wgcna <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/"

# Variables to define
region <- "TL" # "TL" or "DE"
phenotypes <- c("S", "GL")
species <- c("Meeli", "Multifasciatus", "Ornatipinnins")
median_counts <- 2
#samples_to_remove <- c("Me33TL", "Or49TL", "Mu7DE", "Mu40TL", "Or48TL")
samples_to_remove <- c("Me33TL", "Mu40TL")

# 1. Create norm_cts and datTraits (same steps as WGCNA1)
{
  # Read counts and colData dataset 
  cts <- read.csv(file = file_cts, sep = "\t", header = T, row.names = 1)
  meta <- read.csv(file = file_metadata, sep = "\t", header = T)
  rownames(meta) <- meta$id
  
  # Filter samples in meta
  meta_fil <- meta[meta$region == region & meta$phenotype %in% phenotypes & meta$species %in% species, ]
  meta_fil <- meta_fil[,c("phenotype", "species")]
  meta_fil$phenotype <- factor(meta_fil$phenotype, levels = phenotypes) # Categorize phenotype column for later use
  meta_fil$species <- factor(meta_fil$species, levels = species) # Categorize species column for later use
  
  # Apply filter to counts
  cts_fil <- cts[, rownames(meta_fil)]
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = cts_fil,
                                colData = meta_fil,
                                design = ~ phenotype) # phenotype or species
  
  # Remove samples if needed
  dds_rm <- dds[, !colnames(dds) %in% samples_to_remove]
  meta_fil_rm <- meta_fil[!row.names(meta_fil) %in% samples_to_remove, ]
  
  # Run DESeq
  dds_obj <- DESeq(dds_rm)
  
  # Filter by low counts
  keep <- rowMedians(counts(dds_obj)) >= median_counts
  dds_obj_lc <- dds_obj[keep,]
  
  # VST - data needs to be in sample to gene format
  norm_cts <- vst(dds_obj_lc) %>% assay() %>% t()
  norm_cts[] <- sapply(norm_cts, as.numeric)
  
  # Create design matrix
  ## Binarise categorical variables
  datTraits <- binarizeCategoricalColumns(meta_fil_rm, convertColumns = c("phenotype", "species"), includePairwise = FALSE,
                                          includeLevelVsAll = TRUE, minCount = 1)
  rownames(datTraits) <- rownames(meta_fil_rm)
  colnames(datTraits) <- c("S")
}


# 2.Load network
net_name <- "TL_all_signhyb_8_3_20"
load(paste0(path_wgcna, net_name, ".RData"))
net <- TL_all_signhyb_8_3_20# Input network variable loaded
mod <- net$MEs
#rm(net)
length(mod)

# 3. Plot dendrogram
# Plot the dendrogram and the module colors before and after merging underneath 
plotDendroAndColors(net$dendrograms[[1]],  #order
                    cbind(net$unmergedColors, net$colors), 
                    c("unmerged", "merged"),  dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.03, guideHang = 0.05, main = "Title")

# 4. Quantify mdule-trait correlations
moduleTraitCor <- cor(mod, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(norm_cts))
merged_val <- merge(moduleTraitCor,moduleTraitPvalue, by = "row.names", all = TRUE )
colnames(merged_val) <- c("module", "corrS", "corrSPP", "pvalS", "pvalSPP")
merged_ord <- merged_val[order(abs(merged_val$corrS), decreasing = T), ]

write.table(merged_ord,
            paste0(path_wgcna, "corr_", net_name, ".txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)


#Display the correlation values withing a heatmap plot
textMatrix <- paste(base::signif(moduleTraitCor, 2), "\n(",
                    base::signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)

sizeGrWindow(10,6)
par(mar=c(6,8.5,3,3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), 
               yLabels = names(mod), xLabelsAdj = 0.52, ySymbols = names(mod),
               colorLabels = F, colors = colorRampPalette(c("#4472C4", "white", "#BF9000"))(50), zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = F, cex.text = 0,
               main = paste("Module-trait relationships"))


ggsave(module_heatmap, file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Presentations/SPE_2024/fig_wgcna.pdf",
       device = "pdf",
       width = 2,
       height = 8)


#### WGCNA3 #####

{
  library(stringr)
  library(WGCNA)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  library(gprofiler2)
}

path_wgcna <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/"
ds_info_genes <- read.csv(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Ref_genome/gene_info_full.txt",
                          sep = "\t",
                          header = T)

# datTraits and norm_cts files need to be loaded

# Load network
net_name <- "TL_all_sign_18_4_15_010"
#net_name <- "TL_MeSMeGLMuSMuGL_hyb_14_4_10"
load(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", net_name, ".RData"))
mod <- TL_all_sign_18_4_15_010$MEs
net <- TL_all_sign_18_4_15_010

length(mod)

moduleCmodmoduleColors
#Definethe trait of interest
trait_of_int <- as.data.frame(as.numeric(datTraits$S))
#Rearrange module names (remove ME from the begining)
modNames <- substring(names(mod), 3)

# Gene - Module membership - correlation between gene expression and ME (to search for highly connected genes)
geneModuleMembership <- as.data.frame(cor(norm_cts, mod, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(norm_cts)))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)
geneModuleMembership$gene <- rownames(geneModuleMembership) # to join for a module_key table later
geneModuleMembership
MMPvalue

#Correlation value of each gene with trait of interest
geneTraitSignificance <- as.data.frame(cor(norm_cts, trait_of_int, use = "p"))
# pvalues for each gene 
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(norm_cts)))
names(GSPvalue) <- paste0("p.GS", names(trait_of_int))
geneTraitSignificance$gene <- rownames(geneTraitSignificance) # to join for a module_key table later

colnames(geneTraitSignificance) <- c("GScorr", "gene")
geneTraitSignificance <- geneTraitSignificance[,c("gene", "GScorr")]
geneTraitSignificance_ord <- geneTraitSignificance[order(geneTraitSignificance$GScorr, decreasing = T),]

# Write results in a new folder
dir.create(paste0(path_wgcna, net_name))
write.table(geneTraitSignificance_ord, paste0(path_wgcna, net_name, "/gene-trait_sign_", net_name, ".txt"),
            quote = F,
            col.names = T,
            row.names = F)


geneModuleMembership$MMlightcyan1
#Intramodular analysis: identifying genes with hight GS and MM - hub genes
module = "lightcyan1"
column <- match(module, modNames)
moduleGenes <- net$colors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste0("Module Membership in ", module, " module"),
                   ylab = "gene significance for trait of interest",
                   main = paste("MM vs GS\n"), abline = TRUE, 
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2) #, col = module)
table(net$colors)


gene_module_key <- tibble::enframe(net$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))
gene_module_key

merged <- left_join(gene_module_key, geneTraitSignificance, by = "gene")

write.table(merged, file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL_MeSMeGLMuSMuGL_hyb_14_4_10/gene-module.txt",
            sep = "\t",
            col.names = T,
            row.names = F)

# Save the results
#readr::write_tsv(gene_module_key, file = file.path("path", "title.tsv"))

# Try some module
ME_of_interest <- gene_module_key %>% dplyr::filter(module == "MElightcyan1") 
dim(ME_of_interest)


# If there is a moment to do the enrichment, it's here
genes <- ME_of_interest$gene
{
  new_df <- data.frame(genes) # Create df
  colnames(new_df) <- "id"
  new_df$id <- gsub(new_df$id, pattern = "gene-", replacement = "")
  ds_info_genes_sub <- ds_info_genes[,c("id", "symbol_gtf")] # Subset only 2 cols of interest
  new_df_mer <- merge(new_df, ds_info_genes_sub, by = "id", all.x = T) # Merge
  
  ids <- new_df_mer$id
  sum(is.na(ids))
  ids <- na.omit(ids)
  
  gost_res <- gost(query = ids,
                   organism = "oniloticus",
                   #significant = F,
                   sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"))
}

gost_res$result

#Relating the Gene-tait value with the module key data frame
gene_module_key_mer1 <- merge(gene_module_key, geneTraitSignificance, by = "gene")
gene_module_key_fin <- merge(gene_module_key_mer1, geneModuleMembership, by = "gene")
gene_module_key_fin

colnames(new_df) <- "gene"
geneTraitSignificance$gene <- gsub(geneTraitSignificance$gene, pattern = "gene-", replacement = "")
new_df_1 <- merge(new_df, geneTraitSignificance, by = "gene")
colnames(new_df_1) <- c("gene", "corr_with_trait")
write.table(new_df_1, sep = "\t", quote = F, 
            file = paste0(path_wgcna, net_name, "/ME", module, "_genes_traitcorr.txt"))


#Extract only genes that are correlated with the trait of interest >0.5
top_cont_genes <- gene_module_key_fin[abs(gene_module_key_fin$GScorr) > 0.5,]
top_cont_genes_fin <- top_cont_genes[order(top_cont_genes$GScorr, decreasing = T),]  
top_cont_genes_fin



# The darkgreen cluster with 302 genes and cor=0.62
module <- "MEorange"
index_module <- grep(module, colnames(top_cont_genes_fin))
module_of_interest <- top_cont_genes_fin[order(top_cont_genes_fin$GScorr),c(1,2,3,index_module)] %>% dplyr::filter(module == "MEorange")
module_of_interest_fin <- module_of_interest[abs(module_of_interest$MMorange) > 0.75,]  ### 84 genes


# Find file equivalences with genes



module_of_interest_mer <- merge(module_of_interest_fin, allGeneNames, by.x = "gene", by.y = "ENSEMBL")  ### 29 genes 
module_of_interest_mer




#Visualization of Eigengene corelations
MET <- orderMEs(cbind(mod, trait_of_int))

plotEigengeneNetworks(MET,"", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab=0.8, xLabelsAngle = 90)



#Expression plot for each cluster
datTraits$id <- rownames(datTraits)
datTraits$S <- as.factor(datTraits$S)
str(datTraits)

modules_df <- mod %>% tibble::rownames_to_column("accession_code") %>%
  dplyr::inner_join(as.data.frame(datTraits) %>%  dplyr::select(id, S), by = c("accession_code" = "id"))

ggplot(modules_df, aes(x = S, y = MElightcyan1)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) + theme_classic()



#Genes per module heatmap
make_module_heatmap(module_name = "MEdarkgreen", expression_mat = norm_X, metadata_df = colData_to_use,
                    gene_module_key = gene_module_key, module_eigengenes_df = ME_df)

{
  library(stringr)
  library(WGCNA)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(devtools)
  library(gridExtra)
  #install_version("pdp", version = "0.7.0")
  #library(pdp)
  # It allows to run the analysis on multiple threads, reducing greatly the computational time.
  options(stringsAsFactors = FALSE);
  enableWGCNAThreads()
}

# Path
counts_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli_orna_reduced/counts_MeeOrn_datasetR.txt"
file_meta <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli_orna_reduced/metadata_datasetR.txt"

# Variables to define
region <- "TL" # "TL" or "DE"
phenotypes <- c("S", "GL")
species <- c("Meeli")
median_counts <- 5
network_type <- "signed hybrid"
network_short <- "hyb"
picked_power <- 15
dsplit <- 4 # From 1 to 4


# Read counts and colData dataset 
cts <- read.csv(file = counts_file, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file = file_meta, sep = "\t", header = T)
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

###################### Filter from DDS #####################################
# If needed, extra step of filtering on the dds object
#species <- c("Meeli", "Ornatipinnis")
#region <- "DE"
#phenotype <- c("S", "GL")
#filtered_samples <- rownames(colData(dds))[colData(dds)$species %in% species &
#                                             colData(dds)$region == region &
#                                             colData(dds)$phenotype %in% phenotype]
#dds_fil <- dds[, filtered_samples]
############################################################################

# Run DESeq
dds_obj <- DESeq(dds)

# Filter by low counts
keep <- rowMedians(counts(dds_obj)) >= median_counts
dds_obj_lc <- dds_obj[keep,]

# Filter by varance genes. N is the number to keep
#dds_obj <- filter_by_variance(dds_obj, n = 5000)
n <- 2000 # Number of genes to keep
{
  row_var <- rowVars(counts(dds_obj_lc))
  top_n_indices <- order(row_var, decreasing = TRUE)[1:n]
  keep <- rep(FALSE, nrow(counts(dds_obj_lc)))
  keep[top_n_indices] <- TRUE
}

dds_obj_var <- dds_obj_lc[keep,]


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
wgcna_good_samples(dds_obj_var) 


# VST - data needs to be in sample to gene format
norm_cts <- vst(dds_obj) %>% assay() %>% t()


# Create design matrix
## Binarise categorical variables
datTraits <- binarizeCategoricalColumns(meta_fil, convertColumns = c("phenotype", "species"), includePairwise = FALSE,
                                          includeLevelVsAll = TRUE, minCount = 1)
rownames(datTraits) <- rownames(meta_fil)
colnames(datTraits) <- c("S")


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
sft <- pickSoftThreshold(norm_cts, dataIsExpr = TRUE, networkType = network_type, 
                         RsquaredCut = 0.85, powerVector = power)
wgcna_powers_plot(sft$fitIndices)

dim(norm_cts)


# Construct network
# convert matrix to numeric
norm_cts[] <- sapply(norm_cts, as.numeric)

temp_cor <- cor #assign correlation function to function, so WGCNA uses its own function and not the default cor() from R
cor <- WGCNA::cor


### MAKE UNIQUE NAME ON NETWORK
wg_TL_hyb_12_4_2000 = blockwiseModules(norm_cts,  # <= input here
                         # Network construction arguments: correlation options
                         corType="pearson",
                         
                         # == Adjacency Function ==
                         power = picked_power,     # power determined by soft threshold 
                         networkType = network_type,   # Unsigned = negative and positively correlated genes group together, if network is signed, the neg. and pos. correlated genes are stores in separated modules 
                         
                         # == Tree and Block Options ==
                         deepSplit = dsplit, # this goes between 1-4, the biger the value, more refined the clustering more clusters in total
                         pamRespectsDendro = F,
                         detectCutHeight = 0.995, 
                         minModuleSize = 10,       # min number of genes in a module
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
save(wg_TL_hyb_12_4_2000, file = paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli/_",
    region, "_", 
    network_short, "_", 
    picked_power, "_",
    dsplit, ".RData"))

# Load
load("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli_orna_reduced/DE_sig_18_2.RData")
wg_DE_sig_18_2 <- bwnet
load("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli_orna_reduced/TL_hyb_17_2.RData")
wg_TL_hyb_17_2 <- bwnet
load("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/meeli_orna_reduced/DE_hyb_9_2.RData")


load("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/TL_MeSMeGL_hyb_13_4.RData")
net <- TL_MeSMeGL_hyb_13_4
net <- wg_TL_hyb_12_4_2000
mod <- net$MEs

# Plot the dendrogram and the module colors before and after merging underneath 
plotDendroAndColors(net$dendrograms[[1]],  #order
                    cbind(net$unmergedColors, net$colors), 
                    c("unmerged", "merged"),  dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.03, guideHang = 0.05, main = "Title")

#Quantify mdule-trait correlations
moduleTraitCor <- cor(mod, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(norm_cts))

# Write dataset with correlation a p values

net <- TL_MeSMeGL_hyb_13_4
net_name <- "TL_MeSMeGL_hyb_13_4"
mod <- net$MEs

moduleTraitCor <- cor(mod, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(norm_cts))
merged_val <- merge(moduleTraitCor,moduleTraitPvalue, by = "row.names", all = TRUE )
colnames(merged_val) <- c("module", "corr", "pval")
merged_ord <- merged_val[order(abs(merged_val$corr), decreasing = T), ]

write.table(merged_ord,
            paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/group_vs_solitary/06.WGCNA/", net_name, ".txt"),
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
               colorLabels = F, colors = blueWhiteRed(50), zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = F, cex.text = 0.5,
               main = paste("Module-trait relationships"))





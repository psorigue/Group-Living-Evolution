{
  library(WGCNA)
  library(SummarizedExperiment)
  library(DESeq2)
  library(gridExtra)
  library(dplyr)
  library(ggplot2)
  options(stringsAsFactors = FALSE);
  enableWGCNAThreads()
}

# Paths
file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/counts_all_samples_datasetR.txt" 
file_metadata <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/metadata_datasetR.txt"


# Variables to define
region <- "TL" # "TL" or "DE"
samples <- "OrMuMe"
path_wgcna <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/RData/")
median_counts <- 10 # Minimum median counts

# Read counts and colData dataset 
cts <- read.csv(file = file_cts, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file = file_metadata, sep = "\t", header = T)
rownames(meta) <- meta$id

# Filter region in meta and counts
meta_fil <- meta[meta$region == region, ]
meta_fil <- meta_fil[,c("phenotype", "species")]
cts_fil <- cts[, rownames(meta_fil)]

# Categorize variables
meta_fil$phenotype <- factor(meta_fil$phenotype, levels = unique(meta_fil$phenotype)) 
meta_fil$species <- factor(meta_fil$species, levels = unique(meta_fil$species))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts_fil,
                              colData = meta_fil,
                              design = ~ phenotype + species) # phenotype or species

#Filter out samples
samples_to_remove <- "Or49DE"
dds <- dds[, !colnames(dds) %in% samples_to_remove]
meta_fil <- meta_fil[!row.names(meta_fil) %in% samples_to_remove, ]

# Run DESeq
dds_obj <- DESeq(dds)

# Filter out low counts
keep <- rowMedians(counts(dds_obj)) >= median_counts # Defined at the beginning
dds_obj_lc <- dds_obj[keep,]


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
#datTraits <- binarizeCategoricalColumns(meta_fil, convertColumns = c("phenotype", "species"), includePairwise = FALSE,
#                                        includeLevelVsAll = TRUE, minCount = 1)
#rownames(datTraits) <- rownames(meta_fil)
#colnames(datTraits) <- c("GL", "Orna", "Meeli")
#colnames(datTraits) <- c("GL", "Meeli")

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
picked_power <- 18 
network_type <- "signed"
dsplit <- 3
mch <- 0.15

# convert matrix to numeric
norm_cts[] <- sapply(norm_cts, as.numeric)

temp_cor <- cor #assign correlation function to function, so WGCNA uses its own function and not the default cor() from R
cor <- WGCNA::cor


### MAKE UNIQUE NAME ON NETWORK
Sys.time()
OrMuMe_TL = blockwiseModules(norm_cts,  # <= input here
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
                                       mergeCutHeight = mch, #correlation of 0.80 between the clusters to be merged | Default = 0.15
                                       
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
net_name <- "OrMuMe_TL"
save(OrMuMe_TL, file = paste0(path_wgcna, net_name, ".RData"))
Sys.time()
# Save datTraits
comp <- "OrMuMe"
#file_name <- paste0("datTraits_", comp)
#save(datTraits, file = paste0(path_wgcna, file_name, ".RData"))
# Save norm_cts
dataset_name <- paste0("norm_cts_", comp)
save(norm_cts, file = paste0(path_wgcna, dataset_name, ".RData"))




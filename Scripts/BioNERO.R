BiocManager::install("BioNERO")
# Reference: https://www.bioconductor.org/packages/release/bioc/vignettes/BioNERO/inst/doc/vignette_01_GCN_inference.html
library(BioNERO)
library(SummarizedExperiment)

# Paths
counts_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/counts_all_samples_datasetR.txt"
file_meta <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/04.counts/all_samples/metadata_datasetR.txt"

# Variables to define
region <- "TL" # "TL" or "DE"
phenotypes <- c("S", "GL")
species <- c("Meeli", "Multifasciatus")
median_counts <- 2
<<<<<<< HEAD
picked_power <- 18
=======
picked_power <- 12
>>>>>>> 726dee3f553a26110aa956576784fe4282e14b4f
network_type <- "signed hybrid"
network_short <- "hyb"

# Read counts and colData dataset 
cts <- read.csv(file = counts_file, sep = "\t", header = T, row.names = 1)
meta <- read.csv(file = file_meta, sep = "\t", header = T)
rownames(meta) <- meta$id

# Filter samples in meta
meta_fil <- meta[meta$region == region & meta$phenotype %in% phenotypes & meta$species %in% species, ]
meta_fil <- meta_fil[,c("phenotype", "species")]

# If applicable, remove samples
samples_to_remove <- c("Me33TL", "Mu40TL")
meta_fil_rm <- meta_fil[!(rownames(meta_fil) %in% samples_to_remove), ]

# Apply filter to counts
cts_fil <- cts[, rownames(meta_fil)]

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = cts_fil),
  colData = meta_fil
)
dim(se)

# Filter low counts
exp_filt_lc <- remove_nonexp(se, method = "median", min_exp = median_counts)

# Filter by variance
#exp_filt_var <- filter_by_variance(exp_filt_lc, n = 3000)


# Step to filter out outliars
exp_filt_zk <- ZKfiltering(exp_filt_lc, cor_method = "pearson")

# Step to normalise the data to avoid confounding artifacts (https://doi.org/10.1186/s13059-019-1700-9.)
exp_filt_pc <- PC_correction(exp_filt_lc)
#exp_filt_pc_zk <- PC_correction(exp_filt_zk)


# Heatmap
p <- plot_heatmap(exp_filt_pc, type = "samplecor", show_rownames = FALSE)
p


plot_PCA(exp_filt_pc)

sft <- SFT_fit(exp_filt_pc, net_type = "signed", cor_method = "pearson")
sft$plot

# Create network
network_type <- "signed hybrid"
region <- "TL"
network_short <- "MuMe_19samples_pcorrected"

picked_power <- 10

### MAKE UNIQUE NAME ON NETWORK
net <- exp2gcn(
  exp_filt_pc, net_type = network_type, SFTpower = picked_power, 
  cor_method = "pearson"
)
# Save network
save(net, file = paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/07.BioNERO/",
                        region, "_", 
                        network_short, "_", 
                        picked_power, ".RData"))
# Load network
load(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/07.BioNERO/",
            region, "_", 
            network_short, "_", 
            picked_power, ".RData"))

plot_dendro_and_colors(net)


plot_eigengene_network(net)


plot_ngenes_per_module(net)


MEtrait <- module_trait_cor(exp = exp_filt_pc, MEs = net$MEs)


plot_module_trait_cor(MEtrait)

plot_expression_profile(
  exp = exp_filt_pc, 
  net = net, 
  plot_module = TRUE, 
  modulename = "darkgreen",
  metadata_cols = "phenotype"
)


hubs <- get_hubs_gcn(exp_filt_pc, net)

# Get the genes from a module
net$genes_and_modules[net$genes_and_modules$Modules == "darkgreen",]

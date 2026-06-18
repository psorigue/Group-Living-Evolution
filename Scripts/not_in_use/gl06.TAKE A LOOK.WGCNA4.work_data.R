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
  library(dplyr)
}

# Files and data
region <- "TL"
comp <- "MuMe" # "MuMe" or "all"
path <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/", region, "/", comp, "/")
modules_cor <- paste0(path, "high_res_net/corr_module-variables_high_res_net.txt")
gene_sign <- paste0(path, "high_res_net/gene-module_trait_significance_high_res_net.txt")

# Load network
net_name <- "high_res_net"
path_data <- paste0(path, "RData/")
load(paste0(path_data, "datTraits.RData")) # Variable named datTraits
load(paste0(path_data, "norm_cts.RData")) # Variable named norm_cts
load(paste0(path_data, paste0(net_name, ".RData")))
# Change network variable name
net <- high_res_net
mod <- net$MEs
length(mod)


# 1.Take modules with correlation above N
cor_value <- 0.552
modules <- read.csv(modules_cor, sep = "\t", header = T)
sign_modules <- modules[abs(modules$corr_GL) > cor_value,]$module
# To know how many genes in each module
as.data.frame(table(net$colors))

# Are these modules similar? Visualization of Eigengene corelations
trait_of_int <- as.data.frame(datTraits$GL)
colnames(trait_of_int) <- "GL"
MET <- orderMEs(cbind(mod, trait_of_int))

# Visualize in PDF
pdf(paste0(path, "module_similarity.pdf"), width = 35, height = 25)
plot <- plotEigengeneNetworks(MET,"", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                              cex.lab=0.8, xLabelsAngle = 90, excludeGrey = T)
dev.off()

# If I compare the network to mid resolution, do two of the modules merge into one? In this case, do enrichment of all genes
# Obtain genes from the significan modules
genes <- read.csv(gene_sign, sep = "\t", header = T)
genes_in_modules <- genes[genes$module %in% sign_modules,]
# Load genes mid_res_net
gene_sign_mid_res_net <- read.csv(file = paste0(path, "mid_res_net/gene-module_trait_significance_mid_res_net.txt"), sep = "\t", header = T)
# Merge module names of mid_res
genes_both <- merge(genes_in_modules[,c("gene", "module")], gene_sign_mid_res_net[,c("gene", "module")], by = "gene")
colnames(genes_both) <- c("gene", "module_name_high", "module_name_mid")

genes_both_sorted <- genes_both[order(genes_both$module_name_high),]
write.table(genes_both_sorted, file = paste0(path, "sign_modules_genes_high_mid.txt"), col.names = T, sep = "\t", quote = F, row.names = F)


# Is there a group of genes that are significant for MuMe and OrMuMe?
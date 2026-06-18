# This script assigns the phenotype (group vs solitary) to each species 
# and creates a tree with only the species that have phenotype information. 
# This is needed for the trait reconstruction at node analysis. 
# Trait assignment is based on the median Nmax value across all deployments 
# for each species. Species with median Nmax > 1 are classified as group, 
# while those with median Nmax = 1 are classified as solitary.

# Load libraries
library(dplyr) # version 1.1.4
library(ape) # version 5.8-1

# Set paths and files
home <- path.expand("~")
file_grouping <- file.path(home, "01.Phenotypic_characterization", "median_Nmax.csv")
tree_tan_file <- file.path(home, "01.Phenotypic_characterization", "tree_tan.nwk") # Tree from Ronco et al. 2021
file_out <- file.path(home, "01.Phenotypic_characterization", "template_grouping.txt")
tree_out <- file.path(home, "01.Phenotypic_characterization", "tree_grouping.nwk")


# INDEX
# 1. Assign phenotype to each species
# 2. Create tree with species that have phenotype information

# 1. Assign phenotype to each species
# -----------------------------------

# Load file
dt <- read.csv(file_grouping)

# Turn into 2 categories
dt_mod <- dt %>%
  mutate(grouping = case_when(
    median_Nmax > 1 ~ "G",    
    median_Nmax == 1 ~ "S")) %>%
  dplyr::select(spp, grouping)

# Change column names and write file
colnames(dt_mod) <- c("species", "grouping")
write.table(dt_mod, file = file_out, sep = "\t", quote = F, row.names = F)


# 2. Create tree with species that have phenotype information
# -----------------------------------------------------------

# Read tree
tree_tan <- read.nexus(tree_tan_file)

# Get species with phenotype information and present in the tree
spp <- intersect(tree_tan$tip.label, dt_mod$species)

# Prune spp without pheno info
tree_pru <- keep.tip(tree_tan, spp)

# Write pruned tree
write.tree(tree_pru, file = tree_out)

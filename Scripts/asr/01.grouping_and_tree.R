library(dplyr)
library(ape)

file_out <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.1.Phenotypic_characterization/template_grouping.txt"
tree_out <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.1.Phenotypic_characterization/tree_grouping.nwk"

# 1. Assign phenotype to each species

# Load file
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.1.Phenotypic_characterization/median_Nmax.csv"
dt <- read.csv(file)

# Turn into 2 categories
dt_mod <- dt %>%
  mutate(grouping = case_when(
    median_Nmax > 1 ~ "G",    
    median_Nmax == 1 ~ "S")) %>%
  dplyr::select(spp, grouping)

colnames(dt_mod) <- c("species", "grouping")

write.table(dt_mod, file = file_out, sep = "\t", quote = F, row.names = F)



# 2. Create tree with species that have phenotype information
tree_tan <- read.nexus("//files1.igc.gulbenkian.pt/folders/ANB/Pol/map_traits/tree_tan/tree_tan_ed")

# Prune spp without pheno info
spp <- intersect(tree_tan$tip.label, dt_mod$species)
tree_pru <- keep.tip(tree_tan, spp)

# Write tree
write.tree(tree_pru, file = tree_out)

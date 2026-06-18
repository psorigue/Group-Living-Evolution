
# Copied from file in the repository and modified
MBASR.directory <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/ASR_and_traits_correlation/ASR/MBASR/repository/"
setwd(MBASR.directory)
source("MBASR.load.functions.R")

file_tree <- "Group_size/tree.nwk"
file_trait <- "Group_size/2cat/template_gs_2cat.txt"
#file.name.ordered.characters="multitrait/ordered.characters.txt" # Only necessary in multitrait
file_plot_settings <- "plot.settings.txt"


n_samples <- 5000

# Single trait
MBASR(file_tree,file_trait,file_plot_settings,character.type = "unordered", n.samples = n_samples)

# Multi trait
#MBASR.multi.trait(file_tree,file_trait,file_plot_settings, file.name.ordered.characters, n.samples = n_samples)
#Optional – replot without reanalysis… after editing the plot settings file
replot.multi.trait(file.name.tree, file.name.multi.trait.csv, file.name.plot.settings)




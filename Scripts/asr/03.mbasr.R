
# Copied from file in the repository and modified
MBASR.directory <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.2.Node_reconstruction/MBASR/repository/"
setwd(MBASR.directory)
source("MBASR.load.functions.R")

# These files must be inside folder 'input' to run function MBASR
file_tree <- "tree_grouping.nwk"
file_trait <- "template_grouping.txt"
file_plot_settings <- "plot.settings.txt"


n_samples <- 10000

# Single trait
MBASR(file_tree,file_trait,file_plot_settings,character.type = "unordered", n.samples = n_samples)

# Multi trait
#MBASR.multi.trait(file_tree,file_trait,file_plot_settings, file.name.ordered.characters, n.samples = n_samples)
#Optional – replot without reanalysis… after editing the plot settings file
replot.multi.trait(file.name.tree, file.name.multi.trait.csv, file.name.plot.settings)




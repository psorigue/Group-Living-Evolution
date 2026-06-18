# This script performs model selection for ancestral state reconstruction using MBASR
# (Heritage S., 2021; https://doi.org/10.1101/2021.01.10.426107) 
# It uses a bayesian framework to calculate the probabilities of each state at nodes
# using a single-rate model of evolution.

# 1. Set paths and files
home <- path.expand("~")
MBASR.directory <- file.path(home, "02.Trait Reconstruction at Node", "MBASR", "repository")
setwd(MBASR.directory)
source("MBASR.load.functions.R") # This file contains the functions to run MBASR. 
# The source script and additional files needed to run MBASR are not contained 
# in this study's repository, but can be obtained from https://doi.org/10.1101/2021.01.10.426107

# 2. Define files. These files must be inside folder 'input' to run function MBASR
file_tree <- "tree_grouping.nwk"
file_trait <- "template_grouping.txt"
file_plot_settings <- "plot.settings.txt"

# 3. Set number of samples for the MCMC. For binary traits, 5000 samples should be enough
n_samples <- 10000

# 4. Run Single-trait MBASR
MBASR(file_tree,file_trait,file_plot_settings,character.type = "unordered", n.samples = n_samples)

# Results and tree plot are saved in the output folder.

# This script performs model selection for ancestral state reconstruction
# at nodes, using the ace function (phytools) with the single rate (ER)
# and multiple rates (ARD, 2 rates of evolution) models.
# It then reconstructs the ancestral states at nodes using the best-fitting model.
# Finally, it plots the tree with the reconstructed states at nodes.

# Load libraries
library(phytools) # version 2.5-2
library(ape) # version 5.8-1

# Set paths and files
home <- path.expand("~")
tree_file <- file.path(home, "01.Phenotypic characterization", "tree_grouping.nwk")
pheno_data_file <- file.path(home, "01.Phenotypic characterization", "template_grouping.txt")
folder_out <- file.path(home, "02.Trait_Reconstruction at Node/")


# INDEX
# 1. Fit models ER and ADR
# 2. Compare models using AIC and likelihood ratio test
# 3. Obtain marginal ancestral states at nodes
# 4. Plot tree with states probabilities at nodes


# 1. Fit models ER and ADR
# -------------------------

# Read tree and phenotype files
tree <- read.tree(tree_file)
pheno <- read.csv(file = pheno_data_file, sep = "\t", header = T)

# Map states on the tree
vector_states <- setNames(pheno$grouping, pheno$species)
vector_states <- vector_states[tree$tip.label]

# Fit models and reconstruct ancestral states at nodes
fit_ER <- ace(vector_states, tree,model="ER",type="discrete")
fit_ARD <- ace(vector_states, tree,model="ARD",type="discrete")

# Obtain Log-likelihood
logL_ER  <- fit_ER$loglik
logL_ARD <- fit_ARD$loglik


# 2. Compare models using AIC and likelihood ratio test
# -----------------------------------------------------

# Model comparison -> AIC
k_ER  <- 1 # Number of parameters (rates of evolution)
k_ARD <- 2

AIC_ER  <- 2*k_ER  - 2*logL_ER
AIC_ARD <- 2*k_ARD - 2*logL_ARD

delta_AIC <- AIC_ARD - AIC_ER
delta_AIC 

# 4. Model comparison -> Likelihood ratio test
LR <- 2 * (logL_ARD - logL_ER)
pval <- pchisq(LR, df = 1, lower.tail = F)

# Write results
results <- data.frame(
  Model        = c("ER", "ARD"),
  LogLik       = c(fit_ER$loglik, fit_ARD$loglik),
  AIC          = c(AIC_ER, AIC_ARD),
  Rates        = I(list(fit_ER$rates, fit_ARD$rates)) # Rates of evolution
)
# Add model comparison results
results <- rbind(results, c("Model_comparison_pval", pval, "", ""))

# Write results to file
write.table(results, file = paste0(folder_out, "results.txt"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# 3. Obtain Marginal Ancestral States at each node, and write it
# --------------------------------------------------------------
# Because we have support for the simpler model (ΔAIC −0.29, p-val = 0.13; see Results),
# we will use the ER model for ancestral state reconstruction at nodes.

# Get probabilities at nodes
probabilities_at_nodes <- fit_ER$lik.anc
pan <- as.data.frame(probabilities_at_nodes)

# Write probabilities at nodes to file
write.table(pan, file = paste0(folder_out, "probab_nodes_ER.txt"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)


# 4. Plot tree with states probabilities at nodes
# -----------------------------------------------

# Plot tree with node numbers
plotTree(tree,node.numbers=T, fsize = 0.3)

# Plot states at nodes
my_colors <- c("#4472C4", "#BF9000")
mapped_tree <- make.simmap(tree, x=vector_states)
states <- getStates(mapped_tree, "tips")
cols <- setNames(my_colors, sort(unique(states)))

# Tree ER
plotTree(mapped_tree,type="fan",fsize=0.6,ftype="i")
nodelabels(node=1:mapped_tree$Nnode+Ntip(mapped_tree),
           pie=fit_ER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(states,sort(unique(states))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(mapped_tree)),fsize=0.8)

# The divergence between the species of this study
# (N. multifasciatus, N. meeli and L. ornatipinnis)
# occurs at nodes 157 and 158. Therefore, this node
# will be the focus of the trait reconstruction at node analysis. 


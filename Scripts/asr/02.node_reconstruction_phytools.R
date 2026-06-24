library(phytools)
library(ape)


tree_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.1.Phenotypic_characterization/tree_grouping.nwk"
pheno_data_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.1.Phenotypic_characterization/template_grouping.txt"
folder_out <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/00.2.Node_reconstruction/"
dir.create(folder_out, showWarnings = F)

# Read files
tree <- read.tree(tree_file) # tree
pheno <- read.csv(file = pheno_data_file, sep = "\t", header = T)

# Map states on the tree
vector_states <- setNames(pheno$grouping, pheno$species)
vector_states <- vector_states[tree$tip.label]


######## MODELS ###########

# 1. Fit models and reconstruct ancestral states at nodes
fit_ER <- ace(vector_states, tree,model="ER",type="discrete")
fit_ARD <- ace(vector_states, tree,model="ARD",type="discrete")

# 2. Log-likelihood
logL_ER  <- fit_ER$loglik
logL_ARD <- fit_ARD$loglik

# 3. Rates
fit_ER$rates
fit_ARD$rates

# 4. Model comparison -> AIC
k_ER  <- 1 # Number of parameters (rates)
k_ARD <- 2

AIC_ER  <- 2*k_ER  - 2*logL_ER
AIC_ARD <- 2*k_ARD - 2*logL_ARD

delta_AIC <- AIC_ARD - AIC_ER
delta_AIC

# 4. Model comparison -> Likelihood ratio test
LR <- 2 * (logL_ARD - logL_ER)
pval <- pchisq(LR, df = 1, lower.tail = F)

# 5. Write results
results <- data.frame(
  Model        = c("ER", "ARD"),
  LogLik       = c(fit_ER$loglik, fit_ARD$loglik),
  AIC          = c(AIC_ER, AIC_ARD),
  Rates        = I(list(fit_ER$rates, fit_ARD$rates))
)

results <- rbind(results, c("Model_comparison_pval", pval, "", ""))

write.table(results, file = paste0(folder_out, "results.txt"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# 5. Obtain Marginal Ancestral States at each node, and write it
probabilities_at_nodes <- fit_ER$lik.anc
pan <- as.data.frame(probabilities_at_nodes)
write.table(pan, file = paste0(folder_out, "probab_nodes_ER.txt"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)
probabilities_at_nodes <- fit_ARD$lik.anc
pan <- as.data.frame(probabilities_at_nodes)
write.table(pan, file = paste0(folder_out, "probab_nodes_ARD.txt"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)





######## PLOT TREES #############

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

# Tree ARD
plotTree(mapped_tree,type="fan",fsize=0.6,ftype="i")
nodelabels(node=1:mapped_tree$Nnode+Ntip(mapped_tree),
           pie=fit_ARD$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(states,sort(unique(states))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(mapped_tree)),fsize=0.8)







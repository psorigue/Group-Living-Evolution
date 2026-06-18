library(ape)
library(phytools)

tree <- read.tree("//files1.igc.gulbenkian.pt/folders/ANB/Pol/map_traits/tree_tan/tree_tan_ed.nwk")
tree_pruned <- keep.tip(tree, c("Neomul", "Neomee", "Lamorn"))


# Change name of facultative species
tree_pruned$tip.label <- c("Multi_Solitary", "Orna_Solitary", "Meeli_Solitary")
plot(tree_pruned)

# Add nodes

new_tip1 <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "Multi_Group", edge.length = 1.0, Nnode = 1)
new_tip2 <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "Meeli_Group", edge.length = 1.0, Nnode = 2)

class(new_tip1) <- "phylo"
class(new_tip2) <- "phylo"

updated_tree <- bind.tip(tree_pruned, tip.label = "Multi_Group", where = 1, position = 1, edge.length = 1)
updated_tree <- bind.tip(updated_tree, tip.label = "Meeli_Group", where = 4, position = 1, edge.length = 1)


plot(updated_tree)


# Write phylo

write.nexus(updated_tree, file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_vs_Solitary/10.phydget/phylogeny.nex", translate = T)

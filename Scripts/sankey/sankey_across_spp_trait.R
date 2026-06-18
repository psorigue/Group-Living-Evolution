# This script includes only the modules that have correlation with the trait > 0.5 AND corr with
# trait is higher than corr with spp

#### Network with networkD3 package
library("dplyr")
#install.packages("networkD3")
library(networkD3)
library("curl")
library("scales")


original <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/Sankey_all_vs_MuMe/sign_modules_genes_mid_trait_corr.txt",
                     header = T,sep = "\t")

sankey_data <- as.data.frame( original %>%
                                group_by(source = MuMe, target = OrMuMe) %>%
                                summarise(value = n(), .groups = "drop"))
sankey_data$target <- gsub(sankey_data$target, pattern = "ME", replacement = "")

nodes <- data.frame(name=c(as.character(sankey_data$source), as.character(sankey_data$target)) %>% unique(), stringsAsFactors = F)

sankey_data$IDsource <- match(sankey_data$source, nodes$name)-1 
sankey_data$IDtarget <- match(sankey_data$target, nodes$name)-1
sankey_data$value <- as.numeric(sankey_data$value)

# Scale values
#sankey_data$value <- rescale(sankey_data$value, to = c(5, 50))

p <- sankeyNetwork(Links = sankey_data, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", sinksRight = T,
                   width = 700, nodeWidth = 30, fontSize = 12, nodePadding=10, height = 900)
p




library("ggsankey")
library(ggplot2)
library(stringr)

# Transform data set
original$OrMuMe <- gsub(original$OrMuMe, pattern = "ME", replacement = "")
df <- original[c("MuMe", "OrMuMe")] %>%
  make_long(MuMe, OrMuMe)
df

# Order nodes according to the correlation with the trait
## 1st column 
### Take order correlation
MuMe_order <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/MuMe/mid_res_net/corr_module-variables_mid_res_net.txt",
                       sep = "\t",
                       header = T)
MuMe_order_ordered <- MuMe_order[order(-abs(MuMe_order[,2])), 1:2]
### Filter this dataset wiht only the values that appear in the df
levels_MuMe <- unique(original$MuMe)
MuMe_order_ordered_fil <- MuMe_order_ordered[MuMe_order_ordered$module %in% levels_MuMe,]
## 2nd column (mid resolution)
### Take order correlation
OrMuMe_order <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/all/mid_res_net/corr_module-variables_mid_res.txt",
                         sep = "\t",
                         header = T)
OrMuMe_order$module <- gsub(OrMuMe_order$module, pattern = "ME", replacement = "")
OrMuMe_order_ordered <- OrMuMe_order[order(-abs(OrMuMe_order[,2])), 1:2]
### Filter this dataset wiht only the values that appear in the df
levels_OrMuMe <- unique(original$OrMuMe)
OrMuMe_order_ordered_fil <- OrMuMe_order_ordered[OrMuMe_order_ordered$module %in% levels_OrMuMe,]


# Change order nodes. Help: https://github.com/davidsjoberg/ggsankey/issues/16
df_test <- df %>%
  mutate(node = factor(node, levels = rev(c(MuMe_order_ordered_fil$module, "NSa", OrMuMe_order_ordered_fil$module, "NSb")))) %>%
  mutate(color = str_remove(node, "^ME"))


# Help: https://rpubs.com/techanswers88/sankey-with-own-data-in-ggplot
ggplot(df_test, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = factor(node),
                    label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none")


#### Network with networkD3 package

library("dplyr")
#install.packages("networkD3")
library(networkD3)
library("curl")
library("scales")


original <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/MuMe/Sankey/sign_modules_genes_high_mid.txt",
                     header = T,sep = "\t")

sankey_data <- as.data.frame( original %>%
  group_by(source = module_name_high, target = module_name_mid) %>%
  summarise(value = n(), .groups = "drop"))
sankey_data$target <- gsub(sankey_data$target, pattern = "ME", replacement = "")

nodes <- data.frame(name=c(as.character(sankey_data$source), as.character(sankey_data$target)) %>% unique(), stringsAsFactors = F)

sankey_data$IDsource <- match(sankey_data$source, nodes$name)-1 
sankey_data$IDtarget <- match(sankey_data$target, nodes$name)-1
sankey_data$value <- as.numeric(sankey_data$value)

# Scale values
sankey_data$value <- rescale(sankey_data$value, to = c(5, 50))

p <- sankeyNetwork(Links = sankey_data, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", sinksRight = T,
                   width = 700, nodeWidth = 30, fontSize = 12, nodePadding=10, height = 900)
p


# Try changing colors
# Try adding low resolution and/or comparing to OrMuMe


#### Plot ggplot Sankey


# Try ggplot sankey: https://r-charts.com/flow/sankey-diagram-ggplot2/

# Help pages:
  #https://rpubs.com/techanswers88/sankey-with-own-data-in-ggplot
  #https://r-charts.com/flow/sankey-diagram-ggplot2/

#devtools::install_github("davidsjoberg/ggsankey")
library("ggsankey")
library(ggplot2)
library(dplyr)
library(stringr)

# Transform data set
original$module_name_mid <- gsub(original$module_name_mid, pattern = "ME", replacement = "")
df <- original[c("module_name_high", "module_name_mid")] %>%
  make_long(module_name_high, module_name_mid)
df



# Order nodes according to the correlation with the trait
## 1st column (high resolution)
### Take order correlation
high_order <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/MuMe/high_res_net/corr_module-variables_high_res_net.txt",
                       sep = "\t",
                       header = T)
high_order_ordered <- high_order[order(-abs(high_order[,2])), 1:2]
### Filter this dataset wiht only the values that appear in the df
levels_high <- unique(original$module_name_high)
high_order_ordered_fil <- high_order_ordered[high_order_ordered$module %in% levels_high,]
## 2nd column (mid resolution)
### Take order correlation
mid_order <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/group_vs_solitary/06.WGCNA/TL/MuMe/mid_res_net/corr_module-variables_mid_res_net.txt",
                      sep = "\t",
                      header = T)
mid_order$module <- gsub(mid_order$module, pattern = "ME", replacement = "")
mid_order_ordered <- mid_order[order(-abs(mid_order[,2])), 1:2]
### Filter this dataset wiht only the values that appear in the df
levels_mid <- unique(original$module_name_mid)
mid_order_ordered_fil <- mid_order_ordered[mid_order_ordered$module %in% levels_mid,]


# Change order nodes. Help: https://github.com/davidsjoberg/ggsankey/issues/16
df_test <- df %>%
  mutate(node = factor(node, levels = rev(c(high_order_ordered_fil$module, mid_order_ordered_fil$module)))) %>%
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


# If I want the nodes colored
ggplot(df_test, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = color,       # use cleaned color
                    label = node)) +    # still show "MEred" etc.
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_identity() +    # interpret fill values as real colors
  theme_sankey(base_size = 16) +
  theme(legend.position = "none")

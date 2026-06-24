library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)


df <- read.csv(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/Group-Living_Evolution_paper/Figures/Fig3_extras/dataset_Fig3_general.txt"),
               header = T, sep = "\t")

filter <- "KEGG pathways"

df_f <- df[df$X == filter,c("region", "N.multifasciatus", "N.meeli", "Common")]

# Convert to long format
df_long <- df_f %>%
  pivot_longer(
    cols = c(N.multifasciatus, N.meeli, Common),
    names_to = "Category",
    values_to = "GO_terms"
  )

df_long

df_long$Category <- factor(
  df_long$Category,
  levels = c("N.multifasciatus", "N.meeli", "Common")
)
df_long <- df_long %>%
  mutate(region = factor(region, levels = c(
    "Telencephalon",
    "Diencephalon"
  )))

p <- ggplot(df_long,
            aes(x = region, y = GO_terms, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(
    values = c(
      "N.multifasciatus" = "#8B4500",
      "N.meeli" = "#104E8B",
      "Common" = "gold2"
    ),
    labels = c(
      expression(italic("N. multifasciatus")),
      expression(italic("N. meeli")),
      "Common"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none") 



type <- "KEGG"
path_plot <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/Group-Living_Evolution_paper/Figures/Fig3_general_", type, ".pdf")
ggsave(plot = a, filename = path_plot, device = "pdf", width = 4, height = 2)

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

df <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/Group-Living_Evolution_paper/Figures/dataset_Fig4.txt",
               header = T, sep = "\t")

# Order clusters by size
df2 <- df %>%
  arrange(number) %>%
  mutate(name = factor(name, levels = name))

# Convert to long format
df_long <- df2 %>%
  pivot_longer(
    cols = c(N.multifasciatus, N.meeli, Common),
    names_to = "Category",
    values_to = "GO_terms"
  )

# Plot
p <- ggplot(df_long,
            aes(x = name, y = GO_terms, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c(
      "N.multifasciatus" = "#8B4500",
      "N.meeli" = "#104E8B",
      "Common" = "gold2"
    ),
    breaks = c("N.multifasciatus", "N.meeli", "Common")
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    labels = scales::label_number(accuracy = 1)
  ) +
  coord_flip() +
  scale_x_discrete(labels = \(x) str_wrap(x, width = 30)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p
path_plot <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/Group-Living_Evolution_paper/Figures/barplot_Fig4.pdf")
ggsave(plot = p, filename = path_plot, device = "pdf", width = 5, height = 5)


file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/own_counts/counts_clean.csv"

dt <- read.csv(file, header = T)
rownames(dt) <- dt$X

dt$X <- NULL



t <- as.data.frame(t(dt))

t$spp <- substring(rownames(t), 1, 6)

gene_id <- "109201892"
gene <- paste0("X", gene_id)

t[,c("spp", gene)]


# Plot
plot_df <- t[, c("spp", gene)]
colnames(plot_df) <- c("spp", "expr")

plot_df$sample <- rownames(plot_df)

plot_df <- plot_df[order(plot_df$spp), ]
plot_df$sample <- factor(plot_df$sample, levels = plot_df$sample)

library(ggplot2)

ggplot(plot_df, aes(x = sample, y = expr, color = spp)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Samples", y = "Gene expression")




plot_df2 <- aggregate(expr ~ spp, data = plot_df, FUN = mean)
ggplot(plot_df2, aes(x = spp, y = expr, fill = spp)) +
  geom_col() +
  theme_bw()





subset_df <- plot_df[plot_df$spp %in% c("Lamkun", "Neomul", "Lepatt"), ]

ggplot(subset_df, aes(x = spp, y = expr, color = spp)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Species", y = "Gene expression")

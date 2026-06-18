library(dplyr)
library(tidyr)

folder <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Group_size/work_raw_data/"


# 1. Process data

# Read dataset
file <- paste0(folder, "Community2018_10_23.csv")
dt <- read.csv(file = file, header = T)

#Remove cols of all zeros
dt0 <- dt[,colSums(dt) != 0]
# Turn 0s into NAs
dt_na <- dt0
dt_na[dt_na == 0] <- NA

# Filter species of radiation
species <- scan(file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/dataset-array/spp/IDar_treespp.txt", sep = " ", what = "character")
species_to_keep <- intersect(species, names(dt_na))
dt_spp <- dt_na[,c("transect_no", "cam_no", "PicID", species_to_keep)]

# How many species appear?
length(species_to_keep)

# Load function to calculate mode
get_mode <- function(x) {
  x <- x[!is.na(x)]  # remove NAs
  if (length(x) == 0) return(NA)
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}


# 2. Analyse raw column

# 2.1. Max
max_values <- sapply(dt_spp, max, na.rm = TRUE)
vals <- max_values[4:length(max_values)]
hist(vals, breaks = 100,
     main = "Max value",
     labels = T)
names(vals[order(vals, decreasing = T)]) # Returns the species names orderded by max value

# 2.2. Min
min_values <- sapply(dt_spp, min, na.rm = TRUE)
vals <- min_values[4:length(min_values)]

# 2.3. Median
median_values <- sapply(dt_spp, median, na.rm = T)
median_values_red <-  median_values[4:length(median_values)]
hist(median_values_red, breaks = 10, labels = T)
median_values_red[order(median_values_red, decreasing = T)]

# 1.4. Mode
modes <- sapply(dt_spp, get_mode)
modes_red <-  modes[4:length(modes)]
hist(modes_red, breaks = 10)



# 3. Group by deployment (using an Nmax per deployment)
folder_Nmax <- paste0(folder, "Nmax_depl/")

# 3.1. Process data
# Create long data frame and compute Nmax per species, camera, and transect
dt_long <- dt_spp %>%
  select(-PicID) %>%
  pivot_longer(-c(cam_no, transect_no), names_to = "species", values_to = "count") %>%
  group_by(cam_no, transect_no, species) %>%
  summarise(
    max_count = if (all(is.na(count))) NA_integer_ else max(count, na.rm = TRUE),
    .groups = "drop"
  ) 


# 3.2. Write data
# Write max counts file
dt_wide <- dt_long %>%
  pivot_wider(names_from = species, values_from = max_count) 

write.table(dt_wide,
            file = paste0(folder_Nmax, "Nmax_transect-camera.txt"),
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)


# Write max counts file grouped by transect
grouped_by_transect <- dt_wide %>% 
  select(-cam_no) %>%
  pivot_longer(-transect_no, names_to = "species", values_to = "Nmax") %>%
  group_by(transect_no, species) %>%
  summarise(
    max_count = if (all(is.na(Nmax))) NA_integer_ else max(Nmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(species)
grouped_by_transect_wide <- grouped_by_transect %>%
  pivot_wider(names_from = species, values_from = max_count) 

write.table(grouped_by_transect_wide,
            file = paste0(folder_Nmax, "Nmax_grouped_by_transect.txt"),
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

# Stats Nmax per transect
stats_overall_Nmx_transects <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(
    occur_transects = sum(!is.na(max_count)),
    median_Nmax = median(max_count, na.rm = TRUE),
    min_Nmax = min(max_count, na.rm = TRUE),
    max_Nmax = max(max_count, na.rm = TRUE),
    .groups = "drop"
  )

write.table(stats_overall_Nmx_transects,
          file = paste0(folder_Nmax, "stats_overall_Nmx_transects.txt"),
          row.names = F,
          sep = "\t",
          quote = F)

# Dispersion
df_out <- data.frame()
for (spp in species) {
  
  IQR <- quantile(grouped_by_transect[grouped_by_transect$species == spp,]$max_count, 0.75, na.rm = TRUE) - quantile(grouped_by_transect[grouped_by_transect$species == spp,]$max_count, 0.25, na.rm = TRUE)
  sd <- sd(grouped_by_transect[grouped_by_transect$species == spp,]$max_count, na.rm = T)
  vec <- c(spp, IQR, sd)
  df_out <- rbind(df_out, vec)
  
}
colnames(df_out) <- c("species", "IQR", "sd")

write.table(df_out,
            file = paste0(folder_Nmax, "dispersion.txt"),
            sep = "\t",
            row.names = F,
            quote = F)





# 3.2. Plot data
library(ggplot2)
library(dplyr)
library(tidyr)


# Obtain max Nmax, min Nmax and median Nmax per species and transect
folder_Nmax <- paste0(folder, "Nmax_depl/")
df <- read.csv(file = paste0(folder_Nmax, "Nmax_transect-camera.txt"), sep = "\t",
               header = T)
stats_overall_Nmx_transects <- read.csv(file = paste0(folder_Nmax, "stats_overall_Nmx_transects.txt"), 
                                        sep = "\t",
                                        header = T)

grouped_by_transect <- df %>% 
  select(-cam_no) %>%
  pivot_longer(-transect_no, names_to = "species", values_to = "Nmax") %>%
  group_by(transect_no, species) %>%
  summarise(
    max_count = if (all(is.na(Nmax))) NA_integer_ else max(Nmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(species)


# Histogram of maximum Nmax
# One max Nmax per species (across transects)
max_per_species <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(max_Nmax = max(max_count, na.rm = TRUE), .groups = "drop")
# Histogram of max Nmax values
ggplot(max_per_species, aes(x = max_Nmax)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_bw() +
  labs(title = "Histogram of maximum Nmax per species",
       x = "Max Nmax (per species across transects)", 
       y = "Frequency")

# Histogram of median Nmax
# One median Nmax per species (across transects)
median_per_species <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(median_Nmax = median(max_count, na.rm = TRUE), .groups = "drop")
# Histogram of median Nmax values
ggplot(median_per_species, aes(x = median_Nmax)) +
  geom_histogram(bins = 32, fill = "darkgreen", color = "white") +
  theme_bw() +
  labs(title = "Histogram of median Nmax per species",
       x = "Median Nmax (per species across transects)", 
       y = "Frequency")

# Histogram of mean Nmax
# One mean Nmax per species (across transects)
mean_per_species <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(mean_Nmax = mean(max_count, na.rm = TRUE), .groups = "drop")
# Histogram of mean Nmax values
ggplot(mean_per_species, aes(x = mean_Nmax)) +
  geom_histogram(bins = 30, fill = "purple", color = "white") +
  theme_bw() +
  labs(title = "Histogram of mean Nmax per species",
       x = "Mean Nmax (per species across transects)", 
       y = "Frequency")

# Histogram of mode Nmax
get_mode <- function(x) {
  x <- x[!is.na(x)]  # remove NAs
  if (length(x) == 0) return(NA)
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}
# Compute mode Nmax per species
mode_per_species <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(mode_Nmax = get_mode(max_count), .groups = "drop")
# Histogram of mode Nmax values
ggplot(mode_per_species, aes(x = mode_Nmax)) +
  geom_histogram(bins = 50, fill = "darkcyan", color = "white") +
  theme_bw() +
  labs(title = "Histogram of mode Nmax per species",
       x = "Mode Nmax (per species across transects)", 
       y = "Frequency")

# Histogram of min Nmax
# One minimum Nmax per species (across transects)
min_per_species <- grouped_by_transect %>%
  group_by(species) %>%
  summarise(min_Nmax = min(max_count, na.rm = TRUE), .groups = "drop")
# Histogram of minimum Nmax values
ggplot(min_per_species, aes(x = min_Nmax)) +
  geom_histogram(bins = 5, fill = "orange", color = "white") +
  theme_bw() +
  labs(title = "Histogram of minimum Nmax per species",
       x = "Minimum Nmax (per species across transects)", 
       y = "Frequency")


# Gather mode and median per species and write file
dt_mer <- merge(mode_per_species, median_per_species, by = "species")
write.csv(dt_mer, file = paste0(folder_Nmax, "mode_median.csv"),
          quote = F,
          row.names = F)

# Extra plots

# Distribution of Nmax across transects (per species)
ggplot(grouped_by_transect, 
       aes(x = reorder(species, -max_count, FUN = median), y = max_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Nmax distribution across transects per species",
       x = "Species (ordered by median Nmax)", 
       y = "Nmax")

# Summary stats (min, median, max per species)
ggplot(stats_overall_Nmx_transects, 
       aes(x = reorder(species, median_Nmax), y = median_Nmax)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = min_Nmax, ymax = max_Nmax), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(title = "Per-species Nmax (median, min, max across transects)",
       x = "Species", y = "Nmax")

# Heatmap of Nmax (species × transect)
ggplot(grouped_by_transect, aes(x = transect_no, y = species, fill = max_count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(title = "Heatmap of Nmax per species and transect",
       x = "Transect", y = "Species", fill = "Nmax")

# Dispersion (IQR vs SD)
ggplot(df_out, aes(x = as.numeric(IQR), y = as.numeric(sd), label = species)) +
  geom_point(color = "red") +
  ggrepel::geom_text_repel() +
  theme_bw() +
  labs(title = "Dispersion of Nmax per species",
       x = "Interquartile Range (IQR)", y = "Standard Deviation")



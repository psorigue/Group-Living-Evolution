

file_cts <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/cichl_exp_init.csv"
cts <- read.csv(file_cts)
file_gs <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/new_phenotypes/Group_size/Nmax_depl/mode_median.csv"
gs <- read.csv(file_gs)



# 1. The dataset contains a column "phenotype", based on pair-bonding. I will update the info and change the name to "pair_bonding".
# Remove first column called "X" (useless)
cts <- cts[,2:length(cts)]
# Change col name
names(cts)[4] <- "pair_bonding"
# Read updated pair_bonding info
pb_info <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/All_Pathways/correlation_analysis/templates_phenotype/template_pb.txt",
                    sep = "\t")
# Merge phenotype column updated
cts_mer <- merge(cts, pb_info, by = "spp", all.x = T)
# Copy updated column to existing one
cts_mer$pair_bonding <- cts_mer$phenotype
# Remove last column
cts_mer_fin <- cts_mer[,1:(length(cts_mer) - 1)]

# 2. Join info on group size: median and mode Nmax
colnames(gs) <- c("spp", "modeNmax", "medianNmax")
cts_mer_fin_mer <- merge(cts_mer_fin, gs, by = "spp", all.x = T)
# Reorder columns
n <- ncol(cts_mer_fin_mer)
cts_mer_fin_mer_arr <- cts_mer_fin_mer[, c(1:3, 5, 4, (n-1):n, 6:(n-2))]

# Write expression. File name expression_phenotype_cts.csv
write.csv(cts_mer_fin_mer_arr, 
          file = "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Expression/expression_phenotype_cts.csv",
          quote = F,
          row.names = F)



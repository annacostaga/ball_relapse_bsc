## TIGER: Transcription Inference using Gene Expression and Regulatory data

# Load libraries
library(decoupleR)
library(netZooR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr) 

#load data and get collectri
counts <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/input/average_filtered_counts.csv", sep = ",", row.names = 1)

net <- read.csv("GRN_Threshold_1.csv")
#convert GRN to adjency matrix
net_adj <- el2adj(net)

#RUN TIGER

res <- netZooR::tiger(counts,net_adj)

# Save RData
save.image("TIGER/tiger_res_average.RData")

#convert adj matrix to dataframes
# Convert to data frame
GRN <- melt(res$W)


# Rename columns
colnames(GRN) <- c("target", "TF", "weight")

# Remove zero-weight edges 
GRN <- GRN[GRN$weight != 0, ]

hist(GRN$weight,breaks = 1000, xlim = c(-3000,4000),main = "Distribution of weights", col = "lightblue")

# grouped boxplot of important TFs to check variance across samples
TF_activity_sub <- TF_activity[which (TF_activity$TF %in% c("PAX5", "EBF1", "BCL6","IKZF3","IKZF1","PRDM1","IRF4","POU2F2","SPI1","RUNX3","SPIB","TCF3","XBP1")),]

# Create a new column for sample groups (assuming sample names follow "HSC", "HSC_1", etc.)
TF_activity_sub <- TF_activity_sub %>%
  mutate(sample_group = str_extract(sample, "^[^_]+"))  # Extracts common prefix before underscore

# Plot with facet_wrap to separate groups
ggplot(TF_activity_sub, aes(x = TF, y = activity)) + 
  geom_boxplot() + 
  ylim(0,1000) +
  facet_wrap(~sample_group, scales = "free_x") +  # Facet by sample group
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#save GRN and TF activity
write.csv(GRN, file = "TIGER/GRN_B_cells_avg_tiger.csv",quote = F, row.names = F)
write.csv(TF_activity, file = "TIGER/TF_activity_B_cells_avg_tiger.csv",quote = F, row.names = F)


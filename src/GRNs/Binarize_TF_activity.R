########################################
#Import packages
#########################################

library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(decoupleR)
setwd("/home/spere5/Documents/BALL_project/BALL_git_results/")

########################################
#Load RData (Optional)
#########################################

load("/home/spere5/Documents/BALL_project/BALL_git_results/Panda_binarized_TFs.Rdata")

#########################################
#Load RNA data and Ground truth
#########################################

#Load counts
counts <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/input/average_filtered_counts.csv", sep = ",", row.names = 1)

#Ground truth
rel_TFs <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/input/binazired_relevant_TFs.csv", row.names = 1)

##################################################
#Load GRN (or Collectri )
##################################################

net <- read.table(file = "/home/spere5/Documents/BALL_project/BALL_git_results/GRNs/quantile_95/Panda_fil.csv", sep = ",", header = T)

col <- get_collectri()

######################################################
#Run ULM and plot most variable TFs accross samples
#####################################################

#Run ulm
raw_tf_acts <- run_ulm(mat=counts, net=col, .source='source', .target='target', .mor='mor',minsize = 1)

########################################################
# Plot distribution of scores
##################################################

# Filter the data for the selected sources
rel_TFs_list <- c("BCL6", "PRDM1",  "EBF1",  "IKZF1",  "IKZF3",  "IRF4",  "POU2F2",  "PAX5",  "SPI1",  "RUNX3",  "SPIB",  "TCF3",  "XBP1")

filtered_data <- raw_tf_acts[raw_tf_acts$source %in% rel_TFs_list, ]

# Plot the distribution as boxplots
ggplot(filtered_data, aes(x = source, y = score, fill = source)) +
  geom_boxplot(alpha = 0.6) +  # Add boxplots
  theme_minimal() +
  labs(title = "Distribution of Acticity Score for Relevant TFs",
       x = "TFs",
       y = "Activity Score",
       fill = "Source") +
  theme(legend.position = "none",  # Remove legend for boxplot colors
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


##########################################
#Binarize TF activity using quantile
###########################################
# Create a new column with binarized values and set value 1 to scores in top 25% and 0 otherwise
raw_tf_acts <- raw_tf_acts %>%
  group_by(source) %>% 
  mutate(binarized_score = ifelse(score >= quantile(score, 0.75), 1, 0)) %>%
  ungroup()

#################################################
#Function that binarize decoupler output
##################################################
# Define the function
binarize_scores <- function(tf_acts, quantile_threshold = 0.75) {
  tf_acts %>%
    group_by(source) %>%
    mutate(binarized_score = ifelse(score >= quantile(score, quantile_threshold), 1, 0)) %>%
    ungroup()
}

##########################################
#Prepare results to plot
#############################################
# Transform to wide matrix
sample_acts_mat <- raw_tf_acts  %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'binarized_score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Find the common TFs between colnames of sample_acts_mat and rownames of rel_TFs
common_TFs <- intersect(colnames(sample_acts_mat), rownames(rel_TFs))

# Subset sample_acts_mat to include only columns that are common TFs
filtered_sample_acts_mat <- sample_acts_mat[, common_TFs, drop = F]

# Reorder results
filtered_sample_acts_mat_T  <- t(filtered_sample_acts_mat)

roadmap <- c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC")

filtered_sample_acts_mat_T <- filtered_sample_acts_mat_T [, roadmap, drop = FALSE]

##################################
#Plot heatmap
#####################################

# Ensure rel_TFs has the same row order as filtered_sample_acts_mat_T
rel_TFs <- rel_TFs[match(rownames(filtered_sample_acts_mat_T), rownames(rel_TFs)), ]


#Plot results
pheatmap(filtered_sample_acts_mat_T, color = (c("darkred","green")),cluster_rows = F, cluster_cols = F,)

 #############################
#Jaccard index
##############################
# Function to calculate the Jaccard index for each column
jaccard_index <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)  # common 1's
  union <- sum(x == 1 | y == 1)  # either 1 in x or y
  return(intersection / union)
}

# Apply Jaccard index to each column
jaccard_values <- mapply(jaccard_index,as.data.frame(rel_TFs),as.data.frame(filtered_sample_acts_mat_T))

# Show the Jaccard indices for each column
PanDa_jaccard <- jaccard_values

#Collect jaccard results
Panda_jaccard <- as.data.frame(Pandajaccard)
GRNBoost2_jaccard <- as.data.frame(GRNBoostjaccard)
GRaNIE_jaccard <- as.data.frame(GRaNIE_jaccard)
Col_jaccard <- as.data.frame(Collectri_jaccard)

# Combine the dataframes using cbind()
combined_df <- cbind(Panda_jaccard, GRNBoost2_jaccard, GRaNIE_jaccard, Col_jaccard)

# Calculate the column-wise means
average_row <- colMeans(combined_df)

# Convert the result to a dataframe and add it as a new row
combined_df_with_avg <- rbind(combined_df, average_row)

# Optionally, set the row name for the new row
rownames(combined_df_with_avg)[nrow(combined_df_with_avg)] <- "Average_jaccard"

#change colnames
colnames(combined_df_with_avg) <- c("Panda_Jac", "GRNBoost2_Jac", "GRaNIE_Jac", "Collectri_Jac") 
###################################
#save RData
##############################################
write.csv(combined_df_with_avg,"Jaccard_GRNs_vs_Ground_truth.csv",row.names = T,quote = F)

save.image("Collectri_binarized_TFs.Rdata")


########################################
#Import packages
#########################################

library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)


#########################################
#Load RNA data and reorder samples
#########################################

#Load counts
counts <- read.csv("/home/spere5/ball_relapse_bsc/results/rna_seq_primary/rlog_filtered_Bcounts.csv", sep = ",", row.names = 1)

#Define B cell roadmap
roadmap <- c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC")

#Define functions
reorder <- function(mat){
  stage <- sapply(strsplit(colnames(mat), "_"), "[[", 1)
  oo <- match(stage, roadmap)
  mat <- mat[, order(oo)]
  return(mat)
}
get_stage <- function(mat){
  sapply(strsplit(colnames(mat), "_"), "[[", 1)
}

#Reorder counts
counts <- reorder(counts)
stage <- get_stage(counts)

##################################################
#Load GRN
##################################################

net <- read.table(file = "/home/spere5/Documents/BALL_project/BALL_git_results/filtered_GRN.csv", sep = ",", header = T)

######################################################
#Run ULM and plot most variable TFs accross samples
#####################################################

#Run ulm
raw_tf_acts <- run_ulm(mat=counts, net=net, .source='TF', .target='target', .mor='importance',minsize = 1)

#Set number of TFs to plot
n_tfs <- 60

# Transform to wide matrix
sample_acts_mat <- raw_tf_acts  %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- raw_tf_acts  %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)
sample_acts_mat_T <- t(sample_acts_mat)
sample_acts_mat_T  <- reorder(sample_acts_mat_T)

# annots for heatmaps
seps <- which(stage[2:length(stage)]!=stage[1:(length(stage)-1)])
annots <- data.frame(stage, row.names=colnames(counts))

#Define breaks and colors
palette_length = 100
my_color = colorRampPalette(c("blue", "white","red"))(palette_length)
my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1), seq(0.05, 3, length.out=floor(palette_length/2)))

#Plot results
pheatmap(sample_acts_mat_T , cluster_cols = F, border_color = NA, color=my_color, breaks = my_breaks, scale="column", annotation_col = annots, gaps_col = seps, legend = F)

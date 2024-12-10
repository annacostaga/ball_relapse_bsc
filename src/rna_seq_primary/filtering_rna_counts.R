# Import packages and set directory
library(readr)
library(dplyr)
library(stringr)
library(DESeq2)
library(edgeR)
library(reticulate)
use_python("/home/spere5/miniconda3/bin/python3.12", required = TRUE)

setwd("/home/spere5/ball_relapse_bsc/annotation/rna_seq_primary")
dir_out = "/home/spere5/ball_relapse_bsc/results/rna_seq_primary"

# Load counts
counts_data <- read_table("counts_star_B_cell_roadmap.tsv")
# - Counts
counts_data <- as.data.frame(counts_data)

samples_ordered <- c(
                     "HSC", "HSC_1", "HSC_2",
                     "PreProB", "PreProB_1", "PreProB_2",
                     "ProB", "ProB_1", "ProB_2",
                     "PreB", "PreB_1", "PreB_2",
                     "immtransB", "immtransB_1", "immtransB_2",
                     "nB", "nB_1", "nB_2",
                     "GCB", "GCB_1", "GCB_2", "GCB_3",
                     "memB", "memB_1", "memB_2", "memB_3", 
                     "PC", "PC_1", "PC_2", "CMP", "CMP_1", "CMP_2", "CMP_3",
                     "nCD8", "nCD8_1", "nCD8_2", 
                     "Mon", "Mon_1", "Mon_2")

counts_data <- counts_data[,c("Geneid", samples_ordered)]
rownames(counts_data) <- counts_data$Geneid
counts_data <- counts_data[,-1]

#Remove lowly expressed genes
NormByCPM <- cpm(counts_data)
keep <-  rownames(NormByCPM[ rowSums(NormByCPM) > 1, ])
counts_keep <- counts_data[keep,]

# - Rlog normalization by DESeq2
NormByRlog <- rlog(as.matrix(counts_keep) + 1)

#Convert Ensemble IDs to Gene symbols
GSynoX <- import("gsynox")$GSynoX
Builder <- import("gsynox.builder")$Builder
pd <- import("pandas")
g = GSynoX()

symbol <- g$ensembl_gene_to_symbol(rownames(counts_keep))
counts_keep$symbol <- symbol

#Remove duplicates and set symbols as rownames
counts_keep <- counts_keep[!duplicated(counts_keep$symbol), ]
rownames(counts_keep) <- counts_keep$symbol
counts_keep <- counts_keep[, -ncol(counts_keep)]

#save result
write.csv(counts_keep,file = paste0(dir_out, "/rlog_filtered_Bcounts.csv"),row.names = T,quote = F)

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

samples_ordered <- c("CMP_fetal_1", "CMP_fetal_2", "CMP_fetal_3", "CMP_fetal_4",
                     "nCD8_WT_1", "nCD8_WT_2", "nCD8_WT_3", 
                     "Mon_WT_1", "Mon_WT_2", "Mon_WT_3",
                     "HSC_fetal_1", "HSC_fetal_2", "HSC_fetal_3",
                     "PreProB_fetal_1", "PreProB_fetal_2", "PreProB_fetal_3",
                     "ProB_fetal_1", "ProB_fetal_2", "ProB_fetal_3",
                     "PreB_fetal_1", "PreB_fetal_2", "PreB_fetal_3",
                     "immtransB_fetal_1", "immtransB_fetal_2", "immtransB_fetal_3",
                     "nB_WT_1", "nB_WT_2", "nB_WT_3",
                     "GCB_tonsil_1", "GCB_tonsil_3", "GCB_tonsil_4", "GCB_tonsil_5",
                     "memB_WT_1", "memB_WT_2", "memB_WT_3", "memB_WT_4", 
                     "PC_WT_1", "PC_WT_2", "PC_WT_3" )

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

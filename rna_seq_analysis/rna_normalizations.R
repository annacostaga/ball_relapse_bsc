# Packages
library(readr)
library(dplyr)
library(DESeq2)
library(EDASeq)
library(edgeR)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(forcats)

# Load counts
counts_data <- read_table("/home/acost1/BALL_project/data/RNA/counts_star_B_cell_roadmap.tsv")
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
expr_dat <- new("ExpressionSet", exprs = as.matrix(counts_data))

# - Setting the variable GROUP (general, specific)
stage <- colnames(counts_data)
pData(expr_dat)$GROUP_specific <- sub("\\_.*", "", stage)

pData(expr_dat)$GROUP_specific <- factor(pData(expr_dat)$GROUP_specific, 
                                         levels = c("CMP", "nCD8", "Mon", "HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"), 
                                         labels = c("CMP-fetal", "nCD8-WT", "Mon-WT", "HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                                    "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT"))  

# - Rlog normalization by DESeq2
NormByRlog <- rlog(as.matrix(counts_data) + 1)


save(NormByRlog, file = "/home/acost1/BALL_project/data/RNA/counts_rlog_normalized_def.rds")
save(counts_data, file = "/home/acost1/BALL_project/data/RNA/raw_counts.rds")

# - Median of ratios normalization by DESeq2
counts_data_selected <- counts_data[,!pData(expr_dat)$GROUP_specific %in% c("CMP-fetal", "nCD8-WT", "Mon-WT")]
dds <- DESeqDataSetFromMatrix(countData = counts_data_selected, 
                              colData = data.frame(GROUP_specific = pData(expr_dat)[!pData(expr_dat)$GROUP_specific %in% c("CMP-fetal", "nCD8-WT", "Mon-WT"),] ), 
                              design = as.formula(~ GROUP_specific))



dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

save(normalized_counts, file = "/home/acost1/BALL_project/data/RNA/counts_deseq2_normalized.rds")


library(tweeDEseq)
NormByTMM <- normalizeCounts(counts_data, method="TMM")

save(NormByTMM, file = "/home/acost1/BALL_project/data/RNA/counts_tmm_normalized.rds")



MbyT <- log2(NormByTMM[, 1] / NormByTMM[, 2])


# Histogram
par(mfrow=c(1,1))
hist(MbyT, xlab="log2-ratio", main="TMM")
abline(v=0, col="red")



# MA plot
library(edgeR)
maPlot(NormByTMM[,1], NormByTMM[,3], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
grid(col="black")  
title("TMM normalization") 


library(EDASeq)
plotRLE(newSeqExpressionSet(as.matrix(NormByTMM)), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(factor(pData(expr_dat)$GROUP_specific)), 
        main = 'TMM normalized counts') 

apply(NormByTMM, 2, summary)
apply(NormByRlog, 2, summary)





# Normalization for GraphNN

calculate_tpm <- function(counts, gene_lengths) {
  # Normalize by gene length (in kilobases)
  rpk <- apply(counts, 2, 
               function(x) x/(gene_lengths/10^3))
  NormByTPM <- apply(rpk, 2, 
                     function(x) x / sum(as.numeric(x), na.rm=TRUE) * 10^6)
  return(NormByTPM)
}

# Function to calculate CPM from raw counts
calculate_cpm <- function(counts) {
  # Calculate total counts per sample (column)
  total_counts <- colSums(counts)
  # Calculate CPM
  cpm <- t(t(counts) / total_counts) * 1e6
  return(cpm)
}

# Log2 transformation
log2_transform <- function(data) {
  return(log2(data + 1))
}

# Z-score normalization function
zscore_normalize <- function(data) {
  return(t(apply(data, 1, function(x) (x - mean(x)) / sd(x))))
}




# REMOVE LOWLY EXPRESSED GENES (genes that have almost no information in any of the given samples)
##################################

counts_data_selected <- counts_data[,!pData(expr_dat)$GROUP_specific %in% c("CMP-fetal", "nCD8-WT", "Mon-WT")]
NormByCPM <- calculate_cpm(counts_data_selected)
keep <-  rownames(NormByCPM[ rowSums(NormByCPM) > 1, ])
# Subset the rows of countdata to keep the more highly expressed genes
counts_keep <- counts_data_selected[keep,]

# TOTAL # of genes = 60.664
# nrow(counts_data)

# TOTAL # of genes removed = 29.488
# nrow(counts_data) - nrow(counts_keep)


# TOTAL # of genes to be analysed = 31.176
# nrow(counts_keep) 


# NORMALIZATION WITH HIGLY EXPRESSED GENES
################################################

# Normalization CPM
NormByCPM <- calculate_cpm(counts_keep)

# Normalization TPM
gene_lengths <- getGeneLengthAndGCContent(rownames(counts_keep), org = "hsa", mode = "org.db")
gene_lengths <- gene_lengths[,1] 

NormByTPM <- calculate_tpm(counts_keep, gene_lengths)

# Normalization TMM
# make the DGEList:
y <- DGEList(counts_keep)

#/ calculate TMM normalization factors:
y <- calcNormFactors(y)

library_size <- y$samples$lib.size

#/ get the normalized counts:
NormByTMM <- cpm(y, log=FALSE)
# NormByTMM2 <- normalizeCounts(counts_keep)

# DESeq2 normalization
load("/home/acost1/BALL_project/data/RNA/counts_deseq2_normalized.rds")


# rlog
NormByRlog <-rlog(as.matrix(counts_keep) )


# NORMALIZATION WITH ALL GENES
###############################

# Normalization CPM
NormByCPM_all <- calculate_cpm(counts_data_selected)

# Normalization TPM
gene_lengths <- getGeneLengthAndGCContent(rownames(counts_data), org = "hsa", mode = "org.db")
gene_lengths <- gene_lengths[,1] 

NormByTPM_all <- calculate_tpm(counts_data_selected, gene_lengths)

# Normalization TMM
# make the DGEList:
y <- DGEList(counts_data_selected)

#/ calculate TMM normalization factors:
y <- calcNormFactors(y, method = "TMM")

library_size <- y$samples$lib.size

#/ get the normalized counts:
NormByTMM_all <- cpm(y, log=FALSE)
# NormByTMM_all2 <- normalizeCounts(counts_data_selected)


# DESeq2 normalization
load("/home/acost1/BALL_project/data/RNA/counts_deseq2_normalized.rds")


# rlog
NormByRlog_all <-rlog(as.matrix(counts_data_selected) )



# BOXPLOT
###########

cell_types <- pData(expr_dat)$GROUP_specific[!pData(expr_dat)$GROUP_specific %in% c("CMP-fetal", "nCD8-WT", "Mon-WT")]
samples_ordered <- c("HSC_fetal_1", "HSC_fetal_2", "HSC_fetal_3",
                     "PreProB_fetal_1", "PreProB_fetal_2", "PreProB_fetal_3",
                     "ProB_fetal_1", "ProB_fetal_2", "ProB_fetal_3",
                     "PreB_fetal_1", "PreB_fetal_2", "PreB_fetal_3",
                     "immtransB_fetal_1", "immtransB_fetal_2", "immtransB_fetal_3",
                     "nB_WT_1", "nB_WT_2", "nB_WT_3",
                     "GCB_tonsil_1", "GCB_tonsil_3", "GCB_tonsil_4", "GCB_tonsil_5",
                     "memB_WT_1", "memB_WT_2", "memB_WT_3", "memB_WT_4", 
                     "PC_WT_1", "PC_WT_2", "PC_WT_3" )

boxplot_norm <- function(norm_data){
    boxplot_norm <- norm_data %>%
      as.data.frame() %>%
      pivot_longer(cols = colnames(.), names_to = "stage", values_to = "expression") %>%
      mutate(stage_c = sub("\\_.*", "", stage)) %>%
      mutate(stage_c = factor(stage_c ) %>% fct_relevel(c("HSC", "PreProB", "ProB", "PreB", "immtransB","nB", "GCB", "memB", "PC")),
             stage =  fct_relevel(stage, samples_ordered)) %>%
      ggplot(aes(x=stage, y=expression, fill=stage_c)) +
      geom_boxplot() + 
      scale_fill_brewer(palette="Paired") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      xlab("Sample") + ylab("Expression") +  guides(fill=guide_legend(title="Stage")) +
      ylim(-3,1000)
    return(boxplot_norm)
}



# Low expressed genes are removed
boxplot_norm_keep<- ggarrange(boxplot_norm(NormByCPM),
                             boxplot_norm(NormByTPM),
                             boxplot_norm(NormByTMM),
                             boxplot_norm(normalized_counts),
                             boxplot_norm(log2_transform(NormByCPM)),
                             boxplot_norm(log2_transform(NormByTPM)),
                             boxplot_norm(log2_transform(NormByTMM)),
                             boxplot_norm(log2_transform(normalized_counts)),
                          labels = c("CPM", "TPM", "TMM", "DESeq2",
                                     "CPM", "TPM", "TMM", "DESeq2"),
                          common.legend = TRUE, legend = "top",
                          ncol = 4, nrow = 2)

png("/home/acost1/BALL_project/results/gnn/figures/boxplot_norm_keep.png", units="in", width=30, height=15, res=200)
boxplot_norm_keep
dev.off()

png("/home/acost1/BALL_project/results/gnn/figures/boxplot_rlog_keep.png", units="in", width=10, height=5, res=150)
boxplot_norm(NormByRlog)
dev.off()

boxplot_norm_all<- ggarrange(boxplot_norm(NormByCPM_all),
                              boxplot_norm(NormByTPM_all),
                              boxplot_norm(NormByTMM_all),
                              boxplot_norm(normalized_counts),
                              boxplot_norm(log2_transform(NormByCPM_all)),
                              boxplot_norm(log2_transform(NormByTPM_all)),
                              boxplot_norm(log2_transform(NormByTMM_all)),
                              boxplot_norm(log2_transform(normalized_counts)),
                              labels = c("CPM", "TPM", "TMM", "DESeq2",
                                         "CPM", "TPM", "TMM", "DESeq2"),
                              common.legend = TRUE, legend = "top",
                              ncol = 4, nrow = 2)

png("/home/acost1/BALL_project/results/gnn/figures/boxplot_norm_all.png", units="in", width=30, height=15, res=200)
boxplot_norm_all
dev.off()



png("/home/acost1/BALL_project/results/gnn/figures/boxplot_rlog_all.png", units="in", width=10, height=5, res=200)
boxplot_norm(NormByRlog_all)
dev.off()


ggarrange(
boxplot_norm(calculate_tpm(NormByTMM_all, gene_lengths)),
boxplot_norm(NormByTMM_all), nrow = 1 )

boxplot_norm_genes <- function(norm_data){
  boxplot_norm <- melt(norm_data) %>%
    dplyr::rename(gene = Var1, sample = Var2, expression = value) %>%
    ggplot(aes(x=sample, y=expression, fill=gene)) +
    geom_boxplot() + 
    scale_fill_brewer(palette="Paired") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Sample") + ylab("Expression") +  guides(fill=guide_legend(title="Stage")) +
    ylim(-3,20)
  return(boxplot_norm)
}

boxplot_norm_genes(norm_data = NormByCPM)



save(NormByTMM,file = "/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_keep.rds")
save(NormByTMM_all, file = "/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_all.rds")

NormByTMM_all_log <- log2_transform(NormByTMM_all)
save(NormByTMM_all_log, file = "/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_all_log.rds")


NormByTMM_log <- log2_transform(NormByTMM)
save(NormByTMM_log,file = "/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_keep_log.rds")

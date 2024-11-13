# DESeq2 analysis


# Packages
rm(list=ls())  # clear namespace
suppressPackageStartupMessages({
  library(Biobase)
  library(readr)
  library(ggfortify)
  library(EDASeq)
  library(DESeq2)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(optparse)
}) 

# function

deseq <- function(counts_data){
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
  
  # Setting the variable group
  
  # Setting the variable GROUP (general, specific)
  stage <- colnames(counts_data)
  pData(expr_dat)$GROUP_general <- rep(NA, ncol(counts_data))
  pData(expr_dat)$GROUP_specific <- sub("\\_.*", "", stage)
  
  pData(expr_dat)$GROUP_specific <- factor(pData(expr_dat)$GROUP_specific, 
                                           levels = c("CMP", "nCD8", "Mon", "HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"), 
                                           labels = c("CMP-fetal", "nCD8-WT", "Mon-WT", "HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                                      "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT"))   
  
  
  
  # DATA EXPLORATION
  
  # rlog DESeq2 normalization
  
  NormByRlog <- rlog(as.matrix(counts_data))
  
  
  
  # Hierarchical clustering
  
  # Euclidean distance
  dist <- dist(t(NormByRlog) , diag=TRUE)

  # Hierarchical Clustering with hclust
  hc <- hclust(dist)

  # Plot the result
  png(file=paste0(opt$outdir, "dendogram.png"), width=1000, height=650)
  plot(hc)
  dev.off()
  
  cat("      - Hierarchical clustering finished \n")
  
  
  # PCA 

  #compute PCA 
  pcaResults <- prcomp(t(NormByRlog))

  #plot PCA results making use of ggplot2's autoplot function
  autoplot(pcaResults, data = pData(expr_dat), colour = 'GROUP_specific')
  ggsave(paste0(opt$outdir, 'plot_pca.png'))
  
  
  cat("      - PCAfinished \n")
  
  
  # RLE PLOT
  png(file=paste0(opt$outdir, "rle_counts.png"), width=1000, height=650)
  plotRLE(newSeqExpressionSet(as.matrix(counts_data)), 
         outline=FALSE, ylim=c(-2, 2), 
         col = as.numeric(factor(pData(expr_dat)$GROUP_specific)), 
         main = 'Raw counts') 
  dev.off()

  png(file=paste0(opt$outdir, "rle_normalised.png"), width=1000, height=650)
  plotRLE(newSeqExpressionSet(as.matrix(NormByRlog)), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(factor(pData(expr_dat)$GROUP_specific)), 
        main = 'Normalized Counts (rlog - DESeq2)') 
  dev.off()
  
  cat("      - RLE finished \n")
  
  
  # DIFFERENTIAL EXPRESSION ANALYSIS
  # #create a DESeq dataset object from the count matrix and the colData 
  dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                                colData = pData(expr_dat), 
                                design = as.formula(~ GROUP_specific))
  
  
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  
  
  png(file=paste0(opt$outdir, "rle_normalised_dds.png"), width=1000, height=650)
  plotRLE(newSeqExpressionSet(counts(dds, normalized=TRUE)), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(factor(pData(expr_dat)$GROUP_specific)), 
        main = 'Normalized Counts (DESeq2)') 
  dev.off()
  
  
  DEresults_HSC <- dds[ rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("HSC-fetal", "CMP-fetal")], normalized = TRUE ) ) > 2 , ] %>%
    results(., contrast = c("GROUP_specific", 'HSC-fetal', 'CMP-fetal'), altHypothesis = "greaterAbs") %>%
    as.data.frame() 
  
  DEresults_PreProB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("PreProB-fetal", "HSC-fetal")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'PreProB-fetal', 'HSC-fetal')) %>%
    as.data.frame() 
  
  DEresults_ProB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("ProB-fetal", "PreProB-fetal")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'ProB-fetal', 'PreProB-fetal')) %>%
    as.data.frame() 
  
  DEresults_PreB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("PreB-fetal", "ProB-fetal")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'PreB-fetal', 'ProB-fetal')) %>%
    as.data.frame()
  
  DEresults_TransB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("immtransB-fetal", "PreB-fetal")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'immtransB-fetal', 'PreB-fetal')) %>%
    as.data.frame() 
  
  DEresults_nB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("nB-WT", "immtransB-fetal")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'nB-WT', 'immtransB-fetal')) %>%
    as.data.frame()
  
  DEresults_GCB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("GCB-tonsil", "nB-WT")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'GCB-tonsil', 'nB-WT')) %>%
    as.data.frame() 
  
  DEresults_memB <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("memB-WT", "GCB-tonsil")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'memB-WT', 'GCB-tonsil')) %>%
    as.data.frame() 
  
  DEresults_PC <- results(dds[rowSums(counts(dds[,colData(dds)$GROUP_specific %in% c("PC-WT", "GCB-tonsil")], normalized = TRUE)) > 2, ], contrast = c("GROUP_specific", 'PC-WT', 'memB-WT')) %>%
    as.data.frame() 
    
    
  
  # Histogram of p-value without surrogate analysis
  pvalue_prepro <- ggplot(data = as.data.frame(DEresults_PreProB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("PreProB")
  pvalue_pro <- ggplot(data = as.data.frame(DEresults_ProB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("ProB")
  pvalue_pre <- ggplot(data = as.data.frame(DEresults_PreB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("PreB")
  pvalue_trans <- ggplot(data = as.data.frame(DEresults_TransB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("TransB")
  pvalue_nb <- ggplot(data = as.data.frame(DEresults_nB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("nB")
  pvalue_gcb <- ggplot(data = as.data.frame(DEresults_GCB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("GCB")
  pvalue_mem <- ggplot(data = as.data.frame(DEresults_memB), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("memB")
  pvalue_plasma <- ggplot(data = as.data.frame(DEresults_PC), aes(x = pvalue)) +
                  geom_histogram(bins = 100) + ggtitle("Plasma")


  ggarrange(pvalue_prepro, pvalue_pro, pvalue_pre, pvalue_trans,
                    pvalue_nb, pvalue_gcb, pvalue_mem, pvalue_plasma,
                    ncol = 4, nrow = 2)
  ggsave(paste0(opt$outdir, 'histogram_pvalues_deseq2.png'))
  
  
  cat("      - DESeq2 finished \n")
  
  
  de <- list(dds = dds, prepro = DEresults_PreProB, proB = DEresults_ProB, preB = DEresults_PreB, transB = DEresults_TransB, 
                   nB = DEresults_nB, gcb = DEresults_GCB, memB = DEresults_memB, plasma = DEresults_PC)
  
  return(de)
}


# Import data
option_list = list(
  make_option(c("-c", "--counts_file"), type="character", help="Counts file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# mandatory params
if (is.null(opt$counts_file) || opt$counts_file == ''){
  print_help(opt_parser)
  stop("No counts_data file path provided", call.=FALSE)
}
if (is.null(opt$outdir) || opt$outdir == ''){
  print_help(opt_parser)
  stop("No output folder path provided", call.=FALSE)
}


cat("\n\n")
cat("***********************************\n")
cat("*** DESeq2 pipeline ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Countsdata: ", opt$counts_file, "\n")


counts_data <- read_table(opt$counts_file)

# create working directory
dir.create(opt$outdir, recursive = T, showWarnings = F)

results <- deseq(counts_data)

save(results, file=paste0(opt$outdir, 'deseq2_results2.rds'))

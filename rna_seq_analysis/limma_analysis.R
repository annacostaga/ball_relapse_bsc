# limma analysis


# Packages
rm(list=ls())  # clear namespace
suppressPackageStartupMessages({
  library(readr)
  library(Biobase)
  library(limma)
  library(sva)
  library(S4Vectors)
  library(edgeR)
  library(ggpubr)
  library(ggplot2)
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
                                           labels = c("CMPfetal", "nCD8WT", "MonWT", "HSCfetal", "PreProBfetal", "ProBfetal", "PreBfetal",
                                                      "immtransBfetal", "nBWT", "GCBtonsil", "memBWT", "PCWT"))   
  
  
  
  # PREPROCESSING
  
  d0 <- DGEList(counts_data)
  d0 <- calcNormFactors(d0)
  
  cutoff <- 2
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  dim(d) 
  
  colData <- pData(expr_dat)
 
 
  # DIFFERENTIAL EXPRESSION ANALYSIS WITHOUT SURROGATE VARIABLES
  mm <- model.matrix( ~ 0 + GROUP_specific , data=colData)
  
  png(file=paste0(opt$outdir, "voom.png"), width=1000, height=650)
  y <- voom(d, mm, plot = T)
  dev.off()
  
  fit <- lmFit(y, mm)
  head(coef(fit))

  contr <- makeContrasts(GROUP_specificHSCfetal - GROUP_specificCMPfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_HSC <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPreProBfetal - GROUP_specificHSCfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_PreProB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificProBfetal - GROUP_specificPreProBfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_ProB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPreBfetal - GROUP_specificProBfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_PreB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificimmtransBfetal - GROUP_specificPreBfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_TransB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificnBWT - GROUP_specificimmtransBfetal, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_nB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificGCBtonsil - GROUP_specificnBWT, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_GCB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificmemBWT - GROUP_specificGCBtonsil, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_memB <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPCWT - GROUP_specificGCBtonsil, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  DEresults_PC <- topTable(tmp, sort.by = "P", n = Inf)
  
  
  # DIFFERENTIAL EXPRESSION ANALYSIS WITH SURROGATE VARIABLES, Accounting for unwanted source of variability
  
  # Surrogate variables ----
  colData <- pData(expr_dat)
  colData$GROUP_specific <- factor(colData$GROUP_specific)

  mod1 <- model.matrix( ~ GROUP_specific , data=colData)
  mod0 <- model.matrix( ~ 1, data=colData)

  #  compute the surrogate variable
  sv <- svaseq(as.matrix(d), mod = mod1, mod0 = mod0) # 8 variables

  # add them in the colData
  temp <- DataFrame(colData, sv$sv)
  
  
  design.sv <- model.matrix( ~ 0 + GROUP_specific + V1 + V2 + V3 +V4 +V5 +V6 +V7 + V8  , data=temp)

  fit.sv <- lmFit(y, design.sv)
  fit.sv <- eBayes(fit.sv)

  contr <- makeContrasts(GROUP_specificHSCfetal - GROUP_specificCMPfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_HSC_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPreProBfetal - GROUP_specificHSCfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_PreProB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificProBfetal - GROUP_specificPreProBfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_ProB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPreBfetal - GROUP_specificProBfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_PreB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificimmtransBfetal - GROUP_specificPreBfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_TransB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificnBWT - GROUP_specificimmtransBfetal, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_nB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificGCBtonsil - GROUP_specificnBWT, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_GCB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificmemBWT - GROUP_specificGCBtonsil, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_memB_sv <- topTable(tmp, sort.by = "P", n = Inf)

  contr <- makeContrasts(GROUP_specificPCWT - GROUP_specificGCBtonsil, levels = colnames(coef(fit.sv)))
  tmp <- contrasts.fit(fit.sv, contr)
  tmp <- eBayes(tmp)
  DEresults_PC_sv <- topTable(tmp, sort.by = "P", n = Inf)
  
  
  # Histogram of p-value without surrogate analysis
  pvalue_prepro <- ggplot(data = as.data.frame(DEresults_PreProB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("PreProB")
  pvalue_pro <- ggplot(data = as.data.frame(DEresults_ProB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("ProB")
  pvalue_pre <- ggplot(data = as.data.frame(DEresults_PreB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("PreB")
  pvalue_trans <- ggplot(data = as.data.frame(DEresults_TransB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("TransB")
  pvalue_nb <- ggplot(data = as.data.frame(DEresults_nB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("nB")
  pvalue_gcb <- ggplot(data = as.data.frame(DEresults_GCB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("GCB")
  pvalue_mem <- ggplot(data = as.data.frame(DEresults_memB), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("memB")
  pvalue_plasma <- ggplot(data = as.data.frame(DEresults_PC), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("Plasma")


  ggarrange(pvalue_prepro, pvalue_pro, pvalue_pre, pvalue_trans,
                    pvalue_nb, pvalue_gcb, pvalue_mem, pvalue_plasma,
                    ncol = 4, nrow = 2)
  ggsave(paste0(opt$outdir, 'histogram_pvalues_limma_not_sv.png'))


  # Histogram of p-value with surrogate analysis
  pvalue_prepro <- ggplot(data = as.data.frame(DEresults_PreProB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("PreProB - sv")
  pvalue_pro <- ggplot(data = as.data.frame(DEresults_ProB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("ProB - sv")
  pvalue_pre <- ggplot(data = as.data.frame(DEresults_PreB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("PreB - sv")
  pvalue_trans <- ggplot(data = as.data.frame(DEresults_TransB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("TransB - sv")
  pvalue_nb <- ggplot(data = as.data.frame(DEresults_nB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("nB - sv")
  pvalue_gcb <- ggplot(data = as.data.frame(DEresults_GCB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("GCB - sv")
  pvalue_mem <- ggplot(data = as.data.frame(DEresults_memB_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("memB - sv")
  pvalue_plasma <- ggplot(data = as.data.frame(DEresults_PC_sv), aes(x = P.Value)) +
                  geom_histogram(bins = 100) + ggtitle("Plasma - sv")

  ggarrange(pvalue_prepro, pvalue_pro, pvalue_pre, pvalue_trans,
                    pvalue_nb, pvalue_gcb, pvalue_mem, pvalue_plasma,
                    ncol = 4, nrow = 2)
  ggsave(paste0(opt$outdir, 'histogram_pvalues_limma_sv.png'))
                    
                    
  
  de_not_sv <- list(
  	     hsc = DEresults_HSC,
  	     prepro = DEresults_PreProB, proB = DEresults_ProB, preB = DEresults_PreB, transB = DEresults_TransB, 
             nB = DEresults_nB, gcb = DEresults_GCB, memB = DEresults_memB, plasma = DEresults_PC)
             
  de_not_sv <- lapply(de_not_sv, function(x)  
						  dplyr::rename(x,
								log2FoldChange = logFC,
								padj = adj.P.Val
						  ) %>% as.data.frame()
  )          
            
  de_sv <- list(hsc = DEresults_HSC_sv,
  	     prepro = DEresults_PreProB_sv, proB = DEresults_ProB_sv, preB = DEresults_PreB_sv, transB = DEresults_TransB_sv, 
             nB = DEresults_nB_sv, gcb = DEresults_GCB_sv, memB = DEresults_memB_sv, plasma = DEresults_PC_sv) 
             
  
  de_sv <- lapply(de_sv, function(x)  
						  dplyr::rename(x,
								log2FoldChange = logFC,
								padj = adj.P.Val
						  ) %>% as.data.frame() )
  
  return( list( de_not_sv = list(de_not_sv),  de_sv = list(de_sv) )  )
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
cat("*** Limma - voom pipeline ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Countsdata: ", opt$counts_file, "\n")


counts_data <- read_table(opt$counts_file)

# create working directory
dir.create(opt$outdir, recursive = T, showWarnings = F)

results <- deseq(counts_data)
results_not_sv <- results$de_not_sv[[1]]
results_sv <- results$de_sv[[1]]

save(results_not_sv, file=paste0(opt$outdir, 'limma_results_not_sv2.rds'))
save(results_sv, file=paste0(opt$outdir, 'limma_results_sv2.rds'))

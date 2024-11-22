# GSEA analysis


# Packages
rm(list=ls())  # clear namespace
suppressPackageStartupMessages({
  library(Biobase)
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(optparse)
}) 

# function

gsea <- function(results){
   
  # PreProB vs HSC
  
  # GSEA GO
  original_gene_list <- results$prepro$log2FoldChange[results$prepro$padj < 0.01 &  !is.na(results$prepro$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$prepro[results$prepro$padj < 0.01 &  !is.na(results$prepro$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_prepro <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
  
  
  cat("      - GSEA PrePro finished \n")
  
  
  
  
  
  # ProB vs PreProB
  
  # GSEA GO
  original_gene_list <- results$proB$log2FoldChange[results$proB$padj < 0.01 &  !is.na(results$proB$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$proB[results$proB$padj < 0.01 &  !is.na(results$proB$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_proB <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
    
  
  cat("      - GSEA Pro finished \n")
  
  
  
  
  
  # PreB vs ProB
  
  # GSEA GO
  original_gene_list <- results$preB$log2FoldChange[results$preB$padj < 0.01 &  !is.na(results$preB$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$preB[results$preB$padj < 0.01 &  !is.na(results$preB$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_preB <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
               
  
  cat("      - GSEA Pre finished \n")
  
  
  
  
  # TransB vs PreB
  
  # GSEA GO
  original_gene_list <- results$transB$log2FoldChange[results$transB$padj < 0.01 &  !is.na(results$transB$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$transB[results$transB$padj < 0.01 &  !is.na(results$transB$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_transB <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()

               
               
  cat("      - GSEA Trans finished \n")
  
  
  
  
  # nB vs TransB
  
  # GSEA GO
  original_gene_list <- results$nB$log2FoldChange[results$nB$padj < 0.01 &  !is.na(results$nB$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$nB[results$nB$padj < 0.01 &  !is.na(results$nB$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_nB <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
  
  
  cat("      - nB finished \n")
  
  # GCB vs nB
  
  # GSEA GO
  original_gene_list <- results$gcb$log2FoldChange[results$gcb$padj < 0.01 &  !is.na(results$gcb$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$gcb[results$gcb$padj < 0.01 &  !is.na(results$gcb$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_gcb <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
    
 
  cat("      - GSEA GCB finished \n")
  
  # memB vs GSC
  
  # GSEA GO
  original_gene_list <- results$memB$log2FoldChange[results$memB$padj < 0.01 &  !is.na(results$memB$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$memB[results$memB$padj < 0.01 &  !is.na(results$memB$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_memB <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
               
  
  cat("      - GSEA memB finished \n")
  
  # Plasma vs memB
  
  original_gene_list <- results$plasma$log2FoldChange[results$plasma$padj < 0.01 &  !is.na(results$plasma$padj)]
  
  # name the vector
  names(original_gene_list) <- rownames(results$plasma[results$plasma$padj < 0.01 &  !is.na(results$plasma$padj),])
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  go_plasma <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "fdr") %>%
    as.data.frame()
 
               
               
  
  cat("      - GSEA Plasma finished \n")
  
  gsea_all <- list(go_prepro = go_prepro, go_proB = go_proB, go_preB = go_preB, go_transB = go_transB, 
                   go_nB = go_nB, go_gcb = go_gcb, go_memB = go_memB, go_plasma = go_plasma
                   )
  
  return(gsea_all)
}


# Import data
option_list = list(
  make_option(c("-d", "--deseq_file"), type="character", help="DESeq file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# mandatory params
if (is.null(opt$deseq_file) || opt$deseq_file == ''){
  print_help(opt_parser)
  stop("No deseq_data file path provided", call.=FALSE)
}
if (is.null(opt$outdir) || opt$outdir == ''){
  print_help(opt_parser)
  stop("No output folder path provided", call.=FALSE)
}


cat("\n\n")
cat("***********************************\n")
cat("*** GSEA PROCESSING ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("DESeq data: ", opt$deseq_file, "\n")


load(opt$deseq_file)

# create working directory
dir.create(opt$outdir, recursive = T, showWarnings = F)

results_gsea <- gsea(results)

save(results_gsea, file=file_path(opt$outdir, 'gsea_results.rds'), row.names=T, sep='\t', col.names=T, quote=F)

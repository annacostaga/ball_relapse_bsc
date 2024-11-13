# KmlShape analysis

# Install packages


# Packages
rm(list=ls())  # clear namespace
suppressPackageStartupMessages({
  library(kmlShape)
  library(readr)
  library(DESeq2)
  library(edgeR)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(optparse)
}) 


# function

kmlshape_pip <- function(counts_data, results){
  
  ## Input - counts data
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
                                                      
  cat("      - Input counts data finished \n")
                                                      
                                                      
  ## Input - DE Genes

  selected_genes <- unique(c(
			results$prepro %>% filter(padj < 0.01) %>% rownames(),
			results$proB %>% filter(padj < 0.01) %>% rownames(),
			results$preB %>% filter(padj < 0.01) %>% rownames(),
			results$transB %>% filter(padj < 0.01) %>% rownames(),
			results$nB %>% filter(padj < 0.01) %>% rownames(),
			results$gcb %>% filter(padj < 0.01) %>% rownames(),
			results$memB %>% filter(padj < 0.01) %>% rownames(),
			results$plasma %>% filter(padj < 0.01) %>% rownames() ) ) 
  
  cat("      - Input - DE genes finished \n")
  
  
  ## Applying rlog
  
  NormByRlog <- rlog(as.matrix(counts_data))
  NormByRlog_selected <- NormByRlog[selected_genes, !pData(expr_dat)$GROUP_specific %in%
  						     c("CMP-fetal", "Mon-WT", "nCD8-WT")] 
  
  
  cat("      - applying rlog finished \n")
  
  
  ## Preparing data for kmlshape
  
  counts_long <- NormByRlog_selected %>%
               as.data.frame() %>%
               mutate(gene = rownames(NormByRlog_selected) )%>%
               pivot_longer(!gene,names_to = "cell_type", values_to = "count")
  counts_long$group_cell_type <- sub("\\_.*", "", counts_long$cell_type)
  counts_long <- group_by(counts_long,gene, group_cell_type) %>%
               summarise(avg_counts = mean(count)) 
               
  print(head(counts_long) )
  
  cat("\n \n")   
           
  counts_long$group_cell_type <- factor(counts_long$group_cell_type,
                                         levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"),
                                         labels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"))
  counts_long$group_cell_type <- as.numeric(counts_long$group_cell_type)
  counts_long$gene <- as.numeric(as.factor(counts_long$gene))

  print(head(counts_long))

  traj <- rename(counts_long, id = gene, times = group_cell_type, traj = avg_counts)
  traj <- traj[order(traj$id, traj$times), ]
  traj$id <- as.integer(traj$id)
  traj$times <- as.integer(traj$times)

  myClds <- cldsLong(as.data.frame(traj))
  
  cat("      - prepare data for kmlshape finished \n")
  
  #png(file=paste0(opt$outdir, "trajectories.png"), width=1000, height=650)
  #plot(myClds)
  #dev.off()

  cat("      - prepare data for kmlshape finished \n")
   
   
  ## Run kmlShape
  
 # kml <- kmlShape(myClds, nbClusters = 8, timeScale = 0.1, FrechetSumOrMax = "max", 
        # toPlot="both", parAlgo=parKmlShape())
  # ggsave(paste0(opt$outdir, 'kml_plot.png'))    
  
  
 #  cat("      - kmlShape finished \n")  
           
  return(kml)
}


# Import data
option_list = list(
  make_option(c("-c", "--counts_file"), type="character", help="Counts file", metavar="character"),
  make_option(c("-d", "--deseq_file"), type="character", help="DESeq file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# mandatory params
if (is.null(opt$counts_file) || opt$counts_file == ''){
  print_help(opt_parser)
  stop("No counts_data file path provided", call.=FALSE)
}
if (is.null(opt$deseq_file) || opt$deseq_file == ''){
  print_help(opt_parser)
  stop("No deseq_file file path provided", call.=FALSE)
}
if (is.null(opt$outdir) || opt$outdir == ''){
  print_help(opt_parser)
  stop("No output folder path provided", call.=FALSE)
}


cat("\n\n")
cat("***********************************\n")
cat("*** KmlShape pipeline ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Counts data: ", opt$counts_file, "\n")
cat("Deseq data: ", opt$deseq_file, "\n")

counts_data <- read_table(opt$counts_file)
load(opt$deseq_file)

# create working directory
dir.create(opt$outdir, recursive = T, showWarnings = F)

results_kml <- kmlshape_pip(counts_data, results)

save(results_kml, file=paste0(opt$outdir, 'kmlshape_results.rds'))

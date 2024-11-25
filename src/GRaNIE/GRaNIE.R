#FOR MORE INFORMATION AND DETAILS SEE:
## WORKFLOW: https://bioconductor.org/packages/devel/bioc/vignettes/GRaNIE/inst/doc/GRaNIE_workflow.html
## PACKAGE DETAILS: https://bioconductor.org/packages/devel/bioc/vignettes/GRaNIE/inst/doc/GRaNIE_packageDetails.html#methods_peakGene

#########################################
### MAKE CHANGES TO PACKAGE TO RUN OFFLINE (MN5)
################################################

"""
1)

Install GRaNIE package --> Go into GRaNIE-package folder --> R-folder --> helperFunctions.

2)

In helperFunctions.R change diretory of chrSize from

chrSizes = GenomeInfoDb::getChromInfoFromUCSC(genomeAssembly)

To:

 chrSizes = read_csv(file = 'chrSizes')

3)

Install modified GRaNIE-package

"""

########################################
### IMPORT PACKAGES AND SET DIRECTORY
################################

library(readr)
library(GRaNIE)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationHub)
library(biomaRt)
library(ChIPseeker)

setwd("../BALL_git/ball_relapse_bsc/annotation/GRaNIE")

#######################################
### LOAD DATA 
######################################

#Load data ATAC:
file_peaks = read.table("GRaNIE_ATAC_peaks.csv",header = T,sep = ",")

#Load RNAseq and edit
file_RNA = read.table("GRaNIE_RNA_Bcounts.csv",header = T,sep = ",")

# Save the name of the respective ID columns
idColumn_peaks = "peakID"
idColumn_RNA = "ENSEMBL"

#load and edit lichi data
lichi <- read.table(("GRaNIE_lichi.csv"),header = T,sep = ",")


###############################################
### Initialize GRaNIE object 
##################################################

#Choose either hg19, hg38 or mm10 genome. Both enhancers and RNA data must have the same genome assembly
genomeAssembly = "hg38"  

# Optional and arbitrary list with information and metadata that is stored within the GRaNIE object
objectMetadata.l = list(name = paste0("B cell differentiation"))

#choose output directory
dir_output = "../BALL_git/ball_relapse_bsc/results/GRaNIE"

#Create GRN object
GRN = initializeGRN(objectMetadata = objectMetadata.l, outputFolder = dir_output,
                    genomeAssembly = genomeAssembly)

#Normalize and add RNA and ATAC data to GRN object
GRN = addData(GRN, counts_peaks = file_peaks, normalization_peaks = "DESeq2_sizeFactors",
              idColumn_peaks = idColumn_peaks, counts_rna = file_RNA, normalization_rna = "limma_quantile",
              idColumn_RNA = idColumn_RNA, forceRerun = TRUE)

#Quality control 1: PCA plot to check the 500 most variable features in the normalized count data for RNA and peaks
GRN = plotPCA_all(GRN, data = c("rna", "peaks"), topn = 500, type = "normalized")

#Filter RNA and peak counts (optional)
GRN = filterData(GRN, minNormalizedMean_peaks = 5, minNormalizedMeanRNA = 1,
                 maxSize_peaks = 10000, forceRerun = TRUE)

#Add TFs and TFBS data 
motifFolder = tools::file_path_as_absolute("H12INVIVO")

GRN = addTFBS(GRN, motifFolder = motifFolder, TFs = "all", filesTFBSPattern = "_TFBS",
              fileEnding = ".bed.gz", forceRerun = TRUE)

#Get TFBS that overlap with peaks
GRN = overlapPeaksAndTFBS(GRN, nCores = 1, forceRerun = TRUE)

#Save GRN object
GRN_file_outputRDS = paste0(dir_output, "/GRN.rds")
saveRDS(GRN, GRN_file_outputRDS)

##############################################
### ADD PEAK-TF AND PEAK-GENE CONNECTIONS
##############################################

#Add TF and peak connections based on correlation
GRN = addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, connectionTypes = c("expression"),
                             corMethod = "pearson", forceRerun = TRUE)

#Quality control 2: Diagnostic plots for TF-enhancer connections
GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = FALSE)
GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = FALSE)
GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("background"), plotAsPDF = FALSE)

#Run the AR classification and QC to classify TFs as activators or repressors (optional)
GRN = AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 0.05, outputFolder = "plots",
                                plot_minNoTFBS_heatmap = 100, plotDiagnosticPlots = TRUE, forceRerun = TRUE)

#Add peak-gene connections based on correlation, proximity and 3D organization
GRN = addConnections_peak_gene(GRN, corMethod = "pearson",
                               TADs = NULL, knownLinks = lichi_df, nCores = 1, plotDiagnosticPlots = FALSE, plotGeneTypes = list(c("all")),
                               forceRerun = TRUE)

#Quality control 3: Diagnostic plots for enhancer-gene connections
GRN = plotDiagnosticPlots_peakGene(GRN, gene.types = list(c("protein_coding", "lincRNA")),
                                   plotAsPDF = FALSE)

#Combine TF-peak and enhancer-peak connections and filter connections based on FDR
GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.2, peak_gene.fdr.threshold = 0.2, 
peak_gene.fdr.method = "BH", gene.types = c("protein_coding", "lincRNA"), allowMissingTFs = FALSE, allowMissingGenes = FALSE, forceRerun = TRUE)


############################################
####  ADD TF-GENE CONNECTIONS AND BUILD GRAPH
###############################################

#Add TF-gene correlations (optional)
GRN = add_TF_gene_correlation(GRN, corMethod = "pearson", nCores = 1, forceRerun = TRUE)

#Delete intermediate data to reduce memory usage and save GRN object
GRN = deleteIntermediateData(GRN)
saveRDS(GRN, GRN_file_outputRDS)

#Retrieve filtered connections as a data frame from GRN object
GRN_connections.all = getGRNConnections(GRN, type = "all.filtered", include_TF_gene_correlations = TRUE)

#Generate a connection summary for filtered connections
GRN = generateStatsSummary(GRN, TF_peak.connectionTypes = "all", peak_gene.r_range = c(0, 1),
allowMissingGenes = c(FALSE,TRUE), allowMissingTFs = c(FALSE), gene.types = c("protein_coding", "lincRNA"),forceRerun = TRUE)
GRN = plot_stats_connectionSummary(GRN, type = "heatmap", plotAsPDF = FALSE)
GRN = plot_stats_connectionSummary(GRN, type = "boxplot", plotAsPDF = FALSE)

#Construct GRN graph with filtered connections for visualization
GRN = build_eGRN_graph(GRN, forceRerun = TRUE)

#Visualize the GRN
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 100000)

#Filter network for visualization by defining a TF
GRN = filterConnectionsForPlotting(GRN, plotAll = FALSE, TF.ID == "ATF1.H12INVIVO.0.P.B",forceRerun = T)
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, forceRerun = T)

################################################
### ANALYSE GRN AND COMMUNITIES
###################################################

#Network and enrichment analyses for filtered connections
GRN = performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), outputFolder = ".", forceRerun = T)

#General network enrichment
GRN = plotGeneralEnrichment(GRN, plotAsPDF = T)

#
GRN = plotCommunitiesEnrichment(GRN, plotAsPDF = T)

#TF enrichment analysis
GRN = plotTFEnrichment(GRN, plotAsPDF = FALSE, n = 3,forceRerun = T)

#############################################
### SAVE GRN OBJECT
###################################################

GRN = deleteIntermediateData(GRN)
saveRDS(GRN, file = "GRN.rds")

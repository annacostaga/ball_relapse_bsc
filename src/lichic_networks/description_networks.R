library(HiCaptuRe)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(viridis)

#######################
# 1) RAW PEAKMATRIX
#######################


file <- "/home/acost1/BALL_project/data/lichic/lichic_peakmatrix_B_cell_roadmap.txt"
data <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")

data$ID_min <- pmin(data$baitID, data$oeID)
data$ID_max <- pmax(data$baitID, data$oeID)

# 2: Count occurrences of each sorted pair
data$pair <- paste(data$ID_min, data$ID_max, sep = "_")
pair_counts <- table(data$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1])


data_overlap <- filter(data, pair %in% interactions_overlap)
data_overlap %>% group_by(pair) %>% summarise(id1 = length(unique(baitID)),
                                              id2 = length(unique(oeID)))
# We have different chicago score, for different directions of the B-B interaction ---> n = 21.248



###################################################
### 2. UPDATED FUNCTIONS FROM HiCaptuRE package 
####################################################
load_interactions
inter <- load_interactions(file)
inter_annotated <- annotate_interactions(inter, annotation = "/home/acost1/BALL_project/data/lichic/annotation_example.txt")
inter_annotated <- as.data.frame(inter_annotated)


inter_annotated$ID_min <- pmin(inter_annotated$ID_1, inter_annotated$ID_2)
inter_annotated$ID_max <- pmax(inter_annotated$ID_1, inter_annotated$ID_2)

# 2: Count occurrences of each sorted pair
inter_annotated$pair <- paste(inter_annotated$ID_min, inter_annotated$ID_max, sep = "_")
pair_counts <- table(inter_annotated$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1])


data_overlap2 <- filter(inter_annotated, pair %in% interactions_overlap)
data_overlap2 %>% group_by(pair) %>% summarise(id1 = length(unique(ID_1)),
                                                       id2 = length(unique(ID_2)))

# - Number of genes/communities
length(unique(inter_annotated$ID_1)) # n = 21.458

# - Number of genes that interact with other genes
inter_annotated_BB <- filter(inter_annotated, int == "B_B")
length(unique(inter_annotated_BB$ID_1)) # n = 16.514

# - Number of genes that interact with other genes this interaction is bidirectional
inter_annotated_BB$ID_min <- pmin(inter_annotated_BB$ID_1, inter_annotated_BB$ID_2)
inter_annotated_BB$ID_max <- pmax(inter_annotated_BB$ID_1, inter_annotated_BB$ID_2)

inter_annotated_BB$pair <- paste(inter_annotated_BB$ID_min, inter_annotated_BB$ID_max, sep = "_")

inter_annotated_BB$pair 
interactions_overlap

length(unique(inter_annotated_BB_bi$ID_1)) # n = 12.240


# Number of all interactions between baits n = 64.977
length(unique(inter_annotated_BB$pair))

# Number of bidirectional interactions between baits n = 21.348
inter_annotated_BB$pair[inter_annotated_BB$pair %in% interactions_overlap ] %>% unique() %>% length()

# Number of unidirectional interactions between baits n = 43.629
inter_annotated_BB$pair[!inter_annotated_BB$pair %in% interactions_overlap ] %>% unique() %>% length()

####################################################################################################################
### 3. Interactions from updated HiCaptuRE package and annotated with ChromHMM (only significant interactions)
####################################################################################################################

load("/home/acost1/BALL_project/results/lichi_chromhmm/inter_annotated_gene_updated.rds") 


# Find bidirectional interactions inside the edge_file of B_B
data_BB <- dplyr::filter(inter_annotated_gene, int == "B_B") %>% dplyr::select(ID_1, ID_2) %>% dplyr::filter(!duplicated(.))

# 1:
data_BB$ID_min <- pmin(data_BB$ID_1, data_BB$ID_2)
data_BB$ID_max <- pmax(data_BB$ID_1, data_BB$ID_2)

# 2: Count occurrences of each sorted pair
data_BB$pair <- paste(data_BB$ID_min, data_BB$ID_max, sep = "_")
pair_counts <- table(data_BB$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1]) # n = 18230


data_overlap <- filter(data_BB, pair %in% interactions_overlap)
prova <- data_overlap %>% group_by(pair) %>% summarise(id1 = length(unique(ID_1)),
                                                       id2 = length(unique(ID_2)))

table(prova$id1)
table(prova$id2)



# Check how many genes we have, and how many of them have bidirectional interactions with other Baits

# - Number of genes/communities
number_genes <- tapply(inter_annotated_gene$ID_1, inter_annotated_gene$cell_type, function(x) length(unique(x)) )

# - Number of genes that interact with other genes
inter_annotated_BB <- filter(inter_annotated_gene, int == "B_B")
number_genes_BB <- tapply(inter_annotated_BB$ID_1, inter_annotated_BB$cell_type, function(x) length(unique(x)) )

# - Number of genes that interact with other genes this interaction is bidirectional
inter_annotated_BB$ID_min <- pmin(inter_annotated_BB$ID_1, inter_annotated_BB$ID_2)
inter_annotated_BB$ID_max <- pmax(inter_annotated_BB$ID_1, inter_annotated_BB$ID_2)

inter_annotated_BB$pair <- paste(inter_annotated_BB$ID_min, inter_annotated_BB$ID_max, sep = "_")

inter_annotated_BB_bi <- filter(inter_annotated_BB, pair %in% interactions_overlap)

number_genes_BB_bi <- tapply(inter_annotated_BB_bi$ID_1, inter_annotated_BB_bi$cell_type, function(x) length(unique(x)) )


####################################################################################################################
### 4. Interactions from GNN file NOT split
####################################################################################################################

load("/home/acost1/BALL_project/results/gnn/edge_file_notsplit_updated.rds")
load("/home/acost1/BALL_project/results/gnn/seq_bed_file_notsplit_updated.rds")


# Add type of region 1
edge_file <- merge(edge_file, 
      dplyr::select(id_bed, id, cell_type, type),
      by.x = c("ID_1", "cell_type"),
      by.y = c("id", "cell_type"),
      all.x = TRUE) %>% dplyr::rename(type1 = type)


# Add type of region 2
edge_file <- merge(edge_file, 
                   dplyr::select(id_bed, id, cell_type, type),
                   by.x = c("ID_2", "cell_type"),
                   by.y = c("id", "cell_type"),
                   all.x = TRUE) %>% dplyr::rename(type2 = type)

# Build type of the interaction
edge_file$int <- paste(edge_file$type1, edge_file$type2, sep = "-")


# Find bidirectional interactions inside the edge_file of B_B
data_BB <- dplyr::filter(edge_file, int == "B-B") %>% dplyr::select(ID_1, ID_2) %>% dplyr::filter(!duplicated(.))

# 1:
data_BB$ID_min <- pmin(data_BB$ID_1, data_BB$ID_2)
data_BB$ID_max <- pmax(data_BB$ID_1, data_BB$ID_2)

# 2: Count occurrences of each sorted pair
data_BB$pair <- paste(data_BB$ID_min, data_BB$ID_max, sep = "_")
pair_counts <- table(data_BB$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1])


data_overlap <- filter(data_BB, pair %in% interactions_overlap)
prova <- data_overlap %>% group_by(pair) %>% summarise(id1 = length(unique(ID_1)),
                                                       id2 = length(unique(ID_2)))

table(prova$id1)
table(prova$id2)





####################################################################################################################
### 5. Interactions from GNN file WITH split
####################################################################################################################

load("/home/acost1/BALL_project/results/gnn/edge_file_split_updated.rds")
load("/home/acost1/BALL_project/results/gnn/seq_bed_file_split_updated.rds")


# Add type of region 1
edge_file_split <- merge(edge_file_split, 
                   dplyr::select(id_bed_split, id, cell_type, type),
                   by.x = c("ID_1", "cell_type"),
                   by.y = c("id", "cell_type"),
                   all.x = TRUE) %>% dplyr::rename(type1 = type)


# Add type of region 2
edge_file_split <- merge(edge_file_split, 
                   dplyr::select(id_bed_split, id, cell_type, type),
                   by.x = c("ID_2", "cell_type"),
                   by.y = c("id", "cell_type"),
                   all.x = TRUE) %>% dplyr::rename(type2 = type)

# Build type of the interaction
edge_file_split$int <- paste(edge_file_split$type1, edge_file_split$type2, sep = "-")


# Find bidirectional interactions inside the edge_file of B_B
data_BB <- dplyr::filter(edge_file_split, int == "B-B") %>% dplyr::select(ID_1, ID_2) %>% dplyr::filter(!duplicated(.))

# 1:
data_BB$ID_min <- pmin(data_BB$ID_1, data_BB$ID_2)
data_BB$ID_max <- pmax(data_BB$ID_1, data_BB$ID_2)

# 2: Count occurrences of each sorted pair
data_BB$pair <- paste(data_BB$ID_min, data_BB$ID_max, sep = "_")
pair_counts <- table(data_BB$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1])


data_overlap <- filter(data_BB, pair %in% interactions_overlap)
prova <- data_overlap %>% group_by(pair) %>% summarise(id1 = length(unique(ID_1)),
                                                       id2 = length(unique(ID_2)))

table(prova$id1)
table(prova$id2)




#####################################################################################
# 6. Communities and network
#####################################################################################

load("/home/acost1/BALL_project/results/gnn/edge_file_split_updated.rds")
load("/home/acost1/BALL_project/results/gnn/seq_bed_file_split_updated.rds")

id_bed_split$caract <- NA
id_bed_split$caract[id_bed_split$name_id == "EA" ] <- "EA"
id_bed_split$caract[id_bed_split$name_id == "EPo"] <- "EPo"
id_bed_split$caract[id_bed_split$name_id == "EPr"] <- "EPr"
id_bed_split$caract[id_bed_split$name_id == "Sil"] <- "Sil"
id_bed_split$caract[id_bed_split$type == "B"] <- "B"

# Number of genes/communities
tapply(edge_file_split$ID_1, edge_file_split$cell_type, function(x) length(unique(x)) )

# igraph function for a gene
gene_community_igraph <- function(gene_id, edge_file, id_bed_file, cell_type_subset){
  
  # Edge file and id_bed file only for a specific cell_type_subset
  edge_file <- dplyr::filter(edge_file, cell_type == cell_type_subset )%>%
                select(ID_1, ID_2, int)
  
  id_bed_file <- filter(id_bed_file, cell_type == cell_type_subset) %>% select(id, name_id, caract)
  id_bed_file$id <- as.numeric(id_bed_file$id)
  print("Step 1")
  
  # Select only interactions from gene_id of interest
  edges_subset <- edge_file[edge_file$ID_1 %in% gene_id,]
  
  # Select all baits and other ends from the previous interactions
  id_subset <- id_bed_file[id_bed_file$id %in% unique(c(edges_subset$ID_1, edges_subset$ID_2) ),]
  
  
  print("Step 2")
  
  # Check if other baits also interact with our gene_id of interest
  id_otherB <- edges_subset$ID_2[edges_subset$int == "B_B"] %>% unique()
  print(id_otherB)
  
  aux <- edge_file[edge_file$ID_1 %in% id_otherB,]
  print(nrow(aux))
  print("Step 3")
  
  # Add (if any) the other bidirectional interactions
  if(nrow(aux[aux$ID_2 == gene_id,]) != 0){
    edges_subset <- rbind(edges_subset, aux[aux$ID_2 == gene_id,] ) }
  
  graph <- graph_from_data_frame(d = edges_subset, vertices = id_subset, directed = TRUE)
  
  print("Step 4")
  
  # Plot the graph, using 'name_id' for visualization
  coul  <- brewer.pal(5, "Set1")
  
  coul  <- c(
    "#feb9b9", "#FF0000", "#ca1401",
    "#c9eed8", "#66CC00", "#003300",
    "#FFCCFF", "#FF99FF", "#990066",
    "#fef2da", "#f9d744", "#d6af06",
    "#dae7f6", "#0066FF", "#000066"
  )
  names(coul) <- c("BREMOVED", "B", "BADDED", 
                   "EAREMOVED", "EA", "EAADDED",
                   "EPoREMOVED", "EPo", "EPoADDED",
                   "EPrREMOVED", "EPr", "EPrADDED",
                   "SilREMOVED", "Sil", "SilADDED")
  
  my_color <- coul[as.numeric(as.factor(V(graph)$caract))]
  
  
  plot(
    graph,
    vertex.label = V(graph)$name_id, # Use 'name_id' for labels
    edge.arrow.size = 0.5,          # Adjust arrow size
    vertex.label.color="black", vertex.label.dist=1,
    vertex.size=8,
    vertex.color=coul[V(graph)$caract],
    layout = layout_with_kk, # Layout for better visualization
    main = cell_type_subset
  )
  
  
}


# plot rna expression for a gene
plot_rna_seq <- function(gene_ensembl){
  # 10. Load p-value of the gene for each cell type
  load("/home/acost1/BALL_project/results/limma/limma_results_sv.rds")
  gene_pvalues <- lapply(results, function(x) x[gene_ensembl,"padj"]) %>% unlist()
  
  
  # 11. Load expression of the gene for each cell type
  load("/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_all_log.rds")
  
  GROUP_specific <- sub("\\_.*", "", c("HSC_fetal_1", "HSC_fetal_2", "HSC_fetal_3",
                                       "PreProB_fetal_1", "PreProB_fetal_2", "PreProB_fetal_3",
                                       "ProB_fetal_1", "ProB_fetal_2", "ProB_fetal_3",
                                       "PreB_fetal_1", "PreB_fetal_2", "PreB_fetal_3",
                                       "immtransB_fetal_1", "immtransB_fetal_2", "immtransB_fetal_3",
                                       "nB_WT_1", "nB_WT_2", "nB_WT_3",
                                       "GCB_tonsil_1", "GCB_tonsil_3", "GCB_tonsil_4", "GCB_tonsil_5",
                                       "memB_WT_1", "memB_WT_2", "memB_WT_3", "memB_WT_4", 
                                       "PC_WT_1", "PC_WT_2", "PC_WT_3" ))
  GROUP_specific <- factor(GROUP_specific, 
                           levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"), 
                           labels = c("HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                      "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT"))  
  
  
  NormByTMM_selected <- tapply(NormByTMM_all_log[gene_ensembl,], GROUP_specific, median) 
  NormByTMM_selected <- NormByTMM_selected[names(NormByTMM_selected) %in% c("HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                                                               "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT")]
  names(gene_pvalues) <- names(NormByTMM_selected) 
  
  
  print("load results of limma and expression!")
  
  # 12. Data.frame with p-values and expression info
  print(round(gene_pvalues, 3) )
  
  gene_pvalues <- ifelse(gene_pvalues < 0.001, "***", 
                      ifelse(gene_pvalues < 0.01, "**", 
                             ifelse(gene_pvalues < 0.05, "*", "")))
  
  
  rna_info <- data.frame(
    cell_type = names(gene_pvalues),
    Row = 9:1 + 0.45,
    padj = gene_pvalues,
    expr = NormByTMM_selected)
  
  print(" Data.frame with p-values and expression info!")
  
  # 13. Heatmap plot with this info
  rna_info_expr <- data.table::melt(as.matrix(rna_info[,-c(1,2,3)]))
  rna_info_expr$Var1 <- all_celltypes
  print(rna_info_expr)
  
  ggp_rna_expr <- ggplot(rna_info_expr, aes(factor(Var2), factor(Var1))) +
    geom_tile(aes(fill = value), colour = "black") +
    geom_text(aes(label = round(value, 2 )), size = 3, colour = "white") + 
    geom_text(data = rna_info, aes(x = 1, y = Row, label = padj), 
              inherit.aes = FALSE, size = 5, color = "red") +
    scale_fill_viridis(discrete=FALSE, limits = c(0,16)) + 
    xlab("") + ylab("") + labs(fill="Expression") + 
    coord_equal()  +
    scale_y_discrete(limits = all_celltypes[9:1]) +
    theme(legend.position="top", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  ggp_rna_expr
  
}






###############
# BCL6 --> GCB
###############
id_BCL6 <- id_bed_split$id[grep("BCL6", id_bed_split$name_id) ] %>% unique()

all_celltypes <- names(table(edge_file_split$cell_type))


png("/home/acost1/BALL_project/results/networks/BCL6_community.png", width = 50, height = 20, units = "cm", res = 250)
par(mfrow = c(2, 5))
apply(data.frame(all_celltypes), 1, function(x){
  gene_community_igraph(id_BCL6[1], edge_file_split, id_bed_split, x)
})
dev.off()


png("/home/acost1/BALL_project/results/networks/BCL6_rna.png", width = 10, height = 10, units = "cm", res = 250)
plot_rna_seq("ENSG00000113916")
dev.off()



###############################
# IRF4 --> ProB, PreB, GCB, PC
###############################

id_IRF4 <- id_bed_split$id[grep("IRF4", id_bed_split$name_id) ] %>% unique()

png("/home/acost1/BALL_project/results/networks/IRF4_community.png", width = 50, height = 20, units = "cm", res = 250)
par(mfrow = c(2, 5))
apply(data.frame(all_celltypes), 1, function(x){
  gene_community_igraph(id_IRF4[1], edge_file_split, id_bed_split, x)
})
dev.off()


png("/home/acost1/BALL_project/results/networks/IRF4_rna.png", width = 10, height = 10, units = "cm", res = 250)
plot_rna_seq("ENSG00000137265")
dev.off()


############################
# IKAROS - IKZF1 all roadmap
############################


id_IKZF1 <- id_bed_split$id[grep("IKZF1", id_bed_split$name_id) ] %>% unique()

png("/home/acost1/BALL_project/results/networks/IKZF1_community.png", width = 50, height = 20, units = "cm", res = 250)
par(mfrow = c(2, 5))
apply(data.frame(all_celltypes), 1, function(x){
  gene_community_igraph(id_IKZF1[1], edge_file_split, id_bed_split, x)
})
dev.off()


png("/home/acost1/BALL_project/results/networks/IKZF1_rna.png", width = 10, height = 10, units = "cm", res = 250)
plot_rna_seq("ENSG00000185811")
dev.off()



#############################################
# AIOLOS - IKZF3, all roadmap except for HSC
#############################################

id_IKZF3 <- id_bed_split$id[grep("IKZF3", id_bed_split$name_id) ] %>% unique()

png("/home/acost1/BALL_project/results/networks/IKZF3_community.png", width = 50, height = 20, units = "cm", res = 250)
par(mfrow = c(2, 5))
apply(data.frame(all_celltypes), 1, function(x){
  gene_community_igraph(id_IKZF3[1], edge_file_split, id_bed_split, x)
})
dev.off()


png("/home/acost1/BALL_project/results/networks/IKZF3_rna.png", width = 10, height = 10, units = "cm", res = 250)
plot_rna_seq("ENSG00000161405")
dev.off()


#########################################
# PAX5, all roadmap except for GCB and PC
##########################################

id_PAX5 <- id_bed_split$id[grep("PAX5", id_bed_split$name_id) ] %>% unique()

png("/home/acost1/BALL_project/results/networks/PAX5_community.png", width = 50, height = 20, units = "cm", res = 250)
par(mfrow = c(2, 5))
apply(data.frame(all_celltypes), 1, function(x){
  gene_community_igraph(id_PAX5[2], edge_file_split, id_bed_split, x)
})
dev.off()


png("/home/acost1/BALL_project/results/networks/PAX5_rna.png", width = 10, height = 10, units = "cm", res = 250)
plot_rna_seq("ENSG00000196092")
dev.off()






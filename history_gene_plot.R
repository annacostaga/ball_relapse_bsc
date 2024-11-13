# Function to reconstruct the history of a gene

history_gene <- function(inter_annotated, chromm_list, gene_ensembl, gene_symbol){
  
  ## - Interaction of a gene
  inter_annotated_gene <- interactionsByBaits(interactions = inter_annotated,
                                              baits = gene_symbol) %>% 
    as.data.frame() # n = 35
  
  
  if(nrow(inter_annotated_gene) == 0) return(NULL)
  
      ## - Get the coordinates of the gene
  gene_coordinates <- apply(inter_annotated_gene, 1, function(x) {
    if(sum(gene_symbol ==  unlist(str_split(x["B.id1"] ,",") ) ) > 0)
      paste0("chr", x["seqnames1"],":", format(x["start1"], big.mark = ","), "-",
             format(x["end1"], big.mark = ",") ) 
    
  } ) %>% paste() 
  gene_coordinates <- unique(gene_coordinates)
  gene_coordinates <- gene_coordinates[gene_coordinates != "NULL"]
  
  
      ## - Selecting the interactions
  inter_annotated_gene <- apply(inter_annotated_gene, 1, function(x) {
    if(grepl(gene_symbol, x["B.id1"]) == TRUE ){
      aux1 <- c(x["seqnames2"], x["start2"], x["end2"], x["B.id2"], x["node.class2"],
                x["CS_HSC"], x["CS_PreProB"], x["CS_ProB"], x["CS_PreB"], x["CS_immtransB"], 
                x["CS_nB1Mnew"], x["CS_GCB"], x["CS_memB"], x["CS_PC"]) 
      names(aux1)[1:5] <- c("chr", "start", "end", "B.id", "node.class")
      aux1
    } else{
      aux2 <- c(x["seqnames1"], x["start1"], x["end1"], x["B.id1"], x["node.class1"],
                x["CS_HSC"], x["CS_PreProB"], x["CS_ProB"], x["CS_PreB"], x["CS_immtransB"], 
                x["CS_nB1Mnew"], x["CS_GCB"], x["CS_memB"], x["CS_PC"])
      names(aux2)[1:5] <- c("chr", "start", "end", "B.id", "node.class")
      aux2
    }
    
  } ) %>% bind_rows() 
  
  
  counter <- 0
  yes_interactions <- dplyr::select(inter_annotated_gene, CS_HSC, CS_PreProB, CS_ProB, CS_PreB, CS_immtransB, CS_nB1Mnew, CS_GCB, CS_memB, CS_PC) %>%
    apply(1, function(x) {
      counter <<- counter + 1
      if(sum(x >= 5) >= 1) counter }
    ) %>% unlist() %>% as.numeric()
  inter_annotated_gene <- inter_annotated_gene[yes_interactions, ] 
  
  
  print("load interactions!")
  
  # 2. Num. interactions
  inter_annotated_gene$num_interaction <- 1:nrow(inter_annotated_gene)
  print("print num interactions!")
  
  # 3. Chicago score from wide to long
  inter_annotated_gene <- inter_annotated_gene %>%
    pivot_longer(
      cols = starts_with("CS_"),
      names_to = "cell_type",
      names_prefix = "CS_",
      values_to = "chicago_score",
      values_drop_na = TRUE
    )
  print("chicago score long!")
  
  # 4. Load ChromHMM files
  
  hsc_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/HSC05_15_dense_fin.bed",  col_names = FALSE,
                             col_types = cols(X1 = col_character()))
  hsc_chromhmm <- hsc_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr:: rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "HSC")
  
  prepro_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/PreProB05_15_dense_fin.bed",  col_names = FALSE,
                                col_types = cols(X1 = col_character()))
  prepro_chromhmm <- prepro_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr:: rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "PreProB")
  
  pro_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/ProB05_15_dense_fin.bed",  col_names = FALSE,
                             col_types = cols(X1 = col_character()))
  pro_chromhmm <- pro_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr:: rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "ProB")
  
  pre_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/PreB05_15_dense_fin.bed",  col_names = FALSE,
                             col_types = cols(X1 = col_character()))
  pre_chromhmm <- pre_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr:: rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "PreB")
  
  trans_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/immtransB05_15_dense_fin.bed",  col_names = FALSE,
                               col_types = cols(X1 = col_character()))
  trans_chromhmm <- trans_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr:: rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "immtransB")
  
  nb_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/nB05_15_dense_fin.bed",  col_names = FALSE,
                            col_types = cols(X1 = col_character()))
  nb_chromhmm <- nb_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) )  %>%
    mutate(cell_type = "nB1Mnew")
  
  gcb_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/GCB05_15_dense_fin.bed",  col_names = FALSE,
                             col_types = cols(X1 = col_character()))
  gcb_chromhmm <- gcb_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "GCB")
  
  memb_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/memB05_15_dense_fin.bed",  col_names = FALSE,
                              col_types = cols(X1 = col_character()))
  memb_chromhmm <- memb_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "memB")
  
  pc_chromhmm <- read_table("/home/acost1/BALL_project/data/ChromHMM/PC05_15_dense_fin.bed",  col_names = FALSE,
                            col_types = cols(X1 = col_character()))
  pc_chromhmm <- pc_chromhmm %>% dplyr:: select(X1, X2, X3, X4) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, category = X4) %>%
    filter(chr %in% unique(inter_annotated_gene$chr) ) %>%
    mutate(cell_type = "PC")
  
  chromhmm <- rbind(hsc_chromhmm, prepro_chromhmm, pro_chromhmm, pre_chromhmm, trans_chromhmm,
                    nb_chromhmm, gcb_chromhmm, memb_chromhmm, pc_chromhmm)
  
  rm(hsc_chromhmm, prepro_chromhmm, pro_chromhmm, pre_chromhmm, trans_chromhmm,
     nb_chromhmm, gcb_chromhmm, memb_chromhmm, pc_chromhmm)
  
  
  print("load chromhmm files!")
  
  
  # 5. select chromhmm categories within the other end/bait of the interactions
  interaction_region <- inter_annotated_gene %>% dplyr::select(chr, start, end, num_interaction) %>%
    filter(!duplicated(.))
  
  counter <- 0
  chromhmm_interaction <- apply(data.frame(interaction_region), 1, function(x) {
    counter <<- counter + 1;
    
    aux1 <- chromhmm[chromhmm$start == as.numeric(x["start"]) & chromhmm$start == as.numeric(x["end"]),] # same interacting region
    aux2 <- chromhmm[chromhmm$start > as.numeric(x["start"]) & chromhmm$start < as.numeric(x["end"]),] # start position of chromhmm between the start and end position of the interacting region
    aux3 <- chromhmm[chromhmm$end > as.numeric(x["start"]) & chromhmm$end < as.numeric(x["end"]),] # # end position of chromhmm between the start and end position of the interacting region
    aux4 <- chromhmm[chromhmm$start > as.numeric(x["start"]) & chromhmm$end < as.numeric(x["end"]),] # chromhmm within the interacting region
    aux5 <- chromhmm[chromhmm$start < as.numeric(x["start"]) & chromhmm$end > as.numeric(x["end"]),] # start of chromhmm before the start, end of chromhmm after the end of the interacting region
    
    aux <- rbind(aux1, aux2, aux3, aux4, aux5) %>% filter(!duplicated(.))
    
    aux$start[aux$start < as.numeric(x["start"])] <- as.numeric(x["start"])
    aux$end[aux$end > as.numeric(x["end"]) ] <-  as.numeric(x["end"])
    
    
    
    groups <- aux %>% filter(category != "None") %>%
      mutate(size = end - start) %>%
      group_by(cell_type, category) %>%
      summarise(size_n = sum(size))  %>%
      mutate(size = max(size_n) ) %>% filter(size_n == size) %>%
      dplyr::select(-c(size_n, size)) %>%
      group_by(cell_type) %>%
      summarise(category = str_c(category, collapse=" ")) 
    
    groups_def <- aux %>% arrange(cell_type, start) %>%
      group_by(cell_type) %>%
      summarise(category_c = paste0(category, collapse = "-")) %>%
      mutate(num_interaction = counter)  %>%
      merge(groups, by = "cell_type", all.x = TRUE)
   
    groups_def$category[groups_def$category_c == "None"] <- "None"
    
    aux <- NULL
    for(i in 1:length(groups_def$category_c)) aux[i] <- (sum(str_split(groups_def$category_c, "-")[[i]] == "None") ==  length(str_split(groups_def$category_c, "-")[[i]] ) )
    
    groups_def$category[aux] <- "None"
    
    return(groups_def)
    
  } ) %>% bind_rows()
  
  print(chromhmm_interaction)
  
  if(any(chromhmm_interaction$category == "None")){
  chromhmm_interaction$category <- factor(chromhmm_interaction$category)
  chromhmm_interaction$category <- relevel(chromhmm_interaction$category, ref = "None") }
  
  if(any(chromhmm_interaction$category_c == "None")){
  chromhmm_interaction$category_c <- factor(chromhmm_interaction$category_c)
  chromhmm_interaction$category_c <- relevel(chromhmm_interaction$category_c, ref = "None") }
  
  
  print("chromhmm into interactions!")
  
  # 6. put chromHmm categories inside data.frame with Chicago score
  inter_annotated_gene <- merge(inter_annotated_gene, chromhmm_interaction, 
                                by= c("cell_type", "num_interaction"), all.x = TRUE)
  
  inter_annotated_gene$chicago_score_c <- ifelse(inter_annotated_gene$chicago_score >= 5, 1, 0)
  inter_annotated_gene$cell_type <- factor(inter_annotated_gene$cell_type, 
                                           levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB1Mnew", "GCB", "memB", "PC"),
                                           labels = c("HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal", "immtransB-fetal", "nB-WT",
                                                      "GCB-tonsil", "memB-WT", "PC-WT")
  )
  inter_annotated_gene[is.na(inter_annotated_gene)] <- 0
  
  
  print(inter_annotated_gene)
  print("chromhmm into interactions data.frame!")
  
  
  # 7. Heatmap with CHiCAGO scores OK
  ggp <- ggplot(inter_annotated_gene,  aes(factor(num_interaction), cell_type)) + 
    geom_tile(aes(fill = chicago_score)) +
    scale_fill_gradient2(midpoint = 5, mid="#eee8d5", high="#dc322f", low="#268bd2") + 
    geom_text(aes(label = round(chicago_score,2)), size = 3, colour = "#333333") +
    scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1])  + 
    xlab("") + ylab("") + labs(fill="CHiCAGO \n score", title = gene_symbol,
                               subtitle = gene_coordinates ) +
    coord_equal()  + 
    theme(legend.position="left", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # Heatmap with different colors depending on the significance of chicago score
  # ggp <- ggplot(inter_annotated_gene,  aes(factor(num_interaction), cell_type)) + 
  #   geom_tile(aes(fill = chicago_score)) +
  #   scale_fill_manual(limits = chicago_score_c, values = c("grey", "#268bd2")) + 
  #   geom_text(aes(label = round(chicago_score, 2)), size = 3, colour = "#333333") +
  #   scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1])  + 
  #   xlab("") + ylab("") + labs(fill="CHiCAGO \n score", title = gene_symbol,
  #                              subtitle = gene_coordinates ) +
  #   coord_equal()  + 
  #   theme(legend.position="left", 
  #         axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())
  print("heatmap with chicago scores!")
  
  # 8. Heatmap with ChromHMM
  # X axis labels
  x_axis <- apply(interaction_region[,-4], 1, function(x)
    paste0("chr", x["chr"], ":", format(as.numeric(x["start"]), big.mark = ","), "-", 
           format(as.numeric(x["end"]), big.mark = "," ) ) )
  
  # OK option 1: heatmap with colour
  colors_legend <- dplyr::select(data.frame( chromm_list[[1]]), itemRgb, name) %>% filter(!duplicated(.))
  colors_legend$name <- factor(colors_legend$name)
  colors_legend$name <- relevel(colors_legend$name, ref = "None")
  colors_legend <- arrange(colors_legend, name)


  ggp_cat_colour <- ggplot(inter_annotated_gene, aes(factor(num_interaction), cell_type)) +
    geom_tile(aes(fill = factor(category))) + 
    xlab("") + ylab("") + labs(fill="ChromHMM \n categories") +
    coord_equal()  +
    scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1]) +
    scale_fill_manual(limits = colors_legend$name, values = colors_legend$itemRgb) +
    scale_x_discrete(labels=x_axis) + 
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          legend.position="left")
  
  # Option 2: heatmap with the written categories
  inter_annotated_gene$category_2c <- ifelse(inter_annotated_gene$category != "None", 1, 0) %>% as.factor()
  ggp_cat_label <- ggplot(inter_annotated_gene, aes(factor(num_interaction), cell_type)) +
    geom_tile(aes(fill = category_2c)) + 
    geom_text(aes(label = category_c), size = 1.5, colour = "#333333") + 
    xlab("") + ylab("") + 
    coord_equal()  +
    scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1]) +
    scale_fill_manual(values=c("azure2", "lightskyblue2")) + 
    theme(legend.position="none") +
    scale_x_discrete(labels=x_axis) + 
    theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  # ggp, ggp_nodeclass,  ggp_cat_colour, ggp_rna_expr
  
  # option 3: heatmap for each type of ChromHMM category
  # gp_list <- NULL
  # counter <- 0
  # gp_list <- apply(inter_annotated_gene[,10:20], 2, function(x){
  #   counter <<- counter + 1
  #   ggplot(inter_annotated_gene,  aes(factor(num_interaction), cell_type)) + 
  #     geom_tile(aes(fill = factor(x)), colour = "black") +
  #     scale_fill_viridis(discrete=TRUE) + 
  #     scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1])  + 
  #     xlab("") + ylab("")+ labs(fill = colnames(inter_annotated_gene)[10:22][counter]) +
  #     coord_equal() 
  #   #ggsave(paste0("/home/acost1/BALL_project/results/lichi_chromhmm/chromhmm", colnames(inter_annotated_gene)[i], ".png"), 
  #   #       width = 35, height = 15, units = "cm")
  # }
  # )
  # 
  # gp_list <- gp_list[-c(1:9)]
  # ggarrange(plotlist=c(list(ggp), gp_list), ncol = 1)
  
  # 9. Join CHiCAGO and ChromHMM heatmaps
  # ggp_lichic1 <- egg::ggarrange(ggp, ggp_cat_label, ncol = 1)
  
  # ggp_lichic2 <- egg::ggarrange(ggp, ggp_cat_colour, ncol = 1)
    
  print("heatmap with chromhmm!")
  
  # 10. Load p-value of the gene for each cell type
  load("/home/acost1/BALL_project/results/limma/limma_results_sv.rds")
  gene_pvalues <- lapply(results, function(x) -log10(x[gene_ensembl,"padj"])) %>% unlist()
  
  
  # 11. Load expression of the gene for each cell type
  load("/home/acost1/BALL_project/data/RNA/counts_rlog_normalized.rds")
  
  GROUP_specific <- sub("\\_.*", "", c("CMP_fetal_1", "CMP_fetal_2", "CMP_fetal_3", "CMP_fetal_4",
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
                                       "PC_WT_1", "PC_WT_2", "PC_WT_3" ))
  GROUP_specific <- factor(GROUP_specific, 
                            levels = c("CMP", "nCD8", "Mon", "HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC"), 
                            labels = c("CMP-fetal", "nCD8-WT", "Mon-WT", "HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                       "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT"))  
  
  NormByRlog_selected <- tapply(NormByRlog[gene_ensembl,], GROUP_specific, median) 
  NormByRlog_selected <- NormByRlog_selected[names(NormByRlog_selected) %in% c("HSC-fetal", "PreProB-fetal", "ProB-fetal", "PreB-fetal",
                                                                               "immtransB-fetal", "nB-WT", "GCB-tonsil", "memB-WT", "PC-WT")]
  names(gene_pvalues) <- names(NormByRlog_selected) 
  
  
  print("load results of limma and expression!")
  
  # 12. Data.frame with p-values and expression info
  rna_info <- data.frame(
    cell_type = names(gene_pvalues),
    padj = gene_pvalues,
    expr = NormByRlog_selected)
  
  print(" Data.frame with p-values and expression info!")
  
  # 13. Heatmap plot with this info
  rna_info_expr <- data.table::melt(as.matrix(rna_info[,-c(1,2)]))
  rna_info_expr$Var1 <- levels(inter_annotated_gene$cell_type)
  print(rna_info_expr)
  
  ggp_rna_expr <- ggplot(rna_info_expr, aes(factor(Var2), factor(Var1))) +
    geom_tile(aes(fill = value), colour = "black") +
    geom_text(aes(label = round(value, 2 )), size = 3, colour = "white") + 
    scale_fill_viridis(discrete=FALSE, limits = c(1.5,16)) + 
    xlab("") + ylab("") + labs(fill="Expression") + 
    coord_equal()  +
    scale_y_discrete(limits = levels(inter_annotated_gene$cell_type)[9:1]) +
    theme(legend.position="top", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  # OK - MIRAR AIXÒ
  ggp_nodeclass <- inter_annotated_gene %>% filter(cell_type == "GCB-tonsil") %>%
    ggplot(aes(factor(num_interaction), cell_type)) +
    geom_tile(aes(fill = node.class), colour = "black") +
    geom_text(aes(label = node.class), size = 4) +
    scale_fill_manual(values=c("gray", "white")) +
    xlab("") + ylab("") + 
    theme(legend.position="none", axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    coord_equal()  
  
  
   # ggp_lichic1 <- egg::ggarrange(ggp, ggp_nodeclass, ggp_cat_colour, ncol = 1)
  
  
  
  # print(ggp_lichic1)
  
  print("heatmap of rna!")
  
  # Plot a track
  
  hsc_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/HSC05_15_dense_fin.bed")
  prepro_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PreProB05_15_dense_fin.bed")
  pro_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/ProB05_15_dense_fin.bed")
  pre_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PreB05_15_dense_fin.bed")
  trans_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/immtransB05_15_dense_fin.bed")
  nb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/nB05_15_dense_fin.bed")
  gcb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/GCB05_15_dense_fin.bed")
  memb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/memB05_15_dense_fin.bed")
  pc_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PC05_15_dense_fin.bed")


  chromm_list <- list(hsc_chromhmm, prepro_chromhmm, pro_chromhmm, pre_chromhmm, trans_chromhmm,
                      nb_chromhmm, gcb_chromhmm, memb_chromhmm, pc_chromhmm)


  cell_type <- c("HSC", "PProB", "ProB", "PreB", "transB", "nB", "GCB", "memB", "PC")
  counter2 <- 0

  png("/home/acost1/BALL_project/results/lichi_chromhmm/plot_track.png",
      width = 50, height = 40, units = "cm", res = 200)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, max(inter_annotated_gene$num_interaction), heights=c(0.5,0.5))))
  apply(interaction_region, 1, function(y) {
    y["chr"] <- as.numeric(y["chr"])
    counter <- 0
    counter2 <<- counter2 + 1
    annotation_track <- lapply(chromm_list, function(x){
      counter <<- counter + 1
      AnnotationTrack(x , chr = y["chr"],
                      id = x$name, width = 5, , genome= "GRCh38", stacking = "dense",
                      feature = x$itemRgb, name = cell_type[counter], groupAnnotation="feature")
    } )

    if(counter2 == 1){
      pushViewport(viewport(layout.pos.col = counter2, layout.pos.row = 1))
      plotTracks(annotation_track, chromosome = as.numeric(y["chr"]),   from =  as.numeric(y["start"]),  to = as.numeric( y["end"]),
                 "#737272" = "#737272", "#AD6DDE"= "#AD6DDE",  "#3659B3"= "#3659B3",
                 "#8B769C" = "#8B769C", "#947160" = "#947160", "#CC7306"= "#CC7306",
                 "#D44950" = "#D44950", "#59D96A" = "#59D96A", "#D4493C"= "#D4493C",
                 "#D44964" = "#D44964", "#288A35" = "#288A35", "#EDE658" = "#EDE658",
                 "#9C3B7F" = "#9C3B7F", "#C4DEA9" = "#C4DEA9", "#A60C03" = "#A60C03",
                 add = TRUE, title.width = 0.8, main = x_axis[counter2], cex.main=1.2 )
      popViewport(1)
    }

    if(counter2 != 1){
      pushViewport(viewport(layout.pos.col = counter2, layout.pos.row = 1))
      plotTracks(annotation_track, chromosome = as.numeric(y["chr"]),   from =  as.numeric(y["start"]),  to = as.numeric( y["end"]),
                 "#737272" = "#737272", "#AD6DDE"= "#AD6DDE",  "#3659B3"= "#3659B3",
                 "#8B769C" = "#8B769C", "#947160" = "#947160", "#CC7306"= "#CC7306",
                 "#D44950" = "#D44950", "#59D96A" = "#59D96A", "#D4493C"= "#D4493C",
                 "#D44964" = "#D44964", "#288A35" = "#288A35", "#EDE658" = "#EDE658",
                 "#9C3B7F" = "#9C3B7F", "#C4DEA9" = "#C4DEA9", "#A60C03" = "#A60C03",
                 add = TRUE, title.width = 0,  main = x_axis[counter2], cex.main=1.2)
      popViewport(1)
    } } )
  dev.off()

  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend = colors_legend$name, pch=15, pt.cex=1.5, cex=0.75, bty='n',
         col = colors_legend$itemRgb)
  mtext("ChromHmm states", at=0.04, cex=1)

  print("plot track!")
  
  
  if(max(inter_annotated_gene$num_interaction) <= 19){
    heights <-  c(0.1,0.1,0.8)
    widths <- c(5,0.4)
    height <- 3000
    width <- 4000
    res <- 250
    } else if( max(inter_annotated_gene$num_interaction) >= 20 & 
      max(inter_annotated_gene$num_interaction) <= 99
      ){
    heights <-  c(0.1,0.1,0.8)
    widths <- c(5,0.1)
    height <- 3000
    width <- 10000
    res <- 250
    } else {
    heights <-  c(0.1,0.1,0.8)
    widths <- c(5,0.1)
    height <- 4100
    width <- 13200
    res <- 250}
  
  
  final_figure <- egg::ggarrange(ggp, ggp_nodeclass, 
                                 ggp_cat_colour, ggp_rna_expr, nrow = 3, ncol = 2, 
                           byrow = FALSE, 
                           heights = heights,
                           widths = widths)
  
  png(paste0("/home/acost1/BALL_project/results/lichi_chromhmm/figures_biola/",
              gene_symbol, ".png"), width = width, height = height, res = res)
  print(final_figure)
  dev.off()
  
  
  
  
  return(list(final_figure, inter_annotated_gene))
  
  }



# Libraries
library(HiCaptuRe)
library(readr)
library(dplyr)
library(biomaRt)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(anticlust)
library(ggpubr)
library(Biobase)
library(DESeq2)
library(grid)
library(Gviz)




# 1. Load general interaction and for a specific gene
inter <- load_interactions("/home/acost1/BALL_project/data/lichic/lichic_peakmatrix_B_cell_roadmap.txt")

## - Put annotation
inter_annotated <- annotate_interactions(interactions = inter,
                                         annotation = "/home/acost1/BALL_project/data/lichic/annotation_example.txt")

hsc_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/HSC05_15_dense_fin.bed")
prepro_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PreProB05_15_dense_fin.bed")
pro_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/ProB05_15_dense_fin.bed")
pre_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PreB05_15_dense_fin.bed")
trans_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/immtransB05_15_dense_fin.bed")
nb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/nB05_15_dense_fin.bed")
gcb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/GCB05_15_dense_fin.bed")
memb_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/memB05_15_dense_fin.bed")
pc_chromhmm <- import.bed("/home/acost1/BALL_project/data/ChromHMM/PC05_15_dense_fin.bed")

chromm_list <- list(hsc_chromhmm, prepro_chromhmm, pro_chromhmm, pre_chromhmm, trans_chromhmm,
                    nb_chromhmm, gcb_chromhmm, memb_chromhmm, pc_chromhmm)

rm(hsc_chromhmm, prepro_chromhmm, pro_chromhmm, pre_chromhmm, trans_chromhmm,
   nb_chromhmm, gcb_chromhmm, memb_chromhmm, pc_chromhmm)



# 2. Relevant genes for the diferentiation

# Function with all the genes

genes_MUST <- read_csv("BALL_project/data/RNA/genes_MUST.csv")


# Problems with gene symbols
genes_MUST$initial_alias[genes_MUST$initial_alias == "BLIMP1"] <- "PRDM1"
genes_MUST$initial_alias[genes_MUST$initial_alias == "C-KIT"] <- "KIT"
genes_MUST$initial_alias[genes_MUST$initial_alias == "C-MYC"] <- "MYC"
genes_MUST$initial_alias[genes_MUST$initial_alias == "CD10"] <- "MME"
genes_MUST$initial_alias[genes_MUST$initial_alias == "CD20"] <- "MS4A1"
genes_MUST$initial_alias[genes_MUST$initial_alias == "CD21"] <- "CR2"
genes_MUST$initial_alias[genes_MUST$initial_alias == "CXCR7"] <- "ACKR3"
genes_MUST$initial_alias[genes_MUST$initial_alias == "E2A"] <- "TCF3"
genes_MUST$initial_alias[genes_MUST$initial_alias == "GLT"] <- "SLC1A2"
genes_MUST$initial_alias[genes_MUST$initial_alias == "IL-7R"] <- "IL7R"
genes_MUST$initial_alias[genes_MUST$initial_alias == "NF-ΚB1"] <- "NFKB1"
genes_MUST$initial_alias[genes_MUST$initial_alias == "NF-ΚB2"] <- "NFKB2"
genes_MUST$initial_alias[genes_MUST$initial_alias == "OBF1"] <- "POU2AF1"
genes_MUST$initial_alias[genes_MUST$initial_alias == "OCT2"] <- "SLC22A2"
genes_MUST$initial_alias[genes_MUST$initial_alias == "PU.1"] <- "SPI1"
genes_MUST$initial_alias[genes_MUST$initial_alias == "SOX-4"] <- "SOX4"
genes_MUST$initial_alias[genes_MUST$initial_alias == "SPI-B"] <- "SPIB"
genes_MUST$initial_alias[genes_MUST$initial_alias == "TDT"] <- "DNTT"

# AID --> not found
# CD24 --> Not found
# STAT3 -->  not found



apply(genes_MUST, 1,  function(x) history_gene(inter_annotated, chromm_list, x[2], x[1]) )

# history_gene(inter_annotated, chromm_list, "ENSG00000196549", "MME" )



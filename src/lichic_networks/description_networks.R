library(HiCaptuRe)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(viridis)

library(ggraph)
library(graphlayouts)
library(ggforce)
library(scatterpie)

library(patchwork)

library(Gmisc)
library(grid)

library(tidyr)
library(ggpubr)

library(rstatix)


library(GenomicRanges)


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



# Number of B-B and OE interactions
tapply(inter_annotated$pair, inter_annotated$int, function(x) length(unique(x)))


# Number of all interactions between baits n = 64.977
inter_annotated_BB <- filter(inter_annotated, int == "B_B")
length(unique(inter_annotated_BB$pair))

# Number of bidirectional interactions between baits n = 21.348
inter_annotated_BB$pair[inter_annotated_BB$pair %in% interactions_overlap ] %>% unique() %>% length()

# Number of unidirectional interactions between baits n = 43.629
inter_annotated_BB$pair[!inter_annotated_BB$pair %in% interactions_overlap ] %>% length()



# Look for the distribution of the chicago score between directional and unidirectional B-B interactions
#########################################################################################################

# Colname of chicago scores
file <- "/home/acost1/BALL_project/data/lichic/lichic_peakmatrix_B_cell_roadmap.txt"
data <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")

pos_inter <- grep("^CS_", colnames(inter_annotated_BB))
colnames(inter_annotated_BB)[pos_inter] <- paste0("CS_", colnames(data)[12:ncol(data)] )


inter_annotated_BB$direction <- ifelse(inter_annotated_BB$pair %in% interactions_overlap, 1, 0)


# 1) Handle bidirectional chicago scores: 

# - Treat the 2 scores (of each direction) separately 
# - Combine them into single metric with mean
inter_annotated_BB_merged <- inter_annotated_BB %>% group_by(pair, direction) %>%
  summarise(CS_CMP = mean(CS_CMP),
         CS_GCB = mean(CS_GCB),
         CS_HSC = mean(CS_HSC),
         CS_immtransB = mean(CS_immtransB),
         CS_memB = mean(CS_memB),
         CS_PC = mean(CS_PC),
         CS_PreB = mean(CS_PreB),
         CS_PreProB = mean(CS_PreProB),
         CS_nB1Mnew = mean(CS_nB1Mnew),
         CS_nCD8 = mean(CS_nCD8),
         CS_ProB = mean(CS_ProB),
         CS_Mon = mean(CS_Mon),
         CS_nCD4 = mean(CS_nCD4))


# 2) Statistical comparison

# wide to long data
inter_annotated_BB_long <- inter_annotated_BB %>% dplyr::select(pair, direction, starts_with("CS_")) %>%
  pivot_longer(
    cols = starts_with("CS_"),
    names_to = "cell_type",
    names_prefix = "CS_",
    values_to = "chicago_score",
    values_drop_na = TRUE
  )


inter_annotated_BB_long$direction <- factor(inter_annotated_BB_long$direction, levels = c(0, 1),
                                            labels = c("Unidirectional", "Bidirectional"))

inter_annotated_BB_long$cell_type <- factor(inter_annotated_BB_long$cell_type,
       levels = c("HSC", "PreProB", "ProB", "PreB", 
                  "immtransB", "nB1Mnew", "GCB",
                  "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"),
       labels = c("HSC", "PreProB", "ProB", "PreB", 
                  "immtransB", "nB1Mnew", "GCB",
                  "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"))



inter_annotated_BB_merged_long <- inter_annotated_BB_merged  %>%
  pivot_longer(
    cols = starts_with("CS_"),
    names_to = "cell_type",
    names_prefix = "CS_",
    values_to = "chicago_score",
    values_drop_na = TRUE
  )


inter_annotated_BB_merged_long$direction <- factor(inter_annotated_BB_merged_long$direction, levels = c(0, 1),
                                            labels = c("Unidirectional", "Bidirectional"))

inter_annotated_BB_merged_long$cell_type <- factor(inter_annotated_BB_merged_long$cell_type,
                                            levels = c("HSC", "PreProB", "ProB", "PreB", 
                                                       "immtransB", "nB1Mnew", "GCB",
                                                       "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"),
                                            labels = c("HSC", "PreProB", "ProB", "PreB", 
                                                       "immtransB", "nB1Mnew", "GCB",
                                                       "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"))


# Boxplot and Statistical testing 

# - Treating 2 scores from bidirectional interactions separately 
test = inter_annotated_BB_long %>% group_by(cell_type) %>% wilcox_test(chicago_score~direction) %>%
  add_xy_position(x = "cell_type", step.increase = 0.03)

test$p_cat <- ifelse(test$p <= 0.001, "***", ".")

means <- inter_annotated_BB_long %>% group_by(cell_type, direction)%>%summarize(chicago_score = mean(chicago_score))%>%as.data.frame()
means$chicago_score <- round(means$chicago_score,2)


png("/home/acost1/BALL_project/results/networks/chicago_inidividual_biuni_interactions.png", width = 50, height = 30, units = "cm", res = 250)
ggboxplot(inter_annotated_BB_long, x = "cell_type", y = "chicago_score", col = "direction", group = "direction",
           xlab = "Cell type", ylab ="Chicago score (individual)", add = "mean") +
           stat_pvalue_manual(test, label = "p_cat",  tip.length = 0.01)+
           geom_hline(yintercept = 5, size = 1, alpha = 0.5)+
           scale_color_manual(values = c("#00AFBB", "#FC4E07"), name = "B-B interaction") +
  ylim(-1, 150) + 
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2), group = direction), 
               color = "black", size = 4, vjust = -2.5, 
               position = position_dodge(0.8)) 
dev.off()

# chicago score median higher in bidirectional interactions


# We can try to normalize the distribution of chicago scores for each cell type
inter_annotated_BB_long <- inter_annotated_BB_long %>% group_by(cell_type) %>% mutate(chicago_score_c = log10(chicago_score + 1))

ggboxplot(inter_annotated_BB_long, x = "cell_type", y = "chicago_score_c", col = "direction", group = "direction",
          xlab = "Cell type", ylab ="Chicago score", add = "mean")
# still right skewed

# - Merging 2 scores from bidirectional interactions separately 
test = inter_annotated_BB_merged_long %>% group_by(cell_type) %>% wilcox_test(chicago_score~direction) %>%
  add_xy_position(x = "cell_type", step.increase = 0.03)

test$p_cat <- ifelse(test$p <= 0.001, "***", ".")

means <- inter_annotated_BB_merged_long %>% group_by(cell_type, direction)%>%summarize(chicago_score = mean(chicago_score))%>%as.data.frame()
means$chicago_score <- round(means$chicago_score,2)


png("/home/acost1/BALL_project/results/networks/chicago_merged_biuni_interactions.png", width = 50, height = 30, units = "cm", res = 250)
ggboxplot(inter_annotated_BB_merged_long, x = "cell_type", y = "chicago_score", col = "direction", group = "direction",
          xlab = "Cell type", ylab ="Chicago score (merged)", add = "mean") +
  stat_pvalue_manual(test, label = "p_cat",  tip.length = 0.01)+
  geom_hline(yintercept = 5, size = 1, alpha = 0.5)+
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), name = "B-B interaction")  + 
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2), group = direction), 
               color = "black", size = 4, vjust = -2.5, 
               position = position_dodge(0.8)) 
dev.off()





# 3) Compare bidirectional scores internally

# - Prepare data (compute the difference and ratio)
inter_annotated_BB_long2 <- inter_annotated_BB_long %>% 
  filter(direction == "Bidirectional") %>%
  arrange(cell_type, pair, chicago_score) %>%
  group_by(cell_type, pair) %>%
  arrange(chicago_score) %>%
  dplyr::slice(c(1, n())) %>%  # Keep the first and last rows
  summarise(score_1 = chicago_score[1], score_2 = chicago_score[2], .groups = "drop") %>%
  mutate(dif = abs(score_2 - score_1),
         ratio = score_1 / score_2,
         log_ratio = log10(ratio))


# - Scatterplot x --> score_1, y --> score_2
png("/home/acost1/BALL_project/results/networks/bidirection_symmetry.png", width = 50, height = 30, units = "cm", res = 250)
ggplot(inter_annotated_BB_long2, aes(x = log10(score_1 + 1), y = log10(score_2 + 1))) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # Line of symmetry
  theme_minimal() +
  #  xlim(0,30) + ylim(0,30) +
  facet_wrap(~ cell_type) + # Faceted by category
  labs(
    title = "Bidirectional interaction chicago scores by cell type",
    x = "log10(Score in Direction 1 + 1)",
    y = "log10(Score in Direction 2 + 1)"
  )
dev.off()


# - Distribution of the abs(difference) between scores 



medians <- tapply(inter_annotated_BB_long2$dif, inter_annotated_BB_long2$cell_type, function(x) median(x, na.rm = TRUE))

png("/home/acost1/BALL_project/results/networks/bidirectional_dif.png", width = 40, height = 20, units = "cm", res = 250)
boxplot(inter_annotated_BB_long2$dif  ~ inter_annotated_BB_long2$cell_type, ylab = "Difference", xlab = "Cell type", ylim = c(-2,100))
text(x = 1:length(medians), c(-1,-1), labels = round(medians, 2), 
     pos = 1, cex = 1.3, col = "black") 
dev.off()


boxplot(inter_annotated_BB$distance ~ inter_annotated_BB$direction, ylab = "Genomic distance", xlab = "Type of B-B interaction", ylim = c(-1, max(inter_annotated_BB$distance, na.rm = TRUE))) # shows a highly right skewed data



# - Proportions of remaining non significant/ significant and changing between scores
inter_annotated_BB_long2$score_1_c <- ifelse(inter_annotated_BB_long2$score_1 >= 5, 1, 0)
inter_annotated_BB_long2$score_2_c <- ifelse(inter_annotated_BB_long2$score_2 >= 5, 1, 0)

inter_annotated_BB_long2$change <- paste(inter_annotated_BB_long2$score_1_c, inter_annotated_BB_long2$score_2_c, sep = "-")

inter_annotated_BB_long2$change <- factor(inter_annotated_BB_long2$change, levels = c("0-0", "0-1", "1-0", "1-1"),
                                          labels = c("Not significant", "Change", "Change", "Significant"))

paste0(table(inter_annotated_BB_long2$change), " (",
       round(table(inter_annotated_BB_long2$change) %>% prop.table() * 100, 2), "%)" )


# - Distribution of the ratio between scores 
boxplot(inter_annotated_BB_long2$ratio  ~ inter_annotated_BB_long2$cell_type)
tapply(inter_annotated_BB_long2$ratio, inter_annotated_BB_long2$cell_type, summary)


# - Test H_0: Median(ratio) = 1 
inter_annotated_BB_long2 %>% group_by(cell_type) %>% filter(!is.na(ratio)) %>% 
  summarise(
  chicago1 = median(score_1),
  chicago2 = median(score_2),
  ratio = median(ratio),
  wilcox_p_value = wilcox.test(ratio, mu = 1, paired = FALSE, alternative = "two.sided")$p.value,
  wilcox_statistic = wilcox.test(ratio, mu = 1, paired = FALSE, alternative = "two.sided")$statistic,
  .groups = "drop"
)



# 4) Check the distance between promoters

inter_annotated_BB$counts

regions1 <- GRanges(
  seqnames = inter_annotated_BB$seqnames1,
  ranges = IRanges(start = inter_annotated_BB$start1, end = inter_annotated_BB$end1)
)

regions2 <- GRanges(
  seqnames = inter_annotated_BB$seqnames2,
  ranges = IRanges(start = inter_annotated_BB$start2, end = inter_annotated_BB$end2)
)

# Step 2: Calculate distances between corresponding rows
distances <- vector("numeric", nrow(inter_annotated_BB)) # Initialize a vector to store distances

for (i in seq_len(nrow(inter_annotated_BB))) {
  if (seqnames(regions1[i]) %>% as.numeric() != seqnames(regions2[i]) %>% as.numeric()) {
    distances[i] <- NA # Different chromosomes
  } else {
    distances[i] <- distance(regions1[i], regions2[i])
  }
}

# Step 3: Add distances to the original data frame
inter_annotated_BB$distance <- distances


inter_annotated_BB$direction <- factor(inter_annotated_BB$direction, levels = c("0", "1"),
                                       labels = c("Unidirectional", "Bidirectional"))

medians <- tapply(inter_annotated_BB$distance, inter_annotated_BB$direction, function(x) median(x, na.rm = TRUE))

png("/home/acost1/BALL_project/results/networks/distance_bbinteractions.png", width = 20, height = 30, units = "cm", res = 250)
boxplot(inter_annotated_BB$distance ~ inter_annotated_BB$direction, ylab = "Genomic distance", xlab = "Type of B-B interaction", ylim = c(-1, max(inter_annotated_BB$distance, na.rm = TRUE))) # shows a highly right skewed data
text(x = 1:length(medians), y = c(-1,-1), labels = round(medians, 2), 
     pos = 1, cex = 1, col = "black") 
dev.off()



# Add median labels on the boxplot


tapply(inter_annotated_BB$distance, inter_annotated_BB$direction, summary) # median distance: unidirect > bidirectio
wilcox.test(inter_annotated_BB$distance ~ inter_annotated_BB$direction, alternative = "greater") # significant 


# Trying to normalize the data with Box-Cox transformation
par(mfrow = c(1,2))
distance_transformed <- tapply(inter_annotated_BB$distance, inter_annotated_BB$direction, function(data) {
  boxcox_result <- boxcox(data ~ 1)
  
  # Best lambda (transformation parameter)
  best_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
  print(best_lambda)
  
  # Apply the Box-Cox transformation with the best lambda
  transformed_data <- (data^best_lambda - 1) / best_lambda
  
  hist(transformed_data)
  
  return(transformed_data)
  }
  
 
  )
dev.off()
# looks a normal distribution

# Compare skewness value
tapply(inter_annotated_BB$distance, inter_annotated_BB$direction, function(x) skewness(x, na.rm =TRUE))
skewness(distance_transformed$`0`, na.rm =TRUE); skewness(distance_transformed$`1`, na.rm =TRUE)

# boxplot

medians <- sapply(distance_transformed, function(x) median(x, na.rm =TRUE))
png("/home/acost1/BALL_project/results/networks/distance_bbinteractions_boxcox.png", width = 20, height = 30, units = "cm", res = 250)
boxplot(distance_transformed, main = "Genomic distance normalized with Box-Cox transformation")
text(x = 1:length(medians), y = medians, labels = round(medians, 2), 
     pos = 3, cex = 1, col = "black") 
dev.off()

# t.test
t.test(distance_transformed$`0`, distance_transformed$`1`, alternative = "greater") # significant: unidirect > bidirectio




# 5) Biological annotation

# Expression data

bait1 <- dplyr::select(inter_annotated_BB, ID_1, bait_1 ) %>% dplyr::rename(id_bait = ID_1, bait = bait_1)
bait2 <- dplyr::select(inter_annotated_BB, ID_2, bait_2 ) %>% dplyr::rename(id_bait = ID_2, bait = bait_2) 

bait1_id <- unique(bait1$id_bait) 
bait2_id <- unique(bait2$id_bait) 


intersect(bait1_id, bait2_id) %>% length() # both regions: 13.557 Baits
bait1_id[!bait1_id %in% bait2_id] %>% length() # only region 1: 2.957 baits
bait2_id[!bait2_id %in% bait1_id] %>% length() # only region 1: 2.983 baits


length(unique(c(bait1_id, bait2_id))) # 19.497

inter_baits <- rbind( bait1, bait2) %>%
  filter(!duplicated(.)) # only 19.497 regions



# Compute expression for these regions
inter_baits <- inter_baits %>% dplyr::select(id_bait, bait) %>% filter(!duplicated(.))


inter_baits <- apply(data.frame(inter_baits), 1, function(x){
  n <- length(str_split_1(x[["bait"]], ","))
  
  data.frame(id_bait = rep(x[["id_bait"]], n),
             bait = str_split_1(x[["bait"]], ",") 
  ) 
} ) %>% bind_rows()


# See if genes are showed as symbols or ensembl id
inter_baits$gene_type <- ifelse(startsWith(inter_baits$bait, "ENSG") == TRUE, "ensembl", "symbol")

# hgnc_symbol to ensembl_gene_id (those with hgnc_symbol)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_ensembl <- biomaRt::getBM(filters= "hgnc_symbol",
                                attributes = c("hgnc_symbol", "ensembl_gene_id", "start_position", "chromosome_name", "gene_biotype"), 
                                values = unique(inter_baits$bait), 
                                mart = mart) # n = 23.256 --> unique = 20.702

# remove genes that are not annotated in chr 1:23, X, Y (n = 20.757) --> unique = 20.701 (removed n = , unique = 0)
genes_ensembl <- filter(genes_ensembl, chromosome_name %in% names(table(inter_annotated_BB$seqnames1)) ) 

# remove genes that have hgnc_symbol duplicated (n = 20645) --> unique = 20645  (removed n = 56, unique = 56)
hgnc_symbol_dupl <- genes_ensembl$hgnc_symbol[duplicated(genes_ensembl$hgnc_symbol)]
genes_ensembl <- filter(genes_ensembl, !hgnc_symbol %in% hgnc_symbol_dupl )



inter_baits_simple <- merge(inter_baits, genes_ensembl, by.x = "bait", by.y = "hgnc_symbol", all.x = TRUE)

# leave ensembl id the same in the column ensembl_gene_id
inter_baits_simple$ensembl_gene_id[inter_baits_simple$gene_type == "ensembl"] <- inter_baits_simple$bait[inter_baits_simple$gene_type == "ensembl"]


#- Couldn't find 427 gene ensembl id 
# - NOW: Couldn't find 445 gene ensembl id 
rm_bait <- inter_baits_simple$bait[is.na(inter_baits_simple$ensembl_gene_id)]
inter_baits_simple <- inter_baits_simple[!is.na(inter_baits_simple$ensembl_gene_id),] 



# get the expression for these baits
# load("/home/acost1/BALL_project/data/RNA/counts_rlog_normalized.rds")

load("/home/acost1/BALL_project/data/RNA/counts_tmm_normalized_all_log.rds")

# Get the mean expression for combination of cell type and gene

# Pivot long
NormByTMM_long <- NormByTMM_all_log %>% as.data.frame() %>%
  mutate(genes = rownames(NormByTMM_all_log)) %>%
  pivot_longer(!genes, names_to = "cell_type", values_to = "expression")

# Categories of a cell type
NormByTMM_long$cell_type_c <- lapply(str_split(NormByTMM_long$cell_type, "_"),
                                     function(x) x[[1]]) %>% unlist()

# Median expression
NormByTMM_long <- NormByTMM_long %>% group_by(genes, cell_type_c) %>%
  dplyr::summarise(expression = median(expression)) %>%
  filter(!cell_type_c %in% c("CMP", "Mon", "nCD8") ) %>%
  dplyr::rename(B.id1 = genes,
                cell_type = cell_type_c)
NormByTMM_long$cell_type[NormByTMM_long$cell_type == "nB"] = "nB1Mnew"

# NormByTMM_long$cell_type <- factor(NormByTMM_long$cell_type, 
#                                         levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB1Mnew", "GCB", "memB", "PC"),
#                                         labels = c("HSC", "PrePro-B", "Pro-B", "Pre-B", "Trans-B", "nB", "GCB", "Mem-B", "Plasma"))





bait1 <- dplyr::select(inter_annotated_BB, ID_1, bait_1 ) %>% dplyr::rename(id_bait = ID_1, bait = bait_1)
bait2 <- dplyr::select(inter_annotated_BB, ID_2, bait_2 ) %>% dplyr::rename(id_bait = ID_2, bait = bait_2)

bait1_id <- unique(bait1$id_bait) 
bait2_id <- unique(bait2$id_bait) 


intersect(bait1_id, bait2_id) %>% length() # both regions: 14.022 Baits
bait1_id[!bait1_id %in% bait2_id] %>% length() # only region 1: 4.164 baits
bait2_id[!bait2_id %in% bait1_id] %>% length() # only region 2: 554 baits


length(unique(c(bait1_id, bait2_id))) # 18740

inter_baits <- rbind( bait1, bait2) %>%
  filter(!duplicated(.))
unique(inter_baits$id_bait) %>% length() # only 18740 regions



# Compute expression for these regions

inter_baits <- apply(data.frame(inter_baits), 1, function(x){
  n <- length(str_split_1(x[["bait"]], ","))
  
  data.frame(id_bait = rep(x[["id_bait"]], n),
             bait = str_split_1(x[["bait"]], ",")
  ) 
} ) %>% bind_rows()



# 
inter_baits <- merge(inter_baits, inter_baits_simple, by = c("id_bait", "bait"), all.y = TRUE)


# Merge the expression data with the combination of cell type (HSC) and gene (with ensemble id), from BAIT_1
inter_baits <- merge(inter_baits, NormByTMM_long, by.x = c("ensembl_gene_id"), by.y = c("B.id1"), all.x = TRUE)

# NOW n = 5, not found expression
# filter(inter_baits, is.na(expression)) %>% dplyr::select(ensembl_gene_id) %>% filter(!duplicated(.))
rm_expression <- inter_baits[is.na(inter_baits$expression),]
inter_baits_all <- filter(inter_baits, !is.na(expression))

inter_baits <- inter_baits %>% group_by(id_bait, cell_type) %>% dplyr::summarise(expression = median(expression))
inter_baits$id_bait <- as.numeric(inter_baits$id_bait)

inter_baits <- filter(inter_baits, !is.na(expression))


# Add expression
inter_annotated_BB <- merge(inter_annotated_BB, inter_baits, by.x = "ID_1", by.y = "id_bait", all.x = TRUE)
inter_annotated_BB <- merge(inter_annotated_BB, inter_baits, by.x = c("ID_2", "cell_type"), 
                            by.y = c("id_bait", "cell_type"), all.x = TRUE)

# table(inter_annotated_BB$bait_1 == "non-annotated", is.na(inter_annotated_BB$expression.x))
# la majoria de na's en la expression són non-annotated, tant al bait 1 com al bait 2


# Comparison of transcriptional activity between unidirectional and bidirectional interactions
inter_annotated_BB_long <- reshape2::melt(inter_annotated_BB, id.vars = c("pair", "direction", "cell_type"), 
                                measure.vars = c("expression.x", "expression.y"))

inter_annotated_BB_long <- filter(inter_annotated_BB_long, !is.na(value))



inter_annotated_BB_long$cell_type <- factor(inter_annotated_BB_long$cell_type,
                                            levels = c("HSC", "PreProB", "ProB", "PreB", 
                                                       "immtransB", "nB1Mnew", "GCB",
                                                       "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"),
                                            labels = c("HSC", "PreProB", "ProB", "PreB", 
                                                       "immtransB", "nB1Mnew", "GCB",
                                                       "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"))

test = inter_annotated_BB_long %>% group_by(cell_type) %>% wilcox_test(value~direction) %>%
  add_xy_position(x = "cell_type", step.increase = 0.03)

test$p_cat <- ifelse(test$p <= 0.001, "***", ".")

means <- inter_annotated_BB_long %>% group_by(cell_type, direction)%>%summarize(value = mean(value))%>%as.data.frame()
means$value <- round(means$value,2)



png("/home/acost1/BALL_project/results/networks/expression_inidividual_biuni_interactions.png", width = 50, height = 30, units = "cm", res = 250)
ggboxplot(inter_annotated_BB_long, x = "cell_type", y = "value", col = "direction", group = "direction",
          xlab = "Cell type", ylab ="Expression (individual)", add = "mean") +
  stat_pvalue_manual(test, label = "p_cat",  tip.length = 0.01)+
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), name = "B-B interaction")  + 
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2), group = direction), 
               color = "black", size = 4, vjust = -5, 
               position = position_dodge(0.8)) 
dev.off()

# expression median higher in bidirectional interactions



inter_annotated_BB_merged_long <- inter_annotated_BB_long %>% group_by(pair, direction, cell_type, value) %>%
  summarise(value = mean(value))


test = inter_annotated_BB_merged_long %>% group_by(cell_type) %>% wilcox_test(value~direction) %>%
  add_xy_position(x = "cell_type", step.increase = 0.03)
test$p_cat <- ifelse(test$p <= 0.001, "***", ".")

means <- inter_annotated_BB_merged_long %>% group_by(cell_type, direction)%>%summarize(value = mean(value))%>%as.data.frame()
means$value <- round(means$value,2)


png("/home/acost1/BALL_project/results/networks/expression_merged_biuni_interactions.png", width = 50, height = 30, units = "cm", res = 250)
ggboxplot(inter_annotated_BB_merged_long, x = "cell_type", y = "value", col = "direction", group = "direction",
          xlab = "Cell type", ylab ="Expression (merged)", add = "mean") +
  stat_pvalue_manual(test, label = "p_cat",  tip.length = 0.01)+
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), name = "B-B interaction") 
dev.off()





# Correlation analysis

inter_annotated_BB %>% 
  filter(!is.na(cell_type)) %>%
  group_by(direction, cell_type) %>%
  summarise(correlation = cor(expression.x, expression.y, use = "complete.obs"))

inter_annotated_BB$direction <- factor(inter_annotated_BB$direction, levels = c(0, 1),
                                                   labels = c("Unidirectional", "Bidirectional"))

inter_annotated_BB$cell_type <- factor(inter_annotated_BB$cell_type,
                                                   levels = c("HSC", "PreProB", "ProB", "PreB", 
                                                              "immtransB", "nB1Mnew", "GCB",
                                                              "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"),
                                                   labels = c("HSC", "PreProB", "ProB", "PreB", 
                                                              "immtransB", "nB1Mnew", "GCB",
                                                              "memB", "PC", "CMP", "nCD4", "nCD8", "Mon"))


inter_annotated_BB$dif <- abs(inter_annotated_BB$expression.x - inter_annotated_BB$expression.y)
inter_annotated_BB$ratio <- inter_annotated_BB$expression.x / inter_annotated_BB$expression.y


test = inter_annotated_BB %>% filter(!is.na(dif)) %>% group_by(cell_type) %>% wilcox_test(dif~direction) %>%
  add_xy_position(x = "cell_type", step.increase = 0.03)
test$p_cat <- ifelse(test$p <= 0.001, "***", ".")

means <- inter_annotated_BB %>% group_by(cell_type, direction)%>%summarize(dif = mean(dif, na.rm = TRUE))%>%as.data.frame()
means$dif <- round(means$dif,2)

png("/home/acost1/BALL_project/results/networks/expression_dif_biuni_interactions.png", width = 50, height = 30, units = "cm", res = 250)
ggboxplot(inter_annotated_BB %>% filter(!is.na(dif)) %>% group_by(cell_type) , x = "cell_type", y = "dif", col = "direction", group = "direction",
          xlab = "Cell type", ylab ="Difference (individual)", add = "mean") +
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), name = "B-B interaction") 
dev.off()





# Draw flowchart
#################
gp <- gpar(fill = "lightgrey")
width <- .15

# Define the nodes with positions
total_interactions <- boxGrob(
  expression(bold("Total interactions")~"\n"~"N = 563,328"),  # Apply bold to "Total interactions"
  x = 0.5, 
  y = 0.9, 
  just = "center", 
  box_gp = gp,
  width = 0.4
)
b_b_interactions <- boxGrob("B-B interactions\n 64,977 (11.53%)", x = 0.3, y = 0.6, just = "center", box_gp = gp, width = width)
b_oe_interactions <- boxGrob("B-OE interactions \n 498,351 (88.47%)", x = 0.7, y = 0.6, just = "center", box_gp = gp, width = width)
bidirectional <- boxGrob("Bidirectional\n 21,348 (32.85%)", x = 0.2, y = 0.3, just = "center", box_gp = gp, width = width)
unidirectional <- boxGrob("Unidirectional\n 43,629 (67.15%)", x = 0.4, y = 0.3, just = "center",box_gp = gp, width = width)


# Start a new grid page
png("/home/acost1/BALL_project/results/networks/flowchart_interactions.png", width = 40, height = 20, units = "cm", res = 250 )
grid.newpage()

# Draw the boxes
grid.draw(total_interactions)
grid.draw(b_b_interactions)
grid.draw(b_oe_interactions)
grid.draw(bidirectional)
grid.draw(unidirectional)

# Add connections
grid.draw(connectGrob(total_interactions, b_b_interactions, type = "N"))
grid.draw(connectGrob(total_interactions, b_oe_interactions, type = "N"))
grid.draw(connectGrob(b_b_interactions, bidirectional, type = "N"))
grid.draw(connectGrob(b_b_interactions, unidirectional, type = "N"))
dev.off()


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



# Number of total SIGNIFICANT interactions
inter_annotated_gene %>% group_by(cell_type) %>% summarise(n = length(unique(pair)))

# Number of B-B and OE SIGFINICANT interactions
description_interaction <- inter_annotated_gene %>% group_by(cell_type, int) %>% summarise(n_int = length(unique(pair)) ) %>%
                         group_by(cell_type) %>% mutate(n_total = sum(n_int)) %>%
                         mutate(prop_int = round(n_int /n_total* 100, 2) )


# Number of bidirectional and unidirectional B-B SIGNIFICANT interactions
inter_annotated_gene$direction <- 0
inter_annotated_gene$direction[inter_annotated_gene$pair %in% interactions_overlap] <- 1

description_direction <- inter_annotated_gene %>% filter(int == "B_B") %>%
  group_by(cell_type, direction, int) %>% summarise(n_direction = length(unique(pair)) )  %>%
  group_by(cell_type) %>% mutate(n_total_B = sum(n_direction)) %>%
  mutate(prop_direction = round(n_direction /n_total_B* 100, 2) )
  

description_all <- merge(description_interaction, description_direction, by = c("cell_type", "int"), all = TRUE)

description_all$type_interaction <- NA
description_all$type_interaction[description_all$int == "B_B" & description_all$direction == 0] <- "B_B unidirectional"
description_all$type_interaction[description_all$int == "B_B" & description_all$direction == 1] <- "B_B bidirectional"
description_all$type_interaction[description_all$int == "B_OE"] <- "B_OE"

bb_uni <- paste0(description_all$n_direction[description_all$type_interaction ==  "B_B unidirectional"], "-",
      description_all$prop_direction[description_all$type_interaction ==  "B_B unidirectional"], "%")

bb_bi <- paste0(description_all$n_direction[description_all$type_interaction ==  "B_B bidirectional"], "-",
                 description_all$prop_direction[description_all$type_interaction ==  "B_B bidirectional"], "%")

oe <- paste0(description_all$n_int[description_all$type_interaction ==  "B_OE"], "-",
                description_all$prop_int[description_all$type_interaction ==  "B_OE"], "%")

description_all$num_interactions_cat <- NA
description_all$num_interactions_cat[description_all$type_interaction ==  "B_B unidirectional"] <- bb_uni
description_all$num_interactions_cat[description_all$type_interaction ==  "B_B bidirectional"] <- bb_bi
description_all$num_interactions_cat[description_all$type_interaction ==  "B_OE"] <- oe

description_all$num_interactions <- NA
description_all$num_interactions[description_all$type_interaction ==  "B_B unidirectional"] <- description_all$n_direction[description_all$type_interaction ==  "B_B unidirectional"]
description_all$num_interactions[description_all$type_interaction ==  "B_B bidirectional"] <- description_all$n_direction[description_all$type_interaction ==  "B_B bidirectional"]
description_all$num_interactions[description_all$type_interaction ==  "B_OE"] <- description_all$n_int[description_all$type_interaction ==  "B_OE"]



n_inter_celltype <- description_all %>% group_by(cell_type) %>% dplyr::summarise(num_interactions = sum(num_interactions))

png("/home/acost1/BALL_project/results/networks/flowchart_interactions_stage.png", width = 50, height = 30, res = 300, units = "cm")
description_all %>%
  ggplot(aes(x=cell_type, y=num_interactions, fill=type_interaction)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=num_interactions_cat), size = 4, position = position_stack(vjust = 0.5), color="black") +
  scale_fill_brewer(palette="Paired") +
  ylim(0, max(rep(n_inter_celltype$num_interactions, each = 2)) + 10000) + 
  geom_label(aes(y = rep(n_inter_celltype$num_interactions, each = 3), label=rep(n_inter_celltype$num_interactions, each = 3)), vjust=-0.5, color="black", alpha=0) + 
  guides(fill=guide_legend(title="Interaction type")) + 
  ylab("# interactions") + xlab("")
dev.off()


 
# Distribution of the interactions for each promoter (WE NEED TO DIFFER IF THE PROMOTER IS ACTIVE OR NOT)

# - With promoters 
inter_B <- inter_annotated_gene %>% filter(int == "B_B") %>%
  group_by(cell_type, ID_1) %>% dplyr::summarise(n_B = length(unique(ID_2)))


medians <- tapply(inter_B$n_B, inter_B$cell_type, median)
png("/home/acost1/BALL_project/results/networks/bb_interactions_stage.png", width = 50, height = 30, res = 300, units = "cm")
boxplot(inter_B$n_B ~ inter_B$cell_type, xlab = "Cell type", ylab = "Num of Baits")
text(x = 1:length(medians), y = medians, labels = round(medians, 2), 
     pos = 3, cex = 1, col = "black") 
dev.off()




# - With OE
inter_OE <- inter_annotated_gene %>% filter(int == "B_OE") %>%
  group_by(cell_type, ID_1) %>% dplyr::summarise(n_OE = length(unique(ID_2)))


medians <- tapply(inter_OE$n_OE, inter_OE$cell_type, median)
png("/home/acost1/BALL_project/results/networks/boe_interactions_stage.png", width = 50, height = 30, res = 300, units = "cm")
boxplot(inter_OE$n_OE ~ inter_OE$cell_type,  xlab = "Cell type", ylab = "Num of OE")
text(x = 1:length(medians), y = medians, labels = round(medians, 2), 
     pos = 3, cex = 1, col = "black")
dev.off()

tapply(inter_OE$n_OE, inter_OE$cell_type, summary)


# - With OE that contain EA
inter_OE_EA <- inter_annotated_gene %>% filter(int == "B_OE" & EA_c == 1 ) %>%
  group_by(cell_type, ID_1) %>% summarise(n_OE = length(unique(ID_2)))

boxplot(inter_OE_EA$n_OE ~ inter_OE_EA$cell_type)
tapply(inter_OE_EA$n_OE, inter_OE_EA$cell_type, summary)


# - With OE that contain EPo
inter_OE_EPo <- inter_annotated_gene %>% filter(int == "B_OE" & EPo_c == 1 ) %>%
  group_by(cell_type, ID_1) %>% summarise(n_OE = length(unique(ID_2)))

boxplot(inter_OE_EPo$n_OE ~ inter_OE_EPo$cell_type)
tapply(inter_OE_EPo$n_OE, inter_OE_EPo$cell_type, summary)


# - With OE that contain EPr
inter_OE_EPr <- inter_annotated_gene %>% filter(int == "B_OE" & EPr_c == 1 ) %>%
  group_by(cell_type, ID_1) %>% summarise(n_OE = length(unique(ID_2)))

boxplot(inter_OE_EPr$n_OE ~ inter_OE_EPr$cell_type)
tapply(inter_OE_EPr$n_OE, inter_OE_EPr$cell_type, summary)


# - With OE that contain Sil
inter_OE_Sil <- inter_annotated_gene %>% filter(int == "B_OE" & Sil_c == 1 ) %>%
  group_by(cell_type, ID_1) %>% summarise(n_OE = length(unique(ID_2)))

boxplot(inter_OE_Sil$n_OE ~ inter_OE_Sil$cell_type)
tapply(inter_OE_Sil$n_OE, inter_OE_Sil$cell_type, summary)



#####################
#  Search for hubs  #
#####################

# For each Bait (ID_1): 
#######################
# 1) See which other Baits interact with (ID_2)
# 2) See if these other Baits also interact the same B/OE that were interacting with the initial bait
# 3) Compute the proportion of the same interaction regions (OE and B) between 2 baits


compute_jaccard_per_bait_df <- function(df) {
  results <- data.frame(
    Primary_Bait = character(),
    Related_Bait = character(),
    Num_B_Interactions = integer(),
    Num_OE_Interactions = integer(),
    Jaccard_BB = numeric(),
    Jaccard_BOE = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique baits in ID_1
  baits <- unique(df$ID_1)
  
  for (bait in baits) {
    # Step 1: Filter interactions for the bait of interest (S1)
    S1 <- df %>% filter(ID_1 == bait)
    
    # Step 2: Identify related baits (Baits_2)
    Baits_2 <- S1 %>% filter(node.class2 == "B") %>% pull(ID_2) %>% unique()
    
    # Step 3: Compute metrics for each related bait
    for (related_bait in Baits_2) {
      # Interactions for the related bait
      S2 <- df %>% filter(ID_1 == related_bait)
      
      # Shared interactions
      BB_S1 <- S1 %>% filter(node.class2 == "B") %>% pull(ID_2)
      BB_S2 <- S2 %>% filter(node.class2 == "B") %>% pull(ID_2)
      BOE_S1 <- S1 %>% filter(node.class2 == "OE") %>% pull(ID_2)
      BOE_S2 <- S2 %>% filter(node.class2 == "OE") %>% pull(ID_2)
      
      # Compute Jaccard indices
      jaccard_BB <- length(intersect(BB_S1, BB_S2)) / length(union(BB_S1, BB_S2))
      jaccard_BOE <- length(intersect(BOE_S1, BOE_S2)) / length(union(BOE_S1, BOE_S2))
      
      # Count interactions
      num_B_Interactions_B1 <- length(unique(S1 %>% filter(node.class2 == "B") %>% pull(ID_2)))
      num_OE_Interactions_B1 <- length(unique(S1 %>% filter(node.class2 == "OE") %>% pull(ID_2)))
      
      
      num_B_Interactions_B2 <- length(intersect(BB_S1, BB_S2))
      num_OE_Interactions_B2 <- length(intersect(BOE_S1, BOE_S2))
      
      
      # Add row to results
      results <- rbind(
        results,
        data.frame(
          Primary_Bait = bait,
          Related_Bait = related_bait,
          Num_B_Interactions_B1 = num_B_Interactions_B1,
          Num_OE_Interactions_B1 = num_OE_Interactions_B1,
          
          num_B_Interactions_B2,
          num_OE_Interactions_B2,
          
          Jaccard_BB = jaccard_BB,
          Jaccard_BOE = jaccard_BOE,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(results)
}



inter_celltype <- filter(inter_annotated_gene, cell_type == "HSC")

# Get all unique baits in ID_1
baits <- unique(inter_celltype$ID_1)

inter_celltype <- filter(inter_celltype, ID_1 %in% baits) %>% 
  select(ID_1, ID_2, node.class2)

all_hubs <- apply(data.frame(levels(inter_annotated_gene$cell_type) ), 1, 
                  function(x){
                    inter_celltype <- filter(inter_annotated_gene, cell_type == x) %>%
                                        select(ID_1, ID_2, node.class2)
                    return(compute_jaccard_per_bait_df(inter_celltype))
                    
                  } )
  
  
filter(all_hubs[[1]], Num_OE_Interactions_B1 > 10 & Jaccard_BB > 0.7 & Jaccard_BOE > 0.7)

save(all_hubs, file = "/home/acost1/BALL_project/results/networks/all_hubs.rds")


load("/home/acost1/BALL_project/results/networks/all_hubs.rds")


# Distribution of Jaccard B-B
jaccard_B <- data.frame(
  variable = unlist(lapply(all_hubs, function(df) df$Jaccard_BB)),
  cell_type = rep(levels(inter_annotated_gene$cell_type), times = sapply(all_hubs, nrow) )
)

jaccard_B$cell_type <- factor(jaccard_B$cell_type, levels = levels(inter_annotated_gene$cell_type),
                                         labels = levels(inter_annotated_gene$cell_type))

medians <- tapply(jaccard_B$variable * 100, jaccard_B$cell_type, median)

png("/home/acost1/BALL_project/results/networks/jaccard_b_stage.png", width = 50, height = 30, res = 300, units = "cm")
boxplot(variable * 100 ~ cell_type, 
        data = jaccard_B, 
        main = "Shared Baits",
        ylab = "Jaccard index", 
        xlab = "",
        col = viridis(9),  # Different colors for each box
        las = 2, ylim = c(-1, max(jaccard_B$variable * 100))) 
text(x = 1:length(medians), y = c(-1,-1), labels = round(medians, 2), 
     pos = 1, cex = 1, col = "black") 
dev.off()





# Distribution of Jaccard B-OE (only for B of interest-OE interactions > 10  )
jaccard_OE <- data.frame(
  variable = unlist(lapply(all_hubs, function(df) {
    df <- df %>% filter(Num_OE_Interactions_B1 > 4)
    df$Jaccard_BOE 
    } )  ),
  cell_type = rep(levels(inter_annotated_gene$cell_type), times = sapply(all_hubs, function(df) {
    df <- df %>% filter(Num_OE_Interactions_B1 > 4)
    nrow(df)
  } ) )
  
)

jaccard_OE$cell_type <- factor(jaccard_OE$cell_type, levels = levels(inter_annotated_gene$cell_type),
                              labels = levels(inter_annotated_gene$cell_type))


medians <- tapply(jaccard_OE$variable * 100, jaccard_OE$cell_type, median)

png("/home/acost1/BALL_project/results/networks/jaccard_oe_stage.png", width = 50, height = 30, res = 300, units = "cm")
boxplot(variable * 100 ~ cell_type, 
        data = jaccard_OE, 
        main = "Shared OE",
        ylab = "Jaccard index", 
        xlab = "",
        col = viridis(9),  # Different colors for each box
        las = 2, ylim = c(-1, max(jaccard_OE$variable * 100)))
text(x = 1:length(medians), y = c(-1,-1), labels = round(medians, 2), 
     pos = 1, cex = 1, col = "black") 
dev.off()


# Most compacted hubs
compacted_hubs <- lapply(all_hubs, function(x){
  filter(x, Jaccard_BB > 0.5 & Jaccard_BOE > 0.5) %>%
    pull(Primary_Bait) %>% unique() } )

Reduce(intersect, compacted_hubs)

sapply(compacted_hubs, length)
sapply(all_hubs, function(x) length(unique(x$Primary_Bait)))


# Venn diagram between 35108 (TXNIP), 35081 (PIAS3)

s1_b <- filter(inter_annotated_gene, ID_1 == 35108 & int == "B_B" & cell_type == "HSC") %>% pull(bait_2)
s2_b <- filter(inter_annotated_gene, ID_1 == 35081 & int == "B_B" & cell_type == "HSC") %>% pull(bait_2)

s1_oe <- filter(inter_annotated_gene, ID_1 == 35108 & int == "B_OE" & cell_type == "HSC") %>% pull(ID_2)
s2_oe <- filter(inter_annotated_gene, ID_1 == 35081 & int == "B_OE" & cell_type == "HSC") %>% pull(ID_2)


ggVennDiagram(list(TXNIP = s1_b, PIAS3 = s2_b)) + ggtitle("Shared Baits")
ggVennDiagram(list(TXNIP = s1_oe, PIAS3 = s2_oe)) + ggtitle("Shared OE")




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


# ggraph function for a gene
gene_community_ggraph_option1 <- function(gene_id, edge_file, id_bed_file, cell_type_subset){

  # Edge file and id_bed file only for a specific cell_type_subset
  edge_file <- dplyr::filter(edge_file_split, cell_type == cell_type_subset )%>%
    select(ID_1, ID_2, int)
  
  
  # identify the lichi-c region (by using the id_chromhmm)
  load("/home/acost1/BALL_project/results/gnn/chromhmm_interaction_all_updated.rds")

  aux <- merge(id_bed_split, chromhmm_interaction_all %>% select(id_oe, cell_type, id) %>% 
                                                        dplyr::rename(id_chromhmm = id), 
             by.x = c("id","cell_type" ),
             by.y = c("id_oe", "cell_type"),
             all.x = TRUE)
  aux$id_chromhmm[is.na(aux$id_chromhmm)] <- aux$id[is.na(aux$id_chromhmm)]
  id_bed_file <- aux
  print(head(aux))

  id_bed_file <- filter(id_bed_file, cell_type == cell_type_subset) %>% select(id, id_chromhmm, name_id, caract, expression)
  id_bed_file$id <- as.numeric(id_bed_file$id)
  id_bed_file$id_chromhmm <- as.numeric(id_bed_file$id_chromhmm)
  print("Step 1")

  # Select only interactions from gene_id of interest
  edges_subset <- edge_file[edge_file$ID_1 %in% gene_id,]

  # Select all baits and other ends from the previous interactions
  id_subset <- id_bed_file[id_bed_file$id %in% unique(c(edges_subset$ID_1, edges_subset$ID_2) ),]
  print("Step 2")

  # Check if other baits also interact with our gene_id of interest
  id_otherB <- edges_subset$ID_2[edges_subset$int == "B_B"] %>% unique()

  aux <- edge_file[edge_file$ID_1 %in% id_otherB,]
  print(nrow(aux))
  print("Step 3")

  # Add (if any) the other bidirectional interactions
  if(nrow(aux[aux$ID_2 == gene_id,]) != 0){
    edges_subset <- rbind(edges_subset, aux[aux$ID_2 == gene_id,] ) }

  graph <- graph_from_data_frame(d = edges_subset, vertices = id_subset, directed = TRUE)
  print("Step 4")


  #create network with a group structure
  graph <- igraph::simplify(graph)
  # V(graph)$id_chromhmm --> group

  # create a categorical variable: seria com les nostres tipus de reguladors distals + promoters
  # V(graph)$caract

  # contract the network based on the groups
  g_clu <- contract(graph, as.numeric(factor(V(graph)$id_chromhmm)), vertex.attr.comb = "concat")
  E(g_clu)$weight <- 1
  g_clu <- simplify(g_clu, edge.attr.comb = "sum")

  # calculate the number of occurrences of each attribute in each cluster
  V(g_clu)$B <- sapply(V(g_clu)$caract, function(x) sum(x == "B"))
  V(g_clu)$EA <- sapply(V(g_clu)$caract, function(x) sum(x == "EA"))
  V(g_clu)$EPr <- sapply(V(g_clu)$caract, function(x) sum(x == "EPr"))
  V(g_clu)$EPo <- sapply(V(g_clu)$caract, function(x) sum(x == "EPo"))
  V(g_clu)$Sil <- sapply(V(g_clu)$caract, function(x) sum(x == "Sil"))

  # precompute layout
  vertex_name <- V(g_clu)$name %>% unlist() 
  xy <- layout_with_focus(g_clu, v = which(vertex_name == gene_id))
  V(g_clu)$x <- xy$xy[, 1]
  V(g_clu)$y <- xy$xy[, 2]

  elements_colors <- data.frame(
    caract = c("B", "EA", "EPr", "EPo", "Sil"),
    colors = c("#ffffff", "#288A35", "#59D96A", "#EDE658", "#3659B3") )

  elements_colors_subset <- filter(elements_colors, 
                                   caract %in% (V(g_clu)$caract %>% unlist %>% unique() ) )
  
  name_id2 <- V(g_clu)$name_id %>% lapply(., function(x) paste(x, collapse = "-")) %>% unlist()
  caract <- V(g_clu)$caract %>% lapply(., function(x) paste(x, collapse = "-")) %>% unlist()
  
  name_id2[caract != "B"] <- ""

  graph_gene <- ggraph(g_clu, "manual", x = V(g_clu)$x, y = V(g_clu)$y) +
    geom_edge_link() +
    geom_scatterpie(
      cols = elements_colors_subset$caract,
      data = as_data_frame(g_clu, "vertices"),
      colour = "white",
      pie_scale = 3
    ) +
    scale_fill_manual(
      values = elements_colors_subset$colors, 
      name = "Regulatory elements" # Optional legend title for pies
    ) +
    # Highlight "B" vertices with a gradient based on activity
    geom_node_point(data = as_data_frame(g_clu, "vertices") %>% 
                      filter(caract == "B") %>% select(x, y, expression) %>%
                      mutate(
                        x = unlist(x),
                        y = unlist(y),
                        expression = unlist(expression ) ), # Filter for "B" vertices
                    aes(x = x, y = y, alpha = expression, color = expression), # Map activity to fill
                     size = 5, show.legend = TRUE, col = "red") +
    scale_alpha_continuous(limits = c(0, 9.6), oob = scales::squish) +
    scale_color_viridis_c(option = "viridis", limits = c(0, 9.6)) +
    # Vertex labels
    geom_node_text(aes(label = name_id2), repel = TRUE, size = 3) +
    coord_fixed() +
    theme_graph() +
    theme(legend.position = "none") + 
    ggtitle(cell_type_subset)
  
  return(graph_gene)

}

geom_edge_link0() +
  # Highlight "B" vertices with a gradient based on activity
  geom_node_point(data = as_data_frame(g_clu, "vertices") %>% 
                    filter(caract == "B") %>% select(x, y, expression) %>%
                    mutate(
                      x = unlist(x),
                      y = unlist(y),
                      expression = unlist(expression) / 15), # Filter for "B" vertices
                  aes(x = x, y = y, fill = expression), # Map activity to fill
                  size = 5, shape = 21) +
  scale_fill_gradient(
    low = "lightblue", high = "darkblue", 
    na.value = "gray",
    name = "Gene Expression" # Optional legend title for activity
  ) +

ggraph(g_clu, "manual", x = V(g_clu)$x, y = V(g_clu)$y) +
  geom_edge_link0() +
  geom_scatterpie(
    cols = elements_colors_subset$caract,
    data = as_data_frame(g_clu, "vertices"),
    colour = "grey",
    pie_scale = 3)


ggraph(g, "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link0() +
  geom_scatterpie(
    cols = c("A", "B", "C"),
    data = as_data_frame(g, "vertices"),
    colour = NA,
    pie_scale = 2
  ) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = "bottom")

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



graph_gene <- apply(data.frame(all_celltypes), 1, function(x){
  gene_community_ggraph_option1(id_BCL6[1], edge_file_split, id_bed_split, x)
})

png("/home/acost1/BALL_project/results/networks/BCL6_option1.png", width = 55, height = 25, units = "cm", res = 250)
wrap_plots(graph_gene, ncol = 5)
dev.off()

dummy_data <- data.frame(
  x = 1:4,
  y = 1:4,
  group = factor(c( "EA", "EPr", "EPo", "Sil"), levels = c( "EA", "EPr", "EPo", "Sil"),
                 labels = c( "EA", "EPr", "EPo", "Sil"))
)

# Create a plot
plot_with_legend <- ggplot(dummy_data, aes(x = x, y = y, color = group)) +
  geom_point(size = 4) +
  theme_minimal() + # Simple theme
  scale_color_manual(values = c("#288A35", "#59D96A", "#EDE658", "#3659B3")) + 
  labs(color = "Regulatory elements") # Add a legend title

# Extract only the legend
legend_only <- get_legend(plot_with_legend)

# Display the legend
grid::grid.newpage()
grid::grid.draw(legend_only)

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






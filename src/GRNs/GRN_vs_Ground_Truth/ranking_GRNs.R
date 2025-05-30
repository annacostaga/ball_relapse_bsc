##############################################
####### COMBINING and Ranking GRNs
##############################################

###LOAD LIBRARIES
library(dplyr)
library(tidyr)
library(decoupleR)


### LOAD counts
counts <- read.csv("ball_relapse_bsc/results/rna_seq_primary/rlog_filtered_Bcounts.csv")

### LOAD GRNs
grn1 <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/GRNs/GRN_Panda.csv", row.names = 1)
grn2 <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/GRNs/GRN_GRNBoost.csv")
grn3 <- read.csv("/home/spere5/Documents/BALL_project/BALL_git_results/GRNs/GRaNIE_GRN_weights.csv")
grn4 <- read.csv("/gpfs/scratch/bsc08/bsc08981/minet/05_quantile_filtered_GRN.csv")

### Rank within each GRN, if values are tied replace by their highest rank
grn1$rank <- rank(grn1$importance, ties.method = "max")
grn2$rank <- rank(grn2$importance, ties.method = "max")
grn3$rank <- rank(grn3$importance, ties.method = "max")
grn4$rank <- rank(grn4$importance, ties.method = "max")

### Load literature: Col and Xtri
col <- get_collectri()
col <- col %>%
  select(source,target) %>%
  rename(TF = source)

extri <- read.csv("/home/spere5/ranking_GRNs/Extri.csv")
extri_df <- extri %>%
  rename(TF = Transcription.Factor..Associated.Gene.Name., target = Target.Gene..Associated.Gene.Name.) %>%
  select(TF, target ,X.ExTRI..Confidence) %>%
  filter(X.ExTRI..Confidence == "High")
  
### Merge rankings based on TF-gene
ranked_final<- full_join(grn1 %>% select(c(TF,target), rank1 = rank),
                       grn2 %>% select(c(TF,target), rank2 = rank),
                       by = c("TF", "target")) %>%
  full_join(grn3 %>% select(c(TF,target), rank3 = rank), by = c("TF", "target"))

# Merge rank4 into ranked_final
ranked_final<- ranked_final%>%
  full_join(grn4 %>% select(c(TF,target), rank4 = rank), by = c("TF", "target"))



#  Assign TRUE if TF-gene is in xtri, FALSE otherwise
ranked_final<- ranked_final%>%
  mutate(ctri = ifelse(paste(TF, target) %in% paste(col$TF, col$target), T, F))


### Assign TRUE if TF-gene is in xtri, FALSE otherwise
ranked_final <- ranked_final %>%
  mutate(extri = ifelse(paste(TF, target) %in% paste(extri_df$TF, extri_df$target), T, F))

### remove duplicates based on all columns and gene-gene interactions
ranked_final<- ranked_final%>% distinct()

TF_list <- read.csv("/gpfs/scratch/bsc08/bsc08981/pyscenic/allTFs_hg38.txt")

ranked_final <- ranked_final %>%
  filter(TF %in% TF_list$TF)


# Calculate pearson correlation for TF and targets in B-cell expression data
calculate_correlation <- function(TF, gene, counts) {
  if (TF %in% rownames(counts) & gene %in% rownames(counts)) {
    return(cor(as.numeric(counts[TF, ]), as.numeric(counts[gene, ]), method = "pearson"))
  } else {
    return(NA)  # Return NA if TF or gene is not found in the counts matrix
  }
}

# Compute correlation for each TF-gene pair
ranked_final<- ranked_final %>%
  rowwise() %>%  
  mutate(corr = calculate_correlation(TF, gene, counts)) 

### Filter ranked_final as it is very big (especially PANDA and minet)

# Define the 5% quantile threshold for rank1
quantile_threshold_1 <- quantile(ranked_final$rank1, 0.05, na.rm = T)
quantile_threshold_2 <- quantile(ranked_final$rank4, 0.05, na.rm = T)


# Filter the data
filtered_ranked <- ranked_final %>%
  filter(
    (rank1 >= quantile_threshold_1 | 
       (!is.na(rank2) | !is.na(rank3))) &
      (rank4 >= quantile_threshold_2 |
         (!is.na(rank2) | !is.na(rank3))) 
  )

rm(quantile_threshold_1,quantile_threshold_2 )

#Save results
save.image("GRNs_ranked.RData")

###########################################################################
## DEFINE THRESHOLDS FOR FINAL GRN
######################################################################
                                        
### Assign which TF-gene to keep based on set of rules (THRESHOLDS 1)
assign_keep <- function(df) {
  # Calculate rank thresholds for top 1%
  rank1_threshold <- quantile(df$rank1, 0.99, na.rm = T)
  rank2_threshold <- quantile(df$rank2, 0.99, na.rm = T)
  rank3_threshold <- quantile(df$rank3, 0.99, na.rm = T)
  rank4_threshold <- quantile(df$rank4, 0.99, na.rm = T)
  df %>%
    mutate(keep = case_when(
      corr > 0.95 ~ T,  # High correlation alone is enough
      ctri == T | extri == T & corr > 0.9 ~ T,  # TF-gene in Collectri and strong correlation
      (rank1 <= rank1_threshold | rank2 <= rank2_threshold | rank3 <= rank3_threshold | rank4<= rank4_threshold) & corr > 0.9 ~ T,  # Top 1% in at least one rank & moderate correlation
      T ~ F  # Default case (if none of the conditions are met)
    ))
}

filtered_ranked <- assign_keep(filtered_ranked)

### Extract the final GRN from filtered_ranked
GRN <- filtered_ranked %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
  mutate(
    weight = 1  # Assign importance = 1 for all rows
  ) %>%
  select(TF, target, weight)  # Keep only relevant columns

write.csv(GRN, "GRN_first_THRESHOLDS_1.csv", row.names = F, quote =F)

### Assign which TF-gene to keep based on set of rules (THRESHOLDS 2)
assign_keep <- function(df) {
  # Calculate rank thresholds for top 1%
  rank1_threshold <- quantile(df$rank1, 0.99, na.rm = T)
  rank2_threshold <- quantile(df$rank2, 0.99, na.rm = T)
  rank3_threshold <- quantile(df$rank3, 0.99, na.rm = T)
  rank4_threshold <- quantile(df$rank4, 0.99, na.rm = T)
  df %>%
    mutate(keep = case_when(
      corr > 0.5 ~ T,  # High correlation alone is enough
      ctri == T | extri == T & corr > -0.5 ~ T,  # TF-gene in Collectri and strong correlation
      (rank1 <= rank1_threshold | rank2 <= rank2_threshold | rank3 <= rank3_threshold | rank4<= rank4_threshold) & corr > 0.6 ~ T,  # Top 5% in at least one rank & moderate correlation
      T ~ F  # Default case (if none of the conditions are met)
    ))
}

ranked_GRNs <- assign_keep(ranked_GRNs)

### Extract the final GRN from filtered_ranked
GRN <- ranked_GRNs %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
  mutate(
    weight= 1  # Assign importance = 1 for all rows
  ) %>%
  select(TF, target, weight)  # Keep only relevant columns

write.csv(GRN, "GRN_first_THRESHOLDS_2.csv", row.names = F, quote =F)

### Assign which TF-gene to keep based on set of rules (THRESHOLDS 3)
assign_keep <- function(df) {
  # Calculate rank thresholds for top 1%
  rank1_threshold <- quantile(df$rank1, 0.99, na.rm = T)
  rank2_threshold <- quantile(df$rank2, 0.99, na.rm = T)
  rank3_threshold <- quantile(df$rank3, 0.99, na.rm = T)
  rank4_threshold <- quantile(df$rank4, 0.99, na.rm = T)
  df %>%
    mutate(
      # Check if TF-gene is in top 1% in each rank
      top_rank1 = rank1 <= rank1_threshold,
      top_rank2 = rank2 <= rank2_threshold,
      top_rank3 = rank3 <= rank3_threshold,
      top_rank4 = rank4 <= rank4_threshold,
      
      # Count how many ranks are in the top 1%
      top_count = rowSums(across(c(top_rank1, top_rank2, top_rank3,top_rank4))),
      
      # Assign T if any of the criteria are met
      keep = case_when(
      corr > 0.99 ~ T,  # High correlation alone is enough
      ctri == T | extri == T & corr > 0.95 ~ T,  # TF-gene in Collectri and strong correlation
      top_count >= 2 & corr > 0.95 ~ T,
      T ~ F  # Default case (if none of the conditions are met)
    )) %>%
    select(-top_rank1, -top_rank2, -top_rank3, -top_rank4, -top_count)  # Remove helper columns
}

filtered_ranked <- assign_keep(filtered_ranked)

### Extract the final GRN from filtered_ranked
GRN <- filtered_ranked %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
    mutate(weight = 1)%>% # Assign importance = 1 for all rows 
    select(TF, target, weight)  # Keep only relevant columns

write.csv(GRN, "GRN_first_THRESHOLDS_3.csv", row.names = F, quote =F)


### Assign which TF-gene to keep based on set of rules (THRESHOLDS 4)
assign_keep <- function(df) {
  # Calculate rank thresholds for top 1%
  rank1_threshold <- quantile(df$rank1, 0.99, na.rm = T)
  rank2_threshold <- quantile(df$rank2, 0.99, na.rm = T)
  rank3_threshold <- quantile(df$rank3, 0.99, na.rm = T)
  rank4_threshold <- quantile(df$rank4, 0.99, na.rm = T)
  df %>%
    mutate(
      # Check if TF-gene is in top 1% in each rank
      top_rank1 = rank1 <= rank1_threshold,
      top_rank2 = rank2 <= rank2_threshold,
      top_rank3 = rank3 <= rank3_threshold,
      top_rank4 = rank4 <= rank4_threshold,
      
      # Count how many ranks are in the top 1%
      top_count = rowSums(across(c(top_rank1, top_rank2, top_rank3,top_rank4))),
      
      # Assign T if any of the criteria are met
      keep = case_when(
        corr > 0.8 | corr < -0.8 ~ T,  # High correlation alone is enough
        (ctri == T | extri == T) & (corr > 0.5 | corr < -0.5)  ~ T,  # TF-gene in Collectri and strong correlation
        top_count >= 2 & (corr > 0.5 | corr < -0.5) ~ T,
        T ~ F  # Default case (if none of the conditions are met)
      )) %>%
    select(-top_rank1, -top_rank2, -top_rank3, -top_rank4, -top_count)  # Remove helper columns
}

filtered_ranked <- assign_keep(filtered_ranked)

### Extract the final GRN from filtered_ranked
GRN <- filtered_ranked %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
  mutate(weight = ifelse(corr > 0, 1, -1))  %>% #assing signed interactions
  select(TF, target, weight)  # Keep only relevant columns

write.csv(GRN, "GRN_first_THRESHOLDS_4.csv", row.names = F, quote =F)

### Assign which TF-gene to keep based on set of rules (THRESHOLDS GA)

#ranked GRNs
ranked_GRNs <- read.csv("ranked_GRNs.csv", sep = ",")
assign_keep <- function(df) {
  # Calculate rank thresholds for top 1%
  rank1_threshold <- quantile(df$rank1,  0.8109199, na.rm = T)
  rank2_threshold <- quantile(df$rank2,  0.3978305, na.rm = T)
  rank3_threshold <- quantile(df$rank3, 0.3966308, na.rm = T)
  rank4_threshold <- quantile(df$rank4,  0.5512853, na.rm = T)
  df %>%
    mutate(
          keep = case_when(
            (corr > 0.8752730 | corr < -0.6649557) ~ T,  # High correlation alone is enough
          (ctri == T | extri == T) & (corr > 0.7032986 | corr < -0.5437288) ~ T,  # TF-gene in Collectri and strong correlation
          (rank1 <= rank1_threshold | rank2 <= rank2_threshold | rank3 <= rank3_threshold | rank4<= rank4_threshold) 
          & (corr >  0.5167392 | corr < -0.3082980) ~ T,  # Top 1% in at least one rank & moderate correlation
          T ~ F  # Default case (if none of the conditions are met)
        ))
}

ranked_GRNs <- assign_keep(ranked_GRNs)

### Extract the final GRN from filtered_ranked
GRN <- ranked_GRNs  %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
   mutate(weight = ifelse(corr > 0, 1, -1))   %>%#assing signed interactions
  select(TF, target, weight)  # Keep only relevant columns

write.csv(GRN, "GRN_GA_THRESHOLDS_comb_30_may_20205.csv", row.names = F, quote =F)
###############################################
### DEFINED NEW THRESHOLDS
##################################################
# Calculate rank thresholds 
rank1_threshold <- quantile(ranked_GRNs$rank1,0.3936774, na.rm = T)
rank2_threshold <- quantile(ranked_GRNs$rank2,0.6806102, na.rm = T)
rank4_threshold <- quantile(ranked_GRNs$rank4,0.02497388, na.rm = T)

ranked_filtered <- ranked_GRNs  %>%
  mutate(keep = case_when(
    (
      (rank1 <= rank1_threshold | rank2 <= rank2_threshold | rank4 <= rank4_threshold) &
        (
          ((ctri == T | extri == T) & (abs(corr) > abs(0.8411693))) |
            ((ctri == F & extri == F) & (abs(corr) > abs(0.9227458))))) ~ T,
    T ~ F
  ))



GRN <- ranked_filtered  %>%
  filter(keep == T) %>%  # Keep only interactions marked as T
  mutate(weight = ifelse(corr > 0, 1, -1))  %>%
  select(TF, target, weight)  # Keep only relevant columns

#ensure unique interactions
GRN <- unique(GRN)

write.csv(GRN, "GRN_GA_THRESHOLDS_new_defined_thresh.csv", row.names = F, quote =F)

#Save results
save.image("GRNs_ranked.RData")

library(netZooR)


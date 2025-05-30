##############################################
####### COMBINING and Ranking GRNs
##############################################

###LOAD LIBRARIES
library(dplyr)
library(tidyr)


## Setwd
Setwd("/gpfs/projects/bsc08/shared_projects/BALL_RELAPSE/annotations/GRNs/ranking_GRNs")


### LOAD counts
counts <- read.csv("../ball_relapse_bsc/results/rna_seq_primary/rlog_filtered_Bcounts.csv")


### LOAD GRNs
grn1 <- read.csv("GRN_Panda.csv", row.names = 1)
grn2 <- read.csv("GRN_GRNBoost.csv")
grn3 <- read.csv("GRaNIE_GRN_weights.csv")
grn4 <- read.csv("05_quantile_filtered_GRN.csv")

### Rank within each GRN, if values are tied replace by their max
grn1$rank <- rank(grn1$importance, ties.method = "max")
grn2$rank <- rank(grn2$importance, ties.method = "max")
grn3$rank <- rank(grn3$importance, ties.method = "max")
grn4$rank <- rank(grn4$importance, ties.method = "max")

### Merge rankings based on TF-gene
ranked_GRNs <- full_join(grn1 %>% select(c(TF,target), rank1 = rank),
                       grn2 %>% select(c(TF,target), rank2 = rank),
                       by = c("TF", "target")) %>%
  full_join(grn3 %>% select(c(TF,target), rank3 = rank), by = c("TF", "target"))

# Merge rank4 into ranked_GRNs
ranked_GRNs <- ranked_GRNs%>%
  full_join(grn4 %>% select(c(TF,target), rank4 = rank), by = c("TF", "target"))

### remove duplicates based on all columns and gene-gene interactions
ranked_GRNs <- ranked_GRNs%>% distinct()

TF_list <- read.csv("../ball_relapse_bsc/annotation/GRNs/ranking_GRNs/allTFs_hg38.txt")

ranked_GRNs <- ranked_GRNs %>%
  filter(TF %in% TF_list$TF)

### Filter ranked_GRNs as it is very big (especially PANDA and minet)

# Define the 5% quantile threshold 
quantile_threshold_1 <- quantile(ranked_GRNs$rank1, 0.95, na.rm = T)
quantile_threshold_2 <- quantile(ranked_GRNs$rank4, 0.95, na.rm = T)


# Filter the data
filtered_ranked <- ranked_GRNs %>%
  filter(
    (rank1 <= quantile_threshold_1 | 
       (!is.na(rank2) | !is.na(rank3))) &
      (rank4 <= quantile_threshold_2 |
         (!is.na(rank2) | !is.na(rank3))) 
  )

rm(quantile_threshold_1,quantile_threshold_2 )


### Load literature: Col and Xtri
col <- get_collectri()
col <- col %>%
  select(source,target) %>%
  rename(TF = source)

col <- paste(col$TF, col$target)
TF_gene <- paste(ranked_GRNs$TF,ranked_GRNs$target)
  
  # Assign TRUE if TF-gene pair is in col, FALSE otherwise
  ranked_GRNs <- ranked_GRNs %>%
    mutate(col = ifelse(TF_gene %in% col, T, F))

extri <- read.csv("../ball_relapse_bsc/annotation/GRNs/ranking_GRNs/Extri.csv")
extri <- extri %>%
   separate(`TF:TG`, into = c("TF", "target"), sep = ":", remove = T) %>%
    select(TF, target, confidence = `[ExTRI] Confidence`) %>%
    filter(confidence == "High")

extri <- paste(extri$TF, extri$target)

  # Assign TRUE if TF-gene pair is in Extri, FALSE otherwise
  ranked_GRNs <- ranked_GRNs %>%
    mutate(extri = ifelse(TF_gene %in% extri, T, F))

   # Calculate pearson correlation for TF-gene
  calculate_correlation <- function(TF, target, counts) {
    if (TF %in% rownames(counts) & target %in% rownames(counts)) {
      return(cor(as.numeric(counts[TF, ]), as.numeric(counts[target, ]), method = "pearson"))
    } else {
      return(NA)  # Return NA if TF or gene is not found in the counts matrix
    }
  }

# Compute correlation for each TF-gene pair
  ranked_GRNs <- ranked_GRNs %>%
    rowwise() %>%  
    mutate(corr = calculate_correlation(TF,target, counts))

  

#Save results
save.image("GRNs_ranked.RData")


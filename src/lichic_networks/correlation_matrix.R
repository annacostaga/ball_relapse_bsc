
# Packges
library(HiCaptuRe)
library(readr)
library(dplyr)
library(stringr)
library(biomaRt)
library(reshape2)

# Matrix de correlación de los datos de expressión


# - Read counts data
counts_data <- read.table("/home/acost1/BALL_project/data/RNA/counts_star_B_cell_roadmap.tsv", header = TRUE)
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


# - Only cell types from the roadmap (not controls)
counts_data <- counts_data[,-c(1:10)]

# - Matrix de correlación (All genes)
matrix_correlation <- cor(t(counts_data))

save(matrix_correlation, file = "matrix_correlation.rds")



# - Matrix de correlación (Only high expressed)
counts_data_high <- counts_data[ rowSums(counts_data) > 1, ]

matrix_correlation_high <- cor(t(counts_data_high))

save(matrix_correlation_high, file = "matrix_correlation_high.rds")




# Extract genes found in the liCHi-C dataset

# - Load liCHi-C data
file <- "/home/acost1/BALL_project/data/lichic/lichic_peakmatrix_B_cell_roadmap.txt"
inter <- load_interactions(file)
inter_annotated <- annotate_interactions(inter, annotation = "/home/acost1/BALL_project/data/lichic/annotation_example.txt")
inter_annotated <- as.data.frame(inter_annotated)

inter_annotated_BB <- filter(inter_annotated, int == "B_B")
rm(inter_annotated)

# - Extract baits
bait1 <- dplyr::select(inter_annotated_BB, ID_1, bait_1 ) %>% dplyr::rename(id_bait = ID_1, bait = bait_1)
bait2 <- dplyr::select(inter_annotated_BB, ID_2, bait_2 ) %>% dplyr::rename(id_bait = ID_2, bait = bait_2)

bait1_id <- unique(bait1$id_bait) 
bait2_id <- unique(bait2$id_bait) 

# - Number of baits (not genes)
intersect(bait1_id, bait2_id) %>% length() # both regions: 13.557 Baits
bait1_id[!bait1_id %in% bait2_id] %>% length() # only region 1: 2.957 baits
bait2_id[!bait2_id %in% bait1_id] %>% length() # only region 2: 2.983 baits

# - Dataset with the baits
inter_baits <- rbind( bait1, bait2) %>%
  filter(!duplicated(.))
unique(inter_baits$id_bait) %>% length() # only 19.497 regions


# - Dataset with the genes
inter_baits <- apply(data.frame(inter_baits), 1, function(x){
  n <- length(str_split_1(x[["bait"]], ","))
  
  data.frame(id_bait = rep(x[["id_bait"]], n),
             bait = str_split_1(x[["bait"]], ",") 
  ) 
} ) %>% bind_rows()


# - Number of genes
length(unique(inter_baits$bait)) # - 25.817


# - Some of the genes are with Gene symbol, and some others are with ensembl id --> we want all with ensembl id
inter_baits$gene_type <- ifelse(startsWith(inter_baits$bait, "ENSG") == TRUE, "ensembl", "symbol")

# ensembl  symbol 
#   4748   21069 


  # - Gene symbol to ensembl id (n = 21.069)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_ensembl <- biomaRt::getBM(filters= "hgnc_symbol",
                                attributes = c("hgnc_symbol", "ensembl_gene_id", "start_position", "chromosome_name", "gene_biotype"), 
                                values = inter_baits %>% filter(gene_type == "symbol") %>% pull(bait) %>% unique(), 
                                mart = mart) # n genes symbol = 20.702

# remove genes that are not annotated in chr 1:23, X, Y --> n genes = 20.701
genes_ensembl <- filter(genes_ensembl, chromosome_name %in% names(table(inter_annotated_BB$seqnames1)) ) 



# remove genes that have hgnc_symbol duplicated -> n = 56
hgnc_symbol_dupl <- genes_ensembl$hgnc_symbol[duplicated(genes_ensembl$hgnc_symbol)]
genes_ensembl <- filter(genes_ensembl, !hgnc_symbol %in% hgnc_symbol_dupl )


# merge gene ensembl to inter_baits data.frame
inter_baits_simple <- merge(inter_baits, genes_ensembl, by.x = "bait", by.y = "hgnc_symbol", all.x = TRUE)


inter_baits_simple$ensembl_gene_id[inter_baits_simple$gene_type == "ensembl"] <- inter_baits_simple$bait[inter_baits_simple$gene_type == "ensembl"]

# Mirar los duplicados
inter_baits_simple %>% dplyr::select(bait, gene_type, ensembl_gene_id) %>% filter(!duplicated(.)) %>% 
  filter( !is.na(ensembl_gene_id)) %>% filter(duplicated(ensembl_gene_id))



genes_lichic <- unique(inter_baits_simple$ensembl_gene_id)

save(inter_baits_simple, file = "/home/acost1/BALL_project/results/networks/inter_baits_simple.rds")
save(inter_annotated_BB, file = "/home/acost1/BALL_project/results/networks/inter_annotated_BB.rds")
save(genes_lichic, file = "/home/acost1/BALL_project/results/networks/genes_lichic.rds")


# Matrix of correlation only for genes found in the lichic dataset
counts_subset <- counts_data[rownames(counts_data) %in% genes_lichic,]
counts_subset <- counts_subset[rowSums(counts_subset) > 2 , ]
matrix_correlation_long <- melt(cor(t(counts_subset[1:50,])), 
                                varnames = c("ensembl_id1", "ensembl_id2"), value.name = "cor")

s.data.frame(cor(t(counts_subset[1:50,])), row.names = "ensembl_id1") %>%
  rownames_to_column(var = "ensembl_id1") %>%
  pivot_longer(cols = -ensembl_id1, names_to = "ensembl_id2", values_to = "cor")



# remove interactions between same genes
matrix_correlation_long <- matrix_correlation_long[!matrix_correlation_long$ensembl_id1 == matrix_correlation_long$ensembl_id2,]


# put id on baits
matrix_correlation_long <- merge(matrix_correlation_long, inter_baits_simple %>% filter(!is.na(ensembl_gene_id)) %>% dplyr::select(ensembl_gene_id, id_bait) %>% dplyr::rename(ID_1 = id_bait), 
      by.x = "ensembl_id1", by.y = "ensembl_gene_id", all.x = TRUE)

matrix_correlation_long <- merge(matrix_correlation_long, inter_baits_simple %>% filter(!is.na(ensembl_gene_id)) %>% dplyr::select(ensembl_gene_id, id_bait) %>% dplyr::rename(ID_2 = id_bait), 
                                 by.x = "ensembl_id2", by.y = "ensembl_gene_id", all.x = TRUE)

# looks that one gene can be in more than 1 region


# find if there are interactions between genes
matrix_correlation_long <- matrix_correlation_long %>%
  mutate(Significant = ifelse(paste(ID_1, ID_2) %in% paste(inter_annotated_BB$ID_1, inter_annotated_BB$ID_2), 1, 0))
matrix_correlation_long$Significant <- factor(matrix_correlation_long$Significant, levels = c(0,1),
                                              labels = c("Not found in liCHi-c", "Found in liCHi-c"))


png("correlation_distribution.png", width = 50, height = 30, units = "cm", res = 250)
ggplot(matrix_correlation_long, aes(x = cor, fill = Significant)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Correlations distribution between gene interactions",
       x = "Value",
       y = "Frequency") +
  theme_minimal()
dev.off()



library(HiCaptuRe)
library(dplyr)


inter <- load_interactions("/home/acost1/BALL_project/data/lichic/lichic_peakmatrix_B_cell_roadmap.txt")

# Load annotating
inter_annotated <- annotate_interactions(interactions = inter,
                                         annotation = "/home/acost1/BALL_project/data/lichic/annotation_example.txt")
inter_annotated

inter_annotated <- as.data.frame(inter_annotated)

inter_BB <- filter(inter_annotated, int == "B_B")

rm(inter); rm(inter_annotated)

# 1: Normalize each interaction by sorting ID1 and ID2
# Use pmin and pmax to ensure ID1 is always the smaller and ID2 is the larger in each pair
inter_BB$ID_min <- pmin(inter_BB$ID_1, inter_BB$ID_2)
inter_BB$ID_max <- pmax(inter_BB$ID_1, inter_BB$ID_2)

# 2: Count occurrences of each sorted pair
inter_BB$pair <- paste(inter_BB$ID_min, inter_BB$ID_max, sep = "_")
pair_counts <- table(inter_BB$pair)

# 3: Filter pairs that appear more than once (indicating bidirectional interactions)
interactions_overlap <- names(pair_counts[pair_counts > 1])


##############################################################################################################

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
prova <- data_overlap %>% group_by(pair) %>% summarise(id1 = length(unique(baitID)),
                                              id2 = length(unique(oeID)))

table(prova$id1)
table(prova$id2)

# We have different chicago score, for different directions of the interaction


### Annotate interactions

inter <- load_interactions(file)
inter_annotated <- annotate_interactions(inter, annotation = "/home/acost1/BALL_project/data/lichic/annotation_example.txt")


inter_annotated <- as.data.frame(inter_annotated)




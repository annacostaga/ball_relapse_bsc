annotate_interactions <- function (interactions, annotation, ...) {
  annotation_file <- annotation
  if (methods::is(annotation, "character")) {
    annotation <- data.table::fread(annotation, stringsAsFactors = F, 
                                    ...)
  }
  if (any(class(annotation) == "data.frame")) {
    if (ncol(annotation) != 5) {
      stop("File has not 5 columns \nAnnotation must be chr start end fragmentID annotation")
    }
    else {
      cn <- colnames(annotation)
      annotation[[cn[5]]][is.na(annotation[[cn[5]]])] <- "non-annotated"
      annotationGR <- GenomicRanges::makeGRangesFromDataFrame(annotation, 
                                                              seqnames.field = cn[1], start.field = cn[2], 
                                                              end.field = cn[3], keep.extra.columns = T)
      annot <- GenomicRanges::split(annotationGR[, -1], 
                                    as.factor(annotation[[cn[5]]]))
      annotation.features = list(annot = annot)
      suppressMessages(GenomicInteractions::annotateInteractions(interactions, 
                                                                 annotation.features))
      interactions@elementMetadata[, "bait_1"] <- unlist(GenomicInteractions::anchorOne(interactions)@elementMetadata[, 
                                                                                                                      "annot.id"])
      interactions@elementMetadata[, "bait_2"] <- unlist(GenomicInteractions::anchorTwo(interactions)@elementMetadata[, 
                                                                                                                      "annot.id"])
      interactions@elementMetadata[is.na(interactions@elementMetadata[, 
                                                                      "bait_2"]), "bait_2"] <- "."
      interactions@elementMetadata[is.na(interactions@elementMetadata[, 
                                                                      "bait_1"]), "bait_1"] <- "."
      interactions <- annotate_BOE(interactions)
      # cond <- ((interactions$ID_1 > interactions$ID_2) & 
      #            interactions$int == "B_B") | ((interactions$ID_1 < 
      #                                             interactions$ID_2) & interactions$int == "OE_B")
      # a1 <- interactions@anchor1[cond]
      # a2 <- interactions@anchor2[cond]
      # interactions@anchor1[cond] <- a2
      # interactions@anchor2[cond] <- a1
      # cols <- sort(grep("_", colnames(S4Vectors::elementMetadata(interactions[cond]))[1:4], 
      #                   value = T))
      # S4Vectors::elementMetadata(interactions[cond])[cols] <- S4Vectors::elementMetadata(interactions[cond])[cols[c(rbind(seq(2, 
      #                                                                                                                         length(cols), 2), seq(1, length(cols), 2)))]]
      # interactions <- annotate_BOE(interactions)
      # param <- getParameters(interactions)
      # param$annotate <- c(annotation_file = normalizePath(annotation_file))
      # interactions <- setParameters(interactions, param)
      interactions <- interactions[order(interactions$ID_1, 
                                         interactions$ID_2)]
      return(interactions)
    }
  }
}

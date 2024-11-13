#' Digest genome for a specific restriction enzyme
#'
#' This function takes a genome installed and generates its digest for a given restriction enzyme


digest_genome <- function(genome="GRCh38",RE_name="hindIII",motif="AAGCTT",cut_position=1,select_chr=c(1:22,"X","Y"),PAR_mask=T,PAR_file=NULL,...)
{
  genome_name <- genome
  ## Load reference genome from installed ones
  genome <-  tryCatch({
    load_genome(genome = genome)
  }, error = function(e) {
    stop(e)
  })
  ## If Pseudoautosomal Regions (PAR) masked, read provided file for human
  if (PAR_mask)
  {
    if (is.null(PAR_file))
    {
      PAR_file <- system.file("extdata", paste0("PAR_",gsub(" ","_",metadata(genome)$organism),"_coordinates.txt"), package="HiCaptuRe")
      if (file.exists(PAR_file))
      {
        PAR <- read.delim(PAR_file,header = T,...)
      } else
      {
        stop(paste("There is no PAR file provided for",metadata(genome)$organism,"\n It must be a headed file with seqnames,start,end columns"))
      }
    } else
    {
      if (file.exists(PAR_file))
      {
        PAR <- read.delim(PAR_file,header = T,...)
      } else
      {
        stop("The PAR file provided doesn't exist")
      }
    }
  } ## PAR mask
  
  ## Select primary chromosomes
  chrs <- GenomeInfoDb::seqnames(genome)
  
  if (!is.null(select_chr))
  {
    chrs <- chrs[chrs %in% select_chr]
  }
  if (length(chrs) == 0)
  {
    stop("Chromosomes selected in 'select_chr' are not present in this genome. Please try changing 'select_chr' or setting it to 'NULL'")
  }
  
  p <- progressr::progressor(steps = length(chrs))
  
  
  ## Digest genome by chromosomes
  digest <- data.frame()
  for (chr in chrs)
  {
    p(sprintf(paste("Digesting",chr)))
    if (PAR_mask)
    {
      if (grepl(unique(PAR$seqnames),chr))
      {
        chr_seq <- Biostrings::replaceAt(genome[[chr]], IRanges::IRanges(PAR$start[1], PAR$end[1]),
                                         Biostrings::DNAStringSet(strrep("N",length(PAR$start[1]:PAR$end[1]))))
        chr_seq <- Biostrings::replaceAt(chr_seq, IRanges::IRanges(PAR$start[2], PAR$end[2]),
                                         Biostrings::DNAStringSet(strrep("N",length(PAR$start[2]:PAR$end[2]))))
      }
      chr_seq <- genome[[chr]]
    }
    else
    {
      chr_seq <- genome[[chr]]
    }
    
    m <- Biostrings::matchPattern(motif, chr_seq)
    
    correct <- start(m)-1+cut_position
    
    starts <- c(1,correct+1)
    ends <- c(correct,length(chr_seq))
    
    df <- data.frame(seqnames=chr,start=starts,end=ends)
    
    digest <- rbind(digest,df)
  }
  
  if (is.null(PAR_file))
  {
    PAR_file <- "NULL"
  } else
  {
    PAR_file <- normalizePath(PAR_file)
  }
  
  digest$fragment_ID <- 1:nrow(digest)
  output <- list(digest=digest,
                 parameters=c("Genome"=genome_name,
                              "Genome_Package"=genome@pkgname,
                              "Restriction_Enzyme"=RE_name,
                              "Motif"=motif,
                              "Cut_Position"=cut_position,
                              "Selected_Chromosomes"=paste(select_chr,collapse = ","),
                              "PAR_mask"=PAR_mask,
                              "PAR_file"=PAR_file),
                 seqinfo=seqinfo(genome)[chrs])
  return(output)
}



annotate_BOE <- function(interactions){
  a1 <- GenomicInteractions::anchorOne(interactions)
  a1$bait_1 <- interactions@elementMetadata[,"bait_1"]
  a2 <- GenomicInteractions::anchorTwo(interactions)
  a2$bait_1 <- interactions@elementMetadata[,"bait_2"]
  
  a2$bait_1[is.na(a2$bait_1)] <- "."
  
  regions <- unique(c(a1,a2))
  
  regions$bait_1[is.na(regions$bait_1)] <- "non-annotated"
  bait <- regions[(regions$bait_1 != ".")]
  oe <- regions[regions$bait_1 == "."]
  
  baitl <- GenomicRanges::split(bait[, -1], as.factor(bait$bait_1))
  oel <- GenomicRanges::split(oe[, -1], as.factor(oe$bait_1))
  
  annotation.features = list(B = baitl, OE = oel)
  GenomicInteractions::resetAnnotations(interactions)
  suppressMessages(GenomicInteractions::annotateInteractions(interactions, annotation.features))
  GenomicInteractions::annotateRegions(interactions,"fragmentID",unique(sort(c(interactions$ID_1,interactions$ID_2))))
  interactions$int <- paste(GenomicInteractions::anchorOne(interactions)$node.class,GenomicInteractions::anchorTwo(interactions)$node.class, sep = "_")
  
  return(interactions)
}




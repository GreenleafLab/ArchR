##########################################################################################
# Helper Functions for GenomicRanges
##########################################################################################

#' Filters unwanted seqlevels from a Genomic Ranges object or similar object
#'
#' This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
#'
#' @param gr A `GRanges` object or another object containing seqlevels.
#' @param remove A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' @param underscore A boolean value indicating whether to remove all seqlevels whose name contains an underscore (for example "chr11_KI270721v1_random").
#' @param standard A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by `GenomeInfoDb::keepStandardChromosomes()`.
#' @param pruning.mode See ?seqinfo for a description of the pruning modes.
#' @export
filterChrGR <- function(
    gr = NULL, 
    remove = c("chrM"), 
    underscore = TRUE, 
    standard = TRUE, 
    pruning.mode="coarse"
  ){

  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = remove, name = "remove", valid = c("character"))
  .validInput(input = underscore, name = "underscore", valid = c("boolean"))
  .validInput(input = standard, name = "standard", valid = c("boolean"))
  .validInput(input = pruning.mode, name = "pruning.mode", valid = c("character"))

  #first we remove all non standard chromosomes
  if(standard){
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = pruning.mode)
  }
  #Then check for underscores or specified remove
  seqNames <- seqlevels(gr)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  #next we remove all chr specified in remove
  chrRemove <- c(chrRemove, which(seqNames %in% remove))
  if(length(chrRemove) > 0){
    chrKeep <- seqNames[-chrRemove]
  }else{
    chrKeep <- seqNames
  }
  #this function restores seqlevels
  seqlevels(gr, pruning.mode=pruning.mode) <- chrKeep
  
  return(gr)

}

#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved. The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`, `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value that determines whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
	gr = NULL, 
	by = "score", 
	decreasing = TRUE, 
	verbose = FALSE
  ){
  
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = by, name = "by", valid = c("character"))
  .validInput(input = decreasing, name = "decreasing", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))

  stopifnot(by %in% colnames(mcols(gr)))

  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)

    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
    
    }else{
      grAll <- c(grAll, grSelect)
    } 

  }
  message(sprintf("Converged after %s iterations!", i))

  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))

  return(grAll)

}

#' Subset a Genomic Ranges object by the provided seqnames
#'
#' This function returns a subsetted Genomic Ranges object based on a vector of provided seqnames
#'
#' @param gr A `GRanges` object to be subsetted.
#' @param names A character vector containing the `seqnames` to keep from the provided `GRanges` object.
#' @export
subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = names, name = "names", valid = c("character"))
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))
  return(gr)
}

#' Adds seqlength information to the seqnames of a Genomic Ranges object
#'
#' This function adds seqlength information for each of the seqnames in the provided Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param genome The name of a valid genome (for example "hg38", "hg19", or "mm10"). See `ArchR:::validBSgenome()`.
#' @export
addSeqLengthsGR <- function(gr = NULL, genome = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  genome <- .validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <- seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

#' Randomly shuffle a Genomic Ranges object
#'
#' This function randomly shuffles a Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param genome The name of a valid genome (for example "hg38", "hg19", or "mm10"). See `ArchR::validBSgenome()`.
#' @param n The number of permutations to perform during shuffling.
#' @param shuffleChr A boolean value indicating whether to shuffle across chromosomes randomly based on length of chromosomes or use previous knowledge of chromosome distribution.
#' @export
shuffleGR <- function(gr = NULL, genome = NULL, n = 100, shuffleChr = TRUE){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = n, name = "n", valid = c("integer"))
  .validInput(input = shuffleChr, name = "shuffleChr", valid = c("boolean"))
  #adapted from ChIPseeker's shuffle
  cs <- getChromSizes(genome)
  seqL <- seqlengths(cs)
  seqL <- seqL[sort(names(seqL))]
  sub <- subsetSeqnamesGR(gr, names = names(seqL)) #change
  sub <- sub[order(as.character(seqnames(sub)))]
  #stopifnot(identical(unique(as.character(seqnames(sub))), names(seqL)))
  w <- width(sub)
  name <- mcols(sub)$name
  seqN <- as.character(seqnames(sub))
  if(shuffleChr){
    expected <- round(length(sub)*as.numeric(seqL)/sum(as.numeric(seqL)))
    names(expected) <- names(seqL)
    #hackish rounding correction
    diff <- sum(length(sub))-sum(expected)
    r <- sample(length(expected),1)
    expected[r] <- expected[r] + diff
    subPerChr <- expected
  }else{
    subPerChr <- table(seqN)
  }
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)  
  grL <- lapply(seq_len(n), function(x){
    setTxtProgressBar(pb,round(x*100/n,0))
    rand <- sample(length(w))
    ws <- w[rand]
    ns <- name[rand]
    d <- lapply(seq_along(subPerChr), function(i){
      st_i <- sample(seqL[i],subPerChr[i])
      return(data.frame(seq = names(subPerChr)[i], start = st_i))
    }) %>% data.table::rbindlist(.) %>% data.frame
    gr <- GRanges(seqnames=d[,1], ranges=IRanges(d[,2], width=ws), strand="*", name = ns)
    suppressWarnings(seqlengths(gr) <- seqL)
    gr <- trim(gr)
    return(gr)
  })
  grL <- GenomicRangesList(grL)
  return(grL)
}

#' Merge genomic regions within a single Genomic Ranges object
#'
#' This function merges overlapping regions within a single Genomic Ranges object
#'
#' @param gr A `GRanges` object.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
mergeGR <- function(gr, ignore.strand = TRUE){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = ignore.strand, name = "ignore.strand", valid = c("boolean"))
  grR <- reduce(gr,min.gapwidth=0L,ignore.strand = ignore.strand)
  o <- DataFrame(findOverlaps(grR, gr,ignore.strand = ignore.strand))
  o$start <- start(gr[o$subjectHits])
  o$end <- end(gr[o$subjectHits])
  o$chr <- seqnames(gr[o$subjectHits])
  os <- o[order(o$start,decreasing = FALSE),] %>% {.[!duplicated(.$queryHits),c("queryHits", "start")]}
  oe <- o[order(o$end,decreasing = TRUE),] %>% {.[!duplicated(.$queryHits),c("queryHits", "end")]}
  oc <- o[!duplicated(o$queryHits),c("queryHits", "chr")]
  df <- merge(merge(oc, os, by = "queryHits"),oe,  by = "queryHits")
  mGR <- GRanges(df[,2], ranges = IRanges(df[,3], df[,4])) %>% sortSeqlevels %>% sort
  return(mGR)
}

#' Extend regions from a Genomic Ranges object
#'
#' This function extends each region in a Genomic Ranges object by a designated upstream and downstream extension in a strand-aware fashion
#'
#' @param gr A `GRanges` object.
#' @param upstream The number of basepairs upstream (5') to extend each region in `x`. Strand-aware.
#' @param downstream The number of basepairs downstream (3') to extend each region in `x`. Strand-aware.
#' @export
extendGR <-  function(gr = NULL, upstream = NULL, downstream = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  isMinus <- BiocGenerics::which(strand(gr) == "-")
  isOther <- BiocGenerics::which(strand(gr) != "-")
  #Forward
  start(gr)[isOther] <- start(gr)[isOther] - upstream
  end(gr)[isOther] <- end(gr)[isOther] + downstream
  #Reverse
  end(gr)[isMinus] <- end(gr)[isMinus] + upstream
  start(gr)[isMinus] <- start(gr)[isMinus] - downstream
  return(gr)
}

#' Identify the number of bases that overlap two Genomic Ranges objects
#'
#' This function returns a data.frame describing how many basepairs overlap the provided query and subject Genomic Ranges objects
#'
#' @param query A `GRanges` object to be used as the query in `findOverlaps()`.
#' @param subject A `GRanges` object to be used as the subject in `findOverlaps()`.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
nOverlapGR <- function(query = NULL, subject = NULL, ignore.strand = TRUE){
  .validInput(input = query, name = "query", valid = c("GRanges"))
  .validInput(input = subject, name = "subject", valid = c("GRanges"))
  .validInput(input = ignore.strand, name = "ignore.strand", valid = c("boolean"))
  o <- findOverlaps(query, subject, ignore.strand = ignore.strand)
  overlaps <- pintersect(query[queryHits(o)], subject[subjectHits(o)])
  percentOverlap <- width(overlaps) / width(subject[subjectHits(o)])
  l <- unlist(lapply(split(percentOverlap, subjectHits(o)),function(x)sum(x)))
  perOverlap <- sum(l*width(subject[unique(subjectHits(o))]))/sum(width(subject))
  nBP <- perOverlap * sum(width(subject))
  type <- c("queryBP", "sharedBP", "subjectBP")
  nBases <- c(sum(width(query))-nBP, nBP, sum(width(subject))-nBP)
  return(data.frame(type = type, bp = nBases))
}

#' Overlap with many genomic regions
#'
#' This function returns a sparse matrix that describes the overlap with each sub-grouped genomic region as specified in the "by" column.
#'
#' @param query A `GRanges` object to be used as the query in `findOverlaps()`.
#' @param subject A `GRanges` object with a column sub-grouping to be used as the subject in `findOverlaps()`.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be sub-grouped.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
overlapsManyGR <- function(query = NULL, subject = NULL, by = NULL, ignore.strand = TRUE){
  .validInput(input = query, name = "query", valid = c("GRanges"))
  .validInput(input = subject, name = "subject", valid = c("GRanges"))
  .validInput(input = by, name = "by", valid = c("character"))
  .validInput(input = ignore.strand, name = "ignore.strand", valid = c("boolean"))
  o <- DataFrame(findOverlaps(query, subject, ignore.strand = ignore.strand))
  o$name <- mcols(subject)[o$subjectHits, by]
  o$id <- match(o$name, unique(o$name))
  sparse <- Matrix::sparseMatrix(
    i = o[,1],
    j = o[,4],
    x = rep(TRUE,nrow(o)),
    dims = c(length(query),length(unique(o$name)))
  )
  colnames(sparse) <- unique(o$name)
  return(sparse)
}

#' Construct a Genomic Ranges object taking into account strandedness
#'
#' This function creates a Genomic Ranges object accounting for strandedness indicated by the relative orientation of the provided start and end positions
#'
#' @param seqnames A character vector containing the seqnames to be added to the `GRanges` object.
#' @param start A vector of start positions to be added to the `GRanges` object.
#' @param end A vector of end positions to be added to the `GRanges` object.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
constructGR <- function(seqnames = NULL, start = NULL, end = NULL, ignore.strand = TRUE){
  .validInput(input = seqnames, name = "seqnames", valid = c("character"))
  .validInput(input = start, name = "start", valid = c("integer"))
  .validInput(input = end, name = "end", valid = c("integer"))
  .validInput(input = ignore.strand, name = "ignore.strand", valid = c("boolean"))
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignore.strand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-" 
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)
}


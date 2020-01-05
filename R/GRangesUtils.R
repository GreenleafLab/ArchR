##########################################################################################
# Helper Functions for GenomicRanges
##########################################################################################

#' Filters unwanted seqlevels from a Genomic Ranges object or similar object
#'
#' This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
#'
#' @param x A `GRanges` object or another object containing seqlevels.
#' @param remove A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' @param underscore A boolean value indicating whether to remove all seqlevels whose name contains an underscore (for example "chr11_KI270721v1_random").
#' @param standard A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by `GenomeInfoDb::keepStandardChromosomes()`.
#' @export
keepFilteredChromosomes <- function(x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode="coarse"){
  #first we remove all non standard chromosomes
  if(standard){
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  #Then check for underscores or specified remove
  seqNames <- seqlevels(x)
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
  seqlevels(x, pruning.mode=pruning.mode) <- chrKeep
  return(x)
}

#' QQQ Retreive column metadata for overlapping Genomic Ranges regions
#'
#' QQQ This function returns a data.frame containing the overlapping regions and designated metadata columns from mcols(gr)
#'
#' @param query A `GRanges` object to be used as the query in `findOverlaps()`.
#' @param subject A `GRanges` object to be used as the subject in `findOverlaps()`.
#' @param colname The column name from `mcols(gr)` to return for the overlapping regions. Mandatory parameter; cannot be `NULL`.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @param decreasing A boolean value indicating whether the returned `data.frame` should be ordered in decreasing order based on `colname`.
#' @export
columnOverlaps <- function(query, subject, colname = "score", ignore.strand = TRUE, decreasing = TRUE){
  #First get overlaps
  o <- data.frame(findOverlaps(query, subject, ignore.strand = ignore.strand))
  #Then append information
  o$col <- mcols(subject)[[colname]][o[,2]]
  #Order it by the factor to rank
  o <- o[order(o$col, decreasing = decreasing),]
  #Deduplicate
  o <- o[!duplicated(o$queryHits),]
  #Initialize
  val <- rep(0, length(query))
  #Fill Values
  val[o[,1]] <- o$col
  return(val)
}

#' QQQ Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' QQQ This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by QQQ The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved. The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`, `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the lowest "score" in the overlap will be retained.
#' @param decreasing QQQ A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value that determines whether the output should include extra reporting.
#' @export
nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
  
  stopifnot(by %in% colnames(mcols(gr)))
  gr <- .validGRanges(gr)
  
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
#' @param seqNames A character vector containing the `seqnames` to keep from the provided `GRanges` object.
#' @param useNames QQQ A boolean value indicating whether QQQ.
#' @export
subsetSeqnames <- function(gr, seqNames, useNames = FALSE){
  gr <- .validGRanges(gr)
  gr <- gr[which(as.character(seqnames(gr)) %in% seqNames),]
  if(useNames){
    seqlevels(gr) <- seqNames
  }else{
    seqlevels(gr) <- as.character(unique(seqnames(gr)))
  }
  return(gr)
}

#' QQQ Adds seqlength information to the seqnames of a Genomic Ranges object
#'
#' QQQ This function adds seqlength information for each of the seqnames in the provided Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param genome The name of a valid genome (for example "hg38", "hg19", or "mm10"). See `ArchR::validBSgenome()`.
#' @export
addSeqLengths <- function(gr, genome){
  gr <- .validGRanges(gr)
  genome <- validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <- seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

#' QQQ Randomly shuffle a Genomic Ranges object
#'
#' QQQ This function randomly shuffles a Genomic Ranges object.
#'
#' @param subject A `GRanges` object.
#' @param genome The name of a valid genome (for example "hg38", "hg19", or "mm10"). See `ArchR::validBSgenome()`.
#' @param n QQQ The number of permutations to perform during shuffling.
#' @param shuffleChr QQQ A boolean value indicating whether to shuffle across chromosomes randomly or to QQQ (WHAT?) use previous knowledge of chromosome distribution
#' @export
shuffleGRanges <- function(subject, genome, n, shuffleChr=TRUE){
  #adapted from ChIPseeker's shuffle
  cs <- getChromSizes(genome)
  seqL <- seqlengths(cs)
  seqL <- seqL[sort(names(seqL))]
  sub <- subsetSeqnames(subject, seqNames = names(seqL)) #change
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

#' QQQ Merge regions within a single Genomic Ranges object
#'
#' QQQ This function merges overlapping regions within a single Genomic Ranges object
#'
#' @param gr A `GRanges` object.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
mergeGRanges <- function(gr, ignore.strand = TRUE){
  gr <- .validGRanges(gr)
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

#' QQQ Extend regions from a Genomic Ranges object
#'
#' QQQ This function extends each region in a Genomic Ranges object by a designated upstream and downstream extension in a strand-aware fashion
#'
#' @param x A `GRanges` object.
#' @param upstream The number of basepairs upstream (5') to extend each region in `x`. Strand-aware.
#' @param downstream The number of basepairs downstream (3') to extend each region in `x`. Strand-aware.
#' @export
extendGRanges <-  function(x, upstream, downstream){
  #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  isMinus <- BiocGenerics::which(strand(x) == "-")
  isOther <- BiocGenerics::which(strand(x) != "-")
  #Forward
  start(x)[isOther] <- start(x)[isOther] - upstream
  end(x)[isOther] <- end(x)[isOther] + downstream
  #Reverse
  end(x)[isMinus] <- end(x)[isMinus] + upstream
  start(x)[isMinus] <- start(x)[isMinus] - downstream
  return(x)
}

#' QQQ Identify the number of bases that overlap two Genomic Ranges objects
#'
#' QQQ This function returns a data.frame describing how many basepairs overlap the provided query and subject Genomic Ranges objects
#'
#' @param query A `GRanges` object to be used as the query in `findOverlaps()`.
#' @param subject A `GRanges` object to be used as the subject in `findOverlaps()`.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
overlappingBP <- function(query, subject, ignore.strand = TRUE){
  query <- .validGRanges(query)
  subject <- .validGRanges(subject)
  o <- findOverlaps(query, subject, ignore.strand = ignore.strand)
  overlaps <- pintersect(query[queryHits(o)], subject[subjectHits(o)])
  percentOverlap <- width(overlaps) / width(subject[subjectHits(o)])
  l <- unlist(lapply(split(percentOverlap, subjectHits(o)),function(x)sum(x)))
  perOverlap <- sum(l*width(subject[unique(subjectHits(o))]))/sum(width(subject))
  nBP <- perOverlap * sum(width(subject))
  type <- c("queryBP", "sharedBP", "subjectBP")
  nBases <- c(sum(width(query))-nBP, nBP, sum(width(subject))-nBP)
  return(data.frame(type,nBases))
}

#' QQQ PLEASE ADD A DESCRIPTIVE FUNCTION TITLE
#'
#' QQQ This function returns a sparse matrix QQQ PLEASE ADD MORE DETAIL????
#'
#' @param query A `GRanges` object to be used as the query in `findOverlaps()`.
#' @param subject A `GRanges` object to be used as the subject in `findOverlaps()`.
#' @param by QQQ The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved. The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`, `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the lowest "score" in the overlap will be retained.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
overlapsMany <- function(query, subject, by, ignore.strand = TRUE){
  o <- DataFrame(findOverlaps(query, subject, ignore.strand = ignore.strand))
  o$name <- mcols(subject)[o$subjectHits,by]
  o$id <- match(o$name, unique(o$name))
  sparse <- Matrix::sparseMatrix(
    i=o[,1],
    j=o[,4],
    x=rep(TRUE,nrow(o)),
    dims=c(length(query),length(unique(o$name)))
  )
  colnames(sparse) <- unique(o$name)
  return(sparse)
}

#' QQQ Construct a Genomic Ranges object taking into account strandedness
#'
#' QQQ This function creates a Genomic Ranges object accounting for strandedness indicated by the relative orientation of the provided start and end positions
#'
#' @param seqnames A character vector containing the seqnames to be added to the `GRanges` object.
#' @param start A vector of start positions to be added to the `GRanges` object.
#' @param end A vector of end positions to be added to the `GRanges` object.
#' @param ignore.strand A boolean value indicating whether strandedness should be ignored in `findOverlaps()`.
#' @export
constructGRanges <- function(seqnames, start, end, ignore.strand = TRUE){
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


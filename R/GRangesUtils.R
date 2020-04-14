##########################################################################################
# Helper Functions for GenomicRanges
##########################################################################################

#' Filters unwanted seqlevels from a Genomic Ranges object or similar object
#'
#' This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
#'
#' @param gr A `GRanges` object or another object containing `seqlevels`.
#' @param remove A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' If no manual removal is desired, `remove` should be set to `NULL`.
#' @param underscore A boolean value indicating whether to remove all seqlevels whose names contain an underscore (for example "chr11_KI270721v1_random").
#' @param standard A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by
#' `GenomeInfoDb::keepStandardChromosomes()`.
#' @param pruningMode The name of the pruning method to use (from`GenomeInfoDb::seqinfo()`) when seqlevels must be removed from a `GRanges` object.
#' When some of the seqlevels to drop from the given `GRanges` object are in use (i.e. have ranges on them), the ranges on these sequences need
#' to be removed before the seqlevels can be dropped. Four pruning modes are currently defined: "error", "coarse", "fine", and "tidy".
#' @export
filterChrGR <- function(
    gr = NULL, 
    remove = NULL, 
    underscore = TRUE, 
    standard = TRUE, 
    pruningMode="coarse"
  ){

  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = remove, name = "remove", valid = c("character", "null"))
  .validInput(input = underscore, name = "underscore", valid = c("boolean"))
  .validInput(input = standard, name = "standard", valid = c("boolean"))
  .validInput(input = pruningMode, name = "pruningMode", valid = c("character"))

  #first we remove all non standard chromosomes
  if(standard){
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = pruningMode)
  }
  #Then check for underscores or specified remove
  seqNames <- seqlevels(gr)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  #next we remove all chr specified in remove
  if(!is.null(remove)){
    chrRemove <- c(chrRemove, which(seqNames %in% remove))
  }
  if(length(chrRemove) > 0){
    chrKeep <- seqNames[-chrRemove]
  }else{
    chrKeep <- seqNames
  }
  #this function restores seqlevels
  seqlevels(gr, pruning.mode=pruningMode) <- chrKeep
  
  return(gr)

}

#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output should include extra reporting.
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
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE){
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
    grSelect <- .clusterGRanges(
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

#' Extend regions from a Genomic Ranges object
#'
#' This function extends each region in a Genomic Ranges object by a designated upstream and downstream extension in a strand-aware fashion
#'
#' @param gr A `GRanges` object.
#' @param upstream The number of basepairs upstream (5') to extend each region in `gr` in a strand-aware fashion.
#' @param downstream The number of basepairs downstream (3') to extend each region in `gr` in a strand-aware fashion.
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




#--------------------------------------------------------------------------------------------
# Helper Functions for GenomicRanges
#--------------------------------------------------------------------------------------------

#' Filters unwanted chr mainly underscores
#' @param x GRanges or something with seqlevels
#' @param remove remove vector
#' @param underscore remove all underscores?
#' @param standard keep standard chromosomes
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

#' Instead of counting overlaps get columns like max score or etc in query
#' @param query granges query
#' @param subject granges subject
#' @param colname mcols(gr)[[colname]] cannot be null
#' @param decreasing for order
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

#' Instead of counting overlaps get columns like max score or etc in query
#' @param query granges query
#' @param subject granges subject
#' @param colname mcols(gr)[[colname]] cannot be null
#' @param decreasing for order
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

#' Subset by Seqnames
#' @param gr grange
#' @param seqnames seqnames to subset
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

#' Add Seqlengths to genomic ranges
#' @param gr see validGRanges
#' @param genome see validBSgenome
#' @export
addSeqLengths <- function(gr, genome){
  gr <- .validGRanges(gr)
  genome <- validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <- seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

#' Shuffle Genomic Ranges
#' @param subject see validGRanges
#' @param genome see validBSgenome
#' @param n nPermutations
#' @param shuffleChr shuffle across chromosomes randomly vs using previous knowledge of chromosome distribution
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

#' Merge Genomic Ranges
#' @param gr see validGRanges
#' @param ignore.strand ignore strandedness for merging
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

#' Merge Genomic Ranges
#' @param query see validGRanges
#' @param subject see validGRanges
#' @param ignore.strand ignore strandedness for overlaps
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

#' Merge Genomic Ranges
#' @param query see validGRanges
#' @param subject see validGRanges
#' @param ignore.strand ignore strandedness for overlaps
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

#' Overlaps Many includes information from mcols(gr)
#' @param query see validGRanges
#' @param subject see validGRanges
#' @param by column in subject to split overlaps by
#' @param ignore.strand ignore strandedness for overlaps
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

#' Construct GRanges seqnames start end accounting for ends before starts (adding strandedness)
#' @param seqnames seqnames of GRanges
#' @param start start of GRanges
#' @param end end of GRanges
#' @param ignore.strand ignore strandedness for overlaps
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


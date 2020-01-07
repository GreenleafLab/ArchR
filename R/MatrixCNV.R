####################################################################
# Copy Number Variation Methods
####################################################################

#' Add a CNV matrix to ArrowFiles or an ArchRProject
#' 
#' This function for each sample will predict copy number variation from accessibility
#'
#' @param input An `ArchRProject` object or character vector of ArrowFiles.
#' @param chromSizes A named numeric vector containing the chromsome names and lengths. The default behavior is to retrieve this from the `ArchRProject` using `ArchR::getChromSizes()`.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist from calling CNVs. The default behavior is to retrieve this from the `ArchRProject` using `ArchR::getBlacklist()`.
#' @param genome The genome used by the `input`. The default behavior is to retrieve this from the `ArchRProject` using `ArchR::getGenome()`.
#' @param windowSize The size in basepairs for the sliding window used to break up each chromosome to look for CNVs.
#' @param stepSize The size in basepairs for the step used to create sliding window bins across each chromosome.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from CNV analysis.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam QQQ A list of parameters to be passed to QQQ for batch-style parallel computing.
#' @param force A boolean value indicating whether to force the CNV matrix to be overwritten if it already exist for `input`.
#' @export
addCNVMatrix <- function(
  input,
  chromSizes = getChromSizes(input),
  blacklist = getBlacklist(input),
  genome = getGenome(input),
  windowSize = 10e6, 
  stepSize = 2e6,
  excludeChr = c("chrM","chrY"),
  threads = 1,
  parallelParam = NULL,
  force = FALSE,
  ...
  ){

  if(inherits(input, "ArchRProject")){
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
  }else if(inherits(input, "character")){
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
  }else{
    stop("Error Unrecognized Input!")
  }
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  #First we need to make windows for CNV
  windows <- .makeWindows(
    genome = genome,
    chromSizes = chromSizes, 
    blacklist = blacklist, 
    windowSize = windowSize, 
    stepSize = stepSize,
    threads = threads
  )
  windows <- windows[BiocGenerics::which(seqnames(windows) %bcni% excludeChr)]

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addCNVMatrix
  args$windows <- windows
  args$registryDir <- file.path(outDir, "CNVRegistry")

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(inherits(input, "ArchRProject")){
    return(input)
  }else{
    return(unlist(outList))
  }

}

.addCNVMatrix <- function(
  i,
  ArrowFiles, 
  normByNeighbors = TRUE,
  cellNames = NULL, 
  allCells = NULL,
  windows = NULL,
  force = FALSE,
  ...
  ){

  ArrowFile <- ArrowFiles[i]
  o <- h5closeAll()
  
  #Check
  if(!suppressMessages(h5createGroup(file = ArrowFile, "CNVMatrix"))){
    if(force){
      o <- h5delete(file = ArrowFile, name = "CNVMatrix")
      o <- h5createGroup(ArrowFile, "CNVMatrix")
    }else{
      stop(sprintf("%s Already Exists!, set force = TRUE to override!", "CNVMatrix"))
    }
  }
  
  tstart <- Sys.time()
 
  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  ######################################
  # Create Window Matrices Then Summarize to CNV Estimates
  ######################################
  uniqueChr <- as.character(unique(seqnames(windows)@values))

  seWindows <- .safelapply(seq_along(uniqueChr), function(x){

    o <- h5closeAll()
    chr <- uniqueChr[x]
    windowsx <- windows[BiocGenerics::which(seqnames(windows)==chr)]
    rangesx <- ranges(windowsx)
    .messageDiffTime(sprintf("Counting Windows for Chromosome %s of %s!", x, length(uniqueChr)), tstart)    

    #Read in Fragments
    fragments <- .getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)

    #Count Left Insertion
    temp <- IRanges(start = start(fragments), width = 1)
    stopifnot(length(temp) == length(fragments))
    oleft <- findOverlaps(ranges(rangesx), temp)
    oleft <- DataFrame(queryHits=Rle(queryHits(oleft)), subjectHits = subjectHits(oleft))

    #Count Right Insertion
    temp <- IRanges(start = end(fragments), width = 1)
    stopifnot(length(temp) == length(fragments))
    oright <- findOverlaps(ranges(rangesx), temp)
    oright <- DataFrame(queryHits=Rle(queryHits(oright)), subjectHits = subjectHits(oright))
    remove(temp)

    #Correct to RG ID
    oleft$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oleft$subjectHits], cellNames))
    oright$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oright$subjectHits], cellNames))
    remove(fragments)

    #Create Sparse Matrix
    mat <- Matrix::sparseMatrix(
      i = c( oleft$queryHits, oright$queryHits ),
      j = c( oleft$subjectHits, oright$subjectHits ),
      x = rep(1, nrow(oleft) + nrow(oright)),
      dims = c(length(rangesx), length(cellNames))
      )
    colnames(mat) <- cellNames

    #Summarize Split Windows
    windowSummary <- GRanges()
    countSummary <- matrix(nrow=length(unique(mcols(windowsx)$name)), ncol = ncol(mat))
    rownames(countSummary) <- unique(mcols(windowsx)$name)
    colnames(countSummary) <- cellNames

    for(y in seq_len(nrow(countSummary))){
      idx <- which(mcols(windowsx)$name == rownames(countSummary)[y])
      wx <- windowsx[idx]
      wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
      mcols(wo)$name <- mcols(wx)$name[1]
      mcols(wo)$effectiveLength <- sum(width(wx))
      mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
      mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
      mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
      mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
      countSummary[y,] <- Matrix::colSums(mat[idx, ,drop=FALSE]) / sum(width(wx))
      windowSummary <- c(windowSummary, wo)
    }
    seqlevels(windowSummary) <- uniqueChr

    SummarizedExperiment::SummarizedExperiment(assays = SimpleList(counts = countSummary), rowRanges = windowSummary)

  }, threads = 1) %>% Reduce("rbind", .)

  .messageDiffTime("Filtering Low Quality Windows", tstart)

  #Keep only regions with less than 0.1% N
  seWindows <- seWindows[which(mcols(seWindows)$N < 0.001) ]

  #Keep only regions with effective size greater than 90%
  seWindows <- seWindows[rowData(seWindows)$percentEffectiveLength >= 90]

  #Normalize By Nucleotide Content KNN
  .messageDiffTime("Normalizing GC-Bias", tstart)
  k <- 25
  seWindows <- seWindows[order(mcols(seWindows)$GC)]
  assays(seWindows)$log2GCNorm <- log2( (assays(seWindows)$counts + 1e-5) / (apply(assays(seWindows)$counts, 2, function(x) .centerRollMean(x, k)) + 1e-5))

  #Normalize By Cells With Similar Total Counts Genome Wide
  if(normByNeighbors){
    .messageDiffTime("Normalizing Coverage-Bias using Neighbors", tstart)
    totalCounts <- colSums(assays(seWindows)$counts * rowData(seWindows)$effectiveLength)
    seWindows <- seWindows[, order(totalCounts)]
    assays(seWindows)$log2GCNorm <- assays(seWindows)$log2GCNorm - t(apply(assays(seWindows)$log2GCNorm, 1, function(y) .centerRollMean(y, floor(0.025 * ncol(seWindows)))))
    seWindows <- seWindows[, cellNames]
  }

  #Re-Order Ranges and Smooth Across Individual Chromosomes
  .messageDiffTime("Smoothing CNV Scores", tstart)
  seWindows <- sort(sortSeqlevels(seWindows))
  seWindows <- lapply(uniqueChr, function(x){
    seWindowsx <- seWindows[BiocGenerics::which(seqnames(seWindows) %bcin% x)]
    assays(seWindowsx)$log2GCSmooth <- apply(assays(seWindowsx)$log2GCNorm, 2, function(y) .centerRollMean(y, ceiling( (windowSize / stepSize) / 2)))
    seWindowsx
  }) %>% Reduce("rbind", .)

  #Index Windows
  windows <- lapply(split(rowRanges(seWindows), seqnames(seWindows)), function(x) {
    mcols(x)$idx <- seq_along(x)
    x
  }) %>% Reduce("c", .) %>% sortSeqlevels %>% sort
  rowRanges(seWindows) <- windows[rownames(seWindows)]

  #Write To Arrow Files
  .messageDiffTime("Adding to Arrow Files", tstart)
  dfParams <- data.frame(
    windowSize = windowSize, 
    stepSize = stepSize,
    excludeChr = excludeChr,
    stringsAsFactors = FALSE)

  featureDF <- data.frame(
    seqnames = paste0(seqnames(seWindows)), 
    idx = mcols(seWindows)$idx, 
    start = start(seWindows), 
    end = end(seWindows), 
    name = mcols(seWindows)$name,
    GC = mcols(seWindows)$GC,
    effectiveLength = mcols(seWindows)$effectiveLength,
    stringsAsFactors = FALSE)

  ######################################
  # Initialize Mat Group
  ######################################
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = "CNVMatrix",
    Class = "Double",
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = force
  )

  ######################################
  # Write Mat Group
  ######################################
  uniqueChr <- unique(paste0(seqnames(seWindows)))
  seWindows <- seWindows[, cellNames]

  for(x in seq_along(uniqueChr)){

    #Write Matrix to Arrow File!
    o <- .addMatToArrow(
      mat = Matrix::Matrix(assays(seWindows[BiocGenerics::which(seqnames(seWindows) == uniqueChr[x])])$log2GCSmooth, sparse = TRUE), 
      ArrowFile = ArrowFile, 
      Group = paste0("CNVMatrix/", uniqueChr[x]), 
      binarize = FALSE
    ) 

    gc()

  }

  return(ArrowFile)

}

.makeWindows <- function(genome, chromSizes, blacklist, windowSize = 10e6, stepSize = 2e6, threads = 1){
  
  genome <- ArchR:::.validBSgenome(genome)

  #Sliding Windows
  windows <- slidingWindows(x = chromSizes, width = windowSize, step = stepSize) %>% 
    unlist %>% 
      .[which(width(.)==windowSize),]

  #Add OG Info Prior to Breaking Windows Up Etc
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- start(windows)
  mcols(windows)$wEnd <- end(windows)
  
  message("Subtracting Blacklist from Windows...")
  windowsBL <- .safelapply(seq_along(windows), function(x){
      if(x %% 100 == 0){
        message(sprintf("%s of %s", x, length(windows)))
      }
      gr <- GenomicRanges::setdiff(windows[x,], blacklist)
      mcols(gr) <- mcols(windows[x,])
      return(gr)
    }, threads = threads)
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  names(windowsBL) <- paste0("s",seq_along(windowsBL))

  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- .safelapply(seq_along(windowSplit), function(x){
      message(".", appendLF = FALSE)
      chrSeq <- Biostrings::getSeq(genome, chromSizes[BiocGenerics::which(seqnames(chromSizes) == names(windowSplit)[x])])
      grx <- windowSplit[[x]]
      aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
      mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
      mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
      return(grx)
    }, threads = threads) %>% GenomicRangesList %>% unlist %>% sortSeqlevels %>% sort
  message("\n")
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  windowNuc
}

# QQQ THIS SHOULD GET DELETED
# pdf("test.pdf", width = 12, height = 8)
# ArchR:::.ArchRHeatmap(mat = round(t(assays(seWindows)$log2GCSmooth), 3), clusterRows = FALSE, clusterCols = FALSE, limits = c(-2,2))
# dev.off()
# proj2 <- addCNVMatrix(proj, force = TRUE)


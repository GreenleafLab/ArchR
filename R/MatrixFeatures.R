####################################################################
# Peak and Feature Matrix Methods
####################################################################

#' Add a feature matrix to an ArchRProject or a set of ArrowFiles
#' 
#' This function for each sample will independently compute counts for each feature per cell in the provided ArchRProject or set of ArrowFiles.
#'
#' @param input An `ArchRProject` object or character vector of paths to ArrowFiles.
#' @param features A `GRanges` object containing the regions (aka features) to use for counting insertions for each cell.
#' @param matrixName The name to be used for storage of the feature matrix in the provided `ArchRProject` or ArrowFiles.
#' @param ceiling The maximum counts per feature allowed. This is used to prevent large biases in feature counts.
#' @param binarize A boolean value indicating whether the feature matrix should be binarized prior to storage. This can be useful for
#' downstream analyses when working with insertion counts.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exists in the `input`.
#' @export
addFeatureMatrix <- function(
  input = NULL,
  features = NULL,
  matrixName = "FeatureMatrix",
  ceiling = 10^9, 
  binarize = FALSE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE
  ){

  .validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  .validInput(input = features, name = "features", valid = c("GRanges"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  matrixName <- .isProtectedArray(matrixName)

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

  #Add Index To Features
  features <- sort(sortSeqlevels(features), ignore.strand = TRUE)
  features <- split(features, seqnames(features))
  features <- lapply(features, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })
  features <- Reduce("c",features)

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$features <- features
  args$FUN <- .addFeatureMatrix
  args$registryDir <- file.path(outDir, "CountFeaturesRegistry")

  #Remove input from args
  args$input <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(inherits(input, "ArchRProject")){
    return(input)
  }else{
    return(unlist(outList))
  }

}

#' Add a Peak Matrix to the ArrowFiles of an ArchRProject
#' 
#' This function, for each sample, will independently compute counts for each peak
#' per cell in the provided ArchRProject using the "PeakMatrix".
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param ceiling The maximum counts per feature allowed. This is used to prevent large biases in peak counts.
#' @param binarize A boolean value indicating whether the peak matrix should be binarized prior to storage. This can be useful
#' for downstream analyses when working with insertion counts.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the "PeakMatrix" to be overwritten if it already exist in the given `ArchRProject`.
#' @export
addPeakMatrix <- function(
  ArchRProj = NULL,
  ceiling = 4, 
  binarize = FALSE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = ceiling, name = "ceiling", valid = c("numeric"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  if(is.null(ArchRProj@peakSet)){
    stop("No peakSet found in ArchRProject!")
  }

  ArrowFiles <- getArrowFiles(ArchRProj)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)

  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName = "PeakMatrix"
  args$features <- ArchRProj@peakSet
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addFeatureMatrix
  args$registryDir <- file.path(outDir, "CountPeaksRegistry")

  #Remove project from args
  args$ArchRProj <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  readsInPeaks <- lapply(outList, function(x) x$RIP) %>% unlist
  FRIP <- lapply(outList, function(x) x$FRIP) %>% unlist
  ArchRProj <- addCellColData(ArchRProj, data = readsInPeaks, name = "ReadsInPeaks", names(readsInPeaks), force = force)
  ArchRProj <- addCellColData(ArchRProj, data = FRIP, name = "FRIP", names(readsInPeaks), force = force)
  return(ArchRProj)

}

.addFeatureMatrix <- function(
  i = NULL,
  ArrowFiles = NULL, 
  features = NULL,
  cellNames = NULL, 
  allCells = NULL,
  matrixName = "PeakMatrix", 
  ceiling = 4, 
  binarize = FALSE,
  tstart = NULL,
  subThreads = 1,
  force = FALSE
  ){

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)

  o <- h5closeAll()
  
  #Check
  if(!suppressMessages(h5createGroup(ArrowFile, matrixName))){
    if(force){
      o <- h5delete(file = ArrowFile, name = matrixName)
      o <- h5createGroup(ArrowFile, matrixName)
    }else{
      stop(sprintf("%s Already Exists!, set force = TRUE to override!", matrixName))
    }
  }
  
  if(!is.null(tstart)){
    tstart <- Sys.time()
  }
 
  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  dfParams <- data.frame(
    ceiling = ceiling, 
    binarize = binarize * 1,
    stringsAsFactors = FALSE)

  if("name" %in% colnames(mcols(features))){
    featureDF <- data.frame(
      seqnames = paste0(seqnames(features)), 
      idx = mcols(features)$idx, 
      start = start(features), 
      end = end(features), 
      name = mcols(features)$name,
      stringsAsFactors = FALSE)
  }else{
    featureDF <- data.frame(
      seqnames = paste0(seqnames(features)), 
      idx = mcols(features)$idx, 
      start = start(features), 
      end = end(features), 
      stringsAsFactors = FALSE)
  }
  
  ######################################
  # Initialize SP Mat Group
  ######################################
  if(binarize){
    Class <- "binary"
    Units <- "BinarizedCounts"
  }else{
    Class <- "integer"
    Units <- "Counts"
  }
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = Class,
    Units = Units,
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = force
  )

  ######################################
  # Add To SP Mat Group
  ######################################
  strand(features) <- "*" #lets make sure nothing strand related occurs
  uniqueChr <- as.character(unique(seqnames(features)@values))
  insertionsInPeaks <- rep(0, length(cellNames))
  names(insertionsInPeaks) <- cellNames
  totalInsertions <- insertionsInPeaks

  for(z in seq_along(uniqueChr)){

    o <- h5closeAll()
    chr <- uniqueChr[z]
    featurez <- features[BiocGenerics::which(seqnames(features)==chr)]
    .messageDiffTime(sprintf("Adding %s to %s for Chr (%s of %s)!", sampleName, matrixName, z, length(uniqueChr)), tstart)

    #Read in Fragments
    fragments <- .getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)
    tabFrags <- table(mcols(fragments)$RG)

    #Count Left Insertion
    temp <- IRanges(start = start(fragments), width = 1)
    stopifnot(length(temp) == length(fragments))
    oleft <- findOverlaps(ranges(featurez), temp)
    oleft <- DataFrame(queryHits=Rle(queryHits(oleft)), subjectHits = subjectHits(oleft))

    #Count Right Insertion
    temp <- IRanges(start = end(fragments), width = 1)
    stopifnot(length(temp) == length(fragments))
    oright <- findOverlaps(ranges(featurez), temp)
    oright <- DataFrame(queryHits=Rle(queryHits(oright)), subjectHits = subjectHits(oright))
    remove(temp)

    #Feature Idx
    oleft$queryHits@values <- mcols(featurez)$idx[oleft$queryHits@values]
    oright$queryHits@values <- mcols(featurez)$idx[oright$queryHits@values]

    #Correct to RG ID
    oleft$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oleft$subjectHits], cellNames))
    oright$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oright$subjectHits], cellNames))
    remove(fragments)

    #Create Sparse Matrix
    mat <- Matrix::sparseMatrix(
      i = c( oleft$queryHits, oright$queryHits ),
      j = c( oleft$subjectHits, oright$subjectHits ),
      x = rep(1, nrow(oleft) + nrow(oright)),
      dims = c(max(mcols(featurez)$idx), length(cellNames))
      )
    colnames(mat) <- cellNames
    
    #Compute total reads in Peak
    totalInsertions[names(tabFrags)] <- totalInsertions[names(tabFrags)] + 2 * tabFrags
    insertionsInPeaks <- insertionsInPeaks + Matrix::colSums(mat)

    #Ceiling
    if(!is.null(ceiling)){
      mat@x[mat@x > ceiling] <- ceiling
    }
    if(binarize){
      mat@x[mat@x > 0] <- 1
    }

    #Write sparseMatrix to Arrow File!
    o <- .addMatToArrow(
      mat = mat, 
      ArrowFile = ArrowFile, 
      Group = paste0(matrixName,"/", chr), 
      binarize = binarize,
      addColSums = TRUE,
      addRowSums = TRUE
      ) 
    gc()

  }

  out <- list(ArrowFile = ArrowFile, RIP = insertionsInPeaks, FRIP = insertionsInPeaks / totalInsertions)

  return(out)

}





####################################################################
# Tile Matrix Methods
####################################################################

#' Add TileMatrix to Arrows/ArchRProject
#' 
#' This function for each sample will independently compute counts for each tile
#' per cell in the Arrow File
#'
#' @param input An `ArchRProject` object or character vector of ArrowFiles.
#' @param chromSizes A named numeric vector containing the chromsome names and lengths. The default behavior is to retrieve this from the `ArchRProject` using `ArchR::getChromSizes()`.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist from calling CNVs. The default behavior is to retrieve this from the `ArchRProject` using `ArchR::getBlacklist()`.
#' @param tileSize The size of the tiles used for binning counts in the `TileMatrix`.
#' @param binarize QQQ A boolean value indicating whether the `TileMatrix` should be binarized QQQ prior to storage.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from CNV analysis.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam QQQ A list of parameters to be passed to QQQ for batch-style parallel computing.
#' @param force QQQ A boolean value indicating whether to force the `TileMatrix` to be overwritten if it already exist in the given `ArchRProject` or ArrowFiles.
#' @export
addTileMatrix <- function(
  input,
  chromSizes = getChromSizes(input),
  blacklist = getBlacklist(input),
  tileSize = 500, 
  binarize = TRUE, 
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

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addTileMat
  args$chromLengths <- end(chromSizes)
  names(args$chromLengths) <- paste0(seqnames(chromSizes))
  args$registryDir <- file.path(outDir, "CountTilesRegistry")

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(inherits(input, "ArchRProject")){
    return(input)
  }else{
    return(unlist(outList))
  }

}

.addTileMat <- function(
  i,
  ArrowFiles, 
  cellNames = NULL,
  allCells = NULL,
  tileSize = 500, 
  binarize = TRUE, 
  excludeChr = "chrY", 
  blacklist = NULL, 
  chromLengths = NULL, 
  force = FALSE,
  ...){

  ArrowFile <- ArrowFiles[i]

  o <- h5closeAll()
  
  #Check
  if(!suppressMessages(h5createGroup(file = ArrowFile, "TileMatrix"))){
    if(force){
      o <- h5delete(file = ArrowFile, name = "TileMatrix")
      o <- h5createGroup(ArrowFile, "TileMatrix")
    }else{
      stop("TileMatrix Already Exists!, set force = TRUE to override!")
    }
  }


  tstart <- Sys.time()
  if(!is.null(blacklist)){
    blacklist <- split(blacklist, seqnames(blacklist))
  }

  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  chromLengths <- chromLengths[names(chromLengths) %ni% excludeChr]
  if(length(chromLengths)==0){
    stop("Error removed all chromLengths with exclude chr!")
  }

  dfParams <- data.frame(
    seqnames = names(chromLengths), 
    length = as.vector(chromLengths), 
    tileSize = tileSize, 
    binarize = binarize,
    stringsAsFactors=FALSE)

  featureDF <- lapply(seq_along(chromLengths), function(x){
    DataFrame(seqnames = names(chromLengths)[x], idx = seq_len(trunc(chromLengths[x])/tileSize + 1))
  }) %>% Reduce("rbind", .)
  featureDF$start <- (featureDF$idx - 1) * tileSize 

  ######################################
  # Initialize SP Mat Group
  ######################################
  if(binarize){
    Class <- "binary"
  }else{
    Class <- "integer"
  }
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = "TileMatrix",
    Class = Class,
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = force
  )

  ######################################
  # Add To SP Mat Group
  ######################################
  for(z in seq_along(chromLengths)){

    o <- h5closeAll()
    chr <- names(chromLengths)[z]
    .messageDiffTime(sprintf("Adding Tile Matrix for Chromosome %s of %s to Arrow File!", z, length(chromLengths)), tstart)

    #Read in Fragments
    fragments <- .getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)

    #N Tiles
    nTiles <- trunc(chromLengths[z] / tileSize) + 1

    #Create Sparse Matrix
    matchID <- S4Vectors::match(mcols(fragments)$RG, cellNames)
    mat <- Matrix::sparseMatrix(
        i = c(trunc(start(fragments) / tileSize), trunc(end(fragments) / tileSize)) + 1,
        j = c(matchID, matchID),
        x = rep(1,  2*length(fragments)),
        dims = c(nTiles, length(cellNames))
      )
    colnames(mat) <- cellNames
    rm(fragments, matchID)
    gc()
    
    #Binarize
    if(binarize){
      mat@x[mat@x > 0] <- 1
    }

    #Remove Blacklisted Tiles!
    if(!is.null(blacklist)){
      blacklistz <- blacklist[[chr]]
      if(length(blacklistz) > 0){
        tile2 <- floor(tileSize/2)
        blacklistIdx <- unique(trunc(start(unlist(GenomicRanges::slidingWindows(blacklistz,tile2,tile2)))/tileSize) + 1)
        blacklistIdx <- sort(blacklistIdx)
        idxToZero <- which((mat@i + 1) %bcin% blacklistIdx)
        if(length(idxToZero) > 0){
          mat@x[idxToZero] <- 0
          mat <- Matrix::drop0(mat)
        }
      }
    }

    #Write sparseMatrix to Arrow File!
    o <- .addMatToArrow(
      mat = mat, 
      ArrowFile = ArrowFile, 
      Group = paste0("TileMatrix/", chr), 
      binarize = binarize,
      addColSums = TRUE,
      addRowSums = TRUE
      )

    gc()

  }

  return(ArrowFile)

}





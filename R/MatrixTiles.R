####################################################################
# Tile Matrix Methods
####################################################################

#' Add TileMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, for each sample, will independently compute counts for each tile
#'
#' @param input An `ArchRProject` object or character vector of ArrowFiles.
#' @param chromSizes A named numeric vector containing the chromsome names and lengths. The default behavior is to retrieve
#' this from the `ArchRProject` using `getChromSizes()`.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist counting in these tiles. The default behavior
#' is to retrieve this from the `ArchRProject` using `getBlacklist()`.
#' @param tileSize The size of the tiles used for binning counts in the "TileMatrix".
#' @param binarize A boolean value indicating whether the "TileMatrix" should be binarized prior to storage.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from the "TileMatrix".
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the "TileMatrix' to be overwritten if it already exist in the given `input`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addTileMatrix <- function(
  input = NULL,
  chromSizes = if(inherits(input, "ArchRProject")) getChromSizes(input) else NULL,
  blacklist = if(inherits(input, "ArchRProject")) getBlacklist(input) else NULL,
  tileSize = 500, 
  binarize = TRUE, 
  excludeChr = c("chrM", "chrY"),
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = FALSE,
  logFile = createLogFile("addTileMatrix")
  ){

  .validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  .validInput(input = chromSizes, name = "chromSizes", valid = c("GRanges"))
  .validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

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

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addTileMatrix Input-Parameters", logFile = logFile)

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addTileMat
  args$chromLengths <- end(chromSizes)
  names(args$chromLengths) <- paste0(seqnames(chromSizes))
  args$registryDir <- file.path(outDir, "CountTilesRegistry")

  #Remove Input from args
  args$input <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(inherits(input, "ArchRProject")){
    return(input)
  }else{
    return(unlist(outList))
  }

}

.addTileMat <- function(
  i = NULL,
  ArrowFiles = NULL, 
  cellNames = NULL,
  allCells = NULL,
  tileSize = 500, 
  binarize = TRUE, 
  excludeChr = c("chrM", "chrY"), 
  blacklist = NULL, 
  chromLengths = NULL, 
  chromSizes = NULL,
  force = FALSE,
  subThreads = 1,
  tstart = NULL,
  logFile = NULL
  ){

  .validInput(input = i, name = "i", valid = c("integer"))
  .validInput(input = ArrowFiles, name = "ArrowFiles", valid = c("character"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character", "null"))
  .validInput(input = allCells, name = "allCells", valid = c("character", "null"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  .validInput(input = chromLengths, name = "chromLengths", valid = c("numeric"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)
  
  #Check
  o <- h5closeAll()
  o <- .createArrowGroup(ArrowFile = ArrowFile, group = "TileMatrix", force = force, logFile = logFile)

  tstart <- Sys.time()
  if(!is.null(blacklist)){
    if(length(blacklist) > 0){
      blacklist <- split(blacklist, seqnames(blacklist))
    }
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
  .logThis(featureDF, paste0(sampleName, " .addTileMat FeatureDF"), logFile = logFile)

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
    Group = "TileMatrix",
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
  for(z in seq_along(chromLengths)){

    o <- tryCatch({

      o <- h5closeAll()
      chr <- names(chromLengths)[z]
      .logDiffTime(sprintf("Adding TileMatrix to %s for Chr (%s of %s)!", sampleName, z, length(chromLengths)), t1 = tstart, logFile = logFile)

      #Read in Fragments
      fragments <- .getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)

      #N Tiles
      nTiles <- trunc(chromLengths[z] / tileSize) + 1

      #Match Cells
      matchID <- S4Vectors::match(mcols(fragments)$RG, cellNames)

      #Log Info
      .logThis(nTiles, paste0("NTiles_TileMatrix_",z,"_",chr), logFile = logFile)
      .logThis(length(cellNames), paste0("NCells_TileMatrix_",z,"_",chr), logFile = logFile)
      .logThis(trunc(min(start(fragments)) / tileSize) + 1, paste0("MinTile_TileMatrix_",z,"_",chr), logFile = logFile)
      .logThis(trunc(max(end(fragments)) / tileSize) + 1, paste0("MaxTile_TileMatrix_",z,"_",chr), logFile = logFile)
      .logThis(min(matchID), paste0("MinCell_TileMatrix_",z,"_",chr), logFile = logFile)
      .logThis(max(matchID), paste0("MaxCell_TileMatrix_",z,"_",chr), logFile = logFile)

      #Check Fragments for validity in case
      nf1 <- length(fragments)

      #Check 1
      fragmentsBad1 <- fragments[!(start(fragments) >= 1)]
      fragments <- fragments[start(fragments) >= 1]

      #Check 2
      fragmentsBad2 <- fragments[!(end(fragments) <= chromLengths[z])]
      fragments <- fragments[end(fragments) <= chromLengths[z]]

      #Check N
      nf2 <- length(fragments)
      if(nf2 < nf1){
        warning("Skipping over fragments not within chromosome range on Chr:", chr)
        .logThis(fragmentsBad1, "fragmentsBad1", logFile = logFile)
        print("Bad1 (Start not greater than 0): ")
        print(fragmentsBad1)
        print("Bad2 (End greater than chromsome length): ")
        .logThis(fragmentsBad2, "fragmentsBad2", logFile = logFile)
        print(fragmentsBad2)
      }

      #Create Sparse Matrix
      mat <- Matrix::sparseMatrix(
          i = c(trunc(start(fragments) / tileSize), trunc(end(fragments) / tileSize)) + 1,
          j = as.vector(c(matchID, matchID)),
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
        if(length(blacklist) > 0){
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
      }

      #Write sparseMatrix to Arrow File!
      o <- .addMatToArrow(
        mat = mat, 
        ArrowFile = ArrowFile, 
        Group = paste0("TileMatrix/", chr), 
        binarize = binarize,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = TRUE
        )

      gc()

      0

    }, error = function(e){

      errorList <- list(
        ArrowFile = ArrowFile,
        chromLengths = chromLengths,
        blacklist = blacklist,
        chr = chromLengths[z],
        mat = if(exists("mat", inherits = FALSE)) mat else "mat"
      )

      .logError(e, fn = ".addTileMat", info = sampleName, errorList = errorList, logFile = logFile)

    })

  }

  return(ArrowFile)

}





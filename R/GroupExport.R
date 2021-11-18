#' @include HelperUtils.R
NULL

#' @export
exportGroupSE <- function(...){
    .Deprecated("getGroupSE")
    getGroupSE(...)
}

#' Export Group Summarized Experiment 
#' 
#' This function will group, summarize and export a summarized experiment for a assay in a ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the matrix in the ArrowFiles. See getAvailableMatrices to see options
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param divideN A boolean describing whether to divide by the number of cells.
#' @param scaleTo Depth normalize to this value if not NULL.
#' @param threads An integer specifying the number of threads for parallel.
#' @param verbose A boolean specifying to print messages during computation.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getGroupSE <- function(
  ArchRProj = NULL,
  useMatrix = NULL,
  groupBy = "Sample",
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = divideN, name = "divideN", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getGroupSE Input-Parameters", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)
  featureDF <- .getFeatureDF(ArrowFiles, subGroup = useMatrix)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)
  if(!.isDiscrete(Groups)){
    .logStop("groupBy must be a discrete variable!", logFile = logFile)
  }
  Cells <- ArchRProj$cellNames

  .logMessage("Getting Group Matrix", logFile=logFile)
  groupMat <- tryCatch({
    .getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = split(Cells, Groups),
      useIndex = FALSE,
      verbose = verbose
    )
  }, error = function(e){
    errorList <- list(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = split(Cells, Groups),
      useIndex = FALSE,
      verbose = verbose
    )
    .logError(e, fn = ".getGroupMatrix", info = "", errorList = errorList, logFile = logFile)
  })

  if(divideN){
    .logMessage("Normalizing by number of Cells", logFile=logFile)
    nCells <- table(Groups)[colnames(groupMat)]
    groupMat <- t(t(groupMat) / as.vector(nCells))
  }

  #Normalize
  if(!is.null(scaleTo)){
    .logMessage("Depth Normalizing", logFile=logFile)
    groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
  }

  assayList <- SimpleList(
    matrix = groupMat
  )
  names(assayList) <- useMatrix
  rm(groupMat)
  .logThis(assayList, "assayList", logFile = logFile)

  ccd <- getCellColData(ArchRProj)
  .logThis(ccd, "ccd", logFile = logFile)

  idx <- lapply(seq_len(ncol(ccd)), function(x){
    is.numeric(ccd[,x])
  }) %>% unlist %>% which

  if(length(idx) > 0){
    cD <- lapply(seq_along(idx), function(x){
      unlist(lapply(split(ccd[Cells, idx[x]], Groups), median))[colnames(assayList[[1]])]
    }) %>% Reduce("cbind", .) %>% DataFrame
    colnames(cD) <- colnames(ccd)[idx]
    rownames(cD) <- colnames(assayList[[1]])
    cD$nCells <- as.vector(table(Groups)[colnames(assayList[[1]])])
  }else{
    cD <- DataFrame(
      rownames = colnames(assayList[[1]]), 
      nCells = as.vector(table(Groups)[colnames(assayList[[1]])])
    )
  }

  .logThis(cD, "cD", logFile = logFile)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assayList,
    colData = cD,
    rowData = featureDF
  )

  .logThis(se, "se", logFile = logFile)

  .endLogging(logFile = logFile)

  se

}

#' Export Group BigWigs
#' 
#' This function will group, summarize and export a bigwig for each group in an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality. Accepted values are
#' "None", "ReadsInTSS", "nCells", "ReadsInPromoter", or "nFrags".
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param maxCells Maximum number of cells used for each bigwig.
#' @param ceiling Maximum contribution of accessibility per cell in each tile.
#' @param verbose A boolean specifying to print messages during computation.
#' @param threads An integer specifying the number of threads for parallel.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getGroupBW <- function(
  ArchRProj = NULL,
  groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "useMatrix", valid = c("character"))
  .validInput(input = normMethod, name = "groupBy", valid = c("character"))
  .validInput(input = tileSize, name = "divideN", valid = c("integer"))
  .validInput(input = maxCells, name = "scaleTo", valid = c("integer", "null"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()

  normMethods <- c("None", "ReadsInTSS", "nCells", "ReadsInPromoter", "nFrags")
  
  if(tolower(normMethod) %ni% tolower(normMethods)){
    stop(paste0("normMethod (",normMethod,") not in supported normMethods : ", paste0(normMethods, collapse=", ")))
  }

  .startLogging(logFile = logFile)
  .logThis(normMethod, "normMethod", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)

  if(tolower(normMethod) %in% tolower(c("ReadsInTSS", "ReadsInPromoter", "nFrags"))){
    normBy <- getCellColData(ArchRProj = ArchRProj, select = normMethod)
  }else{
    normBy <- NULL
  }

  if(!.isDiscrete(Groups)){
    stop("groupBy must be a discrete variable!")
  }

  Cells <- ArchRProj$cellNames

  cellGroups <- split(Cells, Groups)

  if(!is.null(maxCells)){
  
    gnames <- names(cellGroups)
  
    cellGroups <- lapply(seq_along(cellGroups), function(x){
      if(length(cellGroups[[x]]) > maxCells){
        sample(cellGroups[[x]], maxCells)
      }else{
        cellGroups[[x]]
      }
    })
  
    names(cellGroups) <- gnames
  
  }

  bwDir1 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs")
  bwDir2 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs", groupBy)

  dir.create(bwDir1, showWarnings = FALSE)
  dir.create(bwDir2, showWarnings = FALSE)

  o <- suppressWarnings(file.remove(list.files(bwDir2, full.names = TRUE)))

  #cellsInArrow <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  cellsInArrow <- split(
    rownames(getCellColData(ArchRProj)), 
    stringr::str_split(rownames(getCellColData(ArchRProj)), pattern="\\#", simplify=TRUE)[,1]
  )
  availableChr <- .availableSeqnames(head(getArrowFiles(ArchRProj)))
  chromLengths <- getChromLengths(ArchRProj)
  chromSizes <- getChromSizes(ArchRProj)
  tiles <- unlist(slidingWindows(chromSizes, width = tileSize, step = tileSize))

  if(threads > 1){
    h5disableFileLocking()
  }

  covFiles <- c()

  for(i in seq_along(cellGroups)){

    o <- .createGroupBW(
      i = i, 
      cellGroups = cellGroups,
      ArrowFiles = getArrowFiles(ArchRProj), 
      cellsInArrow = cellsInArrow, 
      availableChr = availableChr,
      chromLengths = chromLengths,
      normMethod = normMethod,
      normBy = normBy,
      ceiling = ceiling,
      tiles = tiles, 
      tileSize = tileSize,
      bwDir = bwDir2,
      tstart = tstart, 
      verbose = verbose,
      logFile = logFile,
      threads = threads
    )

    covFiles <- c(covFiles, o)

  }

  if(threads > 1){
    h5enableFileLocking()
  }

  .endLogging(logFile = logFile)

  covFiles

}

.createGroupBW <- function(
  i = NULL, 
  cellGroups = NULL,
  ArrowFiles = NULL, 
  cellsInArrow = NULL, 
  availableChr = NULL,
  chromLengths = NULL, 
  tiles = NULL,
  ceiling = NULL,
  tileSize = 100,
  normMethod = NULL,
  normBy = NULL,
  bwDir = "bigwigs",
  tstart = NULL, 
  verbose = TRUE,
  logFile = NULL,
  threads = 1
  ){

  .logDiffTime(sprintf("%s (%s of %s) : Creating BigWig for Group", names(cellGroups)[i], i, length(cellGroups)), tstart, logFile = logFile, verbose = verbose)

  #Cells
  cellGroupi <- cellGroups[[i]]
  #print(sum(normBy[cellGroupi, 1]))

  #Bigwig File!
  covFile <- file.path(bwDir, paste0(make.names(names(cellGroups)[i]), "-TileSize-",tileSize,"-normMethod-",normMethod,"-ArchR.bw"))
  rmf <- .suppressAll(file.remove(covFile))

  covList <- .safelapply(seq_along(availableChr), function(k){

    it <- 0

    for(j in seq_along(ArrowFiles)){
      cellsInI <- sum(cellsInArrow[[names(ArrowFiles)[j]]] %in% cellGroupi)
      if(cellsInI > 0){
        it <- it + 1
        if(it == 1){
          fragik <- .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi)
        }else{
          fragik <- c(fragik, .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi))
        }
      }
    }

    tilesk <- tiles[BiocGenerics::which(seqnames(tiles) %bcin% availableChr[k])]

    if(length(fragik) == 0){

      tilesk$reads <- 0

    }else{

      #N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      if (nTiles%%1 != 0) {
          nTiles <- trunc(nTiles) + 1
      }

      #Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$RG, cellGroupi)
      
      mat <- Matrix::sparseMatrix(
          i = c(trunc(start(fragik) / tileSize), trunc(end(fragik) / tileSize)) + 1,
          j = as.vector(c(matchID, matchID)),
          x = rep(1,  2*length(fragik)),
          dims = c(nTiles, length(cellGroupi))
        )

      if(!is.null(ceiling)){
        mat@x[mat@x > ceiling] <- ceiling
      }
      
      mat <- Matrix::rowSums(mat)
      
      rm(fragik, matchID)
         
      tilesk$reads <- mat

      if(tolower(normMethod) %in% c("readsintss", "readsinpromoter", "nfrags")){
        tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
      }else if(tolower(normMethod) %in% c("ncells")){
        tilesk$reads <- tilesk$reads / length(cellGroupi)
      }else if(tolower(normMethod) %in% c("none")){
      }else{
        stop("NormMethod not recognized!")
      }

    }

    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]

    tilesk

  }, threads = threads)

  names(covList) <- availableChr

  covList <- as(covList, "RleList")

  rtracklayer::export.bw(object = covList, con = covFile)

  return(covFile)

}



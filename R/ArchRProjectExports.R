#' Export Group Summarized Experiment JJJ
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
#' @export
exportGroupSE <- function(
  ArchRProj = NULL,
  useMatrix = NULL,
  groupBy = "Sample",
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE
  ){

  ArrowFiles <- getArrowFiles(ArchRProj)
  featureDF <- .getFeatureDF(ArrowFiles, subGroup = useMatrix)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)
  if(!.isDiscrete(Groups)){
    stop("groupBy must be a discrete variable!")
  }
  Cells <- ArchRProj$cellNames

    groupMat <- .getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = split(Cells, Groups),
      useIndex = FALSE,
      verbose = verbose
    )

    if(divideN){
      nCells <- table(Groups)[colnames(groupMat)]
      groupMat <- t(t(groupMat) / as.vector(nCells))
    }

    #Normalize
    if(!is.null(scaleTo)){
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }

  assayList <- SimpleList(
    matrix = groupMat
  )
  names(assayList) <- useMatrix
  rm(groupMat)

  ccd <- getCellColData(ArchRProj)
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

    SummarizedExperiment::SummarizedExperiment(
      assays = assayList,
      colData = cD,
      rowData = featureDF
    )

}

#' Export Group Summarized Experiment JJJ
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
#' @export
exportGroupBW <- function(
  ArchRProj = NULL,
  maxCells = 1000,
  groupBy = "Sample",
  normBy = "ReadsInTSS",
  ceiling = 2,
  threads = getArchRThreads(),
  verbose = TRUE
  ){

  stop("not functional yet!")

  normMethods <- c("None", "ReadsInTSS", "nCells")
  if(normBy %ni% c("none", "readsintss")){
    stop(paste0("normBy not in normMethods : ", paste0(normMethods, collapse=", ")))
  }
  ArrowFiles <- getArrowFiles(ArchRProj)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)
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

}

.createGroupBW <- function(
  i = NULL, 
  cellGroups = NULL,
  ArrowFiles = NULL, 
  cellsInArrow = NULL, 
  availableChr = NULL,
  chromLengths = NULL, 
  tileSize = 100,
  bwDir = "bigwigs",
  tstart = NULL, 
  verboseHeader = TRUE,
  verboseAll = FALSE
  ){

  .messageDiffTime(sprintf("Creating Group BW %s of %s", i, length(cellGroups)), tstart, verbose = verboseHeader)

  #Cells
  cellGroupi <- cellGroups[[i]]
  
  #Bigwig File!
  covFile <- file.path(bwDir, paste0(names(cellGroups)[i], "-ArchR-BW-TileSize-",tileSize,".bw"))
  rmf <- .suppressAll(file.remove(covFile))

  covList <- lapply(seq_along(availableChr), function(k)){

    message(k)

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

    fragik <- coverage(
      resize(IRanges(start = c( start(fragik), end(fragik) ), width = 1), tileSize, "center"), 
      width = chromLengths[availableChr[k]]
    )

    fragik

  }

  return(out)

}




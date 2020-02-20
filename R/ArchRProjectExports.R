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









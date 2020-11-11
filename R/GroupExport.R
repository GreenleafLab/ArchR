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

# #' Export Group Summarized Experiment QQQ
# #' 
# #' This function will group, summarize and export a summarized experiment for a assay in a ArchRProject.
# #'
# #' @param ArchRProj An `ArchRProject` object.
# #' @param useMatrix The name of the matrix in the ArrowFiles. See getAvailableMatrices to see options
# #' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
# #' @param divideN A boolean describing whether to divide by the number of cells.
# #' @param scaleTo Depth normalize to this value if not NULL.
# #' @param threads An integer specifying the number of threads for parallel.
# #' @param verbose A boolean specifying to print messages during computation.
# #' @export
# exportGroupBW <- function(
#   ArchRProj = NULL,
#   maxCells = 1000,
#   groupBy = "Sample",
#   normBy = "ReadsInTSS",
#   ceiling = 2,
#   threads = getArchRThreads(),
#   verbose = TRUE
#   ){

#   stop("not functional yet!")

#   normMethods <- c("None", "ReadsInTSS", "nCells")
#   if(normBy %ni% c("none", "readsintss")){
#     stop(paste0("normBy not in normMethods : ", paste0(normMethods, collapse=", ")))
#   }
#   ArrowFiles <- getArrowFiles(ArchRProj)
#   Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)
#   if(!.isDiscrete(Groups)){
#     stop("groupBy must be a discrete variable!")
#   }
#   Cells <- ArchRProj$cellNames

#   cellGroups <- split(Cells, Groups)

#   if(!is.null(maxCells)){
#     gnames <- names(cellGroups)
#     cellGroups <- lapply(seq_along(cellGroups), function(x){
#       if(length(cellGroups[[x]]) > maxCells){
#         sample(cellGroups[[x]], maxCells)
#       }else{
#         cellGroups[[x]]
#       }
#     })
#     names(cellGroups) <- gnames
#   }

# }

# .createGroupBW <- function(
#   i = NULL, 
#   cellGroups = NULL,
#   ArrowFiles = NULL, 
#   cellsInArrow = NULL, 
#   availableChr = NULL,
#   chromLengths = NULL, 
#   tileSize = 100,
#   bwDir = "bigwigs",
#   tstart = NULL, 
#   verboseHeader = TRUE,
#   verboseAll = FALSE
#   ){

#   .messageDiffTime(sprintf("Creating Group BW %s of %s", i, length(cellGroups)), tstart, verbose = verboseHeader)

#   #Cells
#   cellGroupi <- cellGroups[[i]]
  
#   #Bigwig File!
#   covFile <- file.path(bwDir, paste0(names(cellGroups)[i], "-ArchR-BW-TileSize-",tileSize,".bw"))
#   rmf <- .suppressAll(file.remove(covFile))

#   covList <- lapply(seq_along(availableChr), function(k)){

#     message(k)

#     it <- 0
#     for(j in seq_along(ArrowFiles)){
#       cellsInI <- sum(cellsInArrow[[names(ArrowFiles)[j]]] %in% cellGroupi)
#       if(cellsInI > 0){
#         it <- it + 1
#         if(it == 1){
#           fragik <- .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi)
#         }else{
#           fragik <- c(fragik, .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi))
#         }
#       }
#     }

#     fragik <- coverage(
#       resize(IRanges(start = c( start(fragik), end(fragik) ), width = 1), tileSize, "center"), 
#       width = chromLengths[availableChr[k]]
#     )

#     fragik

#   }

#   return(out)

# }

#' Export to Monocle3
#' 
#' This function will return a monocle3 cell_data_set for a assay in a ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the matrix in the ArrowFiles. See getAvailableMatrices to see options
#' @param threads An integer specifying the number of threads for parallel.
#' @param verbose A boolean specifying to print messages during computation.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @import monocle3
#' @import SingleCellExperiment
#' @export

exportMonocle3 <- function(
  ArchRProj = NULL,
  useMatrix = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  binarize = T,
  logFile = createLogFile("exportMonocle3")
){
  require(SingleCellExperiment)
  require(monocle3)
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportMonocle3 Input-Parameters", logFile = logFile)
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- sapply(ArrowFiles, ArchR:::.validArrow())
  featureDF <- ArchR:::.getFeatureDF(ArrowFiles, subGroup = useMatrix)
  Groups <- getCellColData(ArchRProj = ArchRProj)
  Cells <- ArchRProj$cellNames
  if(is.null(featureDF$end)){
    starts<-featureDF[as.character(featureDF$seqnames) %in% as.character(featureDF$seqnames)[1],]$start
    width<-getmode(starts[-1]-starts[-length(starts)])
    ranges<-GRanges(seqnames = featureDF$seqnames, ranges = IRanges(start = featureDF$start, width = width), names=featureDF$idx)
  }else{
    ranges<-GRanges(seqnames = featureDF$seqnames, ranges = IRanges(start = featureDF$start, end = featureDF$end), names=featureDF$idx)
  }
  ArchR:::.logMessage("Getting Group Matrix", logFile=logFile)
  mat <- tryCatch({
    getMatrixFromArrow(ArrowFiles, useMatrix = useMatrix, binarize = binarize, cellNames = Cells, useSeqnames = as.character(seqnames(ranges)))
  }, error = function(e){
    errorList <- list(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      verbose = verbose
    )
    ArchR:::.logError(e, fn = ".getMatrixFromArrow", info = "", errorList = errorList, logFile = logFile)
  })
  o <- h5read(file = ArrowFiles, name = paste0("/", useMatrix,"/Info/Params"))
  h5closeAll()
  window_size<-getmode(o$tileSize)
  mat@assays@data$counts<-mat@assays@data[[useMatrix]]
  mat@assays@data[[useMatrix]]<-NULL
  rs<-ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix)
  rowRanges(mat)<-ranges
  mat<-mat[, match(rownames(ArchRProj@reducedDims$IterativeLSI$matSVD), colnames(mat))]
  cds<-new_cell_data_set(expression_data = mat@assays@data$counts, cell_metadata = colData(mat))
  rowRanges(cds)<-rowRanges(mat)
  cds<-cds[,rownames(ArchRProj@reducedDims$IterativeLSI$matSVD)]
  reducedDims(cds)<-SimpleList(LSI=ArchRProj@reducedDims$IterativeLSI$matSVD, 
                               UMAP=ArchRProj@embeddings$UMAP$df)
  irlba_rotation = ArchRProj@reducedDims$IterativeLSI$svd$u
  row.names(irlba_rotation) = paste0(ArchRProj@reducedDims$IterativeLSI$LSIFeatures$seqnames, "_", ArchRProj@reducedDims$IterativeLSI$LSIFeatures$idx)
  iLSI<-SimpleList(svd=ArchRProj@reducedDims$IterativeLSI$svd,
                   features=ArchRProj@reducedDims$IterativeLSI$LSIFeatures, 
                   row_sums = ArchRProj@reducedDims$IterativeLSI$rowSm,
                   seed=ArchRProj@reducedDims$IterativeLSI$seed,
                   binarize=ArchRProj@reducedDims$IterativeLSI$binarize,
                   scale_to=ArchRProj@reducedDims$IterativeLSI$scaleTo,
                   num_dim=ArchRProj@reducedDims$IterativeLSI$nDimensions, 
                   resolution=NULL, 
                   granges=ArchRProj@reducedDims$IterativeLSI$LSIFeatures, 
                   LSI_method=ArchRProj@reducedDims$IterativeLSI$LSIMethod, outliers=NULL)
  pp_aux <- SimpleList(iLSI=iLSI, gene_loadings=irlba_rotation)
  cds@preprocess_aux <- pp_aux
  if(is.null(cds@preprocess_aux$iLSI$granges$end)){
    starts<-cds@preprocess_aux$iLSI$granges[as.character(cds@preprocess_aux$iLSI$granges$seqnames) %in% as.character(cds@preprocess_aux$iLSI$granges$seqnames)[1],]$start
    width<-getmode(starts[-1]-starts[-length(starts)])
    ranges<-GRanges(seqnames = cds@preprocess_aux$iLSI$granges$seqnames, ranges = IRanges(start = cds@preprocess_aux$iLSI$granges$start, width = width), names=cds@preprocess_aux$iLSI$granges$idx)
  }else{
    ranges<-GRanges(seqnames = cds@preprocess_aux$iLSI$granges$seqnames, ranges = IRanges(start = cds@preprocess_aux$iLSI$granges$start, end = cds@preprocess_aux$iLSI$granges$end), names=cds@preprocess_aux$iLSI$granges$idx)
  }
  cds@preprocess_aux$iLSI$granges<-ranges
  cds@clusters[["UMAP"]]$clusters[colnames(exprs(cds))]<-as.character(ArchRProj@cellColData[colnames(exprs(cds)),]$Clusters)
  cds@reduce_dim_aux<-SimpleList(UMAP=SimpleList(scale_info=NULL, model_file=ArchRProj@embeddings$UMAP$params$uwotModel, num_dim=cds@preprocess_aux$iLSI$num_dim))
  cds
  
}

#' @export
getmode<-function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


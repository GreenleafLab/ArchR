####################################################################
# Trajectory Analysis Methods
####################################################################

#' Add a Supervised Trajectory to an ArchR Project
#' 
#' This function will fit a supervised trajectory in a lower dimensional space that 
#' can then be used for downstream analyses.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory to be added in `cellColData`.
#' @param trajectory The order of cell groups to be used for constraining the initial supervised fitting procedure.
#' For example, to get a trajectory from Cluster1 to Cluster2 to Cluster3, input should be c("Cluster1", "Cluster2", "Cluster3").
#' Cells will then be used from these 3 groups to constrain an initial fit in the group order.
#' @param groupBy A string indicating the column name from `cellColData` that contains the cell group definitions used in
#' `trajectory` to constrain the initial supervised fitting procedure.
#' @param reducedDims A string indicating the name of the `reducedDims` object from the `ArchRProject` that should be used for distance computation.
#' @param embedding A string indicating the name of the `embedding` object from the `ArchRProject` that should be used for distance computation.
#' @param preFilterQuantile Prior to the initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping
#' center is above the provided quantile will be excluded.
#' @param postFilterQuantile After initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping center
#' is above the provided quantile will be excluded.
#' @param useAll A boolean describing whether to use cells outside of trajectory groups for post-fitting procedure.
#' @param dof The number of degrees of freedom to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param spar The sparsity to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exists in the given `ArchRProject`.
#' @param seed A number to be used as the seed for random number generation for trajectory creation.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  trajectory = NULL, 
  groupBy = "Clusters",
  reducedDims = "IterativeLSI",
  embedding = NULL,
  preFilterQuantile = 0.9, 
  postFilterQuantile = 0.9,
  useAll = FALSE, 
  dof = 250,
  spar = 1,
  force = FALSE,
  seed = 1,
  logFile = createLogFile("addTrajectory")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = trajectory, name = "trajectory", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = embedding, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = preFilterQuantile, name = "preFilterQuantile", valid = c("numeric"))
  .validInput(input = postFilterQuantile, name = "postFilterQuantile", valid = c("numeric"))
  .validInput(input = useAll, name = "useAll", valid = c("boolean"))
  .validInput(input = dof, name = "dof", valid = c("integer"))
  .validInput(input = spar, name = "spar", valid = c("numeric"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  if(!is.null(seed)) set.seed(seed)

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addTrajectory Input-Parameters", logFile=logFile)
  
  groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
  groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]

  if(sum(unique(groupDF[,1]) %in% trajectory)==0){
    .logStop("trajectory does not span any groups in groupBy! Are you sure your input is correct?", logFile = logFile)
  }

  if(sum(unique(groupDF[,1]) %in% trajectory) < 3){
    if(!force){
      .logStop("trajectory must span at least 3 groups in groupBy!", logFile = logFile)
    }
  }

  if(is.null(embedding)){
    mat <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims)
  }else{
    mat <- getEmbedding(ArchRProj = ArchRProj, embedding = embedding)
  }
  mat <- mat[rownames(groupDF),,drop = FALSE]

  ######################################################
  #Filter Outliers
  ######################################################
  .logMessage("Filtering outliers", logFile = logFile)
  filterObj <- lapply(seq_along(trajectory), function(x){
      
      #Subset
      groupsx <- rownames(groupDF)[groupDF[,1]==trajectory[x]]
      matx <- mat[groupsx,,drop = FALSE]

      #Filter Distance
      matMeanx <- colMeans(matx)
      diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
      idxKeep <- which(diffx <= quantile(diffx, preFilterQuantile))
      
      #Filter
      list(mat = matx[idxKeep,,drop=FALSE], groups = groupsx[idxKeep])

  })

  matFilter <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat) %>% Reduce("rbind", .)
  groupsFilter <- groupDF[lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups) %>% Reduce("c", .),,drop=FALSE]

  ######################################################
  #Now Initial Alignment
  ######################################################
  .logMessage("Initial Alignment Before Spline Fit", logFile = logFile)
  initialTime <- lapply(seq_along(trajectory), function(x){
      
      groupsx <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x]]
      matx <- matFilter[groupsx,,drop = FALSE]
      
      #Get Differences
      if(x != length(trajectory)){
          groupsxp1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x + 1]]
          meanx <- colMeans(matFilter[groupsxp1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- (1 - .getQuantiles(diffx)) + x
      }else{
          groupsxm1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x - 1]]
          meanx <- colMeans(matFilter[groupsxm1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- .getQuantiles(diffx) + x
      }
      
      timex

  }) %>% unlist

  ######################################################
  #Fit Cubic Splines
  ######################################################
  .logMessage("Spline Fit", logFile = logFile)
  matSpline <- lapply(seq_len(ncol(matFilter)), function(x){
    tryCatch({
      stats::smooth.spline(
          x = initialTime, 
          y = matFilter[names(initialTime), x], 
          df = dof, 
          spar = spar
      )[[2]]
    }, error = function(e){
      errorList <- list(
        it = x,
        x = initialTime, 
        y = matFilter[names(initialTime), x], 
        df = dof, 
        spar = spar
      )
      .logError(e, fn = "smooth.spline", info = "", errorList = errorList, logFile = logFile)      
    })
  }) %>% Reduce("cbind",.) %>% data.frame()

  ######################################################
  # 1. KNN Fit vs Actual
  ######################################################
  .logMessage("KNN to Spline", logFile = logFile)
  knnObj <- nabor::knn(
      data =  matSpline,
      query = mat, 
      k = 3
  )

  #Estimate place along trajectory
  knnIdx <- knnObj[[1]]
  knnDist <- knnObj[[2]]
  knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
  knnDistQ <- .getQuantiles(knnDist[,1])

  #Filter Outlier Cells to Trajectory for High Resolution
  idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], postFilterQuantile))
  dfTrajectory <- DataFrame(
      row.names = rownames(mat),
      Distance = knnDist[, 1],
      DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
  )[idxKeep, , drop = FALSE]

  ######################################################
  # 2. Fit cells not in trajectory clusters
  ######################################################
  if(useAll){
    .logMessage("Aligning cells not in trajectory", logFile = logFile)
    
    if(is.null(embedding)){
      mat2 <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims)
    }else{
      mat2 <- getEmbedding(ArchRProj = ArchRProj, embedding = embedding)
    }
    
    groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
    groupDF <- groupDF[groupDF[,1] %ni% trajectory,,drop=FALSE]
    mat2 <- mat2[rownames(groupDF),,drop = FALSE]

    #Nearest Neighbors
    knnObj2 <- nabor::knn(
        data =  matSpline,
        query = mat2, 
        k = 3
    )

    #Estimate place along trajectory
    knnIdx2 <- knnObj2[[1]]
    knnDist2 <- knnObj2[[2]]
    knnDiff2 <- ifelse(knnIdx2[,2] > knnIdx2[,3], 1, -1)
    knnDistQ2 <- .getQuantiles(knnDist2[,1])

    #Keep Cells that are within the maximum distance of a cluster
    idxKeep <- which(knnDist2[,1] < max(dfTrajectory[,1]))
    dfTrajectory2 <- DataFrame(
        row.names = rownames(mat2),
        Distance = knnDist2[, 1],
        DistanceIdx = knnIdx2[, 1] + knnDiff2 * knnDistQ2
    )[idxKeep, , drop = FALSE]

    #Final Output
    dfTrajectory3 <- rbind(dfTrajectory, dfTrajectory2)
  }else{
    dfTrajectory3 <- dfTrajectory
  }
  
  dfTrajectory3$Trajectory <- 100 * .getQuantiles(dfTrajectory3[,2])
  
  #Add To ArchR Project
  ArchRProj <- addCellColData(
      ArchRProj = ArchRProj,
      data = dfTrajectory3$Trajectory,
      name = name,
      cells = rownames(dfTrajectory3),
      force = force
  )

  .endLogging(logFile = logFile)

  ArchRProj

}


#' Get Supervised Trajectory from an ArchR Project
#' 
#' This function will get a supervised trajectory from an `ArchRProject` (see `addTrajectory`), get data
#' from a desired matrix, and smooth each value across the input trajectory.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory in `cellColData` to retrieve from the given `ArchRProject`.
#' @param useMatrix The name of the data matrix from the `ArrowFiles` to get numerical values for each cell from. Recommended
#' matrices are "GeneScoreMatrix", "PeakMatrix", or "MotifMatrix".
#' @param groupEvery The number of sequential percentiles to group together when generating a trajectory. This is similar to smoothing
#' via a non-overlapping sliding window across pseudo-time. If `groupEvery = 2`, the values for percentiles [1 and 2], [3 and 4],
#' [5 and 6], etc. will be grouped together.
#' @param log2Norm A boolean value that indicates whether the summarized trajectory matrix should be log2 transformed. If you are using
#' a "MotifMatrix" set to FALSE.
#' @param scaleTo Once the sequential trajectory matrix is created, each column in that matrix will be normalized to a column sum
#' indicated by `scaleTo`. Setting this to `NULL` will prevent any normalization and should be done in certain circumstances
#' (for ex. if you are using a "MotifMatrix").
#' @param smoothWindow An integer value indicating the smoothing window in size (relaive to `groupEvery`) for the sequential 
#' trajectory matrix to better reveal temporal dynamics.
#' @param threads The number of threads to be used for parallel computing.
#' @export
getTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  groupEvery = 1,
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = 11,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = groupEvery, name = "groupEvery", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric", "null"))
  .validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  trajectory <- getCellColData(ArchRProj, name)
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(0, 100, groupEvery)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }

  groupList <- lapply(seq_along(breaks), function(x){
      if(x == 1){
          NULL
      }else{
          rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
      }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])

  featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))

  message("Creating Trajectory Group Matrix..")
  groupMat <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = featureDF,
      groupList = groupList, 
      threads = threads, 
      verbose = FALSE, 
      useMatrix = useMatrix
  )

  #Scale
  if(!is.null(scaleTo)){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }else{
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }

  if(log2Norm){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }else{
      groupMat <- log2(groupMat + 1)
    }
  }

  if(!is.null(smoothWindow)){
    
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) .centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          smoothMat = as.matrix(smoothGroupMat), 
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }else{

    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }

  metadata(seTrajectory)$Params <- list(
    useMatrix = useMatrix, 
    matrixClass = matrixClass,
    scaleTo = scaleTo, 
    log2Norm = log2Norm, 
    smoothWindow = smoothWindow, 
    date = Sys.Date()
  )

  seTrajectory

}

#' @export
trajectoryHeatmap <- function(...){
    .Deprecated("plotTrajectoryHeatmap")
    plotTrajectoryHeatmap(...)
}

#' Plot a Heatmap of Features across a Trajectory
#' 
#' This function will plot a heatmap of the results from getTrajectory
#' 
#' @param seTrajectory A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param varCutOff The "Variance Quantile Cutoff" to be used for identifying the top variable features across the given trajectory.
#' Only features with a variance above the provided quantile will be retained.
#' @param maxFeatures The maximum number of features, ordered by variance, to consider from `useMatrix` when generating a trajectory.
#' This prevents smoothing a large number number of features which can be very time consuming.
#' @param scaleRows A boolean value that indicates whether row-wise z-scores should be computed on the matrix provided by `seTrajectory`.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies
#' rownames from `seTrajectory` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the `rownames` of `seTrajectory` that should be labeled on the side of the heatmap.
#' @param labelTop A number indicating how many of the top N features, based on variance, in `seTrajectory` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param rowOrder If wanting to set the order of rows to be plotted, the indices (integer or character correpsonding 
#' to rownmaes) can be provided here.
#' @param useSeqnames A character vector that indicates which `seqnames` should be plotted in the heatmap. Features from
#' `seqnames` that are not listed will be ignored. In the context of a `Sparse.Assays.Matrix`, such as a matrix containing chromVAR
#' deviations, the `seqnames` do not correspond to chromosomes, rather they correspond to the sub-portions of the matrix, for example
#' raw deviations ("deviations") or deviation z-scores ("z") for a chromVAR deviations matrix.
#' @param returnMatrix A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param force If useSeqnames is longer than 1 if matrixClass is "Sparse.Assays.Matrix" to continue. This is not recommended because these matrices
#' can be in different units.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotTrajectoryHeatmap <- function(
  seTrajectory = NULL,
  varCutOff = 0.9,
  maxFeatures = 25000,
  scaleRows = TRUE,
  limits = c(-1.5, 1.5),
  grepExclude = NULL,
  pal = NULL,
  labelMarkers = NULL,
  labelTop = 50,
  labelRows = FALSE,
  rowOrder = NULL, 
  useSeqnames = NULL,
  returnMatrix = FALSE,
  force = FALSE,
  logFile = createLogFile("plotTrajectoryHeatmap")
  ){

  .validInput(input = seTrajectory, name = "seTrajectory", valid = c("SummarizedExperiment"))
  .validInput(input = varCutOff, name = "varCutOff", valid = c("numeric", "null"))
  .validInput(input = maxFeatures, name = "maxFeatures", valid = c("integer", "null"))
  .validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
  .validInput(input = limits, name = "limits", valid = c("numeric"))
  .validInput(input = grepExclude, name = "grepExclude", valid = c("character", "null"))
  .validInput(input = pal, name = "pal", valid = c("palette", "null"))
  .validInput(input = labelMarkers, name = "labelMarkers", valid = c("character", "null"))
  .validInput(input = labelTop, name = "labelTop", valid = c("integer"))
  .validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
  .validInput(input = rowOrder, name = "rowOrder", valid = c("vector", "null"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character", "null"))
  .validInput(input = returnMatrix, name = "returnMatrix", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotTrajectoryHeatmap Input-Parameters", logFile = logFile)

  if(metadata(seTrajectory)$Params$matrixClass == "Sparse.Assays.Matrix"){
    if(is.null(useSeqnames) || length(useSeqnames) > 1){
      .logMessage("useSeqnames is NULL or greater than 1 with a Sparse.Assays.Matrix trajectory input.", verbose = TRUE, logFile = logFile)
      if(force){
        .logMessage("force=TRUE thus continuing", verbose = verbose, logFile = logFile)
      }else{
        useSeqnames <- rev(unique(rowData(seTrajectory)$seqnames))[1]
        .logMessage(paste0("force=FALSE thus continuing with subsetting useSeqnames = ", useSeqnames) , verbose = TRUE, logFile = logFile)
      }
    }
  }

  if(!is.null(useSeqnames)){
    seTrajectory <- seTrajectory[paste0(rowData(seTrajectory)$seqnames) %in% paste0(useSeqnames), ]
  }

  if(nrow(seTrajectory) == 0){
    .logStop("No features left in seTrajectory, please check input!", logFile = logFile)
  }

  mat <- assay(seTrajectory)

  if(!is.null(grepExclude)){
    idxExclude <- grep(grepExclude, rownames(mat))
    if(length(idxExclude) > 0){
      mat <- mat[-grep(grepExclude, rownames(mat)), , drop = FALSE]
    }
  }
  
  #Rows with NA
  rSNA <- rowSums(is.na(mat))
  if(sum(rSNA > 0) > 0){
    .logMessage("Removing rows with NA values...", verbose = TRUE, logFile = logFile)
    mat <- mat[rSNA == 0, ]#Remove NA Rows
  }
  .logThis(mat, "mat-pre", logFile = logFile)
  varQ <- .getQuantiles(matrixStats::rowVars(mat))
  .logThis(varQ, "varQ", logFile = logFile)
  orderedVar <- FALSE
  if(is.null(rowOrder)){
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if(is.null(varCutOff) & is.null(maxFeatures)){
      n <- nrow(mat)
    }else if(is.null(varCutOff)){
      n <- maxFeatures
    }else if(is.null(maxFeatures)){
      n <- (1-varCutOff) * nrow(mat)
    }else{
      n <- min((1-varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n ),]
  }
  .logThis(mat, "mat-post", logFile = logFile)

  #rownames(mat) <- rowData(seTrajectory)$name
  
  if(!is.null(labelTop)){
    if(orderedVar){
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    }else{
      idxLabel <- rownames(mat)[order(varQ,decreasing=TRUE)][seq_len(labelTop)]
    }
  }else{
    idxLabel <- NULL
  }
  .logThis(idxLabel, "idxLabel", logFile = logFile)

  if(!is.null(labelMarkers)){
    idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
    idxLabel2 <- idxLabel2[idxLabel2 > 0]
  }else{
    idxLabel2 <- NULL
  }
  .logThis(idxLabel2, "idxLabel2", logFile = logFile)

  idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])
  .logThis(idxLabel, "idxLabel", logFile = logFile)

  if(scaleRows){
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    .logThis(mat, "mat-zscores", logFile = logFile)
  }

  if(nrow(mat) == 0){
    stop("No Features Remaining!")
  }

  if(is.null(pal)){
    if(is.null(metadata(seTrajectory)$Params$useMatrix)){
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }else if(tolower(metadata(seTrajectory)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blueYellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }
  }

  if(!is.null(rowOrder)){
    idx <- rowOrder
  }else{
    idx <- order(apply(mat, 1, which.max))
  }
  .logThis(idx, "idx", logFile = logFile)

  ht <- tryCatch({
    
    .ArchRHeatmap(
      mat = mat[idx, ],
      scale = FALSE,
      limits = c(min(mat), max(mat)),
      color = pal, 
      clusterCols = FALSE, 
      clusterRows = FALSE,
      labelRows = labelRows,
      labelCols = FALSE,
      customRowLabel = match(idxLabel, rownames(mat[idx,])),
      showColDendrogram = TRUE,
      name = metadata(seTrajectory)$Params$useMatrix,
      draw = FALSE
    )
  
  }, error = function(e){
    
    errorList = list(
      mat = mat[idx, ],
      scale = FALSE,
      limits = c(min(mat), max(mat)),
      color = pal, 
      clusterCols = FALSE, 
      clusterRows = FALSE,
      labelRows = labelRows,
      labelCols = FALSE,
      customRowLabel = match(idxLabel, rownames(mat[idx,])),
      showColDendrogram = TRUE,
      name = metadata(seTrajectory)$Params$useMatrix,
      draw = FALSE
    )
  
    .logError(e, fn = ".ArchRHeatmap", info = "", errorList = errorList, logFile = logFile)      
  
  })

  .endLogging(logFile = logFile)

  if(returnMatrix){
    return(mat[idx, ])
  }else{
    return(ht)
  }

}

#' Visualize a Trajectory from ArchR Project
#' 
#' This function will plot a trajectory that was created onto an embedding.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the embedding to use to visualize the given `trajectory`. See `addEmbedding()` for more information.
#' @param trajectory The column name in `cellColData` that refers the trajectory to be plotted. See `addTrajectory()` for more information.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData")
#' or by a data matrix in the associated ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName/rowname of the data matrix to be used for plotting. 
#' For example if colorBy is "cellColData" then `name` refers to a column name in the cellcoldata (see `getCellcoldata()`). If `colorBy`
#' is "GeneScoreMatrix" then `name` refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values from `colorBy`.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells'
#' values. See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for coloring cells.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels,
#' just the internal portions of the plot.
#' @param quantCut If this is not `NULL`, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values. 
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the lower threshold and y is 
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the 2.5th percentile and 97.5 percentile of values and
#' set values below/above to the value of the 2.5th and 97.5th percentile values respectively.
#' @param quantHex The numeric xth quantile of all dots within each individual hexagon will determine the numerical value for
#' coloring to be displayed. This occurs when (i) `plotAs` is set to "hex" or (ii) `plotAs` is set to `NULL` and the values of `colorBy` are numeric.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a discrete color set is desired.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a continuous color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster
#' being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed arrow in the embedding based on the aligned trajectory.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted. By default
#' if `colorBy` is numeric, then `plotAs` is set to "hex".
#' @param smoothWindow An integer value indicating the smoothing window for creating inferred Arrow overlay on to embedding.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' @export
plotTrajectory <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  trajectory = "Trajectory",
  colorBy = "colData",
  name = "Trajectory",
  log2Norm = NULL,
  imputeWeights = if(!grepl("coldata",tolower(colorBy[1]))) getImputeWeights(ArchRProj),
  pal = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5,
  logFile = createLogFile("plotTrajectory"),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = embedding, name = "reducedDims", valid = c("character"))
  .validInput(input = trajectory, name = "trajectory", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = quantCut, name = "quantCut", valid = c("numeric", "null"))
  .validInput(input = quantHex, name = "quantHex", valid = c("numeric"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character", "null"))
  .validInput(input = continuousSet, name = "continuousSet", valid = c("character", "null"))
  .validInput(input = randomize, name = "randomize", valid = c("boolean"))
  .validInput(input = keepAxis, name = "keepAxis", valid = c("boolean"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = addArrow, name = "addArrow", valid = c("boolean"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character", "null"))
  .validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .requirePackage("ggplot2", source = "cran")

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)

  if(is.null(quantCut)){
    quantCut <- c(0, 1)
  }

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }
  allColorBy <-  c("colData", "cellColData", .availableArrays(getArrowFiles(ArchRProj)))
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy (",colorBy,") must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]

  ##############################
  # Plot Helpers
  ##############################
  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }

  ##############################
  # Get Trajectory
  ##############################
  dfT <- getCellColData(ArchRProj, select = trajectory)
  idxRemove <- which(is.na(dfT[,1]))
  .logThis(dfT, "dfT", logFile = logFile)

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  .logThis(df, "embedding", logFile = logFile)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("x", "y", "PseudoTime")

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    
    plotParams$color <- as.vector(getCellColData(ArchRProj, select = name, drop = FALSE)[rownames(df), 1])
    plotParams$discrete <- .isDiscrete(plotParams$color)
    plotParams$continuousSet <- "horizonExtra"
    plotParams$discreteSet <- "stallion"
    plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }

  }else{

    units <- tryCatch({
        .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
      },error=function(e){
        "values"
    })

    plotParams$continuousSet <- "solarExtra"

    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
      plotParams$continuousSet <- "horizonExtra"
    }

    if(is.null(log2Norm)){
      log2Norm <- FALSE
    }

    plotParams$color <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = FALSE)

    if(!all(rownames(df) %in% colnames(plotParams$color))){
      .logMessage("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.", logFile = logFile)
      stop("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.")
    }

    plotParams$color <- plotParams$color[, rownames(df), drop = FALSE]

    .logThis(plotParams$color, "colorMat-Before-Impute", logFile = logFile)

    plotParams$discrete <- FALSE
    plotParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }
    if(plotAs=="hexplot"){
      plotParams$fun <- .summarizeHex
    }

  }

  #Additional Params!
  plotParams$xlabel <- gsub("_", " ",stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,2])
  plotParams$ylabel <- gsub("_", " ",stringr::str_split(colnames(df)[2],pattern="#",simplify=TRUE)[,2])

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  
  if(plotParams$discrete){
    plotParams$color <- paste0(plotParams$color)
  }

  if(!plotParams$discrete){
   
    if(!is.null(imputeWeights)){
      message("Imputing Matrix")
      mat <- matrix(as.vector(plotParams$color), nrow = 1)
      colnames(mat) <- rownames(df)
      plotParams$color <- imputeMatrix(mat = mat, imputeWeights = imputeWeights, logFile = logFile)
      .logThis(plotParams$color, "colorMat-After-Impute", logFile = logFile)
    }

    plotParams$color <- as.vector(plotParams$color)

    if(name != trajectory){
      plotParams$color <- .quantileCut(plotParams$color, min(quantCut), max(quantCut))
    }

    if(!is.null(log2Norm)){
      if(log2Norm){
        plotParams$color <- log2(plotParams$color + 1)
        plotParams$colorTitle <- paste0("Log2(",units," + 1)")
      }else{
        plotParams$colorTitle <- units
      }
    }
    
    plotParams$color[idxRemove] <- NA
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)
    
    if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
      plotParams$addPoints <- TRUE
      if(is.null(plotParams$bins)){
        plotParams$bins <- 150
      }
      message("Plotting")
      .logThis(plotParams, name = "PlotParams", logFile = logFile)
      out <- do.call(ggHex, plotParams)
    }else{
      message("Plotting")
      .logThis(plotParams, name = "PlotParams", logFile = logFile)
      out <- do.call(ggPoint, plotParams)
    }
  
  }else{
    message("Plotting")
    .logThis(plotParams, name = "PlotParams", logFile = logFile)
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  .logMessage("Plotting Trajectory", logFile = logFile)

  #Prep Trajectory Vector
  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  .logThis(dfT, "TrajectoryDF", logFile = logFile)

  #Plot Pseudo-Time
  out2 <- ggPoint(
    x = dfT$PseudoTime, 
    y = dfT$value, 
    color = dfT$PseudoTime, 
    discrete = FALSE,
    xlabel = "PseudoTime", 
    ylabel = name, 
    pal = plotParams$pal,
    ratioYX = 0.5, 
    rastr = TRUE
  ) + geom_smooth(color = "black")

  attr(out2, "ratioYX") <- 0.5


  if(addArrow){

    .logMessage("Adding Inferred Arrow Trajectory to Plot", logFile = logFile)
 
    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>% 
      lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame
    dfArrow$x <- .centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- .centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + geom_path(
            data = dfArrow, aes(x, y, color=NULL), size= 1, 
            arrow = arrow(type = "open", length = unit(0.1, "inches"))
          )
  }

  .endLogging(logFile = logFile)

  list(out, out2)

}


###################################################################
# New Trajectory Analyses
#
# - Support for Monocle3 based trajectory analysis
# - Support for Slingshot based trajectory analysis
#
###################################################################

###################################################################
# Monocle3
###################################################################

#' Get a Monocle CDS of Trajectories that can be added to an ArchRProject #NEW
#' 
#' This function will use monocle3 to find trajectories and then returns a monocle CDS object that can be used as
#' input for `addMonocleTrajectory`.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column
#' in `cellColData`. This limits the groups used to identify trajectories.
#' @param principalGroup The principal group which represents the group that will be the starting point for all trajectories.
#' @param groupBy A string indicating the column name from `cellColData` that contains the cell group definitions used in
#' `useGroups` to constrain trajectory analysis.
#' @param embedding A string indicating the name of the `embedding` object from the `ArchRProject` that should be used for trajectory analysis.
#' @param clusterParams A list of parameters to be added when clustering cells for monocle3 with `monocle3::cluster_cells`.
#' @param graphParams A list of parameters to be added when learning graphs for monocle3 with `monocle3::learn_graph`.
#' @param seed A number to be used as the seed for random number generation for trajectory creation.
#' @export
getMonocleTrajectories <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useGroups = NULL,
  principalGroup = NULL,
  groupBy = NULL,
  embedding = NULL,
  clusterParams = list(),
  graphParams = list(),
  seed = 1
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character"))
  .validInput(input = principalGroup, name = "principalGroup", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = clusterParams, name = "clusterParams", valid = c("list"))
  .validInput(input = graphParams, name = "graphParams", valid = c("list"))
  .validInput(input = seed, name = "seed", valid = c("numeric"))

  .requirePackage("monocle3")

  set.seed(seed)

  message("Running Monocole3 Trajectory Infrastructure!")

  #Create CDS
  sce <- SingleCellExperiment(
    assays = SimpleList(
      counts = as(matrix(rnorm(nCells(ArchRProj) * 3), ncol = nCells(ArchRProj), nrow = 3), "dgCMatrix")
    ),
    colData = getCellColData(ArchRProj)
  )

  cds <- methods::new(
    "cell_data_set", 
    assays = SummarizedExperiment::Assays(list(counts = methods::as(assay(sce), "dgCMatrix"))), 
    colData = colData(sce), 
    int_elementMetadata = int_elementMetadata(sce), 
      int_colData = int_colData(sce), 
      int_metadata = int_metadata(sce), 
      metadata = metadata(sce), 
      NAMES = NULL, 
      elementMetadata = elementMetadata(sce)[, 0], 
      rowRanges = rowRanges(sce)
  )
  metadata(cds)$cds_version <- Biobase::package.version("monocle3")

  rm(sce)

  #Add Embedding
  message("Adding Embedding")
  reducedDims(cds)$UMAP <- getEmbedding(ArchRProj, embedding = embedding)

  if(!is.null(useGroups)){
    cds <- cds[, which(colData(cds)[, groupBy] %in% useGroups)]
  }

  #Check principalGroup
  pCells <- which(colData(cds)[, groupBy] == principalGroup)
  if(length(pCells) == 0){
    stop("No Cells in groupBy are equal to principalGroup")
  }

  #Run Clustering on Embedding and LearnGraph
  message("Clustering Embedding")
  clusterParams$cds <- cds
  cds <- do.call(monocle3::cluster_cells, clusterParams)
  rm(clusterParams)
  gc()

  message("Learning Graphs")
  graphParams$cds <- cds
  cds <- do.call(monocle3::learn_graph, graphParams)
  rm(graphParams)
  gc()

  #Get Prinicipal Node
  message("Getting Principal Node")
  closestVertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closestVertex <- as.matrix(closestVertex[colnames(cds), ])
  rootNodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closestVertex[pCells,]))))]

  #Order Cells
  message("Ordering Cells")
  cds <- order_cells(cds, root_pr_nodes = rootNodes)

  #Get Pseudotime
  cds@principal_graph_aux[[1]]$pseudotime <- ArchR:::.getQuantiles(cds@principal_graph_aux[[1]]$pseudotime) * 100

  #Plot Results
  canRaster <- requireNamespace("ggrastr", quietly = TRUE)

  p1 <- plot_cells(cds,
         color_cells_by = groupBy,
         rasterize = canRaster,
         label_groups_by_cluster=FALSE,
         label_leaves=FALSE,
         label_branch_points=FALSE) + 
    scale_colour_manual(values = paletteDiscrete(values = colData(cds)[,groupBy])) + theme_ArchR() + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p2 <- plot_cells(cds,
         color_cells_by = "pseudotime",
         label_cell_groups=FALSE,
         rasterize = canRaster,
         label_leaves=FALSE,
         label_branch_points=FALSE,
         graph_label_size=1.5) + theme_ArchR() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank())


  path <- file.path(getOutputDirectory(ArchRProj), "Monocole3", paste0("Plot-Results-", name, ".pdf"))

  message("Plotting Results - ", path)
  pdf(path, width = 6, height = 6, useDingbats = FALSE)
  ArchR:::.fixPlotSize(p1)
  ArchR:::.fixPlotSize(p2, newPage = TRUE)
  dev.off()

  cds

}

#' Add a Monocle Trajectory to an ArchR Project #NEW
#' 
#' This function will add a trajectory from a monocle CDS created from `getMonocleTrajectories` to an
#' ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory to be added in `cellColData`.
#' @param useGroups The cell groups to be used for creating trajectory analysis.
#' @param groupBy A string indicating the column name from `cellColData` that contains the cell group definitions used in
#' `useGroups` to constrain trajectory analysis.
#' @param monocleCDS A monocle CDS object created from `getMonocleTrajectories`.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exists in the given `ArchRProject`.
#' @export
addMonocleTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useGroups = NULL, 
  groupBy = "Clusters",
  monocleCDS = NULL, 
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))

  .requirePackage("monocle3")

  groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
  groupDF <- groupDF[groupDF[,1] %in% useGroups,,drop=FALSE]

  if(sum(unique(groupDF[,1]) %in% useGroups)==0){
    stop("useGroups does not span any groups in groupBy! Are you sure your input is correct?")
  }

  monoclePT <- pseudotime(monocleCDS)
  monoclePT <- monoclePT[rownames(groupDF)]
  monoclePT <- ArchR:::.getQuantiles(monoclePT) * 100

  #Add To ArchR Project
  ArchRProj <- addCellColData(
      ArchRProj = ArchRProj,
      data = as.vector(monoclePT),
      name = name,
      cells = names(monoclePT),
      force = force
  )

  ArchRProj

}

###################################################################
# Slingshot
###################################################################


#' Add a Slingshot Trajectories to an ArchR Project #NEW
#' 
#' This function will fit a supervised trajectory in a lower dimensional space that 
#' can then be used for downstream analyses.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory to be added in `cellColData`.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column
#' in `cellColData`. This limits the groups used to identify trajectories.
#' @param principalGroup The principal group which represents the group that will be the starting point for all trajectories.
#' @param groupBy A string indicating the column name from `cellColData` that contains the cell group definitions used in
#' `useGroups` to constrain trajectory analysis.
#' @param embedding A string indicating the name of the `embedding` object from the `ArchRProject` that should be used for trajectory analysis.
#' @param reducedDims A string indicating the name of the `reducedDims` object from the `ArchRProject` that should be used for trajectory analysis. `embedding` must equal NULL to use.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exists in the given `ArchRProject`.
#' @param seed A number to be used as the seed for random number generation for trajectory creation.
#' @export
addSlingShotTrajectories <- function(
  ArchRProj = NULL,
  name = "SlingShot",
  useGroups = NULL,
  principalGroup = NULL,
  groupBy = NULL,
  embedding = NULL,
  reducedDims = NULL,
  force = FALSE,
  seed = 1
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character"))
  .validInput(input = principalGroup, name = "principalGroup", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character", "null"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("numeric"))

  .requirePackage("slingshot")

  set.seed(seed)

  if(!is.null(embedding)){
    rD <- getEmbedding(ArchRProj, embedding = embedding)
  }else{
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims)
  }

  groups <- getCellColData(ArchRProj, groupBy)

  if(!is.null(useGroups)){
    idx <- which(groups[,1] %in% useGroups)
    rD <- rD[idx, , drop = FALSE]
    groups <- groups[idx, , drop = FALSE]
  }

  sds <- slingshot(
    data = rD, 
    clusterLabels = groups[rownames(rD), ], 
      start.clus = principalGroup
  )

  #Get PseudoTimes
  pt <- slingPseudotime(sds)
  colnames(pt) <- paste0(name, ".Curve", seq_len(ncol(pt)))

  #Scale
  ptn <- apply(pt, 2, ArchR:::.getQuantiles) * 100

  for(i in seq_len(ncol(ptn))){

    ArchRProj <- addCellColData(
      ArchRProj = ArchRProj,
      data = as.vector(ptn[, i]),
      name = colnames(ptn)[i],
      cells = rownames(ptn),
      force = force
    )

  }

  ArchRProj

}

###################################################################
# STREAM
###################################################################

#' Get a PeakMatrix stored in an ArchRProject and write out for STREAM
#' 
#' This function gets a PeakMatrix from an `ArchRProject` and writes it to a set of files for STREAM (https://github.com/pinellolab/STREAM)
#'
#' @param ArchRProj An `ArchRProject` object to get data matrix from.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportPeakMatrixForSTREAM <- function(
  ArchRProj = NULL,
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("exportMatrixForSTREAM")
  ){

  mat <- getMatrixFromProject(
    ArchRProj = ArchRProj, 
    useMatrix = "PeakMatrix", 
    useSeqnames = useSeqnames, 
    verbose = verbose, 
    binarize = binarize,
    threads = threads,
    logFile = logFile
  )

  featureDF <- ArchR:::.getFeatureDF(getArrowFiles(ArchRProj)[1], "PeakMatrix")

  stopifnot(all(featureDF$idx == rowData(mat)$idx))

  countsDF <- Matrix::summary(assay(mat))
  peaksDF <- data.frame(as.vector(featureDF[,1]), featureDF[,3], featureDF[,4])
  cellsDF <- data.frame(colnames(mat))

  data.table::fwrite(countsDF, file = "STREAM_Counts.tsv.gz", sep = "\t", row.names = FALSE, col.names = FALSE)
  data.table::fwrite(peaksDF, file = "STREAM_Regions.tsv.gz", sep = "\t", row.names = FALSE, col.names = FALSE)
  data.table::fwrite(cellsDF, file = "STREAM_Sample.tsv.gz", sep = "\t", row.names = FALSE, col.names = FALSE)

  return(0)

}































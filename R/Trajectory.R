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
#' @param dof The number of degrees of freedom to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param spar The sparsity to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exists in the given `ArchRProject`.
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
  useAll = TRUE, 
  dof = 250,
  spar = 1,
  force = FALSE,
  logFile = createLogFile("addTrajectory"),
  seed = 1
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = trajectory, name = "trajectory", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = embedding, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = preFilterQuantile, name = "preFilterQuantile", valid = c("numeric"))
  .validInput(input = postFilterQuantile, name = "postFilterQuantile", valid = c("numeric"))
  .validInput(input = dof, name = "dof", valid = c("integer"))
  .validInput(input = spar, name = "spar", valid = c("numeric"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  if(!is.null(seed)) set.seed(seed)

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addTrajectory Input-Parameters", logFile=logFile)
  
  groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
  groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]

  if(sum(unique(groupDF[,1]) %in% trajectory)==0){
    .logMessage("trajectory does not span any groups in groupBy! Are you sure your input is correct?", logFile = logFile)
    stop("trajectory does not span any groups in groupBy! Are you sure your input is correct?")
  }

  if(sum(unique(groupDF[,1]) %in% trajectory) < 3){
    .logMessage("trajectory must span at least 3 groups in groupBy!", logFile = logFile)
    stop("trajectory must span at least 3 groups in groupBy!")
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
#' @param varCutOff The "Variance Quantile Cutoff" to be used for identifying the top variable features across the given trajectory.
#' Only features with a variance above the provided quantile will be retained.
#' @param maxFeatures The maximum number of features, ordered by variance, to consider from `useMatrix` when generating a trajectory.
#' This prevents smoothing a large number number of features which can be very time consuming.
#' @param groupEvery The number of sequential percentiles to group together when generating a trajectory. This is similar to smoothing
#' via a non-overlapping sliding window across pseudo-time. If `groupEvery = 2`, the values for percentiles [1 and 2], [3 and 4],
#' [5 and 6], etc. will be grouped together.
#' @param threads The number of threads to be used for parallel computing.
#' @param log2Norm A boolean value that indicates whether the summarized trajectory matrix should be log2 transformed. If you are using
#' a "MotifMatrix" set to FALSE.
#' @param scaleTo Once the sequential trajectory matrix is created, each column in that matrix will be normalized to a column sum
#' indicated by `scaleTo`. Setting this to `NULL` will prevent any normalization and should be done in certain circumstances
#' (for ex. if you are using a "MotifMatrix").
#' @param smoothWindow A boolean value indicating whether the sequential trajectory matrix should be furthered smooth to better reveal temporal dynamics.
#' @export
getTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  groupEvery = 1,
  threads = getArchRThreads(),
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = 3
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = groupEvery, name = "groupEvery", valid = c("numeric"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))

  trajectory <- getCellColData(ArchRProj, name)
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(0, 100, groupEvery)

  groupList <- lapply(seq_along(breaks), function(x){
      if(x == 1){
          NULL
      }else{
          rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
      }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])

  featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)

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
    
    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          smoothMat = smoothGroupMat, 
          mat = groupMat
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }else{

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          mat = groupMat
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
    scaleTo = scaleTo, 
    log2Norm = log2Norm, 
    smoothWindow = smoothWindow, 
    date = Sys.Date()
  )

  seTrajectory

}


#' @export
correlateTrajectories <- function(
  seTrajectory1 = NULL,
  seTrajectory2 = NULL,
  corCutOff = 0.5,
  varCutOff1 = 0.8,
  varCutOff2 = 0.8,
  removeFromName1 = c("underscore", "dash"),
  removeFromName2 = c("underscore", "dash"),
  useRanges = FALSE,
  fix1 = "center",
  fix2 = "start",
  maxDist = 250000,
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  force = FALSE,
  logFile = createLogFile("correlateTrajectories")
  ){

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "correlateTrajectories Input-Parameters", logFile=logFile)

  featureDF1 <- rowData(seTrajectory1)
  featureDF2 <- rowData(seTrajectory2)

  .logThis(featureDF1, "featureDF1", logFile = logFile)
  .logThis(featureDF2, "featureDF2", logFile = logFile)

  if("name" %in% colnames(featureDF1)){
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
  }else{
    if(!useRanges){
      .logMessage("seTrajectory1 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
      stop("seTrajectory1 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!")
    }
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
  }

  if("name" %in% colnames(featureDF2)){
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
  }else{
    if(!useRanges){
      .logMessage("seTrajectory2 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
      stop("seTrajectory2 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!")
    }
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
  }

  .logThis(rownames(featureDF1), "rownames(featureDF1)", logFile = logFile)
  .logThis(rownames(featureDF2), "rownames(featureDF2)", logFile = logFile)

  if(useRanges){

    if("start" %ni% colnames(featureDF1)){
      .logMessage("start is not in seTrajectory1, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
      stop("start is not in seTrajectory1, this is not a ranges object. Please set useRanges = FALSE")
    }

    if("start" %ni% colnames(featureDF2)){
      .logMessage("start is not in seTrajectory2, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
      stop("start is not in seTrajectory2, this is not a ranges object. Please set useRanges = FALSE")
    }

    if("strand" %in% colnames(featureDF1)){
      ranges1 <- GRanges(
        seqnames = featureDF1$seqnames, 
        IRanges(
          ifelse(featureDF1$strand == 2, featureDF1$end, featureDF1$start),
          ifelse(featureDF1$strand == 2, featureDF1$start, featureDF1$end)
        ),
        strand = ifelse(featureDF1$strand == 2, "-", "+")
      )
    }else{
      ranges1 <- GRanges(featureDF1$seqnames, IRanges(featureDF1$start, featureDF1$end))
    }
    #mcols(ranges1) <- featureDF1
    names(ranges1) <- rownames(featureDF1)
    rowRanges(seTrajectory1) <- ranges1
    rm(ranges1)

    if("strand" %in% colnames(featureDF2)){
      ranges2 <- GRanges(
        seqnames = featureDF2$seqnames, 
        IRanges(
          ifelse(featureDF2$strand == 2, featureDF2$end, featureDF2$start),
          ifelse(featureDF2$strand == 2, featureDF2$start, featureDF2$end)
        ),
        strand = ifelse(featureDF2$strand == 2, "-", "+")
      )
    }else{
      ranges2 <- GRanges(featureDF2$seqnames, IRanges(featureDF2$start, featureDF2$end))
    }
    #mcols(ranges2) <- featureDF2
    names(ranges2) <- rownames(featureDF2)
    rowRanges(seTrajectory2) <- ranges2
    rm(ranges2)

    .logThis(ranges1, "ranges1", logFile = logFile)
    .logThis(ranges2, "ranges2", logFile = logFile)

    #Find Associations to test
    isStranded1 <- any(as.integer(strand(seTrajectory1)) == 2)
    isStranded2 <- any(as.integer(strand(seTrajectory2)) == 2)

    if(fix1 == "center" & isStranded1){
      if(!force){
        .logMessage("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
        stop("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.")
      }else{
        message("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE")
      }         
    }

    if(fix1 != "center" & !isStranded1){
      if(!force){
        .logMessage("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
        stop("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.")
      }else{
        message("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE") 
      }         
    }

    if(fix2 == "center" & isStranded2){
      if(!force){
        .logMessage("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
        stop("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.")
      }else{
        message("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE")
      }         
    }

    if(fix2 != "center" & !isStranded2){
      if(!force){
        .logMessage("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
        stop("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.")
      }else{
        message("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE") 
      }         
    }

    #Overlaps
    mappingDF <- DataFrame(
      findOverlaps( 
        resize(rowRanges(seTrajectory1), 1, fix1), 
        .suppressAll(resize(resize(seTrajectory2, 1, fix2), 2 * maxDist + 1, "center")),
        ignore.strand = TRUE
      )
    )

    #Get Distance 
    mappingDF$distance <- distance(
      x = ranges(rowRanges(seTrajectory1)[mappingDF[,1]]), 
      y = ranges(rowRanges(resize(seTrajectory2, 1, fix2))[mappingDF[,2]])
    )

  }else{

    #Create Match Names
    featureDF1$matchName <- featureDF1$name
    if("underscore" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\_.*","",featureDF1$matchName)
    }
    if("dash" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\-.*","",featureDF1$matchName)
    }
    if("numeric" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("[0-9]+","",featureDF1$matchName)
    }

    featureDF2$matchName <- featureDF2$name
    if("underscore" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\_.*","",featureDF2$matchName)
    }
    if("dash" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\-.*","",featureDF2$matchName)
    }
    if("numeric" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("[0-9]+","",featureDF2$matchName)
    }
    
    #Now Lets see how many matched pairings
    matchP1 <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
    matchP2 <- sum(featureDF2$matchName %in% featureDF1$matchName) / nrow(featureDF2)
    matchP <- max(matchP1, matchP2)

    .logThis(featureDF1$matchName, "featureDF1$matchName", logFile)
    .logThis(featureDF2$matchName, "featureDF2$matchName", logFile)

    if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
      .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!", logFile = logFile)
      stop("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!")
    }
    if(matchP < 0.05){
      if(!force){
        .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Set force = TRUE to continue!", logFile = logFile)
        stop("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Set force = TRUE to continue!")
      }else{
        message("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Continuing since force = TRUE.")
      }
    }

    #Create Mappings
    mappingDF <- lapply(seq_len(nrow(featureDF1)), function(x){
      idx <- which(paste0(featureDF2$matchName) %in% paste0(featureDF1$matchName[x]))
      if(length(idx) > 0){
        expand.grid(x, idx)
      }else{
        NULL
      }
    }) %>% Reduce("rbind", .)

  }

  colnames(mappingDF)[1:2] <- c("idx1", "idx2")
  mappingDF <- DataFrame(mappingDF)

  .logThis(mappingDF, "mappingDF", logFile = logFile)

  if(!useRanges){
    mappingDF$matchname1 <- featureDF1$matchName[mappingDF$idx1]
    mappingDF$matchname2 <- featureDF2$matchName[mappingDF$idx2]
    mappingDF$name1 <- rownames(featureDF1)[mappingDF$idx1]
    mappingDF$name2 <- rownames(featureDF2)[mappingDF$idx2]
  }

  mappingDF$Correlation <- rowCorCpp(
    idxX = as.integer(mappingDF[,1]), 
    idxY = as.integer(mappingDF[,2]), 
    X = assays(seTrajectory1)[["mat"]], 
    Y = assays(seTrajectory2)[["mat"]]
  )
  mappingDF$VarAssay1 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory1)[["mat"]]))[as.integer(mappingDF[,1])]
  mappingDF$VarAssay2 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory2)[["mat"]]))[as.integer(mappingDF[,2])]
  mappingDF$TStat <- (mappingDF$Correlation / sqrt((1-mappingDF$Correlation^2)/(ncol(seTrajectory1)-2))) #T-statistic P-value
  mappingDF$Pval <- 2 * pt(-abs(mappingDF$TStat), ncol(seTrajectory1) - 2)
  mappingDF$FDR <- p.adjust(mappingDF$Pval, method = "fdr")

  idxPF <- which(mappingDF$Correlation > corCutOff & mappingDF$VarAssay1 > varCutOff1 & mappingDF$VarAssay2 > varCutOff2)
  message("Found ", length(idxPF), " Correlated Pairings!")

  .logThis(mappingDF[idxPF,], "mappingDF-PF", logFile = logFile)

  out <- SimpleList(
    correlatedMappings = mappingDF[idxPF,],
    allMappings = mappingDF,
    seTrajectory1 = seTrajectory1,
    seTrajectory2 = seTrajectory2
  )

  out

}


#' Plot a Heatmap of Features across a Trajectory
#' 
#' This function will plot a heatmap of the results from getTrajectory
#' 
#' @param seTrajectory A `SummarizedExperiment` object that results from calling `markerFeatures()`.
#' @param scaleRows A boolean value that indicates whether row-wise z-scores should be computed on the matrix provided by `seTrajectory`.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies
#' rownames from `seTrajectory` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the `rownames` of `seTrajectory` that should be labeled on the side of the heatmap.
#' @param labelTop A number indicating how many of the top N features, based on variance, in `seTrajectory` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMat A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @export
trajectoryHeatmap <- function(
  seTrajectory = NULL,
  varCutOff = 0.9,
  maxFeatures = 25000,
  scaleRows = TRUE,
  limits = c(-2,2),
  grepExclude = NULL,
  pal = NULL,
  labelMarkers = NULL,
  labelTop = 50,
  labelRows = FALSE,
  rowOrder = NULL, 
  returnMat = FALSE
  ){

  .validInput(input = seTrajectory, name = "seTrajectory", valid = c("SummarizedExperiment"))
  .validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
  .validInput(input = limits, name = "limits", valid = c("numeric"))
  .validInput(input = grepExclude, name = "grepExclude", valid = c("character", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = labelMarkers, name = "labelMarkers", valid = c("character", "null"))
  .validInput(input = labelTop, name = "labelTop", valid = c("integer"))
  .validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
  .validInput(input = returnMat, name = "returnMat", valid = c("boolean"))

  mat <- assay(seTrajectory)
  varQ <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory)[["mat"]]))
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

  if(!is.null(labelMarkers)){
    idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
    idxLabel2 <- idxLabel2[idxLabel2 > 0]
  }else{
    idxLabel2 <- NULL
  }

  idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])

  if(scaleRows){
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
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

  ht <- .ArchRHeatmap(
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

  if(returnMat){
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
  .validInput(input = quantCut, name = "quantCut", valid = c("numeric"))
  .validInput(input = quantHex, name = "quantHex", valid = c("numeric"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character", "null"))
  .validInput(input = continuousSet, name = "continuousSet", valid = c("character", "null"))
  .validInput(input = randomize, name = "randomize", valid = c("boolean"))
  .validInput(input = keepAxis, name = "keepAxis", valid = c("boolean"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = addArrow, name = "addArrow", valid = c("boolean"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character", "null"))

  .requirePackage("ggplot2", source = "cran")

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)

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

































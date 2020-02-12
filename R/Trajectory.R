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
#' @param trajectory The order of cell groups to be used for constraining the initial supervised fitting procedure. For example, to get a trajectory from Cluster1 to Cluster2 to Cluster3, input should be c("Cluster1", "Cluster2", "Cluster3"). Cells will then be used from these 3 groups to constrain an initial fit in the group order.
#' @param groupBy A string indicating the column name from `cellColData` that contains the cell group definitions used in `trajectory` to constrain the initial supervised fitting procedure.
#' @param reducedDims A string indicating the name of the `reducedDims` object from the `ArchRProject` that should be used for distance computation.
#' @param embedding A string indicating the name of the `embedding` object from the `ArchRProject` that should be used for distance computation.
#' @param preFilterQ Prior to the initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping center is above the provided quantile will be excluded.
#' @param postFilterQ After initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping center is above the provided quantile will be excluded.
#' @param dof The number of degrees of freedom to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param spar The sparsity to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exist in the given `ArchRProject`.
#' @export
addTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  trajectory = NULL, 
  groupBy = "Clusters",
  reducedDims = "IterativeLSI",
  embedding = NULL,
  preFilterQ = 0.9, 
  postFilterQ = 0.9, 
  dof = 250,
  spar = 1,
  force = FALSE
  ){

    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = name, name = "name", valid = c("character"))
    .validInput(input = trajectory, name = "trajectory", valid = c("character"))
    .validInput(input = groupBy, name = "groupBy", valid = c("character"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character", "null"))
    .validInput(input = embedding, name = "reducedDims", valid = c("character", "null"))
    .validInput(input = preFilterQ, name = "preFilterQ", valid = c("numeric"))
    .validInput(input = postFilterQ, name = "postFilterQ", valid = c("numeric"))
    .validInput(input = dof, name = "dof", valid = c("integer"))
    .validInput(input = spar, name = "spar", valid = c("numeric"))
    .validInput(input = force, name = "force", valid = c("boolean"))

    set.seed(1)

    #Knn Method
    if(requireNamespace("nabor", quietly = TRUE)){
        knnMethod <- nabor::knn
    }else if(requireNamespace("RANN", quietly = TRUE)){
        knnMethod <- RANN::nn2
    }else if(requireNamespace("FNN", quietly = TRUE)){
        knnMethod <- FNN::get.knnx
    }else{
        stop("Computing KNN requires package nabor, RANN or FNN")
    }
    
    groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
    groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]

    if(is.null(embedding)){
      mat <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims)
    }else{
      mat <- getEmbedding(ArchRProj = ArchRProj, embedding = embedding)
    }
    mat <- mat[rownames(groupDF),,drop = FALSE]

    ######################################################
    #Filter Outliers
    ######################################################
    filterObj <- lapply(seq_along(trajectory), function(x){
        
        #Subset
        groupsx <- rownames(groupDF)[groupDF[,1]==trajectory[x]]
        matx <- mat[groupsx,,drop = FALSE]

        #Filter Distance
        matMeanx <- colMeans(matx)
        diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
        idxKeep <- which(diffx <= quantile(diffx, preFilterQ))
        
        #Filter
        list(mat = matx[idxKeep,,drop=FALSE], groups = groupsx[idxKeep])

    })

    matFilter <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat) %>% Reduce("rbind", .)
    groupsFilter <- groupDF[lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups) %>% Reduce("c", .),,drop=FALSE]

    ######################################################
    #Now Initial Alignment
    ######################################################
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
    matSpline <- lapply(seq_len(ncol(matFilter)), function(x){
        smooth.spline(
            x = initialTime, 
            y = matFilter[names(initialTime), x], 
            df = dof, 
            spar = spar
        )[[2]]
    }) %>% Reduce("cbind",.) %>% data.frame()

    ######################################################
    # 1. KNN Fit vs Actual
    ######################################################
    knnObj <- knnMethod(
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
    idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], postFilterQ))
    dfTrajectory <- DataFrame(
        row.names = rownames(mat),
        Distance = knnDist[, 1],
        DistanceIdx = knnIdx[, 1] + knnDistQ
    )[idxKeep, , drop = FALSE]

    ######################################################
    # 2. Fit cells not in trajectory clusters
    ######################################################
    mat2 <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims)
    groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
    groupDF <- groupDF[groupDF[,1] %ni% trajectory,,drop=FALSE]
    mat2 <- mat2[rownames(groupDF),,drop = FALSE]

    #Nearest Neighbors
    knnObj2 <- knnMethod(
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
        DistanceIdx = knnIdx2[, 1] + knnDistQ2
    )[idxKeep, , drop = FALSE]

    #Final Output
    dfTrajectory3 <- rbind(dfTrajectory, dfTrajectory2)
    dfTrajectory3$Trajectory <- 100 * .getQuantiles(dfTrajectory3[,2])
    
    #Add To ArchR Project
    ArchRProj <- addCellColData(
        ArchRProj = ArchRProj,
        data = dfTrajectory3$Trajectory,
        name = name,
        cells = rownames(dfTrajectory3),
        force = force
    )

    ArchRProj

}

#' Get Supervised Trajectory from an ArchR Project
#' 
#' This function will get a supervised trajectory from an `ArchRProject` (see `addTrajectory`), get data
#' from a desired matrix, and smooth each value across the input trajectory.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory in `cellColData` to retrieve from the given `ArchRProject`.
#' @param useMatrix The name of the data matrix from the `ArrowFiles` to get numerical values for each cell from. Recommended matrices are "GeneScoreMatrix", "PeakMatrix", or "MotifMatrix".
#' @param varCutOff The "Variance Quantile Cutoff" to be used for identifying the top variable features across the given trajectory. Only features with a variance above the provided quantile will be retained.
#' @param maxFeatures The maximum number of features, ordered by variance, to consider from `useMatrix` when generating a trajectory. This prevents smoothing a large number number of features which can be very time consuming.
#' @param groupEvery The number of sequential percentiles to group together when generating a trajectory. This is similar to smoothing via a non-overlapping sliding window across pseudo-time. If `groupEvery = 2`, the values for percentiles [1 and 2], [3 and 4], [5 and 6], etc. will be grouped together.
#' @param threads The number of threads to be used for parallel computing.
#' @param log2Norm A boolean value that indicates whether the summarized trajectory matrix should be log2 transformed. If you are using a "MotifMatrix" set to FALSE.
#' @param scaleTo Once the sequential trajectory matrix is created, each column in that matrix will be normalized to a column sum indicated by `scaleTo`. Setting this to `NULL` will prevent any normalization and should be done in certain circumstances (for ex. if you are using a "MotifMatrix").
#' @param smooth A boolean value indicating whether the sequential trajectory matrix should be furthered smooth to better reveal temporal dynamics.
#' @param smoothFormula The smoothing formula to use in the generalized additive model. See the `formula` parameter in `mgcv::gam()` for additional information.
#' @export
getTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  varCutOff = 0.9,
  maxFeatures = 25000,
  groupEvery = 2,
  threads = getArchRThreads(),
  log2Norm = TRUE,
  scaleTo = 10000,
  smooth = TRUE,
  smoothFormula = "y ~ s(x, bs = 'cs')"
  ){

    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = name, name = "name", valid = c("character"))
    .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
    .validInput(input = varCutOff, name = "varCutOff", valid = c("numeric"))
    .validInput(input = maxFeatures, name = "maxFeatures", valid = c("numeric"))
    .validInput(input = groupEvery, name = "groupEvery", valid = c("numeric"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
    .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    .validInput(input = smooth, name = "smooth", valid = c("boolean"))
    .validInput(input = smoothFormula, name = "smoothFormula", valid = c("character"))

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
      if(any(groupMat) < 0){
        message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
      }else{
        groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
      }
    }

    if(log2Norm){
      if(any(groupMat) < 0){
        message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
      }else{
        groupMat <- log2(groupMat + 1)
      }
    }

    if(!is.null(varCutOff)){
      rV <- matrixStats::rowVars(groupMat)
      idx <- head(order(rV, decreasing = TRUE), nrow(groupMat) * (1-varCutOff))
      groupMat <- groupMat[idx, ,drop=FALSE]
    }

    if(nrow(groupMat) > maxFeatures){
      rV <- matrixStats::rowVars(groupMat)
      idx2 <- head(order(rV, decreasing = TRUE), nrow(groupMat) * (1-varCutOff))
      idx <- idx[idx2]
      groupMat <- groupMat[idx2, ,drop=FALSE]
    }

    if(smooth){
      
      message("Smoothing with mgcv::gam formula : ", smoothFormula)

      t <- breaks[-length(breaks)] + 0.5 * groupEvery
      groupMat2 <- matrix(NA, nrow=nrow(groupMat),ncol=ncol(groupMat))
      colnames(groupMat2) <- colnames(groupMat)
      rownames(groupMat2) <- rownames(groupMat)
      n <- ncol(groupMat2)

      pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
      for(x in seq_len(nrow(groupMat))){
        setTxtProgressBar(pb,round(x*100/nrow(groupMat),0))
        groupMat2[x, ] <- stats::predict(
            mgcv::gam(formula = eval(parse(text=smoothFormula)), data = data.frame(x = t, y = groupMat[x, ])), 
            newdata = data.frame(x = t), 
            se.fit = FALSE,
            level = 0.95, 
            interval = "confidence"
          )
      }

      #Rename and remove
      groupMat <- groupMat2
      rm(groupMat2)
      gc()

      message("\n")

    }

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(mat = groupMat), 
        rowData = featureDF[idx, ,drop = FALSE]
    )
    metadata(seTrajectory)$Params <- list(useMatrix = useMatrix, 
      scaleTo = scaleTo, log2Norm = log2Norm, smooth = smooth, smoothFormula = smoothFormula, date = Sys.Date())

    seTrajectory

}

#' Plot a Heatmap of Features across a Trajectory
#' 
#' This function will plot a heatmap of the results from getTrajectory
#' 
#' @param seTrajectory A `SummarizedExperiment` object that results from calling `markerFeatures()`.
#' @param scaleRows A boolean value that indicates whether row-wise z-scores should be computed on the matrix provided by `seTrajectory`.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies rownames from `seTrajectory` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the `rownames` of `seTrajectory` that should be labeled on the side of the heatmap.
#' @param labelTop A number indicating how many of the top N features, based on variance, in `seTrajectory` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMat A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @export
trajectoryHeatmap <- function(
  seTrajectory = NULL,
  scaleRows = TRUE,
  limits = c(-2,2),
  grepExclude = NULL,
  pal = NULL,
  labelMarkers = NULL,
  labelTop = 50,
  labelRows = FALSE,
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
  rownames(mat) <- rowData(seTrajectory)$name
  
  if(!is.null(labelTop)){
    idxLabel <- rownames(mat)[seq_len(labelTop)]
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

  idx <- order(apply(mat, 1, which.max))

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
    draw = FALSE
  )

  if(returnMat){
    return(mat)
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
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData") or by a data matrix in the associated ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName/rowname of the data matrix to be used for plotting. 
#' For example if colorBy is `cellColData` then name refers to a column name in the cellcoldata (see `getCellcoldata()`), if colorBy is `GeneScoreMatrix` then name refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values from `colorBy`.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells' values. See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for coloring cells.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param quantCut If this is not `NULL`, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values. 
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the lower threshold and y is 
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the 2.5th percentile and 97.5 percentile of values and set values below/above to the value of 
#' the 2.5th and 97.5th percentile values respectively.
#' @param quantHex JJJ The numeric xth quantile of all dots within each individual hexagon will determine the the numerical value for coloring to be displayed. This occurs when `plotAs` = "hex" or `NULL` (if numerical values by default).
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a discrete color set is desired.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a continuous color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed arrow in the embedding based on the aligned trajectory.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted. By default if `colorBy` is numeric this is "hex".
#' @param plotParams Additional parameters to pass to `ggPoint()` or `ggHex()`.
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
  size = 0.5,
  rastr = TRUE,
  quantCut = c(0.05, 0.95),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  plotParams = list()
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
  .validInput(input = plotParams, name = "plotParams", valid = c("list"))

  .requirePackage("ggplot2", source = "cran")

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
  .quantileCut0 <- function (x = NULL, lo = 0, hi = 0.975, rm0 = TRUE){
    q <- quantile(x, probs = c(lo, hi), na.rm = TRUE)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(x)
  }

  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex)
  }

  ##############################
  # Get Trajectory
  ##############################
  dfT <- getCellColData(ArchRProj, select = trajectory)
  idxRemove <- which(is.na(dfT[,1]))

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("x", "y", "PseudoTime")

  #Parameters
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    
    plotParams$color <- as.vector(getCellColData(ArchRProj)[,name])
    plotParams$discrete <- .isDiscrete(plotParams$color)
    plotParams$continuousSet <- "solarExtra"
    plotParams$discreteSet <- "stallion"
    plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }

  }else{
    if (tolower(colorBy) == "genescorematrix"){
      if(is.null(log2Norm)){
        log2Norm <- TRUE
      }
      plotParams$continuousSet <- "comet"
    }else{
      plotParams$continuousSet <- "solarExtra"
    }
    plotParams$color <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)[1,]
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

  plotParams$color[idxRemove] <- NA

  if(!plotParams$discrete){
    if(!is.null(imputeWeights)){
      imputeWeights <- imputeWeights$Weights[rownames(df), rownames(df)]
      plotParams$color <- (imputeWeights %*% as(as.matrix(plotParams$color), "dgCMatrix"))[,1] 
    }else{
      plotParams$color <- .quantileCut0(plotParams$color, min(quantCut), max(quantCut))
    }
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)
    if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
      plotParams$addPoints <- TRUE
      out <- do.call(ggHex, plotParams)
    }else{
      out <- do.call(ggPoint, plotParams)
    }
  }else{
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  #Prep Trajectory Vector
  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  #Plot Pseudo-Time
  out2 <- ggPoint(dfT$PseudoTime, dfT$value, dfT$PseudoTime, 
    discrete = FALSE, xlabel = "PseudoTime", ylabel = name, ratioYX = 0.5, rastr = TRUE) +
    geom_smooth(color = "black")

  attr(out2, "ratioYX") <- 0.5

  if(addArrow){
    dfArrow <- .splitEvery(dfT, floor(nrow(dfT) / 10)) %>% 
      lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame
    out <- out + geom_path(
            data = data.frame(dfArrow), aes(x, y, color=NULL), size= 1, 
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches"))
          )
  }

  list(out, out2)

}


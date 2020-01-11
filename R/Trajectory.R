####################################################################
# Trajectory Analysis Methods
####################################################################

#' Add Supervised Trajectory to an ArchR Project
#' 
#' This function will fit a supervised trajectory in a lower dimensional space that 
#' can then be used for downstream analyses.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory to be added in `cellColData`.
#' @param trajectory Supervised order of groups to constrain supervised fitting (ie c("Cluster1", "Cluster2", "Cluster3") )
#' @param groupBy initial group column in cellColData to constrain supervised fit
#' @param reducedDims A string indicating the name of the `reducedDims` object from the `ArchRProject` that should be used for distance computation.
#' @param preQuantile Prior to supervised trajectory fitting, the quantile for filtering cells that are far (by euclidean distance) from cluster centers.
#' @param postQuantile Post supervised trajectory fitting, the quantile for determining the cutoff for cells not in the groups to be aligned to the trajectory.
#' @param dof The number of degrees of freedom to be used in spline fit (see smooth.spline)
#' @param spar The sparsity to be used in spline fit (see smooth.spline)
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exist in the given `ArchRProject`.
#' @param ... additional args
#' @export
addTrajectory <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  trajectory = NULL, 
  groupBy = "Clusters",
  reducedDims = "IterativeLSI",
  preQuantile = 0.1, 
  postQuantile = 0.1, 
  dof = 250,
  spar = 1,
  force = FALSE,
  ...
  ){

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

    mat <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims)
    groupDF <- getCellColData(ArchRProj = ArchRProj, select = groupBy)
    groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]
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
        idxKeep <- which(diffx <= quantile(diffx, 1 - preQuantile))
        
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
    idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], 1 - postQuantile))
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
        force = force,
        ...
    )

    ArchRProj

}

#' Get Supervised Trajectory from an ArchR Project
#' 
#' This function will fit get a supervised trajectory from an ArchRProject and aggregate signal
#' from a matrix and smooth across the trajectory
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name A string indicating the name of the fitted trajectory in `cellColData` to retrieve from the given `ArchRProject`.
#' @param useMatrix The name of the data matrix from the `ArrowFiles` to get numerical values for each cell from.
#' @param varCutOff Variance Quantile Cutoff for identifying top variable features across trajectory.
#' @param maxFeatures The maximum number of features (ordered by variance) to consider from `useMatrix` when generating a trajectory.
#' @param groupEvery The number of sequential percentiles (0-100 every x percentile) to group together when generating a trajectory (similar to smoothing).
#' @param threads The number of threads to be used for parallel computing.
#' @param scaleTo Scale group data matrix to this value for normalizaiton.
#' @param log2Norm A boolean value that indicates whether the summarized trajectory matrix should be log2 transformed.
#' @param smooth A boolean value indicating whether the summarized trajectory matrix should be smoothed in a row-wise fashion.
#' @param smoothFormula The smoothing formula to use in the generalized additive model. See the `formula` parameter in `mgcv::gam()`.
#' @param ... additional args
#' @export
getTrajectory <- function(
  ArchRProj,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  varCutOff = 0.1,
  maxFeatures = 25000,
  groupEvery = 2,
  threads = 1,
  scaleTo = 10000,
  log2Norm = TRUE,
  smooth = TRUE,
  smoothFormula = "y ~ s(x, bs = 'cs')",
  ...
  ){

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

    message("Creating Trajectory Group Matrix...")
    groupMat <- .getGroupMatrix(
        ArrowFiles = getArrowFiles(ArchRProj), 
        featureDF = featureDF,
        groupList = groupList, 
        threads = threads, 
        verbose = FALSE, 
        useMatrix = useMatrix
    )

    #Scale
    groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo

    if(log2Norm){
      groupMat <- log2(groupMat + 1)
    }

    if(!is.null(varCutOff)){
      rV <- matrixStats::rowVars(groupMat)
      idx <- head(order(rV, decreasing = TRUE), nrow(groupMat) * varCutOff)
      groupMat <- groupMat[idx, ,drop=FALSE]
    }

    if(nrow(groupMat) > maxFeatures){
      rV <- matrixStats::rowVars(groupMat)
      idx2 <- head(order(rV, decreasing = TRUE), nrow(groupMat) * varCutOff)
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
#' @param seTrajectory A `SummarizedExperiment` object that results from `markerFeatures()`.
#' @param scaleRows A boolean value that indicates whether row-wise z-scores should be computed on matrix in `seTrajectory`.
#' @param limits A numeric vector of two numbers that represent the lower and upper color limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies rownames from `seTrajectory` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see paletteContinuous) used to override the continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the `rownames` of `seTrajectory` that should be labeled on the side of the heatmap.
#' @param labelTop A boolean value that indicates whether the top features for each column in `seTrajectory` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMat A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param ... additional args
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
  returnMat = FALSE,
  ...
  ){

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

  idxLabel <- c(idxLabel, idxLabel2)

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
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }else if(tolower(metadata(seTrajectory)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blue_yellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }
  }

  idx <- order(apply(mat, 1, which.max))

  ht <- ArchR:::.ArchRHeatmap(
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
    draw = FALSE,
    ...
  )

  if(returnMat){
    return(mat)
  }else{
    return(ht)
  }

}

#' Visualize a Trajectory from ArchR Project
#' 
#' This function will plot a trajectory that was created from addTrajectory onto an embedding.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the embedding to use to visualize the given `trajectory`. See `ArchR::addEmbedding()` for more information.
#' @param trajectory The name of the trajectory as a column in `cellColData` to plot.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in cellColData ("cellColData") or by a data matrix in the ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrx", "PeakMatrix").
#' @param name The name of the trajectory as a column in `cellColData` to plot.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values prior to plotting.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for coloring cells.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param quantCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values. 
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the upper threshold and y is 
#' the lower threshold. For example, quantileCut = c(0.975,0.025) will take the top and bottom 2.5% of values and set them to the value of 
#' the 97.5th and 2.5th percentile values respectively.
#' @param quantHex The quantile to 
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing colorBy in the embedding.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing colorBy in the embedding.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed arrow in the embedding based on the aligned trajectory.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted.
#' @param plotParams Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' @param ... additional args
#' @export
plotTrajectory <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  trajectory = "Trajectory",
  colorBy = "colData",
  name = "Trajectory",
  log2Norm = NULL,
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
  plotParams = list(),
  ...
  ){

  .requirePackage("ggplot2")

  ##############################
  # Plot Helpers
  ##############################
  .quantileCut0 <- function (x, lo = 0, hi = 0.975, rm0 = TRUE){
    q <- quantile(x, probs = c(lo, hi), na.rm = TRUE)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(x)
  }

  .summarizeHex <- function(x){
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
  df <- getEmbedding(ArchRProj, embedding = embedding, return = "df")
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
    plotParams$continuousSet <- "solar_extra"
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
      plotParams$continuousSet <- "white_blue_purple"
    }else{
      plotParams$continuousSet <- "solar_extra"
    }
    plotParams$color <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)
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
    plotParams$color <- .quantileCut0(plotParams$color, min(quantCut), max(quantCut))
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


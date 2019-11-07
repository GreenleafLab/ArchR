#' Add Supervised Trajectory to an ArchR Project
#' 
#' This function will fit a supervised trajectory in a lower dimensional space that 
#' can then be used for downstream analyses.
#'
#' @param ArchRProj ArchRProject
#' @param name name of fitted trajectory to be added in cellColData
#' @param trajectory trajectory of groups to constrain supervised fitting (in order)
#' @param groupBy initial group column in cellColData to constrain supervised fit
#' @param name name of column in cellColData or Feature in Array in Arrows
#' @param reducedDims name of reduced dimensions used for distance computation
#' @param preFilter pre filtering quantile for supervised trajectory fit
#' @param postFilter post filtering quantile for supervised trajectory fit
#' @param dof degrees of freedom
#' @param spar sparsity
#' @param force force addition into ArchRProject if the column name exists
#' @param ... additional args
#' @export
addTrajectory <- function(
  ArchRProj,
  name = "Trajectory",
  trajectory = NULL, 
  groupBy = "Clusters",
  reducedDims = "IterativeLSI",
  preFilter = 0.1, 
  postFilter = 0.1, 
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
        idxKeep <- which(diffx <= quantile(diffx, 1 - preFilter))
        
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
    idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], 1 - postFilter))
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
#' @param ArchRProj ArchRProject
#' @param name name of fitted trajectory in cellColData
#' @param useMatrix matrix to summarize across the trajectory
#' @param groupEvery group cells every x quantile
#' @param scaleTo scale summarized matrix to
#' @param log2Norm log2 normalize the scaled summarized matrix
#' @param threads number of threads to use
#' @param ... additional args
#' @export
getTrajectory <- function(
  ArchRProj,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  groupEvery = 2,
  threads = 1,
  scaleTo = 10000,
  log2Norm = TRUE,
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

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(mat = groupMat), 
        rowData = featureDF
    )

    seTrajectory

}


#' Visualize Embedding from ArchR Project
#' 
#' This function will plot an embedding that was created from
#' computeEmbedding
#'
#' @param ArchRProj ArchRProject
#' @param embedding embedding to visualize (see computeEmbedding)
#' @param colorBy colorBy cellColData or Arrays in Arrows (ie GeneScoreMatrix)
#' @param name name of column in cellColData or Feature in Array in Arrows
#' @param log2Norm log2 Normalize features if they are continuous
#' @param pal custom palette to use for plotting
#' @param size size of points in plot
#' @param rastr rastr points in plot
#' @param quantCut quantile cut of continuous features
#' @param quantHex quantile evaluation for each hex in geom_hex
#' @param discreteSet discrete palette for visualizing embedding
#' @param continuousSet continuous palette for visualizing embedding
#' @param randomize randomize points prior to plotting
#' @param keepAxis keep x and y axis for plot
#' @param baseSize base size for text in plot
#' @param plotAs how to plot (points vs hex)
#' @param plotParams additional params to pass to ggPoint/ggHex
#' @param ... additional args
#' @export
VisualizeTrajectory <- function(
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
  .quantileCut <- function (x, lo = 0, hi = 0.975, rm0 = TRUE){
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
    plotParams$color <- .quantileCut(plotParams$color, min(quantCut), max(quantCut))
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

.splitEvery <- function(x, n){
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if(is.atomic(x)){
    split(x, ceiling(seq_along(x) / n))
  }else{
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}






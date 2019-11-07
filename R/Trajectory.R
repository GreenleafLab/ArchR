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








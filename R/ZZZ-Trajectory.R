#' Get Optimal Plotting Windows
#'
#' @param rle rle
#' @export
alignCellsToTrajectory <- function(mat, groups, trajectory, filterInitial = 0.1, finterFinal = 0.25, dof = 250, spar = 1){
    
    #Filter by Distance
    matFilter <- lapply(seq_along(trajectory), function(x){
        #Subset
        matx <- mat[groups==trajectory[x],]
        groupsx <- groups[groups==trajectory[x]]
        #Filter Distance
        matMeanx <- colMeans(matx)
        diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
        idxKeep <- which(diffx <= quantile(diffx, 1 - filterInitial))
        #Filter
        matx <- matx[idxKeep,]
        matx$Groups <- groupsx[idxKeep]
        matx
    }) %>% Reduce("rbind",.)

    #Now Initial Alignment
    initialTime <- lapply(seq_along(trajectory), function(x){
        #Subset
        matx <- matFilter[matFilter$Groups==trajectory[x], which(colnames(matFilter) %ni% "Groups")]
        #Get Differences
        if(x!=length(trajectory)){
            matMeanxp1 <- colMeans(matFilter[matFilter$Groups==trajectory[x+1], which(colnames(matFilter) %ni% "Groups")])
            diffx1 <- sqrt(colSums((t(matx) - matMeanxp1)^2))
            timex <- (1 - getQuantiles(diffx1)) + x
        }else{
            matMeanxm1 <- colMeans(matFilter[matFilter$Groups==trajectory[x-1], which(colnames(matFilter) %ni% "Groups")])
            diffx1 <- sqrt(colSums((t(matx) - matMeanxm1)^2))
            timex <- getQuantiles(diffx1) + x
        }
        timex
    }) %>% unlist

    #Fit Splines
    dfSplineFit <- lapply(seq_len(ncol(matFilter[, which(colnames(matFilter) %ni% "Groups")])), function(x){
        smooth.spline(initialTime, matFilter[, which(colnames(matFilter) %ni% "Groups")][,x], df = dof, spar = spar)[[2]]
    }) %>% Reduce("cbind",.) %>% data.frame()

    #Trajectories
    dfTrajectory <- mat[groups %in% trajectory,]
    dfTrajectory$Groups <- groups[groups %in% trajectory]

    #Nearest Neighbors
    knnMat <- get.knnx(
        data = dfSplineFit[,which(colnames(matFilter) %ni% "Groups")],
        query = dfTrajectory[,which(colnames(matFilter) %ni% "Groups")], 
        k = 3)

    #Lets Create Pseudotime
    knn_index <- knnMat[[1]]
    knn_dist <- knnMat[[2]]
    knn_diff <- ifelse(knn_index[,2] > knn_index[,3], 1, -1)
    knn_distq <- getQuantiles(knn_dist[,1])

    #Filter
    idxKeep <- which(knn_dist[,1] < quantile(knn_dist[,1], 1 - filterFinal))
    finalTrajectory <- data.frame(
        Cells = rownames(dfTrajectory)[idxKeep],
        Distance = knn_dist[idxKeep, 1],
        DistanceIdx = knn_index[idxKeep, 1] + knn_distq[idxKeep]
        )

    #Lets Align the non Clusters
    dfTrajectory <- mat[groups %ni% trajectory,]
    dfTrajectory$Groups <- groups[groups %ni% trajectory]

    #Nearest Neighbors
    knnMat <- get.knnx(
        data = dfSplineFit[,which(colnames(matFilter) %ni% "Groups")],
        query = dfTrajectory[,which(colnames(matFilter) %ni% "Groups")], 
        k = 3)

    #Lets Create Pseudotime
    knn_index <- knnMat[[1]]
    knn_dist <- knnMat[[2]]
    knn_diff <- ifelse(knn_index[,2] > knn_index[,3], 1, -1)
    knn_distq <- getQuantiles(knn_dist[,1])
    idxKeep <- which(knn_dist[,1] < max(finalTrajectory[,2]))
    finalTrajectoryAdditional <- data.frame(
        Cells = rownames(dfTrajectory)[idxKeep],
        Distance = knn_dist[idxKeep, 1],
        DistanceIdx = knn_index[idxKeep, 1] + knn_distq[idxKeep]
        )

    #Final Matrix
    finalTrajectory <- rbind(finalTrajectory, finalTrajectoryAdditional)
    finalTrajectory$PseudoTime <- 100*getQuantiles(finalTrajectory[,3])
    finalTrajectory
}













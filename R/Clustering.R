##########################################################################################
# Clustering Methods
##########################################################################################

#' Add cluster information to an ArchRProject
#' 
#' This function will identify clusters from a reduced dimensions object in an ArchRProject or from a supplied reduced dimensions matrix.
#' 
#' @param input Either (i) an `ArchRProject` object containing the dimensionality reduction matrix passed by `reducedDims` or (ii) a dimensionality reduction matrix. This object will be used for cluster identification.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`. Not required if input is a matrix.
#' @param name The column name of the cluster label column to be added to `cellColData` if `input` is an `ArchRProject` object.
#' @param sampleCells QQQ  UNCLEAR. FIX AND ADD LEAVE JJJ FOR ME TO CHECK. An integer specifying the number of cells to subset perform clustering and assign the remainder cells by euclidean distance.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param method A string indicating the clustering method to be used. Supported methods are "Seurat" and "Scran".
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param knnAssign The number of nearest neighbors to be used during clustering for assignment of outliers (clusters with less than nOutlier cells).
#' @param nOutlier QQQ IS THERE EVER A TIME WHEN A CELL DOESNT GET ASSIGNED TO ANY CLUSTER? IF SO, UPDATE THIS. The minimum number of cells required for a group of cells to be called as a cluster. If a group of cells does not reach this threshold, then the cells will be considered outliers and assigned to nearby clusters.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param tstart QQQ THIS IS ODD. WHAT ARE YOU ALLOWED TO MANUALLY SUPPLY A TIME? WHY IS THIS USEFUL? The time at which the function run was started. Useful for keeping track of how long clustering takes relative to a start time. 
#' @param force A boolean value that indicates whether or not to overwrite data in a given column when the value passed to `name` already exists as a column name in `cellColData`.
#' @param ... Additional arguments to be provided to Seurat::FindClusters or scran::buildSNNGraph (for example, knn = 50, jaccard = TRUE)
#' @export
#'
addClusters <- function(
    input = NULL, 
    reducedDims = "IterativeLSI",
    name = "Clusters",
    sampleCells = NULL,
    seed = 1, 
    method = "Seurat", 
    dimsToUse = NULL, 
    corCutOff = 0.75,
    knnAssign = 10, 
    nOutlier = 20, 
    verbose = TRUE,
    tstart = NULL,
    force = FALSE,
    ...
    ){

    .validInput(input = input, name = "input", valid = c("ArchRProj", "matrix"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    .validInput(input = name, name = "name", valid = c("character"))
    .validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
    .validInput(input = seed, name = "seed", valid = c("integer"))
    .validInput(input = method, name = "method", valid = c("character"))
    .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
    .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
    .validInput(input = knnAssign, name = "knnAssign", valid = c("integer"))
    .validInput(input = nOutlier, name = "nOutlier", valid = c("integer"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = tstart, name = "tstart", valid = c("timestamp","null"))
    .validInput(input = force, name = "force", valid = c("boolean"))

    if(is.null(tstart)){
        tstart <- Sys.time()
    }

    if(inherits(input, "ArchRProject")){
        #Check
        input <- addCellColData(ArchRProj = input, data = rep(NA, nCells(input)), name = name, force = force)
        if(reducedDims %ni% names(input@reducedDims)){
            stop("Error reducedDims not available!")
        }
        matDR <- getReducedDims(ArchRProj = input, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    }else if(inherits(input, "matrix")){
        matDR <- input
    }else{
        stop("Input an ArchRProject or Cell by Reduced Dims Matrix!")
    }

    #Subset Matrix
    set.seed(seed)
    nr <- nrow(matDR)

    if(!is.null(sampleCells)){
        if(sampleCells < nrow(matDR)){
            .messageDiffTime("Estimating Clusters by Sampling", tstart, verbose = verbose)
            estimatingClusters <- 1
            idx <- sample(seq_len(nrow(matDR)), sampleCells)
            matDRAll <- matDR
            matDR <- matDR[idx,]
        }else{
            estimatingClusters <- 0
        }
    }else{
        estimatingClusters <- 0
    }


    #################################################################################
    # Decide on which clustering setup to use
    #################################################################################
    if(grepl("seurat",tolower(method))){

        clustParams <- list(...)
        clustParams$verbose <- verbose
        clustParams$tstart <- tstart
        clust <- .clustSeurat(mat = matDR, clustParams = clustParams)

    }else if(grepl("scran",tolower(method))){

        clustParams <- list(...)
        clustParams$x <- matDR
        clustParams$d <- ncol(matDR)
        clustParams$k <- ifelse(!is.null(...$k), ...$k, 25)
        clust <- .clustScran(clustParams)

    }else if(grepl("louvainjaccard",tolower(method))){

        stop("LouvainJaccard method not currently functional!")
        clust <- .clustLouvain(matDR, ...)

    }else{

        stop("Clustering Method Not Recognized!")

    }

    #################################################################################
    # If estimating clsuters we will assign to nearest neighbor cluster
    #################################################################################
    if(estimatingClusters == 1){
        .messageDiffTime("Finding Nearest Clusters", tstart, verbose = verbose)
        knnAssigni <- .computeKNN(matDR, matDRAll[-idx,], knnAssign)
        clustUnique <- unique(clust)
        clustMatch <- match(clust, clustUnique)
        knnAssigni <- apply(knnAssigni, 2, function(x) clustMatch[x])

        .messageDiffTime("Assigning Nearest Clusters", tstart, verbose = verbose)
        clustAssign <- lapply(seq_along(clustUnique), function(x){
            rowSums(knnAssigni == x)
        }) %>% Reduce("cbind", .) %>% apply(., 1, which.max)
        clustOld <- clust
        clust <- rep(NA, nr)
        clust[idx] <- clustOld
        clust[-idx] <- clustUnique[clustAssign]
        matDR <- matDRAll
        remove(matDRAll)
        gc()
    }

    #################################################################################
    # Test if clusters are outliers identified as cells with fewer than nOutlier
    #################################################################################
    .messageDiffTime("Testing Outlier Clusters", tstart, verbose = verbose)
    tabClust <- table(clust)
    clustAssign <- which(tabClust < nOutlier)
    if(length(clustAssign) > 0){
        .messageDiffTime(sprintf("Assigning Outlier Clusters (n = %s, nOutlier < %s cells) to NN", length(clustAssign), nOutlier), tstart, verbose = verbose)
        for(i in seq_along(clustAssign)){
            clusti <- names(clustAssign[i])
            idxi <- which(clust==clusti)
            knni <- .computeKNN(matDR[-idxi,], matDR[idxi,], knnAssign)
            clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
            clust[idxi] <- clustf
        }
    }

    #################################################################################
    # Renaming Clusters based on Proximity in Reduced Dimensions
    #################################################################################
    .reLabel <- function(labels, oldLabels, newLabels){
        labels <- paste0(labels)
        oldLabels <- paste0(oldLabels)
        newLabels <- paste0(newLabels)
        labelsNew <- labels
        for(i in seq_along(oldLabels)){
            labelsNew[labels == oldLabels[i]] <- newLabels[i]
        }
        paste0(labelsNew)
    }
    .messageDiffTime(sprintf("Assigning Cluster Names to %s Clusters", length(unique(clust))), tstart, verbose = verbose)
    meanSVD <- t(.groupMeans(t(matDR), clust))
    meanKNN <- .computeKNN(meanSVD, meanSVD, nrow(meanSVD))
    idx <- sample(seq_len(nrow(meanSVD)), 1)
    clustOld <- c()
    clustNew <- c()
    for(i in seq_len(nrow(meanSVD))){
        clustOld[i] <- rownames(meanSVD)[idx]
        clustNew[i] <- paste0("Cluster", i)
        if(i != nrow(meanSVD)){
            idx <- meanKNN[idx, ][which(rownames(meanSVD)[meanKNN[idx, ]] %ni% clustOld)][1]
        }
    }
    out <- .reLabel(clust, oldLabels = clustOld, newLabels = clustNew)

    if(inherits(input, "ArchRProject")){
        input <- .suppressAll(addCellColData(
                input, 
                data = out, 
                name = name, 
                cells = rownames(matDR),
                force = TRUE
            ))
    }else if(!inherits(input, "ArchRProject")){
        return(out)
    }

}

#Simply a wrapper on Seurats FindClusters
.clustSeurat <- function(mat, clustParams){

    .requirePackage("Seurat")
    .messageDiffTime("Running Seurats FindClusters (Stuart et al. Cell 2019)", clustParams$tstart, verbose=clustParams$verbose)
    set.seed(1)

    #Arxiv Seurat 2.3.4 method
    tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
    colnames(tmp) <- rownames(mat)
    rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))

    obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
    obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
    clustParams$object <- obj
    clustParams$reduction <- "pca"

    obj <- suppressWarnings(do.call(Seurat::FindNeighbors, clustParams))
    clustParams$object <- obj

    cS <- Matrix::colSums(obj@graphs$RNA_snn)

    if(cS[length(cS)] == 1){

        #Error Handling with Singletons
        idxSingles <- which(cS == 1)
        idxNonSingles <- which(cS != 1)

        rn <- rownames(mat) #original order
        mat <- mat[c(idxSingles, idxNonSingles), ,drop = FALSE]

        set.seed(1)

        tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
        colnames(tmp) <- rownames(mat)
        rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))

        obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
        obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
        clustParams$object <- obj
        clustParams$reduction <- "pca"

        obj <- .suppressAll(do.call(Seurat::FindNeighbors, clustParams))
        clustParams$object <- obj

        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))

        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)
        clust <- clust[rn]

    }else{

        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))

        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)

    }

    clust

}

.clustScran <- function(clustParams){
    .requirePackage("scran")
    .requirePackage("igraph")
    #See Scran Vignette!
    set.seed(1)
    #clustParams$x <- matDR
    snn <- do.call(scran::buildSNNGraph, clustParams)
    cluster <- igraph::cluster_walktrap(snn)$membership
    paste0("Cluster", cluster)
}

# #Need to work on making this work
# .clustLouvain <- function(matDR, knn = 50, jaccard = TRUE){

#     getEdges <- function(X, knn, jaccard) {
#         nearest <- RANN::nn2(X, X, k = knn + 1, treetype = "bd", searchtype = "priority")
#         nearest$nn.idx <- nearest$nn.idx[, -1]
#         nearest$nn.dists <- nearest$nn.dists[, -1]
#         nearest$nn.sim <- 1 * (nearest$nn.dists >= 0)
#         edges <- reshape2::melt(t(nearest$nn.idx))
#         colnames(edges) = c("B", "A", "C")
#         edges = edges[, c("A", "B", "C")]
#         edges$B <- edges$C
#         edges$C <- 1
#         edges <- unique(transform(edges, A = pmin(A, B), B = pmax(A, B)))
#         if (jaccard) {
#           message("Calculating Jaccard Distance...")
#             a <- Matrix::tcrossprod(nearest$nn.idx[edges[,1],], nearest$nn.idx[edges[,2],])
#             bi <- Matrix::rowSums(nearest$nn.idx[edges[,1],])
#             bj <- Matrix::rowSums(nearest$nn.idx[edges[,2],])
#             jaccardDist <- a / (rep(ncol(a), bi) + t(rep(nrow(a), bj)) - a)
#             pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
#             #RCPPP?
#             # jaccardDist <- unlist(lapply(seq_len(nrow(edges)), function(x){
#             #     setTxtProgressBar(pb,round(x*100/nrow(edges),0))
#             #     jInt   <- intersect(nearest$nn.idx[edges[x,1],], nearest$nn.idx[edges[x,2],])
#             #     jUnion <- union(nearest$nn.idx[edges[x,1],], nearest$nn.idx[edges[x,2],])
#             #     length(jInt) / length(jUnion)
#             # }))
#             edges$C <- jaccardDist
#             edges <- subset(edges, C != 0)
#             edges$C <- edges$C/max(edges$C)
#         }
#         edges <- Matrix::sparseMatrix(i = edges$A, j = edges$B, x = edges$C, dims = c(nrow(X),nrow(X)), symmetric = TRUE)
#         return(edges)
#     }

#     assignClusters <- function(edges, jaccard) {
#         if (jaccard) {
#             weights <- TRUE
#         }else {
#             weights <- NULL
#         }
#         g <- igraph::graph.adjacency(edges, mode = "undirected", weighted = weights)
#         graphOut <- igraph::cluster_louvain(g)
#         clustAssign <- factor(graphOut$membership, levels = sort(unique(graphOut$membership)))
#         names(clustAssign) <- graphOut$names
#         k = order(table(clustAssign), decreasing = TRUE)
#         newLevels <- rep(1, length(unique(graphOut$membership)))
#         newLevels[k] <- seq_len(length(unique(graphOut$membership)))
#         levels(clustAssign) <- newLevels
#         clustAssign <- factor(clustAssign, levels = seq_len(length(unique(graphOut$membership))))
#         return(paste0("Cluster", clustAssign))
#     }

#     require(RANN)
#     require(cluster)
#     require(igraph)
#     require(Matrix)
#     message("Running Louvian Jaccard Graph Clustering...")
#     message("Adapted from Comprehensive Classification of Retinal Bipolar Neurons by Single-Cell Transcriptomics. Cell 2016.")
    
#     message("Calculating Edges...")
#     edges <- getEdges(X = matDR, knn = knn, jaccard = jaccard)
#     message("\nAssigning Clusters...")
#     clustAssign <- assignClusters(edges = edges, jaccard = jaccard)
    
#     return(clustAssign)

# }

#' @export
.computeKNN <- function(data = NULL, query = NULL, k = 50, method = NULL, includeSelf = FALSE, ...){

  .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = method, name = "method", valid = c("character", "null"))
  .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))

  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }

  if(is.null(method)){
    if(requireNamespace("nabor", quietly = TRUE)){
        method <- "nabor"
    }else if(requireNamespace("RANN", quietly = TRUE)){
        method <- "RANN"
    }else if(requireNamespace("FNN", quietly = TRUE)){
        method <- "FNN"
    }else{
        stop("Computing KNN requires package nabor, RANN or FNN")
    }
  }

  if(tolower(method)=="nabor"){
    
    .requirePackage("nabor")
    if(searchSelf & !includeSelf){
      knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
      knnIdx <- knnIdx[,-1]
    }else{
      knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
    }
  
  }else if(tolower(method)=="rann"){
    
    .requirePackage("RANN")
    if(searchSelf & !includeSelf){
      knnIdx <- RANN::nn2(data = data, query = query, k = k + 1, ...)$nn.idx
      knnIdx <- knnIdx[,-1]
    }else{
      knnIdx <- RANN::nn2(data = data, query = query, k = k, ...)$nn.idx
    }

  }else if(tolower(method)=="fnn"){

    .requirePackage("FNN")
    if(searchSelf & !includeSelf){
      knnIdx <- FNN::get.knnx(data = data, query = query, k = k + 1, ...)$nn.index
      knnIdx <- knnIdx[,-1]
    }else{
      knnIdx <- FNN::get.knnx(data = data, query = query, k = k, ...)$nn.index
    }

  }else{

    stop(sprintf("KNN Method %s not Recognized!", method))

  }

  knnIdx

}


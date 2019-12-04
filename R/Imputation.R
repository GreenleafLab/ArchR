#' Add TileMatrix to Arrows/ArchRProject
#' 
#' This function for each sample will independently compute counts for each tile
#' per cell in the Arrow File
#'
#' @param input ArchRProject or ArrowFiles
#' @param chromSizes chromomosome sizes used for identifying number of tiles to count
#' @param windowSize size for each window to break up each chromosome
#' @param binarize save as a Sparse.Binary.Matrix or Sparse.Integer.Matrix
#' @param excludeChr exclude chromosomes from this analysis
#' @param threads number of threads
#' @param parallelParam parallel parameters for batch style execution
#' @param force force overwriting previous TileMatrix in ArrowFile
#' @export
addMagicWeights <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  dimsToUse = NULL, 
  td = 3,
  ka = 4,
  sampleCells = max(5000, floor(nCells(ArchRProj) / 10)),
  k = 15,
  weighted = TRUE
  ){

  #Adapted From
  #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R

  set.seed(1)

  tstart <- Sys.time()

  #Get Reduced Dims
  if(!is.null(dimsToUse)){
      matDR <- getReducedDims(ArchRProj, reducedDims = reducedDims)[, dimsToUse, drop = FALSE]
  }else{
      matDR <- getReducedDims(ArchRProj, reducedDims = reducedDims)
  }
  N <- nrow(matDR)
  rn <- rownames(matDR)

  idx <- sample(seq_len(nrow(matDR)), nrow(matDR))

  if(is.null(sampleCells)){
    sampleCells <- N
  }

  cutoffs <- lapply(seq_len(1000), function(x){
    N / x
  }) %>% unlist
  binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)

  groups <- split(idx, ceiling(seq_along(idx)/binSize))

  Wt <- lapply(seq_along(groups), function(x){

    .messageDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)", x, length(groups)), tstart)

    ix <- groups[[x]]
    Nx <- length(ix)

    #Compute KNN
    if(requireNamespace("nabor", quietly = TRUE)) {
        knnObj <- nabor::knn(data = matDR[ix,], query = matDR[ix, ], k = k)
        knnIdx <- knnObj$nn.idx
        knnDist <- knnObj$nn.dists
    }else if(requireNamespace("RANN", quietly = TRUE)) {
        knnObj <- RANN::nn2(data = matDR[ix,], query = matDR[ix, ], k = k)
        knnIdx <- knnObj$nn.idx
        knnDist <- knnObj$nn.dists
    }else if(requireNamespace("FNN", quietly = TRUE)) {
        knnObj <- FNN::get.knnx(data = matDR[ix,], query = matDR[ix, ], k = k)
        knnIdx <- knnObj$nn.index
        knnDist <-knnObj$nn.dist
    }else{
        stop("Computing KNN requires package nabor, RANN or FNN")
    }
    rm(knnObj)

    if(ka > 0){
      knnObj$dist <- knnObj$dist / knnObj$dist[,ka]
    }

    if (weighted) {
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
    } else {
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=1, dims = c(Nx, Nx)) # unweighted kNN graph
    }

    W <- W + t(W)

    #Compute Kernel
    if(epsilon > 0){
      W@x <- exp(-(W@x / epsilon^2))
    }

    #Markov normalization
    W <- W / Matrix::rowSums(W) 

    #Initialize Matrix
    Wt <- W

    #Computing Diffusion Matrix
    for(i in seq_len(td)){
      #message(i, " ", appendLF = FALSE)
        Wt <- Wt %*% W
    }
    #message("\n", appendLF = FALSE)

    Wt <- Matrix::summary(Wt)
    Wt[,1] <- ix[Wt[,1]]
    Wt[,2] <- ix[Wt[,2]]

    rm(knnObj)
    rm(W)
    gc()

    Wt

  }) %>% Reduce("rbind", .) %>% {Matrix::sparseMatrix(i = .[,1], j = .[,2], x = .[,3], dims = c(N, N))}

  .messageDiffTime(sprintf("Completed Getting Magic Weights (Size = %s GB)!", round(object.size(Wt) / 10^9, 3)), tstart)

  rownames(Wt) <- rownames(matDR)
  colnames(Wt) <- rownames(matDR)

  #ArchRProj@imputeWeights <- Wt
  
  #ArchRProj

  Wt

}

#' Get outputDirectory in ArchRProject
#' 
#' This function gets outputDirectory from ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getImputeWeights <- function(ArchRProj, ...){
  
  ArchRProj <- .validArchRProject(ArchRProj)
  iW <- ArchRProj@imputeWeights
  
  if(!is.null(iW)){
    iW <- iW[rownames(getCellColData(ArchRProj)),rownames(getCellColData(ArchRProj))]
  }else{
    iW <- NULL
  
  }
  
  iW

}














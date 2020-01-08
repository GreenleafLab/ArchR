##########################################################################################
# Imputation Methods
##########################################################################################

#' QQQ
#' 
#' This function QQQ
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims QQQ The name of the `reducedDims` object to retrieve from the designated `ArchRProject`. Options include QQQ. QQQ Not required if input is a matrix.
#' @param dimsToUse QQQ A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param td diffusion time (number of iterations) for Magic
#' @param ka kNN autotune parameter for Magic
#' @param sampleCells number of cells to sample per block of estimated imputation matrix
#' @param k number of nearest neighbors to use for Magic
#' @param epsilon a value for the standard deviation of the kernel for Magic
#' @param ... additional params
#' @export
addImputeWeights <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  dimsToUse = NULL, 
  td = 3,
  ka = 4,
  sampleCells = max(5000, floor(nCells(ArchRProj) / 10)),
  k = 15,
  epsilon = 1,
  ...
  ){

  #Adapted From
  #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R

  set.seed(1)

  tstart <- Sys.time()
  .messageDiffTime("Computing Impute Weights Using Magic (Cell 2018)", tstart)

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

  blocks <- split(idx, ceiling(seq_along(idx)/binSize))

  Wt <- lapply(seq_along(blocks), function(x){

    .messageDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)", x, length(blocks)), tstart)

    ix <- blocks[[x]]
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
      knnDist <- knnDist / knnDist[,ka]
    }

    if(epsilon > 0){
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
    } else {
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=1, dims = c(Nx, Nx)) # unweighted kNN graph
    }
    W <- W + Matrix::t(W)

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

    rm(knnIdx)
    rm(knnDist)
    rm(W)
    gc()

    Wt

  }) %>% Reduce("rbind", .) %>% {Matrix::sparseMatrix(i = .[,1], j = .[,2], x = .[,3], dims = c(N, N))}

  .messageDiffTime(sprintf("Completed Getting Magic Weights (Size = %s GB)!", round(object.size(Wt) / 10^9, 3)), tstart)

  rownames(Wt) <- rownames(matDR)
  colnames(Wt) <- rownames(matDR)

  ArchRProj@imputeWeights <- SimpleList(
    Weights = Wt, 
    Blocks = blocks, 
    Params = 
      list(
        reducedDims = reducedDims, 
        td = td, 
        k = k, 
        ka = ka,
        epsilon = epsilon,
        weighted = weighted
        )
      )
  
  ArchRProj

}

#' QQQ
#' 
#' This function QQQ
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param ... additional args
#' @export
getImputeWeights <- function(ArchRProj, ...){  
  ArchRProj <- .validArchRProject(ArchRProj)
  ArchRProj@imputeWeights
}



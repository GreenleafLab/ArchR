##########################################################################################
# Imputation Methods
##########################################################################################

#' Add Imputation Weights to an ArchRProject
#' 
#' This function computes imputations weights that describe each cell as a linear combination of many cells based on a MAGIC diffusion matrix.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded.
#' @param td QQQ IN GENERAL, I FIND THE MAGIC PARAM EXPLANATIONS USELESS. NOBODY COULD INTERPRET THESE WITHOUT GOING TO MAGIC FIRST. The diffusion time (number of iterations) for MAGIC.
#' @param ka The k-nearest neighbors autotune parameter for MAGIC.
#' @param sampleCells QQQ UNCLEAR WHAT A BLOCK IS. The number of cells to sample per block of the estimated imputation matrix.
#' @param k The number of nearest neighbors to use for MAGIC.
#' @param epsilon QQQ UNCLEAR. IS THIS LIKE A MINIMUM STD DEV FOR A FILTER? The value for the standard deviation of the kernel for MAGIC.
#' @param ... additional params
#' @export
addImputeWeights <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  dimsToUse = NULL, 
  corCutOff = 0.75, 
  td = 3,
  ka = 4,
  sampleCells = max(5000, floor(nCells(ArchRProj) / 10)),
  k = 15,
  epsilon = 1,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = td, name = "td", valid = c("integer"))
  .validInput(input = ka, name = "ka", valid = c("integer"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("integer"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = epsilon, name = "epsilon", valid = c("numeric"))

  #Adapted From
  #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R

  set.seed(1)

  tstart <- Sys.time()
  .messageDiffTime("Computing Impute Weights Using Magic (Cell 2018)", tstart)

  #Get Reduced Dims
  matDR <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
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
        epsilon = epsilon
        )
      )
  
  ArchRProj

}

#' Get Imputation Weights from ArchRProject
#' 
#' This function gets imputation weights from an ArchRProject to impute numeric values.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param ... additional args
#' @export
getImputeWeights <- function(ArchRProj = NULL, ...){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchRProj@imputeWeights
}



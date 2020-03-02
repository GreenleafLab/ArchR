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
#' @param scaleDims A boolean that indicates whether to z-score the reduced dimensions for each cell. This is useful forminimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded.
#' @param td The diffusion time parameter determines the number of smoothing iterations to be performed (see MAGIC from van Dijk et al Cell 2018).
#' @param ka The k-nearest neighbors autotune parameter to equalize the effective number of neighbors for each cell, thereby diminishing
#' the effect of differences in density. (see MAGIC from van Dijk et al Cell 2018).
#' @param sampleCells The number of cells to sub-sample to compute an imputation block. An imputation block is a cell x cell matrix that
#' describes the linear combination for imputation for numerical values within these cells. ArchR creates many blocks to keep this
#' cell x cell matrix sparse for memory concerns.
#' @param k The number of nearest neighbors for smoothing to use for MAGIC (see MAGIC from van Dijk et al Cell 2018).
#' @param epsilon The value for the standard deviation of the kernel for MAGIC (see MAGIC from van Dijk et al Cell 2018).
#' @export
addImputeWeights <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  dimsToUse = NULL, 
  scaleDims = NULL,
  corCutOff = 0.75, 
  td = 3,
  ka = 4,
  sampleCells = max(5000, floor(nCells(ArchRProj) / 50)),
  k = 15,
  epsilon = 1
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
#' @export
getImputeWeights <- function(ArchRProj = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  if(length(ArchRProj@imputeWeights) == 0){
    return(NULL)
  }
  ArchRProj@imputeWeights
}

#' @export
addImputeWeights2 <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  dimsToUse = NULL, 
  scaleDims = NULL,
  corCutOff = 0.75, 
  td = 3,
  ka = 4,
  sampleCells = 5000,
  nRep = 2,
  k = 15,
  epsilon = 1,
  useHdf5 = TRUE,
  randomSuffix = FALSE
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

  if(is.null(sampleCells)){
    sampleCells <- N
  }

  cutoffs <- lapply(seq_len(1000), function(x){
    N / x
  }) %>% unlist
  binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)

  if(useHdf5){
    dir.create(file.path(getOutputDirectory(ArchRProj), "ImputeWeights"), showWarnings = FALSE)
    if(randomSuffix){
      randomString <- ArchR:::.randomStr(n = 1)
      randomString <- paste0(randomString, "-", seq_len(nRep))
      weightFiles <- file.path(getOutputDirectory(ArchRProj), "ImputeWeights", paste0("Impute-Weights-Rep-", randomString))
    }else{
      weightFiles <- file.path(getOutputDirectory(ArchRProj), "ImputeWeights", paste0("Impute-Weights-Rep-", seq_len(nRep)))
    }
  }

  weightList <- ArchR:::.safelapply(seq_len(nRep), function(y){

    .messageDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)", y, nRep), tstart)

    idx <- sample(seq_len(nrow(matDR)), nrow(matDR))

    blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx)/binSize))

    weightFile <- weightFiles[y]

    if(useHdf5){
      o <- h5createFile(weightFile)
    }

    blockList <- lapply(seq_along(blocks), function(x){

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
          Wt <- Wt %*% W
      }
      rownames(Wt) <- rownames(matDR)[ix]
      colnames(Wt) <- rownames(matDR)[ix]

      rm(knnIdx)
      rm(knnDist)
      rm(W)
      gc()

      if(useHdf5){
        o <- ArchR:::.suppressAll(h5createGroup(file = weightFile, paste0("block", x)))      
        o <- ArchR:::.suppressAll(h5write(obj = ix, file = weightFile, name = paste0("block", x, "/Names"), level = 0))
        o <- ArchR:::.suppressAll(h5write(obj = as.matrix(Wt), file = weightFile, name = paste0("block", x, "/Weights"), level = 0))
        return(weightFile)
      }else{
        Wt
      }

    }) %>% SimpleList

    names(blockList) <- paste0("b",seq_along(blockList))

    blockList

  }) %>% SimpleList
  names(weightList) <- paste0("w",seq_along(weightList))

  .messageDiffTime(sprintf("Completed Getting Magic Weights!", round(object.size(weightList) / 10^9, 3)), tstart)

  ArchRProj@imputeWeights <- SimpleList(
    Weights = weightList, 
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

#' @export
imputeMatrix <- function(
  mat = NULL, 
  imputeWeights = NULL,
  threads = getArchRThreads()
  ){

  weightList <- imputeWeights$Weights

  if(inherits(weightLis, "dgCMatrix")){
    as.matrix(mat) %*% as.matrix(weightList)
  }

  if(!inherits(weightList, "SimpleList") & !inherits(weightList, "list")){
    stop("Weights are not a list, Please re-run addImputeWeights (update)!")
  }

  start <- Sys.time()
  
  imputeMat <- lapply(seq_along(weightList), function(x){
    
    ArchR:::.messageDiffTime(sprintf("Imputing Matrix (%s of %s)", x, length(weightList)), start)

    if(is.character(weightList[[x]])){
      if(!file.exists(weightList[[x]])){
        message("Weight File Does Not Exist! Please re-run addImputeWeights!")
      }

      h5df <- h5ls(weightList[[x]])
      blocks <- gtools::mixedsort(grep("block",h5df$name,value=TRUE))
      matx <- ArchR:::.safelapply(seq_along(blocks), function(y){

        message(y, " ", appendLF = FALSE)

        #Read In Weights and Names
        bn <- h5read(weightList[[x]], paste0(blocks[y], "/Names"))
        by <- h5read(weightList[[x]], paste0(blocks[y], "/Weights"))
        colnames(by) <- bn

        #Multiply
        as.matrix(mat[, paste0(bn)]) %*% by
      
      }, threads = threads) %>% Reduce("cbind", .)

      message("")

    }else{

      matx <- ArchR:::.safelapply(seq_along(weightList[[x]]), function(y){
        
        message(y, " ", appendLF = FALSE)

        as.matrix(mat[, colnames(weightList[[x]][[y]])]) %*% as.matrix(weightList[[x]][[y]])

      }, threads = threads) %>% Reduce("cbind", .)

      message("")

    }

    matx[, colnames(mat)] #Return Ordered
  
  }) %>% Reduce("+", .)

  Sys.time() - start

  #Compute Average
  imputeMat <- imputeMat / length(weightList)

  ArchR:::.messageDiffTime("Finished Imputing Matrix", start)

  imputeMat

}





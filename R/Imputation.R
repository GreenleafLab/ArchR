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
#' @param nRep An integer representing the number of imputation replicates to create when downsampling extremely low.
#' @param k The number of nearest neighbors for smoothing to use for MAGIC (see MAGIC from van Dijk et al Cell 2018).
#' @param epsilon The value for the standard deviation of the kernel for MAGIC (see MAGIC from van Dijk et al Cell 2018).
#' @param useHdf5 A boolean value that indicates whether HDF5 format should be used to store the impute weights.
#' @param randomSuffix A boolean value that indicates whether a random suffix should be appended to the saved imputation weights hdf5 files.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addImputeWeights <- function(
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
  randomSuffix = FALSE,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  logFile = createLogFile("addImputeWeights")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = td, name = "td", valid = c("integer"))
  .validInput(input = ka, name = "ka", valid = c("integer"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
  .validInput(input = nRep, name = "nRep", valid = c("integer"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = epsilon, name = "epsilon", valid = c("numeric"))
  .validInput(input = useHdf5, name = "useHdf5", valid = c("boolean"))
  .validInput(input = randomSuffix, name = "randomSuffix", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addImputeWeights Input-Parameters", logFile = logFile)  

  #Adapted From
  #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R

  set.seed(seed)

  tstart <- Sys.time()
  .logDiffTime("Computing Impute Weights Using Magic (Cell 2018)", t1 = tstart, verbose = verbose, logFile = logFile)

  #Get Reduced Dims
  matDR <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
  N <- nrow(matDR)
  rn <- rownames(matDR)

  if(!is.null(sampleCells)){
    if(sampleCells > nrow(matDR)){
      sampleCells <- NULL
    }
  }
  if(is.null(sampleCells)){
    binSize <- N
    nRep <- 1
  }else{
    cutoffs <- lapply(seq_len(1000), function(x){
      N / x
    }) %>% unlist
    binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)
  }

  if(useHdf5){
    dir.create(file.path(getOutputDirectory(ArchRProj), "ImputeWeights"), showWarnings = FALSE)
    if(randomSuffix){
      weightFiles <- .tempfile("Impute-Weights", tmpdir = file.path(gsub(paste0(getwd(),"/"),"",getOutputDirectory(ArchRProj)), "ImputeWeights"))
      weightFiles <- paste0(weightFiles, "-Rep-", seq_len(nRep))
    }else{
      weightFiles <- file.path(getOutputDirectory(ArchRProj), "ImputeWeights", paste0("Impute-Weights-Rep-", seq_len(nRep)))
    }
    o <- suppressWarnings(file.remove(weightFiles))
  }

  weightList <- .safelapply(seq_len(nRep), function(y){

    .logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)", y, nRep), t1 = tstart, verbose = FALSE, logFile = logFile)

    if(!is.null(sampleCells)){
      idx <- sample(seq_len(nrow(matDR)), nrow(matDR))
      blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx)/binSize))
    }else{
      blocks <- list(rownames(matDR)) 
    }

    if(useHdf5){
      weightFile <- weightFiles[y]
      o <- h5createFile(weightFile)
    }

    blockList <- lapply(seq_along(blocks), function(x){

      if(x %% 10 == 0){
        .logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s, Iteration %s of %s)", y, nRep, x, length(blocks)), 
          t1 = tstart, verbose = FALSE, logFile = logFile)
      }

      ix <- blocks[[x]]
      Nx <- length(ix)

      #Compute KNN
      knnObj <- nabor::knn(data = matDR[ix,], query = matDR[ix, ], k = k)
      knnIdx <- knnObj$nn.idx
      knnDist <- knnObj$nn.dists
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
      rownames(Wt) <- ix
      colnames(Wt) <- ix

      rm(knnIdx)
      rm(knnDist)
      rm(W)
      gc()

      if(useHdf5){
        o <- .suppressAll(h5createGroup(file = weightFile, paste0("block", x)))      
        o <- .suppressAll(h5write(obj = ix, file = weightFile, name = paste0("block", x, "/Names"), level = 0))
        o <- .suppressAll(h5write(obj = as.matrix(Wt), file = weightFile, name = paste0("block", x, "/Weights"), level = 0))
        return(weightFile)
      }else{
        Wt
      }

    }) %>% SimpleList

    if(useHdf5){
      return(weightFile)
    }else{
      names(blockList) <- paste0("b",seq_along(blockList))
      blockList
    }
  }, threads = threads) %>% SimpleList
  names(weightList) <- paste0("w",seq_along(weightList))

  .logDiffTime(sprintf("Completed Getting Magic Weights! Object size - %s.", round(object.size(weightList) / 10^9, 3)), 
    t1 = tstart, verbose = FALSE, logFile = logFile)

  ArchRProj@imputeWeights <- SimpleList(
    Weights = weightList, 
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
  message("Getting ImputeWeights")
  if(length(ArchRProj@imputeWeights) == 0){
    message("No imputeWeights found, returning NULL")
    return(NULL)
  }
  ArchRProj@imputeWeights
}

#' Impute a matrix with impute weights
#' 
#' This function gets imputation weights from an ArchRProject to impute a numerical matrix
#' 
#' @param mat A matrix or sparseMatrix of class dgCMatrix to be imputed.
#' @param imputeWeights An R object containing impute weights as returned by `getImputeWeights(ArchRProj)`. See `addImputeWeights()` for more details.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
imputeMatrix <- function(
  mat = NULL, 
  imputeWeights = NULL,
  threads = getArchRThreads(),
  verbose = FALSE,
  logFile = createLogFile("imputeMatrix")
  ){

  .validInput(input = mat, name = "mat", valid = c("matrix", "sparseMatrix"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  if(!inherits(imputeWeights$Weights, "SimpleList") & !inherits(imputeWeights$Weights, "list")){
    .logMessage("Weights are not a list, Please re-run addImputeWeights (update)!", logFile = logFile)
    stop("Weights are not a list, Please re-run addImputeWeights (update)!")
  }

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "imputeMatrix Input-Parameters", logFile = logFile)  

  weightList <- imputeWeights$Weights
  .logThis(mat, "mat", logFile = logFile)
  .logThis(weightList, "weightList", logFile = logFile)

  tstart <- Sys.time()
  
  imputeMat <- lapply(seq_along(weightList), function(x){
    
    .logDiffTime(sprintf("Imputing Matrix (%s of %s)", x, length(weightList)), tstart, verbose = verbose, logFile = logFile)

    if(is.character(weightList[[x]])){

      .logMessage("Using weights on disk", logFile = logFile)

      if(!file.exists(weightList[[x]])){
        .logMessage("Weight File Does Not Exist! Please re-run addImputeWeights!", logFile = logFile)
        stop("Weight File Does Not Exist! Please re-run addImputeWeights!")
      }

      h5df <- h5ls(weightList[[x]])
      blocks <- gtools::mixedsort(grep("block",h5df$name,value=TRUE))
      matx <- .safelapply(seq_along(blocks), function(y){

        if(verbose) message(y, " ", appendLF = FALSE)
        .logMessage(paste0(y, " of ", length(blocks)),  logFile = logFile)

        #Read In Weights and Names
        bn <- h5read(weightList[[x]], paste0(blocks[y], "/Names"))
        by <- h5read(weightList[[x]], paste0(blocks[y], "/Weights"))
        colnames(by) <- bn
        rownames(by) <- bn

        #Multiply
        if(!all(paste0(bn) %in% colnames(mat))){
          .logThis(paste0(bn), "Block cellNames", logFile = logFile)
          .logThis(colnames(mat), "Matrix cellNames", logFile = logFile)
          .logThis(paste0(bn)[paste0(bn) %ni% colnames(mat)], "Block cellNames not in matrix", logFile = logFile)
          .logMessage("Not all cellNames from imputeWeights are present. If you subsetted cells from the original imputation, please re-run with addImputeWeights!")
          stop("Not all cellNames from imputeWeights are present. If you subsetted cells from the original imputation, please re-run with addImputeWeights!")
        }

        Matrix::t(by %*% Matrix::t(mat[, paste0(bn), drop = FALSE]))
      
      }, threads = threads) %>% Reduce("cbind", .)

      if(verbose) message("")

    }else{

      .logMessage("Using weights in memory", logFile = logFile)

      matx <- .safelapply(seq_along(weightList[[x]]), function(y){
        
        if(verbose) message(y, " ", appendLF = FALSE)
        .logMessage(paste0(y, " of ", length(weightList[[x]])), logFile = logFile)
        Matrix::t(as.matrix(weightList[[x]][[y]]) %*% Matrix::t(mat[, paste0(colnames(weightList[[x]][[y]])), drop = FALSE]))

      }, threads = threads) %>% Reduce("cbind", .)

      if(verbose) message("")

    }

    matx[, colnames(mat)] #Return Ordered
  
  }) %>% Reduce("+", .)

  #Compute Average
  imputeMat <- imputeMat / length(weightList)

  .logDiffTime("Finished Imputing Matrix", tstart, verbose = verbose, logFile = logFile)

  imputeMat

}





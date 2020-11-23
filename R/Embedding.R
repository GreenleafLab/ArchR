##########################################################################################
# Embedding Methods
##########################################################################################

#' Add a UMAP embedding of a reduced dimensions object to an ArchRProject
#' 
#' This function will compute a UMAP embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param name The name for the UMAP embedding to store in the given `ArchRProject` object.
#' @param nNeighbors An integer describing the number of nearest neighbors to compute a UMAP. This argument is passed to `n_neighbors` in `uwot::umap()`.
#' @param minDist A number that determines how tightly the UMAP is allowed to pack points together. This argument is passed to `min_dist` in
#' `uwot::umap()`. For more info on this see https://jlmelville.github.io/uwot/abparams.html.
#' @param metric A number that determines how distance is computed in the `reducedDims` to compute a UMAP. This argument is passed to `metric` in `uwot::umap()`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param sampleCells An integer specifying the number of cells to subsample and perform UMAP Embedding on. The remaining cells
#' that were not subsampled will be re-projected using uwot::umap_transform to the UMAP Embedding. This enables a decrease in run time
#' and memory but can lower the overal quality of the UMAP Embedding. Only recommended for extremely large number of cells.
#' @param outlierQuantile A numeric (0 to 1) describing the distance quantile in the subsampled cels (see `sampleCells`) to use to filter poor quality re-projections.
#' This is necessary because there are lots of outliers if undersampled significantly.
#' @param saveModel A boolean value indicating whether or not to save the UMAP model in an RDS file for downstream usage such as projection of data into the UMAP embedding.
#' @param verbose A boolean value that indicates whether printing UMAP output.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param force A boolean value that indicates whether to overwrite the relevant data in the `ArchRProject` object if the embedding indicated by
#' `name` already exists.
#' @param threads The number of threads to be used for parallel computing. Default set to 1 because if set to high can cause C stack usage errors.
#' @param ... Additional parameters to pass to `uwot::umap()`
#' @export
addUMAP <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  sampleCells = NULL,
  outlierQuantile = 0.9,
  saveModel = TRUE,
  verbose = TRUE,
  seed = 1,
  force = FALSE,
  threads = 1,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  .validInput(input = nNeighbors, name = "nNeighbors", valid = c("integer", "null"))
  .validInput(input = minDist, name = "minDist", valid = c("numeric", "null"))
  .validInput(input = metric, name = "metric", valid = c("character", "null"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
  .validInput(input = outlierQuantile, name = "outlierQuantile", valid = c("numeric"))
  .validInput(input = saveModel, name = "saveModel", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  if(name %in% names(ArchRProj@embeddings)){
    if(!force){
      stop("Embedding Already Exists! Either set force = TRUE or use a different name!")
    }
  }

  #############################################################################################
  # Default Parameters for Input Embeddings!
  #############################################################################################

  #Merge Parameters
  embeddingParams <- list(...)
  embeddingParams$X <- getReducedDims(
      ArchRProj = ArchRProj, 
      reducedDims = reducedDims, 
      dimsToUse = dimsToUse, 
      corCutOff = corCutOff, 
      scaleDims = scaleDims
  )
  embeddingParams$n_neighbors <- nNeighbors
  embeddingParams$min_dist <- minDist
  embeddingParams$verbose <- verbose
  embeddingParams$metric <- metric

  estimateUMAP <- FALSE
  projectDF <- DataFrame(row.names = rownames(embeddingParams$X), projected = rep(0, nrow(embeddingParams$X))) #Projection ID
  if(!is.null(sampleCells)){
    if(sampleCells < nrow(embeddingParams$X)){
      message("Creating an Estimated UMAP by sub-sampling cells N = ", sampleCells, "!")
      saveModel <- TRUE
      idx <- sample(seq_len(nrow(embeddingParams$X)), sampleCells)
      cellNames <- rownames(embeddingParams$X)
      saveX <- embeddingParams$X[-idx, , drop = FALSE]
      embeddingParams$X <- embeddingParams$X[idx, , drop = FALSE]
      estimateUMAP <- TRUE
      projectDF[idx, 1] <- 1
    }
  }

  if(saveModel){
    embeddingParams$ret_nn <- TRUE
    embeddingParams$ret_model <- TRUE 
  }else{
    embeddingParams$ret_nn <- FALSE
    embeddingParams$ret_model <- FALSE      
  }

  #############################################################################################
  # Run Embedding
  #############################################################################################
  #Seed
  set.seed(seed)
  uwot_umap <- do.call(uwot::umap, embeddingParams)

  if(estimateUMAP){
    uwot_umap2 <- uwot::umap_transform(X = saveX, model = uwot_umap, n_threads = as.integer(threads), verbose = verbose)
    #We should check the distances
    knnRef <- as.vector(nabor::knn(data = uwot_umap[[1]], query = uwot_umap[[1]], k = 2)$nn.dists[,-1])
    knnProj <- as.vector(nabor::knn(data = uwot_umap[[1]], query = uwot_umap2, k = 1)$nn.dists)
    idxExclude <- which(knnProj >= quantile(knnRef, outlierQuantile))
    uwot_umap2[idxExclude, ] <- NA
  }

  #############################################################################################
  # Add Embedding to Project
  #############################################################################################
  nc <- ncol(embeddingParams$X)
  nr <- nrow(embeddingParams$X)

  if(saveModel){
    dir.create(file.path(getOutputDirectory(ArchRProj), "Embeddings"), showWarnings = FALSE)
    modelFile <- .tempfile(
      pattern = paste0("Save-Uwot-UMAP-Params-",reducedDims), 
      tmpdir = file.path(getOutputDirectory(ArchRProj), "Embeddings"),
      fileext = ".tar",
      addDOC = TRUE
    )
    #file.path(getOutputDirectory(ArchRProj), "Embeddings", paste0("Save-Uwot-UMAP-Params-",reducedDims,"-",.randomStr(),".tar"))
    saveModelTmp <- .saveUWOT(uwot_umap, modelFile)
    if(!file.exists(modelFile)){
      warning("Model was not saved properly, continuing without saving model!")
      modelFile <- NA
    }
    dfEmbedding <- data.frame(uwot_umap[[1]])
    colnames(dfEmbedding) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)
    embeddingParams$X <- NULL

    if(estimateUMAP){
      dfEmbedding2 <- data.frame(uwot_umap2)
      colnames(dfEmbedding2) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding2)))
      rownames(dfEmbedding2) <- rownames(saveX)
      rm(uwot_umap2)
      dfEmbedding <- rbind(dfEmbedding, dfEmbedding2)
      dfEmbedding <- dfEmbedding[cellNames,,drop=FALSE]
    }

    ArchRProj@embeddings[[name]] <- SimpleList(
      df = dfEmbedding, 
      params = c(
        embeddingParams,
        dimsToUse = dimsToUse,
        scaleDims = scaleDims,
        corCutOff = corCutOff,
        nr=nr,
        nc=nc,
        uwotModel = modelFile,
        estimateUMAP = estimateUMAP,
        projectID = projectDF
      )
    )
  }else{
    dfEmbedding <- data.frame(uwot_umap)    
    colnames(dfEmbedding) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)
    embeddingParams$X <- NULL
    ArchRProj@embeddings[[name]] <- SimpleList(
      df = dfEmbedding, 
      params = c(
        embeddingParams,
        dimsToUse = dimsToUse,
        scaleDims = scaleDims,
        corCutOff = corCutOff,
        nr=nr,
        nc=nc,
        uwotModel = NA,
        estimateUMAP = estimateUMAP,
        projectID = projectDF
      )
    )
  }   

  return(ArchRProj)

}

#New Save UWOT
.saveUWOT <- function(model, file){
  tryCatch({
    uwot::save_uwot(model = model, file = file, verbose = TRUE)
  }, error = function(e){
    .saveUWOT_Deprecated(model = model, file = file) #backwards to previous version
  })
}

#save_uwot does not work because tarring doesnt work for some reason on Stanford's compute server
#Adapted from save_uwot
.saveUWOT_Deprecated <- function(model, file){
  file <- file.path(normalizePath(dirname(file)), basename(file))
  wd <- getwd()
  mod_dir <- tempfile(pattern = "dir")
  dir.create(mod_dir)
  uwot_dir <- file.path(mod_dir, "uwot")
  dir.create(uwot_dir)
  model_tmpfname <- file.path(uwot_dir, "model")
  .safeSaveRDS(model, file = model_tmpfname)
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in seq_len(n_metrics)) {
      nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
      if (n_metrics == 1) {
          model$nn_index$save(nn_tmpfname)
          model$nn_index$unload()
          model$nn_index$load(nn_tmpfname)
      }
      else {
          model$nn_index[[i]]$save(nn_tmpfname)
          model$nn_index[[i]]$unload()
          model$nn_index[[i]]$load(nn_tmpfname)
      }
  }
  setwd(mod_dir)
  system2("tar", "-cvf uwot.tar uwot", stdout = NULL, stderr = NULL)
  o <- .fileRename("uwot.tar", file)
  setwd(wd)
  if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
  }
  return(o)
}

#New Save UWOT
.loadUWOT <- function(file){
  tryCatch({
    uwot::load_uwot(file = file, verbose = TRUE)
  }, error = function(e){
    .loadUWOT_Deprecated(file = file, nDim = nDim) #backwards to previous version
  })
}

#Adapted from load_uwot
.loadUWOT_Deprecated <- function(file, nDim = NULL){
    model <- NULL
    tryCatch({
        mod_dir <- tempfile(pattern = "dir")
        dir.create(mod_dir)
        utils::untar(file, exdir = mod_dir)
        model_fname <- file.path(mod_dir, "uwot/model")
        if (!file.exists(model_fname)) {
            stop("Can't find model in ", file)
        }
        model <- readRDS(file = model_fname)
        metrics <- names(model$metric)
        n_metrics <- length(metrics)
        for (i in seq_len(n_metrics)){
            nn_fname <- file.path(mod_dir, paste0("uwot/nn", i))
            if (!file.exists(nn_fname)) {
                stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
            }
            metric <- metrics[[i]]
            if(length(model$metric[[i]]) == 0){
              if(!is.null(nDim)){
                nDim2 <- nDim
              }else{
                nDim2 <- length(model$metric[[i]])
              }
            }
            if(!is.null(nDim)){
              nDim2 <- nDim
            }
            ann <- uwot:::create_ann(metric, ndim = nDim2)
            ann$load(nn_fname)
            if (n_metrics == 1) {
                model$nn_index <- ann
            }else{
                model$nn_index[[i]] <- ann
            }
        }
    }, finally = {
        if (file.exists(mod_dir)) {
            unlink(mod_dir, recursive = TRUE)
        }
    })
    model 
}

#' Add a TSNE embedding of a reduced dimensions object to an ArchRProject
#' 
#' This function will compute a TSNE embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param method The method for computing a TSNE embedding to add to the `ArchRProject` object. Possible options
#' are "RTSNE", which uses `Rtsne::Rtsne()`, and "FFRTSNE", which uses `Seurat::RunTSNE()`.
#' @param name The name for the TSNE embedding to store in the given `ArchRProject` object.
#' @param perplexity An integer describing the number of nearest neighbors to compute an `Rtsne`. This argument is passed to `perplexity` in `Rtsne::Rtsne()`.
#' @param maxIterations An integer describing the maximum number of iterations when computing a TSNE. This argument is passed to `max_iter` in `Rtsne::Rtsne()`.
#' @param learningRate An integer controlling how much the weights are adjusted at each iteration. This argument is passed to `eta` in `Rtsne::Rtsne()`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing
#' depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param verbose A boolean value that indicates whether printing TSNE output.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param force A boolean value that indicates whether to overwrite the relevant data in the `ArchRProject` object if the embedding indicated by
#' `name` already exists.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... Additional parameters for computing the TSNE embedding to pass to `Rtsne::Rtsne()` (when `method = "RTSNE"`) or to `Seurat::RunTSNE()` (when method = "FFRTSNE"). 
#' @export
addTSNE <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  method = "RTSNE",
  name = "TSNE",
  perplexity = 50,
  maxIterations = 1000,
  learningRate = 200,
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  saveModel = FALSE,
  verbose = TRUE,
  seed = 1,
  force = FALSE,
  threads = max(floor(getArchRThreads() / 2), 1),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = method, name = "method", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  .validInput(input = perplexity, name = "perplexity", valid = c("integer"))
  .validInput(input = maxIterations, name = "maxIterations", valid = c("integer"))
  .validInput(input = learningRate, name = "learningRate", valid = c("integer"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))  

  #TSNE Params
  .validInput(input = perplexity, name = "perplexity", valid = c("integer", "null"))

  if(name %in% names(ArchRProj@embeddings)){
    if(!force){
      stop("Embedding Already Exists! Either set force = TRUE or use a different name!")
    }
  }

  #############################################################################################
  # Default Parameters for Input Embeddings!
  #############################################################################################

  #Merge Parameters
  embeddingParams <- list(...)
  embeddingParams$X <- getReducedDims(
      ArchRProj = ArchRProj, 
      reducedDims = reducedDims, 
      dimsToUse = dimsToUse, 
      corCutOff = corCutOff, 
      scaleDims = scaleDims
  )

  if(tolower(method)!="rtsne"){
    message("Methods other than Rtsne not currently supported, defaulting to Rtsne!")
    method <- "rtsne" 
  }

  if(tolower(method)=="rtsne"){

    .requirePackage("Rtsne", source = "cran")
    
    embeddingParams$X <- getReducedDims(
        ArchRProj = ArchRProj, 
        reducedDims = reducedDims, 
        dimsToUse = dimsToUse, 
        corCutOff = corCutOff, 
        scaleDims = scaleDims
    )

    embeddingParams$pca <- FALSE
    embeddingParams$verbose <- verbose
    embeddingParams$num_threads <- threads
    embeddingParams$max_iter <- maxIterations
    embeddingParams$eta <- learningRate

    Rtsne_tsne <- do.call(Rtsne::Rtsne, embeddingParams)
    dfEmbedding <- data.frame(Rtsne_tsne$Y)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#TSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)

  }else if(tolower(method)=="fftrtsne" | tolower(method)=="fit-tsne"){

    .requirePackage("Seurat", source = "cran")
    
    embeddingParams$X <- getReducedDims(
        ArchRProj = ArchRProj, 
        reducedDims = reducedDims, 
        dimsToUse = dimsToUse, 
        corCutOff = corCutOff, 
        scaleDims = scaleDims
    )
    
    embeddingParams$assay <- NULL
    embeddingParams$dim.embed <- 2
    embeddingParams$tsne.method <- "FIt-SNE"
    fftrtsne_tsne <- do.call(Seurat::RunTSNE, embeddingParams)

    dfEmbedding <- data.frame(fftrtsne_tsne@cell.embeddings)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#TSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)

  }else{

    stop("TSNE Method Not Currently Supported!")
  
  }

  #############################################################################################
  # Add Embedding to Project
  #############################################################################################

  embeddingParams$X <- NULL
  ArchRProj@embeddings[[name]] <- SimpleList(df = dfEmbedding, params = embeddingParams)
  
  return(ArchRProj)

}



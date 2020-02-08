##########################################################################################
# Embedding Methods
##########################################################################################

#' Add a UMAP embedding of a reduced dimensions object to an ArchRProject JJJ
#' 
#' This function will compute a UMAP embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param name The name for the UMAP embedding to be stored as in the given `ArchRProject` object.
#' @param nNeighbors An integer describing the number of nearest neighbors to compute a `umap`. This argument is passed to `n_neighbors` in `uwot::umap`.
#' @param minDist A numeric describing how tightly the `umap` is allowed to pack points together. This argument is passed to `min_dist` in `uwot::umap`. For more info on this see https://jlmelville.github.io/uwot/abparams.html.
#' @param metric A numeric describing how distance is computed in the `reducedDims` to compute a `umap`. This argument is passed to `metric` in `uwot::umap`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean describing whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution of strong biases 
#' (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. 
#' If `NULL` this will scale the dimensions depending on if this were set true when the `reducedDims` were created by the dimensionality reduction method.
#' This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param saveModel A boolean value indicating whether to save the UMAP model for downstream usage such as projection of data into the UMAP embedding.
#' @param verbose A boolean value that indicates whether printing UMAP output.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param force A boolean value that indicates whether to overwrite the relevant data in the `ArchRProject` object if the embedding named by `name` already exists.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... Additional Params to add to `uwot::umap`
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
  saveModel = FALSE,
  verbose = TRUE,
  seed = 1,
  force = FALSE,
  threads = getArchRThreads(),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = saveModel, name = "saveModel", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  #Umap Params
  .validInput(input = nNeighbors, name = "nNeighbors", valid = c("integer", "null"))
  .validInput(input = minDist, name = "minDist", valid = c("numeric", "null"))
  .validInput(input = metric, name = "metric", valid = c("character", "null"))

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

  if(saveModel){
    message("Saving UMAP model not currently supported (will be shortly), running without model!")
    #embeddingParams$ret_nn <- TRUE
    #embeddingParams$ret_model <- TRUE
    embeddingParams$ret_nn <- FALSE
    embeddingParams$ret_model <- FALSE    
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

  #############################################################################################
  # Add Embedding to Project
  #############################################################################################

  if(saveModel){
    dfEmbedding <- data.frame(uwot_umap[[1]])
  }else{
    dfEmbedding <- data.frame(uwot_umap)    
  }   
  colnames(dfEmbedding) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
  rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  embeddingParams$X <- NULL
  ArchRProj@embeddings[[name]] <- SimpleList(df = dfEmbedding, params = embeddingParams)
  
  return(ArchRProj)

}

#' Add a TSNE embedding of a reduced dimensions object to an ArchRProject JJJ
#' 
#' This function will compute a TSNE embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param method The method for computing TSNE embedding to add to `ArchRProject` object. Possible options are "RTSNE", and "FFRTSNE".
#' @param name The name for the TSNE embedding to be stored as in the given `ArchRProject` object.
#' @param perplexity An integer describing the number of nearest neighbors to compute a `Rtsne`. This argument is passed to `perplexity` in `Rtsne::Rtsne`.
#' @param maxIterations An integer describing the maximum number of iterations when computing a TSNE. This argument is passed to `max_iter` in `Rtsne::Rtsne`.
#' @param learningRate An integer controlling how much the weights are adjusted at each iteration. This argument is passed to `eta` in `Rtsne::Rtsne`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean describing whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution of strong biases 
#' (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. 
#' If `NULL` this will scale the dimensions depending on if this were set true when the `reducedDims` were created by the dimensionality reduction method.
#' This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param verbose A boolean value that indicates whether printing TSNE output.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param force A boolean value that indicates whether to overwrite the relevant data in the `ArchRProject` object if the embedding named by `name` already exists.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... Additional Params to add to for metohd = RTSNE (`Rtsne::Rtsne`) or method = FFRTSNE (`Seurat::RunTSNE`) for computing the TSNE embedding. 
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
  threads = getArchRThreads(),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
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
  .validInput(input = method, name = "method", valid = c("character"))

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

    .requirePackage("Rtsne")
    
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
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(method)=="fftrtsne" | tolower(method)=="fit-tsne"){

    .requirePackage("Seurat")
    
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
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

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



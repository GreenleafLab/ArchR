##########################################################################################
# Embedding Methods
##########################################################################################

#' Add an embedding of a reduced dimensions object to an ArchRProject
#' 
#' This function will compute an embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims QQQ The name of the `reducedDims` object to use from the designated `ArchRProject`. Possible options include "IterativeLSI", QQQ.
#' @param embedding QQQ The name of the embedding to add to the `ArchRProject` object. Possible options include "UMAP", "TUMAP", "RTSNE", and "FFRTSNE".
#' @param embeddingOut QQQ
#' @param dimsToUse QQQ A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param corCutOff QQQ A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is QQQ greater than the corCutOff, it will be excluded from analysis.
#' @param saveModel QQQ A boolean value indicating whether QQQ.
#' @param seed QQQ A number to be used as the seed for random number generation required in cluster determination. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param force A boolean value that indicates whether or not to overwrite the relevant data in the `ArchRProject` object if the given `embedding` already exists.
#' @param threads The number of threads to use for embedding generation computations.
#' @param embeddingParams A list of extra parameters to pass to the designated `embedding` function.
#' @param ... additional args
#' @export
addEmbedding <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  embedding = "UMAP",
  embeddingOut = NULL,
  dimsToUse = NULL,
  corCutOff = 0.75,
  saveModel = TRUE,
  seed = 1,
  force = FALSE,
  threads = floor(detectCores()/2),
  embeddingParams = list(),
  ...
  ){

  if(is.null(embeddingOut)){
    embeddingOut <- embedding
  }

  if(embeddingOut %in% names(ArchRProj@embeddings)){
    if(!force){
      stop("Embedding Already Exists! Either set force = TRUE or use a different name!")
    }
  }

  #############################################################################################
  # Default Parameters for Input Embeddings!
  #############################################################################################
  if(tolower(embedding)=="umap"){
    defaultEmbeddingParams <- list(
      n_neighbors = 40,
      min_dist = 0.4,
      metric = "euclidean",
      n_threads = threads, 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="tumap"){
    defaultEmbeddingParams <- list(
      n_neighbors = 40,
      min_dist = 0.4,
      metric = "euclidean",
      n_threads = threads, 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="rtsne"){
    defaultEmbeddingParams <- list(
      perplexity = 50,
      pca = FALSE,
      num_threads = threads, 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="fit-tsne" | toupper(embedding)=="fftrtsne"){
    defaultEmbeddingParams <- list(
      perplexity = 50,
      pca = FALSE,
      num_threads = threads, 
      verbose = TRUE
    )
  }else{
    defaultEmbeddingParams <- list()
  }

  #Merge Parameters
  embeddingParams <- .mergeParams(embeddingParams, defaultEmbeddingParams)

  #############################################################################################
  # Run Embedding
  #############################################################################################
  #Seed
  set.seed(seed)

  if(tolower(embedding)=="umap"){

    .requirePackage("uwot")
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    if(saveModel){
      embeddingParams$ret_nn <- TRUE
      embeddingParams$ret_model <- TRUE
    }else{
      embeddingParams$ret_nn <- FALSE
      embeddingParams$ret_model <- FALSE      
    }
    uwot_umap <- do.call(uwot::umap, embeddingParams)
    if(saveModel){
      dfEmbedding <- data.frame(uwot_umap[[1]])
    }else{
      dfEmbedding <- data.frame(uwot_umap)    
    }   
    colnames(dfEmbedding) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="tumap"){
    
    .requirePackage("uwot")
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    if(saveModel){
      embeddingParams$ret_nn <- TRUE
      embeddingParams$ret_model <- TRUE
    }else{
      embeddingParams$ret_nn <- FALSE
      embeddingParams$ret_model <- FALSE      
    }
    uwot_umap <- do.call(uwot::umap, embeddingParams)
    if(saveModel){
      dfEmbedding <- data.frame(uwot_umap[[1]])
    }else{
      dfEmbedding <- data.frame(uwot_umap)    
    }   
    colnames(dfEmbedding) <- paste0(reducedDims,"#TUMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="rtsne"){

    .requirePackage("Rtsne")
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    embeddingParams$pca <- FALSE
    Rtsne_tsne <- do.call(Rtsne::Rtsne, embeddingParams)
    dfEmbedding <- data.frame(Rtsne_tsne$Y)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#RTSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="fftrtsne" | tolower(embedding)=="fit-tsne"){

    #Seurat::RunTSNE(dims = , tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    embeddingParams$pca <- FALSE
    fftrtsne_tsne <- do.call(.fftRtsne, embeddingParams)
    dfEmbedding <- data.frame(fftrtsne_tsne)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#FITTSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else{

    stop("Embedding Method Not Currently Supported!")
  
  }

  #############################################################################################
  # Add Embedding to Project
  #############################################################################################
  embeddingParams$X <- NULL
  ArchRProj@embeddings[[embeddingOut]] <- SimpleList(df = dfEmbedding, params = embeddingParams)
  return(ArchRProj)
  
}


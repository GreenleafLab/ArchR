#' Compute Embedding from Reduced Dimensions in ArchR Project
#' 
#' This function will plot an embedding that was created from
#' computeEmbedding
#'
#' @param ArchRProj ArchRProject
#' @param reducedDims reduced dimensions to use
#' @param embedding embedding type (umap, tumap, rtsne, fftrtsne)
#' @param colorBy colorBy cellColData or Arrays in Arrows (ie GeneScoreMatrix)
#' @param name name of column in cellColData or Feature in Array in Arrows
#' @param log2Norm log2 Normalize features if they are continuous
#' @param pal custom palette to use for plotting
#' @param size size of points in plot
#' @param rastr rastr points in plot
#' @param quantCut quantile cut of continuous features
#' @param quantHex quantile evaluation for each hex in geom_hex
#' @param discreteSet discrete palette for visualizing embedding
#' @param continuousSet continuous palette for visualizing embedding
#' @param randomize randomize points prior to plotting
#' @param keepAxis keep x and y axis for plot
#' @param baseSize base size for text in plot
#' @param plotContinuous how to plot continuous features (points and hex)
#' @param plotParams additional params to pass to ggPoint/ggHex
#' @param plotWidth plot width used for creating a consistent plot independent of legend size
#' @param plotHeight plot height used for creating a consistent plot independent of legend size
#' @param ... additional args
#' @export
addEmbedding <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  embedding = "UMAP",
  embeddingOut = NULL,
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
      num_threads = threads, 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="fit-tsne" | toupper(embedding)=="fftrtsne"){
    defaultEmbeddingParams <- list(
      perplexity = 50,
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
    embeddingParams$X <- ArchRProj@reducedDims[[reducedDims]][[1]]
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
    embeddingParams$X <- ArchRProj@reducedDims[[reducedDims]][[1]]
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
    embeddingParams$X <- ArchRProj@reducedDims[[reducedDims]][[1]]
    embeddingParams$pca <- FALSE
    Rtsne_tsne <- do.call(Rtsne::Rtsne, embeddingParams)
    dfEmbedding <- data.frame(Rtsne_tsne$Y)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#RTSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="fftrtsne" | tolower(embedding)=="fit-tsne"){

    #Seurat::RunTSNE(dims = , tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
    embeddingParams$X <- ArchRProj@reducedDims[[reducedDims]][[1]]
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











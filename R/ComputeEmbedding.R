#' Add embedding of a reduced dimensions object in an ArchRProject
#' 
#' This function will compute an embedding and add to an ArchRProject.
#'
#' @param ArchRProj An ArchRProject object.
#' @param reducedDims QQQ The name of the reducedDims object to use. Possible options include "IterativeLSI", QQQ.
#' @param embedding QQQ The name of the embedding to create. Possible options include "UMAP", "TUMAP", "RTSNE", and "FFRTSNE".
#' @param colorBy QQQ colorBy cellColData or Arrays in Arrows (ie GeneScoreMatrix)
#' @param name QQQ name of column in cellColData or Feature in Array in Arrows
#' @param log2Norm A boolean value that indicates whether log2 Normalization should be performed on the features if they are continuous.
#' @param pal The name or numeric index of a custom palette from ArchR_palettes to use for plotting the individual points of the embedding visualization.
#' @param size The numeric size of points to plot.
#' @param rastr A boolean valut that indicates that the plot should be rasterized. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param quantCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution. This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the upper threshold and y is the lower threshold. For example, quantileCut = c(0.975,0.025) will take the top and bottom 2.5% of values and set them to the value of the 97.5th and 2.5th percentile values respectively.
#' @param quantHex QQQ quantile evaluation for each hex in geom_hex
#' @param discreteSet QQQ The name or numeric index of a discrete palette from ArchR_palettesdiscrete to use for plotting QQQ.
#' @param continuousSet QQQcontinuous palette for visualizing embedding
#' @param randomize A boolean value that determines whether to randomly order the plotting of points to avoid (for ex.) all points from a single cluster being plotted as the top-most layer of the plot.
#' @param keepAxis QQQ keep x and y axis for plot
#' @param baseSize QQQ The numeric font size to be used in the plot. This applies to all plot labels.
#' @param plotContinuous QQQ how to plot continuous features (points and hex)
#' @param plotParams QQQ additional params to pass to ggPoint/ggHex
#' @param plotWidth QQQ plot width used for creating a consistent plot independent of legend size
#' @param plotHeight QQQ plot height used for creating a consistent plot independent of legend size
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











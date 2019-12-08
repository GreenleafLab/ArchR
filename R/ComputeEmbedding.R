#' Compute the embedding from a reduced dimensions object in an ArchRProject
#' 
#' This function will compute an embedding for a given ArchRProject.
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
ComputeEmbedding <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  embedding = "UMAP",
  embeddingOut = NULL,
  saveModel = TRUE,
  seed = 1,
  force = FALSE,
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
      n_threads = floor(detectCores()/2), 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="tumap"){
    defaultEmbeddingParams <- list(
      n_neighbors = 40,
      min_dist = 0.4,
      metric = "euclidean",
      n_threads = floor(detectCores()/2), 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="rtsne"){
    defaultEmbeddingParams <- list(
      perplexity = 50,
      num_threads = floor(detectCores()/2), 
      verbose = TRUE
    )
  }else if(tolower(embedding)=="fit-tsne" | toupper(embedding)=="fftrtsne"){
    defaultEmbeddingParams <- list(
      perplexity = 50,
      num_threads = floor(detectCores()/2), 
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

.fftRtsne <- function(
  X, 
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  max_iter = 1000,
  fft_not_bh = TRUE,
  ann_not_vptree = TRUE,
  stop_early_exag_iter = 250,
  exaggeration_factor = 12.0,
  no_momentum_during_exag = FALSE,
  start_late_exag_iter = -1.0,
  late_exag_coeff = 1.0,
  mom_switch_iter = 250,
  momentum = 0.5,
  final_momentum = 0.8,
  learning_rate = 200,
  n_trees = 50, 
  search_k = -1, 
  rand_seed = -1,
  nterms = 3, 
  intervals_per_integer = 1,
  min_num_intervals = 50, 
  K = -1,
  sigma = -30,
  initialization = NULL,
  data_path = NULL,
  result_path = NULL,
  load_affinities = NULL,
  fast_tsne_path = NULL,
  nthreads = 0,
  perplexity_list = NULL, 
  get_costs = FALSE,
  df = 1.0
  ){
  
  tstart <- Sys.time()
  .messageDiffTime("Running FIt-SNE version 1.1.0 from https://github.com/KlugerLab/FIt-SNE/", tstart, addHeader = TRUE)
  #version_number <- '1.1.0'

  .messageDiffTime("Checking Input", tstart)
  if (is.null(fast_tsne_path)) {
    if (.Platform$OS.type == "unix") {
      fast_tsne_path <- file.path(FAST_TSNE_SCRIPT_DIR, "bin", "fast_tsne")
    } else {
      fast_tsne_path <- file.path(FAST_TSNE_SCRIPT_DIR, "bin", "FItSNE.exe")
    }
  }

  if (is.null(data_path)) {
    data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
  }
  if (is.null(result_path)) {
    result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
  }
  if (is.null(fast_tsne_path)) {
    fast_tsne_path <- system2('which', 'fast_tsne', stdout = TRUE)
  }
  fast_tsne_path <- normalizePath(fast_tsne_path)
  if (!file_test('-x', fast_tsne_path)) {
    stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
  }

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0) ) {
    stop("Incorrect theta.")
  }
  if (nrow(X) - 1 < 3 * perplexity){
    stop("Perplexity is too large.")
  }
  if (!is.matrix(X)) {
    stop("Input X is not a matrix")
  }
  if (!(max_iter > 0)) {
    stop("Incorrect number of iterations.")
  }
  if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter < 0) {
    stop("stop_early_exag_iter should be a positive integer")
  }
  if (!is.numeric(exaggeration_factor)) {
    stop("exaggeration_factor should be numeric")
  }
  if (!is.numeric(df)) {
    stop("df should be numeric")
  }
  if (!is.wholenumber(dims) || dims <= 0) {
    stop("Incorrect dimensionality.")
  }

  if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else { 
      search_k <- n_trees * K
    }
  }

  if (fft_not_bh) {
    nbody_algo <- 2
  } else {
    nbody_algo <- 1
  }

  if (is.null(load_affinities)) {
    load_affinities <- 0
  } else {
    if (load_affinities == 'load') {
      load_affinities <- 1
    } else if (load_affinities == 'save') {
      load_affinities <- 2
    } else {
      load_affinities <- 0
    }
  }

  if (ann_not_vptree) {
    knn_algo <- 1
  } else {
    knn_algo <- 2
  }
  tX <- as.numeric(t(X))

  .messageDiffTime("Writing Data for FIt-SNE", tstart)
  f <- file(data_path, "wb")
  n <- nrow(X)
  D <- ncol(X)
  writeBin(as.integer(n), f, size = 4)
  writeBin(as.integer(D), f, size = 4)
  writeBin(as.numeric(theta), f, size = 8) #theta
  writeBin(as.numeric(perplexity), f, size = 8)

  if (perplexity == 0) {
    writeBin(as.integer(length(perplexity_list)), f, size = 4)
    writeBin(perplexity_list, f) 
  }

  writeBin(as.integer(dims), f, size = 4)
  writeBin(as.integer(max_iter), f, size = 4)
  writeBin(as.integer(stop_early_exag_iter), f, size = 4)
  writeBin(as.integer(mom_switch_iter), f, size = 4)
  writeBin(as.numeric(momentum), f, size = 8)
  writeBin(as.numeric(final_momentum), f, size = 8)
  writeBin(as.numeric(learning_rate), f, size = 8)
  writeBin(as.integer(K), f, size = 4) #K
  writeBin(as.numeric(sigma), f, size = 8) #sigma
  writeBin(as.integer(nbody_algo), f, size = 4)  #not barnes hut
  writeBin(as.integer(knn_algo), f, size = 4) 
  writeBin(as.numeric(exaggeration_factor), f, size = 8) #compexag
  writeBin(as.integer(no_momentum_during_exag), f, size = 4) 
  writeBin(as.integer(n_trees), f, size = 4) 
  writeBin(as.integer(search_k), f, size = 4) 
  writeBin(as.integer(start_late_exag_iter), f, size = 4) 
  writeBin(as.numeric(late_exag_coeff), f, size = 8) 

  writeBin(as.integer(nterms), f, size = 4) 
  writeBin(as.numeric(intervals_per_integer), f, size = 8) 
  writeBin(as.integer(min_num_intervals), f, size = 4) 
  writeBin(tX, f) 
  writeBin(as.integer(rand_seed), f, size = 4) 
  writeBin(as.numeric(df), f, size = 8)
  writeBin(as.integer(load_affinities), f, size = 4) 
  if (!is.null(initialization)) { writeBin( c(t(initialization)), f) }    
  close(f) 

  .messageDiffTime("Executing system FIt-SNE", tstart)
  flag <- system2(command = fast_tsne_path, 
                  args = c(version_number, data_path, result_path, nthreads))
  if (flag != 0) {
    stop('tsne call failed')
  }
  f <- file(result_path, "rb")
  n <- readBin(f, integer(), n = 1, size = 4)
  d <- readBin(f, integer(), n = 1, size = 4)
  Y <- readBin(f, numeric(), n = n * d)
  Y <- t(matrix(Y, nrow = d))
  if (get_costs) {
    readBin(f, integer(), n = 1, size = 4)
    costs <- readBin(f, numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(f)
  file.remove(data_path)
  file.remove(result_path)
  
  .messageDiffTime("Successfully finished FIt-SNE", tstart)

  return(Yout)

}


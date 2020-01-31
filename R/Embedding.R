##########################################################################################
# Embedding Methods
##########################################################################################

#' Add an embedding of a reduced dimensions object to an ArchRProject
#' 
#' This function will compute an embedding and add it to an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param embedding The type of embedding to add to the `ArchRProject` object. Possible options include "UMAP", "TUMAP", "RTSNE", and "FFRTSNE".
#' @param embeddingOut The name for the embedding to be stored in the given `ArchRProject` object.
#' @param embeddingParams A list of extra parameters to pass to the designated `embedding` function.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean describing whether to rescale the total variance for each principal component. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. If `NULL` this will scale the dimensions depending on if this were set true when the `reducedDims` were created by `addIterativeLSI`.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param saveModel A boolean value indicating whether to save the UMAP model for downstream usage such as projection of data into the UMAP embedding. Only relevant if `embedding` is set to "UMAP".
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param force A boolean value that indicates whether to overwrite the relevant data in the `ArchRProject` object if the embedding named by `embeddingOut` already exists.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... additional args
#' @export
addEmbedding <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI", 
  embedding = "UMAP",
  embeddingOut = NULL,
  embeddingParams = list(),
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  saveModel = FALSE,
  seed = 1,
  force = FALSE,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = embeddingOut, name = "embeddingOut", valid = c("character", "null"))
  .validInput(input = embeddingParams, name = "embeddingParams", valid = c("list"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = saveModel, name = "saveModel", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

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
      metric = "cosine",
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
      perplexity = 50
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
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, 
      corCutOff = corCutOff, scaleDims = scaleDims)

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
    uwot_umap <- do.call(uwot::umap, embeddingParams)
    if(saveModel){
      dfEmbedding <- data.frame(uwot_umap[[1]])
    }else{
      dfEmbedding <- data.frame(uwot_umap)    
    }   
    colnames(dfEmbedding) <- paste0(reducedDims,"#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="rtsne"){

    .requirePackage("Rtsne")
    embeddingParams$X <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, 
      corCutOff = corCutOff, scaleDims = scaleDims)

    embeddingParams$pca <- FALSE
    Rtsne_tsne <- do.call(Rtsne::Rtsne, embeddingParams)
    dfEmbedding <- data.frame(Rtsne_tsne$Y)   
    colnames(dfEmbedding) <- paste0(reducedDims,"#RTSNE_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(ArchRProj@reducedDims[[reducedDims]][[1]])

  }else if(tolower(embedding)=="fftrtsne" | tolower(embedding)=="fit-tsne"){

    .requirePackage("Seurat")
    embeddingParams$object <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, 
      corCutOff = corCutOff, scaleDims = scaleDims)
    
    embeddingParams$assay <- NULL
    embeddingParams$dim.embed <- 2
    embeddingParams$tsne.method <- "FIt-SNE"
    fftrtsne_tsne <- do.call(RunTSNE, embeddingParams)

    dfEmbedding <- data.frame(fftrtsne_tsne@cell.embeddings)   
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


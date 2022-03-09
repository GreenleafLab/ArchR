####################################################################
# Import Multi-Modal Data
####################################################################

#' Import Feature Matrix from 10x Feature HDF5 file.
#' 
#' This function will import the feature matrix from a 10x feature hdf5 file.
#'
#' @param input A character of paths to 10x feature hdf5 file(s). These will traditionally have a suffix similar to "filtered_feature_bc_matrix.h5".
#' @param names A character of sample names associated with each input file.
#' @param featureType The name of the feature to extract from the 10x feature file. 
#' See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices for more information.
#' @export
import10xFeatureMatrix <- function(
  input = NULL, 
  names = NULL, 
  featureType = "Gene Expression"
  ){

  .validInput(input = input, name = "input", valid = c("character"))
  .validInput(input = names, name = "names", valid = c("character"))
  .validInput(input = featureType, name = "featureType", valid = c("character"))
  
  if (!all(file.exists(input))) {
    stop("Not all input file paths exist!")
  }
  featureMats <- lapply(seq_along(input), function(y) {
    message("Importing Feature Matrix ", y, " of ", length(input))
    .importFM(featureMatrix = input[y], featureType = featureType, 
              name = names[y])
  })

  #if more than one filtered feature barcode matrix is supplied, then merge the RSE objects
    if (length(featureMats) > 1) {
    rse_final <- featureMats[[1]]
    for (i in 2:length(featureMats)) {
      print(paste0("featureMats[[",i,"]]"))
      if (!identical(rownames(rse_final), rownames(featureMats[[i]]))) {
        stop("Error - rownames (genes) of individual RNA objects are not equivalent.")
      }
      if (!identical(colnames(rowData(rse_final)), colnames(rowData(featureMats[[i]])))) {
        stop("Error - rowData (gene metadata) of individual RNA objects have different columns. This is highly unusual and merging has been aborted.")
      }
      if (!identical(names(assays(rse_final)), names(assays(featureMats[[i]])))) {
        stop("Error - available assays of individual RNA objects are not equivalent. Each object is expected to only have one assay named 'counts'.")
      }
      
      #check each column in rowData to check for mismatches that should be thrown as warnings
      #occasionally, it seems like 10x is annotating different ensembl IDs to the same gene which seems like a bad way to go
      #this is a bit heavy-handed but it seems like the safest thing to do is report any mismatch rather than merge blindly
      mismatchWarning <- TRUE
      for (x in 1:ncol(rowData(rse_final))) {
        if (!identical(rowData(rse_final)[,x], rowData(featureMats[[i]])[,x])) {
          if(mismatchWarning) {
            message(sprintf("Warning! Some values within column \"%s\" the rowData of your objects do not precisely match!", colnames(rowData(rse_final))[x]))
            message("This is often caused by slight variations in Ensembl IDs used by cellranger. ArchR will ignore these mismatches and allow merging to proceed but you should check to make sure that these are ok for your data.\n")
            mismatchWarning <- FALSE
          }
          
          mismatch <- which(rowData(rse_final)[,x] != rowData(featureMats[[i]])[,x])
          for (y in 1:length(mismatch)) {
            message(sprintf("Mismatch in column \"%s\" row %s for %s: %s does not exactly match %s!", colnames(rowData(rse_final))[x], mismatch[y], names[i], rowData(rse_final)[mismatch[y],x], rowData(featureMats[[i]])[mismatch[y],x]))
          }
        }
      }
      
      rse_final <- SummarizedExperiment::cbind(rse_final, featureMats[[i]])
    }
    return(rse_final)
  }
  else {
    return(featureMats)
  }
}

.importFM <- function(featureMatrix = NULL, featureType = NULL, name = NULL){

  o <- h5closeAll()
  barcodes <- h5read(featureMatrix, "/matrix/barcodes")
  data <- h5read(featureMatrix, "/matrix/data")
  indices <- h5read(featureMatrix, "/matrix/indices")
  indptr <- h5read(featureMatrix, "/matrix/indptr")
  shape <- h5read(featureMatrix, "/matrix/shape")

  spMat <- sparseMatrix(
    i = indices, 
    p = indptr, 
    x = data, 
    dims = shape,
    index1 = FALSE
  )

  colnames(spMat) <- paste0(name, "#", barcodes)

  features <- h5read(featureMatrix, "/matrix/features")
  features <- lapply(seq_along(features), function(x){
    if(length(features[[x]]) == nrow(spMat)){
      if(object.size(features[[x]]) > object.size(Rle(features[[x]]))){
        df <- DataFrame(x = Rle(features[[x]]))
      }else{
        df <- DataFrame(x = features[[x]])
      }
      colnames(df) <- names(features)[x]
      df
    }else{
      NULL
    }
  })
  features <- Reduce("cbind",features[!unlist(lapply(features,is.null))])

  se <- SummarizedExperiment(assays = SimpleList(counts = spMat), rowData = features)

  rownames(se) <- features$name

  if("feature_type" %in% colnames(rowData(se))){
    if(!is.null(featureType)){
      idx <- BiocGenerics::which(rowData(se)$feature_type %bcin% featureType)
      if(length(idx) == 0){
        stop("Error featureType not within provided features!")
      }
      se <- se[idx]
    }
  }

  if("interval" %in% colnames(rowData(se))){
    idxNA <- which(rowData(se)$interval=="NA")
    if(length(idxNA) > 0){
      se <- se[-idxNA, ]
    }
    rr <- GRanges(paste0(rowData(se)$interval))
    mcols(rr) <- rowData(se)
    se <- SummarizedExperiment(assays = SimpleList(counts = assay(se)), rowRanges = rr)
  }

  idxDup <- which(rownames(se) %in% rownames(se[duplicated(rownames(se))]))
  names(idxDup) <- rownames(se)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(se[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    se <- se[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  gc()

  se

}

####################################################################
# Combined Modalities
####################################################################

#' Combine two or more modalities dimensionality reductions.
#' 
#' This function will combine two or more modalities dimensionality reductions into a single reduction.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name for the combinedDims to be stored as.
#' @param reducedDims The name of the `reducedDims` objects (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param dimWeights A vector of weights to be used to weight each dimensionality reduction when combining.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` objects to use.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @export
addCombinedDims <- function(
  ArchRProj = NULL, 
  name = "CombinedDims",
  reducedDims = NULL, 
  dimWeights = NULL,
  dimsToUse = NULL, 
  scaleDims = NULL,
  corCutOff = 0.75
  ){

  if(is.null(dimWeights)){
    dimWeights <- rep(1, length(reducedDims))
  }

  combinedDims <- lapply(seq_along(reducedDims), function(x){
    rD <- getReducedDims(
      ArchRProj = ArchRProj,
      reducedDims = reducedDims[x],
      dimsToUse = dimsToUse,
      scaleDims = scaleDims,
      corCutOff = corCutOff
    )
    cV <- apply(as.matrix(rD), 2, var) 
    normV <- 1 / sqrt(sum(cV))
    rD * normV * dimWeights[x]
  }) %>% Reduce("cbind", .)

  ArchRProj@reducedDims[[name]] <- SimpleList(
    matRD = combinedDims, 
    scaleDims = NA, 
    corToDepth = list(scaled = rep(0, ncol(combinedDims)), none = rep(0, ncol(combinedDims)))
  )

  ArchRProj

}













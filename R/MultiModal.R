####################################################################
# Import Multi-Modal Data
####################################################################

#' Import Feature Matrix from 10x Feature HDF5 file.
#' 
#' This function will import the feature matrix from a 10x feature hdf5 file.
#'
#' @param input A character of paths to 10x feature hdf5 file(s). These will traditionally have a suffix similar to "filtered_feature_bc_matrix.h5".
#' @param names A character of sample names associated with each input file.
#' @param strictMatch Only relevant when multiple input files are used. A boolean that indictes whether rows (genes) that do not match perfectly in the matrices
#' should be removed (`strictMatch = TRUE`) or coerced (`strictMatch = FALSE`). Cell Ranger seems to occassionally use different ensembl ids for the same gene across
#' different samples. If you are comfortable tolerating such mismatches, you can coerce all matrices to fit together, in which case the gene metadata present in
#' the first listed sample will be applied to all matrices for that particular gene entry. Sometimes, the h5 files from Cell Ranger have mismatches in the actual gene names.
#' See the `force` paramter for handling mismatched gene names.
#' @param verbose Only relevant when multiple input files are used. A boolean that indicates whether messaging about mismatches should be verbose (`TRUE`) or minimal (`FALSE`)
#' @param featureType The name of the feature to extract from the 10x feature file. 
#' See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices for more information.
#' @param features A genomic ranges object containing a "name" column to help fill missing 10x intervals for RSE.
#' For example, in hg38 features provided could be using Bioconductors `genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)`.
#' @param genesToExclude A character vector of gene names to be excluded from the returned `SummarizedExperiment` object. This can be useful when dealing with very large
#' datasets because R has an upper limit for how large a given object can be. If you have >1.5 million cells, you cannot create a matrix in R with all ~30,000 genes so you
#' can use this parameter to remove genes upfront.
#' @export
import10xFeatureMatrix <- function(
  input = NULL, 
  names = NULL,
  strictMatch = TRUE,
  verbose = TRUE,
  featureType = "Gene Expression",
  features = NULL,
  genesToExclude = NULL,
  force = FALSE
){
  
  .validInput(input = input, name = "input", valid = c("character"))
  .validInput(input = names, name = "names", valid = c("character"))
  .validInput(input = strictMatch, name = "strictMatch", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = featureType, name = "featureType", valid = c("character"))
  .validInput(input = features, name = "features", valid = c("GRanges", "NULL"))
  .validInput(input = genesToExclude, name = "genesToExclude", valid = c("character","null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  
  if (!all(file.exists(input))) {
    stop("Not all input file paths exist!")
  }
  
  message("Importing Feature Matrix ", 1, " of ", length(input))
  rse_final <- .import10xToSE(
    h5 = input[1], 
    type10x = featureType, 
    name = names[1],
    ranges = features
  )
  #convert the gene names to exclude to gene indices to exclude. These must be removed from each imported SE upfront
  #currently, this assumes that all other samples will follow the same indexing of the first sample which may be a risky
  #assumption but I dont want to match each sample individually based on gene names as we already know that some samples
  #use different gene names
  if(!is.null(genesToExclude)) {
    genesToExclude <- which(rownames(rse_final) %in% genesToExclude)
  }
  
  if(length(genesToExclude) > 0) {
    rse_final <- rse_final[-genesToExclude,]
  }
  
  rowsToRemove <- c() #row indices that need to be removed from rse_final because they do not match across samples being merged
    
  #if more than one filtered feature barcode matrix is supplied, then merge the RSE objects
  if (length(input) > 1) {
    message("Merging individual RNA objects...")
    
    #merge all others into 1st -- the first sample is the template sample
    
    #for each additional feature matrix (starting with the second), look for mismatches with rse_final and merge accordingly
    for (i in seq(2, length(input))){
      
      message(sprintf("\nMerging %s", names[i]))
      message("Importing Feature Matrix ", i, " of ", length(input))
      
      #Read RSE
      rse_i <- .import10xToSE(
        h5 = input[i], 
        type10x = featureType, 
        name = names[i],
        ranges = features
      )
      
      if(length(genesToExclude) > 0) {
        rse_i <- rse_i[-genesToExclude,]
      }
      
      mismatchWarning <- TRUE #a boolean to prevent output of the warning message many times and only output it once
      
      if (!identical(rownames(rse_final), rownames(rse_i))) {
        if(force) {
          #remove rows from rse_i that dont exist in rse_final
          iNotFinal <- which(rownames(rse_i) %ni% rownames(rse_final))
          geneMismatch <- rownames(rse_i)[iNotFinal]
          if(length(iNotFinal) > 0) {
            rse_i <- rse_i[-iNotFinal,]
          }
          
          #earmark mismatched rows for downstream removal that dont exist in rse_i but do exist in rse_final
          finalNotI <- which(rownames(rse_final) %ni% rownames(rse_i))
          rowsToRemove <- unique(c(rowsToRemove, finalNotI))
          geneMismatch <- unique(c(geneMismatch, rownames(rse_final)[finalNotI]))
          #add dummy rows to rse_i with these gene names and rowData. this allows for proper alignment but these rows will eventually be removed
          dummySE <- rse_i[1:length(finalNotI)]
          rowData(dummySE) <- rowData(rse_final[rownames(rse_final)[finalNotI],])
          rowRanges(dummySE) <- rowRanges(rse_final[rownames(rse_final)[finalNotI],])
          rownames(dummySE) <- rownames(rse_final)[finalNotI]
          rse_i <- SummarizedExperiment::rbind(rse_i, dummySE)
          rse_i <- ArchR:::.sortRSE(rse_i)
          
          #check to ensure that the rownames now exactly match. if not something is wrong
          if(!identical(rownames(rse_final), rownames(rse_i))) {
            stop("Error - Something has gone wrong with attempting to align gene names across your individual h5 files. Please report to GitHub Issues.")
          }
          
          message(sprintf("Warning! Some gene names do not match between sample %s and the reference sample (sample #1). \"force = TRUE\" so the following genes will be removed:\n%s\n", i, paste(geneMismatch, collapse = ",")))
          
        } else {
          stop("Error - rownames (genes) of individual RNA objects are not equivalent. To bypass this, set \"force = TRUE\" and the offending genes will be removed.")
        }
      }
      if (!identical(colnames(rowData(rse_final)), colnames(rowData(rse_i)))) {
        stop("Error - rowData (gene metadata) of individual RNA objects have different columns. This is highly unusual and merging has been aborted. Please check your h5 input files..")
      }
      if (!identical(names(assays(rse_final)), names(assays(rse_i)))) {
        stop("Error - available assays of individual RNA objects are not equivalent. Each object is expected to only have one assay named 'counts'. Merging has been aborted. Please check your h5 input files.")
      }
      
      #check each column in rowData to check for mismatches that should be thrown as warnings
      #occasionally, it seems like 10x is annotating different ensembl IDs to the same gene which seems like a bad way to go.
      #this is a bit heavy-handed but it seems like the safest thing to do is report any mismatch rather than merge blindly
      
      for (x in seq_len(ncol(rowData(rse_final)))){
        if (!identical(rowData(rse_final)[,x], rowData(rse_i)[,x])) {
          if(mismatchWarning) {
            message(sprintf("Warning! Some values within column \"%s\" of the rowData (gene metadata) of your objects do not precisely match!", colnames(rowData(rse_final))[x]))
            message("This is often caused by slight variations in Ensembl IDs and gene locations used by cellranger across different samples. ArchR will ignore these mismatches and allow merging to proceed but you should check to make sure that these are ok for your data. See the strictMatch parameter for how these mismatches will be handled.\n")
            mismatchWarning <- FALSE
          }
          
          #detect all of the mismatches betwenn rse_final and the current featureMat
          mismatch <- which(rowData(rse_final)[,x] != rowData(rse_i)[,x])
          #for each detected mismatch, handle the mismatch according to the value of strictMatch
          for (y in 1:length(mismatch)) {
            if (verbose) {
              message(sprintf("Mismatch in column \"%s\" row %s for %s: %s does not exactly match %s!", colnames(rowData(rse_final))[x], mismatch[y], names[i], rowData(rse_final)[mismatch[y],x], rowData(featureMats[[i]])[mismatch[y],x]))
            }
            if (strictMatch) {
              if (verbose) {
                message("strictMatch = TRUE so the corresponding gene entry with mismatching information will be removed.")
              }
              rowsToRemove <- unique(c(rowsToRemove, mismatch[y]))
              #temporarily force the data to match so that merging can occur easily. Mismatched rows will be removed later
              rowData(rse_i)[mismatch[y],] <- rowData(rse_final)[mismatch[y],]
              rowRanges(rse_i)[mismatch[y]] <- rowRanges(rse_final)[mismatch[y]]
            } else {
              if (verbose) {
                message("strictMatch = FALSE so mismatching information will be coerced to match the first sample provided.")
              }
              rowData(rse_i)[mismatch[y],] <- rowData(rse_final)[mismatch[y],]
              rowRanges(rse_i)[mismatch[y]] <- rowRanges(rse_final)[mismatch[y]]
            }
          }
        }
      }
      
      #merge rse_i into res_final and garbage collect
      rse_final <- SummarizedExperiment::cbind(rse_final, rse_i)
      remove(rse_i)
      gc()
    }
    
  }
  
  if(length(rowsToRemove) > 0){
    rse_final <- rse_final[-rowsToRemove,]
  }

  if("rowRanges" %in% slotNames(rse_final)){
    idxNA <- which(seqnames(rse_final)=="chrNA")
    if(length(idxNA) > 0){
      namesNA <- paste0(rowRanges(rse_final)$name[idxNA], collapse=",")
      warning(
        "These features had no interval present and were given a fake GRanges location on `chrNA`:\n\n",
        namesNA,
        "\n\nEither modify the output SE or see `features` input to use a reference GRanges to add to these features!"
      )
    }
  }

  rse_final
  
}

.import10xToSE <- function(
  h5 = NULL, 
  type10x = NULL, 
  name = NULL,
  ranges = NULL
){
  
  #Shape
  shape <- h5read(h5, "/matrix/shape")
  
  #Read features10x
  features10x <- h5read(h5, "/matrix/features")
  features10x <- lapply(seq_along(features10x), function(x){
    if(length(features10x[[x]]) == shape[1]){
      if(object.size(features10x[[x]]) > object.size(Rle(features10x[[x]]))){
        df <- DataFrame(x = Rle(as.vector(features10x[[x]])))
      }else{
        df <- DataFrame(x = as.vector(features10x[[x]]))
      }
      colnames(df) <- names(as.vector(features10x))[x]
      df
    }else{
      NULL
    }
  })
  features10x <- Reduce("cbind",features10x[!unlist(lapply(features10x,is.null))])
  
  #Determine Idx
  if(!is.null(type10x)){
    idx <- which(paste0(features10x$feature_type) %in% type10x)
  }else{
    idx <- seq_len(nrow(features10x))
  }
  if(length(idx)==0){
    stop(
      paste0(
        h5,
        "\nMissing `type10x`! Feature Types in h5:\n",
        "\t", paste0(unique(features10x$feature_type),collapse="; ")
      )
    )
  }
  
  #Subset
  features10x <- features10x[idx, , drop=FALSE]
  
  #Interval
  if("interval" %in% colnames(features10x)){
    
    idxNA <- which(features10x$interval=="NA")
    
    if(length(idxNA) > 0){
      
      #if ranges = NULL, then interval is ignored(?)
      if(!is.null(ranges) & is(ranges, "GRanges")){
        #Fix ranges
        idx1 <- paste0(seqnames(ranges)) %in% c(1:22, "X", "Y", "MT")
        if(length(idx1) > 0){
          ranges2 <- GRanges(
            seqnames = ifelse(idx1, paste0("chr",seqnames(ranges)), paste0(seqnames(ranges))),
            ranges = GenomicRanges::ranges(ranges)
          )
          mcols(ranges2) <- mcols(ranges)
        }
        
        #Try To Use Ranges To Match
        features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
        namesNA <- features10xNA$name
        idxFix <- match(namesNA, mcols(ranges2)[, grep("name", colnames(mcols(ranges)), ignore.case=TRUE)])
        if(length(idxFix[!is.na(idxFix)]) > 0){
          message("Correcting missing intervals...")
          idx2 <- which(!is.na(idxFix))
          rangesFix <- ranges2[idxFix[idx2]]
          strand(rangesFix) <- "*"
          features10xNA$interval[idx2] <- paste0(rangesFix)
          features10x[which(features10x$interval=="NA"), ] <- features10xNA
        }
        
        #NA add Fake Chromosome
        features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
        if(nrow(features10xNA) > 0){
          features10xNA$interval <- paste0("chrNA:1-1")
          features10x[which(features10x$interval=="NA"), ] <- features10xNA
        }
        
        features10x$ranges <- GRanges(paste0(features10x$interval))
        features10x$interval <- NULL

      }else{

        #NA add Fake Chromosome
        features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
        if(nrow(features10xNA) > 0){
          features10xNA$interval <- paste0("chrNA:1-1")
          features10x[which(features10x$interval=="NA"), ] <- features10xNA
        }

        features10x$ranges <- GRanges(paste0(features10x$interval))
        features10x$interval <- NULL

      }
      
    }else {
      features10x$ranges <- GRanges(paste0(features10x$interval))
      features10x$interval <- NULL
    }
    
  }
  
  #Read Matrix
  mat <- sparseMatrix(
    i = h5read(h5, "/matrix/indices"), 
    p = h5read(h5, "/matrix/indptr"), 
    x = h5read(h5, "/matrix/data"), 
    dims = shape,
    index1 = FALSE
  )
  barcodes <- h5read(h5, "/matrix/barcodes")
  if(!is.null(name)){
    colnames(mat) <- paste0(name, "#", barcodes)
  }else{
    colnames(mat) <- barcodes
  }
  
  #Subset
  mat <- mat[idx, , drop = FALSE]
  gc()
  
  #Summarized Experiment
  if("ranges" %in% colnames(features10x)){
    mat <- SummarizedExperiment(
      assays = SimpleList(
        data = mat
      ),
      rowRanges = features10x$ranges
    )
    rowData(mat) <- features10x
    rowData(mat)$ranges <- NULL
  }else{
    mat <- SummarizedExperiment(
      assays = SimpleList(
        data = mat
      ),
      rowData = features10x
    )    
  }
  
  rownames(mat) <- rowData(mat)$name
  .sortRSE(mat)
  
}

.sortRSE <- function(rse){
  if(!is.null(rowRanges(rse))){
    sort.GenomicRanges(sortSeqlevels(rse), ignore.strand = TRUE)
  }else{
    rse[order(rowData(rse)$name)]
  }
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













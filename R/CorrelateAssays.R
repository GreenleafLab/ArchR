##########################################################################################
# Assay Correlation Methods
##########################################################################################

#' Correlate Matrices within an ArchRProject
#' 
#' This function will correlate 2 matrices within an ArchRProject by name matching.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix1 The name of a valid matrix in the `ArchRProject` to use as the first matrix in the correlation comparison.
#' See `getAvailableMatrices()` for valid options.
#' @param useMatrix2 The name of a valid matrix in the `ArchRProject` to use as the second matrix in the correlation comparison.
#' See `getAvailableMatrices()` for valid options.
#' @param useSeqnames1 JJJ An optional character vector of the subset of `seqnames` (i.e. chromosomes) to use from the matrix identified by `useMatrix1`.
#' @param useSeqnames2 JJJ An optional character vector of the subset of `seqnames` (i.e. chromosomes) to use from the matrix identified by `useMatrix2`.
#' @param removeFromName1 JJJ A character vector describing how column names from matrix 1 should be adjusted.  Valid options include "underscore", "dash", and "numeric".
#' For example, `removeFromName1 = c("dash","underscore")` indicates that all dashes ("-") and underscores ("_") should be removed from the column names of the matrix.
#' The remainder of each column name will be left untouched. This allows for better matching of column names across diverse input matrices.
#' @param removeFromName2 JJJ A character vector describing how column names from matrix 2 should be adjusted.  Valid options include "underscore", "dash", and "numeric".
#' For example, `removeFromName1 = c("dash","underscore")` indicates that all dashes ("-") and underscores ("_") should be removed from the column names of the matrix.
#' The remainder of each column name will be left untouched. This allows for better matching of column names across diverse input matrices.
#' @param log2Norm1 A boolean value that indicates whether or not to log2 normalize matrix 1.
#' @param log2Norm2 A boolean value that indicates whether or not to log2 normalize matrix 2.
#' @param reducedDims JJJ WHAT IS THIS USED FOR? The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param dimsToUse JJJ WHY IS AN EMBEDDING BEING COMPUTED? A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims JJJ NOT USED IN FUNCTION? A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param knnMethod JJJ The method to be used for k-nearest neighbor computations. Options are "nabor", "RANN", and "FNN" and the corresponding package is required.
#' @param threads The number of threads to be used for parallel computing.
#' @export
correlateMatrices <- function(
  ArchRProj = NULL,
  useMatrix1 = NULL,
  useMatrix2 = NULL,
  useSeqnames1 = NULL,
  useSeqnames2 = NULL,
  removeFromName1 = c("underscore", "dash"),
  removeFromName2 = c("underscore", "dash"),
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,#JJJ should this have a default value that isnt NULL?
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  seed = 1, 
  knnMethod = NULL,
  threads = getArchRThreads()
  ){

  #JJJ .validInput(input = input, name = "input", valid = c("ArchRProj", "matrix"))
  #JJJ .validInput(input = useMatrix1, name = "useMatrix1", valid = c("character"))
  #JJJ .validInput(input = useMatrix2, name = "useMatrix2", valid = c("character"))
  #JJJ .validInput(input = useSeqnames1, name = "useSeqnames1", valid = c("character", "null"))
  #JJJ .validInput(input = useSeqnames2, name = "useSeqnames2", valid = c("character","null"))
  #JJJ .validInput(input = removeFromName1, name = "removeFromName1", valid = c("character","null"))
  #JJJ .validInput(input = removeFromName2, name = "removeFromName2", valid = c("character","null"))
  #JJJ .validInput(input = log2Norm1, name = "log2Norm1", valid = c("boolean"))
  #JJJ .validInput(input = log2Norm2, name = "log2Norm2", valid = c("boolean"))
  #JJJ .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  #JJJ .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  #JJJ .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean"))
  #JJJ .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  #JJJ .validInput(input = k, name = "k", valid = c("integer"))
  #JJJ .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  #JJJ .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  #JJJ .validInput(input = seed, name = "seed", valid = c("numeric"))
  #JJJ .validInput(input = knnMethod, name = "knnMethod", valid = c("character"))
  #JJJ .validInput(input = threads, name = "threads", valid = c("integer"))

  tstart <- Sys.time()

  set.seed(seed)

  #Get Available Matrices
  matrixNames <- getAvailableMatrices(proj)

  if(useMatrix1 %ni% matrixNames){
    stop(paste0("useMatrix1 (",useMatrix1,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")))
  }

  if(useMatrix2 %ni% matrixNames){
    stop(paste0("useMatrix2 (",useMatrix2,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")))
  }

  #Get Matrix Classes
  matrixClass1 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix1, "/Info/Class")))
  matrixClass2 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix2, "/Info/Class")))

  #Get Feature DFs
  featureDF1 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix1)
  featureDF2 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix2)

  #Check Seqnames
  featureDF1 <- .checkSeqnames(featureDF1, useMatrix1, useSeqnames1, matrixClass1)
  featureDF2 <- .checkSeqnames(featureDF2, useMatrix2, useSeqnames2, matrixClass2)

  #Create Match Names
  featureDF1$matchName <- featureDF1$name
  if("underscore" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\_.*","",featureDF1$matchName)
  }
  if("dash" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\-.*","",featureDF1$matchName)
  }
  if("numeric" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("[0-9]+","",featureDF1$matchName)
  }

  featureDF2$matchName <- featureDF2$name
  if("underscore" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\_.*","",featureDF2$matchName)
  }
  if("dash" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\-.*","",featureDF2$matchName)
  }
  if("numeric" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("[0-9]+","",featureDF2$matchName)
  }
  
  #Now Lets see how many matched pairings
  matchP <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
  if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
    stop("Matching of useMatrix1 and useMatrix2 resulted in no mappings!")
  }
  if(matchP < 0.05){
    if(force){
      stop("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Set force = TRUE to continue!")
    }else{
      message("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Continuing since force = TRUE.")
    }
  }
  matchedNames <- intersect(featureDF1$matchName, featureDF2$matchName)
  featureDF1m <- featureDF1[featureDF1$matchName %in% matchedNames, ]
  featureDF2m <- featureDF2[featureDF2$matchName %in% matchedNames, ]

  #Create Mappings
  mappingDF <- lapply(seq_len(nrow(featureDF1m)), function(x){
    expand.grid(x, which(paste0(featureDF2m$matchName) %in% paste0(featureDF1m$matchName[x])))
  }) %>% Reduce("rbind", .)

  #Test Mappings
  message("Testing ", nrow(mappingDF), " Mappings!")

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .messageDiffTime("Computing KNN", tstart)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k, method = knnMethod)

  #Determin Overlap
  .messageDiffTime("Identifying Non-Overlapping KNN pairs", tstart)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .messageDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), tstart)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Get Group Matrices
  .messageDiffTime("Getting Group Matrix 1", tstart)
  groupMat1 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF1m, 
    groupList = knnObj, 
    useMatrix = useMatrix1,
    threads = threads,
    verbose = FALSE
  )
 
  .messageDiffTime("Getting Group Matrix 2", tstart)
  groupMat2 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF2m, 
    groupList = knnObj, 
    useMatrix = useMatrix2,
    threads = threads,
    verbose = FALSE
  )

  #We need to divide by number of cells for the mean
  groupMat1 <- t(t(groupMat1) / k)
  groupMat2 <- t(t(groupMat2) / k)

  #Now we can normalize
  if(log2Norm1){
    if(any(groupMat1 < 0)){
      message("Some entries in groupMat1 are less than 0 continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.")
    }else{
      groupMat1 <- log2(groupMat1 + 1)
    }
  }
  if(log2Norm2){
    if(any(groupMat2 < 0)){
      message("Some entries in groupMat2 are less than 0 continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.")
    }else{
      groupMat2 <- log2(groupMat2 + 1)
    }
  }
  
  #Row Correlate
  rowTest <- .rowCorTest(
    X = groupMat1,
    Y = groupMat2,
    idxX = mappingDF[,1],
    idxY = mappingDF[,2],
    padjMethod = "bonferroni"
  )

  #Output DF
  colnames(featureDF1m) <- paste0(useMatrix1, "_", colnames(featureDF1m))
  colnames(featureDF2m) <- paste0(useMatrix2, "_", colnames(featureDF2m))

  df <- DataFrame(
    cbind(rowTest, featureDF1m[mappingDF[,1],,drop=FALSE], featureDF2m[mappingDF[,2],,drop=FALSE])
  )

  frontOrder <- c(paste0(useMatrix1, "_name"), paste0(useMatrix2, "_name"), "cor", "padj", "pval")
  df <- df[, c(frontOrder, colnames(df)[colnames(df) %ni% frontOrder])]


  return(df)

}

.rowCorTest <- function(
  X = NULL, 
  Y = NULL, 
  idxX = seq_row(X), 
  idxY = seq_row(Y), 
  padjMethod = "BH", 
  min = 10, 
  use="complete", 
  threads = 1
  ){
  message("Getting Correlations...")
  corTestList <- .safelapply(seq_along(idxX), function(i){
    if(i %% 1000 == 0){
      message("Computing Correlation (",i," of ", length(idxX), ")")
    }
    if(length(which(!is.na(X[idxX[i],]))) > min && length(which(!is.na(Y[idxY[i],]))) > min){
      corx <- .suppressAll(cor.test(X[idxX[i],],Y[idxY[i],],use=use))
      return(list(cor=corx$estimate[1], pval=corx$p.value[1]))
    }else{
      return(list(cor=NA, pval=NA))
    }
  }, threads = threads)
  corTest <- data.frame(
    cor = lapply(corTestList, function(x) x[[1]]) %>% unlist,
    pval = lapply(corTestList, function(x) x[[2]]) %>% unlist
  )
  rm(corTestList)
  corTest$padj <- p.adjust(corTest$pval, method=padjMethod)
  return(corTest)
}

.checkSeqnames <- function(featureDF = NULL, useMatrix = NULL, useSeqnames = NULL, matrixClass = NULL){

  seqnames <- unique(as.vector(featureDF$seqnames))
  useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
  if(length(useSeqnames)==0){
    useSeqnames <- NULL
  }

  if(!is.null(useSeqnames)){
    if(length(useSeqnames) == 1){
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }else{
      message("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
        "Continuing with first seqname '", seqnames[1], "'!\n",
        "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!")
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }else{
    if(matrixClass == "Sparse.Assays.Matrix"){
      if(all(seqnames %in% c("deviations", "z"))){
        seqnames <- c("z", "deviations")
      }
      message("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
        "Continuing with first seqname '", seqnames[1], "'!\n",
        "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!")
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }
  if(!(nrow(featureDF) > 1)){
    stop("Less than 1 feature is remaining in featureDF please check input!")
  }

  featureDF

}

##########################################################################################
# Co-accessibility Methods
##########################################################################################

#' Add Peak Co-Accessibility to an ArchRProject
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param knnMethod The method to be used for k-nearest neighbor computations. Options are "nabor", "RANN", and "FNN" and the corresponding package is required.
#' @param threads The number of threads to be used for parallel computing.
#' @export
addCoAccessibility <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 100000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1, 
  knnMethod = NULL,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = knnMethod, name = "knnMethod", valid = c("character", "null"))
  #JJJ .validInput(input = seed, name = "seed", valid = c("numeric"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  tstart <- Sys.time()

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .messageDiffTime("Computing KNN", tstart)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k, method = knnMethod)

  #Determin Overlap
  .messageDiffTime("Identifying Non-Overlapping KNN pairs", tstart)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .messageDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), tstart)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
  stopifnot(identical(chri,chrj))

  #Create Ranges
  peakSummits <- resize(peakSet, 1, "center")
  peakWindows <- resize(peakSummits, maxDist, "center")

  #Create Pairwise Things to Test
  o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
  o <- o[o[,1] != o[,2],]
  o$seqnames <- seqnames(peakSet)[o[,1]]
  o$idx1 <- peakSet$idx[o[,1]]
  o$idx2 <- peakSet$idx[o[,2]]
  o$correlation <- NA

  #Peak Matrix ColSums
  cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

  for(x in seq_along(chri)){
  
    .messageDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), tstart)

    #Features
    featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]),]
    featureDF$seqnames <- chri[x]

    #Group Matrix
    groupMat <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = featureDF, 
      groupList = knnObj, 
      useMatrix = "PeakMatrix",
      threads = threads,
      verbose = FALSE
    )
    
    #Scale
    groupMat <- t(t(groupMat) / gS) * scaleTo

    if(log2Norm){
      groupMat <- log2(groupMat + 1)
    }

    #Correlations
    idx <- BiocGenerics::which(o$seqnames==chri[x])
    o[idx,]$correlation <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = groupMat, Y = groupMat)

  }
  
  o$idx1 <- NULL
  o$idx2 <- NULL
  o <- o[!is.na(o$correlation),]
  mcols(peakSet) <- NULL
  o@metadata$peakSet <- peakSet

  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  ArchRProj

}

#' Get the peak co-accessibility from an ArchRProject
#' 
#' This function obtains co-accessibility data from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the co-accessibility signal as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export
getCoAccessibility <- function(
  ArchRProj = NULL, 
  corCutOff = 0.5, 
  resolution = 1000, 
  returnLoops = TRUE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){

    return(NULL)
  }

  
  if(is.null(metadata(ArchRProj@peakSet)$CoAccessibility)){
  
    return(NULL)
  
  }else{
   
    coA <- metadata(ArchRProj@peakSet)$CoAccessibility
    coA <- coA[coA$correlation >= corCutOff,,drop=FALSE]

    if(returnLoops){
      
      peakSummits <- resize(metadata(coA)$peakSet, 1, "center")

      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
      }
    
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[coA[,1]]),
        start = summitTiles[coA[,1]],
        end = summitTiles[coA[,2]]
      )
      mcols(loops)$value <- coA$correlation

      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))

      loops <- GenomicRangesList(CoAccessibility = loops)

      return(loops)

    }else{

      return(coA)

    }

  }

}

.constructGR <- function(
  seqnames = NULL, 
  start = NULL, 
  end = NULL, 
  ignoreStrand = TRUE
  ){
  .validInput(input = seqnames, name = "seqnames", valid = c("character", "rleCharacter"))
  .validInput(input = start, name = "start", valid = c("integer"))
  .validInput(input = end, name = "end", valid = c("integer"))
  .validInput(input = ignoreStrand, name = "ignoreStrand", valid = c("boolean"))
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignoreStrand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-" 
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)
}

##########################################################################################
# Peak2Gene Links Methods
##########################################################################################

ArchRProj = proj
reducedDims = "IterativeLSI"
dimsToUse = 1:30
scaleDims = NULL
corCutOff = 0.75
k = 100
knnIteration = 500
overlapCutoff = 0.8
maxDist = 250000
scaleTo = 10^4
log2Norm = TRUE
predictionCutoff = 0.4
seed = 1
knnMethod = NULL
threads = getArchRThreads()


#' Add Peak2GeneLinks to an ArchRProject JJJ
#' 
#' This function will add peak-to-gene links to a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param knnMethod The method to be used for k-nearest neighbor computations. Options are "nabor", "RANN", and "FNN" and the corresponding package is required.
#' @param threads The number of threads to be used for parallel computing.
#' @export
addPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  seed = 1, 
  knnMethod = NULL,
  threads = max(floor(getArchRThreads() / 2), 1)
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  #JJJ .validInput(input = predictionCutoff, name = "predictionCutoff", valid = c("numeric"))
  #JJJ .validInput(input = seed, name = "seed", valid = c("numeric"))
  .validInput(input = knnMethod, name = "knnMethod", valid = c("character", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  tstart <- Sys.time()
  .messageDiffTime("Getting Available Matrices", tstart)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)

  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }

  if("GeneIntegrationMatrix" %ni% AvailableMatrices){
    stop("GeneIntegrationMatrix not in AvailableMatrices")
  }

  ArrowFiles <- getArrowFiles(ArchRProj)

  tstart <- Sys.time()

  dfAll <- .safelapply(seq_along(ArrowFiles), function(x){
    DataFrame(
      cellNames = paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], "GeneIntegrationMatrix/Info/CellNames")),
      predictionScore = h5read(ArrowFiles[x], "GeneIntegrationMatrix/Info/predictionScore")
    )
  }, threads = threads) %>% Reduce("rbind", .)

  .messageDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
    sum(dfAll[,2] < predictionCutoff), 
    nrow(dfAll), 
    round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), tstart)

  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)

  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, "GeneIntegrationMatrix", threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .messageDiffTime("Computing KNN", tstart)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k, method = knnMethod)

  #Determin Overlap
  .messageDiffTime("Identifying Non-Overlapping KNN pairs", tstart)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .messageDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), tstart)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)

  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)

  #Group Matrix RNA
  .messageDiffTime("Getting Group RNA Matrix", tstart)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = "GeneIntegrationMatrix",
    threads = threads,
    verbose = TRUE
  )

  #Group Matrix ATAC
  .messageDiffTime("Getting Group ATAC Matrix", tstart)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = TRUE
  )

  .messageDiffTime("Normalizing Group Matrices", tstart)

  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }

  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj

  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC), 
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj

  rm(groupMatRNA, groupMatATAC)
  gc()

  #Overlaps
  .messageDiffTime("Finding Peak Gene Pairings", tstart)
  o <- DataFrame(
    findOverlaps(
      .suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )

  #Get Distance from Fixed point A B 
  o$distance <- distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")

  #Null Correlations
  #.messageDiffTime("Computing Background Correlations", tstart)
  #nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)

  .messageDiffTime("Computing Correlations", tstart)
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((1-o$Correlation^2)/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  #o$EmpPval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  #o$EmpFDR <- p.adjust(o$EmpPval, method = "fdr")
  
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart

  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  saveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  saveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA

  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out

  .messageDiffTime("Completed Peak2Gene Correlations!", tstart)

  ArchRProj

}

.getNullCorrelations <- function(seA, seB, o, n){

  o$seq <- seqnames(seA)[o$A]

  nullCor <- lapply(seq_along(unique(o$seq)), function(i){

    #Get chr from olist
    chri <- unique(o$seq)[i]
    #message(chri, " ", appendLF = FALSE)

    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)

    return(out)

  }) %>% SimpleList
  #message("")

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

}

#' Get the peak-to-gene links from an ArchRProject
#' 
#' This function obtains peak-to-gene links from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-gene correlation to return.
#' @param FDRCutOff A numeric describing the maximum numeric peak-to-gene false discovery rate to return.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the peak-to-gene links as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export
getPeak2GeneLinks <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  resolution = 1000, 
  returnLoops = TRUE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  #JJJ .validInput(input = FDRCutOff, name = "FDRCutOff", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    return(NULL)
  }

  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
  
    return(NULL)
  
  }else{
   
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]

    if(returnLoops){
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")

      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
    
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR

      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))

      loops <- GRangesList(Peak2GeneLinks = loops)

      return(loops)

    }else{

      return(p2g)

    }

  }

}

#' Plot Peak2Gene Heatmap from an ArchRProject
#' 
#' This function plots side by side heatmaps of linked ATAC and Gene regions from `addPeak2GeneLinks`.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-gene correlation used to define the set of peak-to-gene links to return.
#' @param FDRCutOff A numeric describing the maximum numeric peak-to-gene false discovery rate used to define the set of peak-to-gene links to return.
#' @param k An integer describing the number of k-means clusters to create when grouping peak-to-gene links prior to plotting heatmaps.
#' @param nPlot An integer describing the maximum number of peak-to-gene links to plot in the heatmap.
#' @param limitsATAC JJJ A numeric vector that determines the upper and lower limits to be used for the color scale for ATAC-seq data. Values below or above these
#' limits will be thresholded to the value of the specified limit.
#' @param limitsRNA JJJ A numeric vector that determines the upper and lower limits to be used for the color scale for RNA-seq data. Values below or above these
#' limits will be thresholded to the value of the specified limit.
#' @param groupBy JJJ THIS DOESNT MAKE SENSE TO ME. THE K GROUPS ARE GROUPS OF PEAKS RIGHT? The name of the column in `cellColData` to use for labeling KNN groupings. The maximum group appeared in the KNN groupings is used.
#' @param palGroup A color palette describing the colors in `groupBy`. For example, if `groupBy = "Clusters"` try `paletteDiscrete(ArchRProj$Clusters)` for a color palette.
#' @param palATAC A color palette describing the colors to be used for the ATAC heatmap. For example, `paletteContinuous("solarExtra")`.
#' @param palRNA A color palette describing the colors to be used for the RNA heatmap. For example, `paletteContinuous("blueYellow")`.
#' @export
peak2GeneHeatmap <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  k = 25,
  nPlot = 25000,
  limitsATAC = c(-2, 2),
  limitsRNA = c(-2, 2),
  groupBy = "Clusters",
  palGroup = NULL,
  palATAC = paletteContinuous("solarExtra"),
  palRNA = paletteContinuous("blueYellow")
  ){

  #JJJ .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  #JJJ .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  #JJJ .validInput(input = FDRCutOff, name = "FDRCutOff", valid = "numeric")
  #JJJ .validInput(input = k, name = "k", valid = c("integer"))
  #JJJ .validInput(input = nPlot, name = "nPlot", valid = "integer")
  #JJJ .validInput(input = limitsATAC, name = "limitsATAC", valid = "numeric")
  #JJJ .validInput(input = limitsRNA, name = "limitsRNA", valid = "numeric")
  #JJJ .validInput(input = groupBy, name = "groupBy", valid = "character")
  #JJJ .validInput(input = palGroup, name = "palGroup", valid = "JJJ")
  #JJJ .validInput(input = palATAC, name = "palATAC", valid = "JJJ")
  #JJJ .validInput(input = palRNA, name = "palRNA", valid = "JJJ")

  tstart <- Sys.time()

  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    stop("No Peak2GeneLinks Found! Try addPeak2GeneLinks!")
  }
  
  #########################################
  # Get Inputs
  #########################################
  ccd <- getCellColData(ArchRProj, select = groupBy)
  p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]

  mATAC <- assay(readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ])
  mRNA <- assay(readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ])
  gc()

  #########################################
  # Determine Groups from KNN
  #########################################
  .messageDiffTime("Determining KNN Groups!", tstart)
  KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
  KNNGroups <- lapply(seq_along(KNNList), function(x){
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
  }) %>% unlist
  cD <- DataFrame(row.names=paste0("K", seq_len(ncol(mATAC))), groupBy = KNNGroups)
  pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
  if(!is.null(palGroup)){
    pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
  }
  colorMap <- list(groupBy = pal)
  attr(colorMap[[1]], "discrete") <- TRUE

  #########################################
  # Organize Matrices
  #########################################
  mATAC <- .rowZscores(mATAC)
  mRNA <- .rowZscores(mRNA)
  rownames(mATAC) <- NULL
  rownames(mRNA) <- NULL
  colnames(mATAC) <- paste0("K", seq_len(ncol(mATAC)))
  colnames(mRNA) <- paste0("K", seq_len(ncol(mRNA)))
  rownames(mATAC) <- paste0("P2G", seq_len(nrow(mATAC)))
  rownames(mRNA) <- paste0("P2G", seq_len(nrow(mRNA)))

  .messageDiffTime("Ordering Peak2Gene Links!", tstart)
  k1 <- kmeans(mATAC, k)
  if(nrow(mATAC) > nPlot){
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      DataFrame(k = x, idx = sample(splitK[[x]], floor(nPK[x])))
    }) %>% Reduce("rbind", .)
  }else{
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
  }
  bS <- .binarySort(t(.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
  rowOrder <- rownames(bS[[1]])
  colOrder <- colnames(bS[[1]])
  kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))

  #########################################
  # Plot Heatmaps
  #########################################
  .messageDiffTime("Constructing ATAC Heatmap!", tstart)
  htATAC <- .ArchRHeatmap(
    mat = mATAC[kDF[,2],colOrder],#[rowOrder, colOrder],
    scale = FALSE,
    limits = limitsATAC,
    color = palATAC, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
  )

  .messageDiffTime("Constructing RNA Heatmap!", tstart)
  htRNA <- .ArchRHeatmap(
    mat = mRNA[kDF[,2],colOrder], #mRNA[rowOrder, colOrder],
    scale = FALSE,
    limits = limitsRNA,
    color = palRNA, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
  )

  htATAC + htRNA

}











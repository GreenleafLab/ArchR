##########################################################################################
# Assay Correlation Methods
##########################################################################################

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
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  seed = 1, 
  knnMethod = NULL,
  threads = getArchRThreads()
  ){

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
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
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
#' the `ArchRBrowser()` or as an `ArchRRegionTrack()`.
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
      loops <- loops[width(loops) >= resolution]
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

#' Add Peak2GeneLinks to an ArchRProject JJJ
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#' QQQ IM NOT CLEAR ON HOW THIS IS DIFFERENT FROM addCoAccessibility? THE FUNCTION DESCRIPTION IS IDENTICAL
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims QQQ DOUBLE CHECK A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
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
  maxDist = 150000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
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
  .validInput(input = threads, name = "threads", valid = c("integer"))

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

  b <- .getFeatureDF(ArrowFiles, "GeneScoreMatrix", threads = threads)

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

  #Create Ranges
  peakSummits <- resize(peakSet, 1, "center")
  geneWindows <- .suppressAll(resize(geneStart, maxDist, "center"))

  #Add IDs
  mcols(peakSummits)$uniqueID <- seq_along(peakSummits)
  mcols(geneWindows)$uniqueID <- seq_along(geneWindows)

  #Create Pairwise Things to Test
  o <- DataFrame(findOverlaps(peakSummits, geneWindows, ignore.strand = TRUE))
  o$seqnames <- seqnames(peakSet)[o[,1]]
  o$correlation <- NA
  o$idx1 <- mcols(peakSummits)$idx[o[,1]]
  o$idx2 <- mcols(geneWindows)$idx[o[,2]]

  #Get Null Correlations
  splitPeaks <- split(peakSummits, seqnames(peakSummits))
  splitGenes <- split(geneWindows, seqnames(geneWindows))

  oNull <- lapply(seq_along(chrij), function(x){
    chrx <- chrij[x]
    genesx <- splitGenes[[chrx]]
    peaksx <- splitPeaks[[chrx]]
    idx1 <- sample(seq_along(genesx), 50000, replace = TRUE)
    idx2 <- sample(seq_along(peaksx), 50000, replace = TRUE)
    idx3 <- which(distance(genesx[idx1], peaksx[idx1]) > maxDist)
    DataFrame(peakIdx = peaksx$uniqueID[idx2[idx3]], geneIdx = genesx$uniqueID[idx1[idx3]], seqnames = Rle(chrx), correlation = NA)
  }) %>% Reduce("rbind", .) 
  oNull <- oNull[order(oNull[,1], oNull[,2]), ]
  oNull$idx1 <- mcols(peakSummits)$idx[oNull[,1]]
  oNull$idx2 <- mcols(geneWindows)$idx[oNull[,2]]

  #Peak Matrix ColSums
  cS <- .getColSums(getArrowFiles(ArchRProj), chrij, verbose = FALSE, useMatrix = "PeakMatrix", threads = threads)
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))


  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- chrij[x]
  peakDF$seqnames <- chrij[x]

  #Group Matrix RNA
  message("RNA ", appendLF = FALSE)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = "GeneIntegrationMatrix",
    threads = threads,
    verbose = TRUE
  )

  #Group Matrix ATAC
  message("ATAC ", appendLF = FALSE)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = TRUE
  )

  for(x in seq_along(chrij)){
  
    .messageDiffTime(sprintf("Computing Peak-2-Gene Links %s (%s of %s)", chri[x], x, length(chri)), tstart)

    #Features
    geneDF <- mcols(geneStart)[BiocGenerics::which(seqnames(geneStart) == chrij[x]),]
    peakDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chrij[x]),]
    geneDF$seqnames <- chrij[x]
    peakDF$seqnames <- chrij[x]

    #Group Matrix RNA
    message("RNA ", appendLF = FALSE)
    groupMatRNA <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = geneDF, 
      groupList = knnObj, 
      useMatrix = "GeneIntegrationMatrix",
      threads = threads,
      verbose = FALSE
    )

    #Group Matrix ATAC
    message("ATAC ", appendLF = FALSE)
    groupMatATAC <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = peakDF, 
      groupList = knnObj, 
      useMatrix = "PeakMatrix",
      threads = threads,
      verbose = FALSE
    )
    
    #Scale
    groupMatATAC <- t(t(groupMatATAC) / gS) * scaleTo

    if(log2Norm){
      groupMatRNA  <- log2(groupMatRNA + 1)
      groupMatATAC <- log2(groupMatATAC + 1)
    }

    #Correlations
    message("Correlating", appendLF = FALSE)
    idx <- BiocGenerics::which(o$seqnames==chrij[x])
    o[idx,]$correlation <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = groupMatATAC, Y = groupMatRNA)

    idx <- BiocGenerics::which(oNull$seqnames==chrij[x])
    oNull[idx,]$correlation <- ArchR:::rowCorCpp(idxX = oNull[idx,]$idx1, idxY = oNull[idx,]$idx2, X = groupMatATAC, Y = groupMatRNA)

    message("")

    rm(groupMatATAC, groupMatRNA)
    gc()

  }
  
  o$idx1 <- NULL
  o$idx2 <- NULL
  o <- o[!is.na(o$correlation),]
  mcols(peakSet) <- NULL
  o@metadata$peakSet <- peakSet

  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  ArchRProj

}

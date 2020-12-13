##########################################################################################
# Assay Correlation Methods
##########################################################################################

#' Correlate Matrices within an ArchRProject
#' 
#' This function will correlate 2 matrices within an ArchRProject by name matching.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix1 A character describing the first matrix to use. See `getAvailableMatrices` for valid options.
#' @param useMatrix2 A character describing the second matrix to use. See `getAvailableMatrices` for valid options.
#' @param useSeqnames1 A character vector describing which seqnames to use in matrix 1.
#' @param useSeqnames2 A character vector describing which seqnames to use in matrix 2.
#' @param removeFromName1 A character vector describing how to filter names in matrix 1. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param removeFromName2 A character vector describing how to filter names in matrix 2. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param log2Norm1 A boolean describing whether to log2 normalize matrix 1.
#' @param log2Norm2 A boolean describing whether to log2 normalize matrix 2.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
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
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
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
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  seed = 1, 
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("correlateMatrices")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix1, name = "useMatrix1", valid = c("character"))
  .validInput(input = useMatrix2, name = "useMatrix2", valid = c("character"))
  .validInput(input = useSeqnames1, name = "useSeqnames1", valid = c("character", "null"))
  .validInput(input = useSeqnames2, name = "useSeqnames2", valid = c("character", "null"))
  .validInput(input = removeFromName1, name = "removeFromName1", valid = c("character", "null"))
  .validInput(input = removeFromName2, name = "removeFromName2", valid = c("character", "null"))
  .validInput(input = log2Norm1, name = "log2Norm1", valid = c("boolean"))
  .validInput(input = log2Norm2, name = "log2Norm2", valid = c("boolean"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "correlateMatrices Input-Parameters", logFile = logFile)

  set.seed(seed)

  #Get Available Matrices
  matrixNames <- getAvailableMatrices(ArchRProj)

  if(useMatrix1 %ni% matrixNames){
    .logStop(paste0("useMatrix1 (",useMatrix1,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")), logFile = logFile)
  }

  if(useMatrix2 %ni% matrixNames){
    .logStop(paste0("useMatrix2 (",useMatrix2,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")), logFile = logFile)
  }

  #Get Matrix Classes
  matrixClass1 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix1, "/Info/Class")))
  matrixClass2 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix2, "/Info/Class")))
  .logThis(matrixClass1, name = "matrixClass1", logFile = logFile)
  .logThis(matrixClass2, name = "matrixClass2", logFile = logFile)

  #Get Feature DFs
  featureDF1 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix1)
  featureDF2 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix2)
  .logThis(featureDF1, name = "featureDF1", logFile = logFile)
  .logThis(featureDF2, name = "featureDF2", logFile = logFile)

  #Check Seqnames
  featureDF1 <- .checkSeqnames(featureDF1, useMatrix1, useSeqnames1, matrixClass1, logFile)
  featureDF2 <- .checkSeqnames(featureDF2, useMatrix2, useSeqnames2, matrixClass2, logFile)

  #Create Match Names
  featureDF1$matchName <- toupper(featureDF1$name)
  if("underscore" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\_.*","",featureDF1$matchName)
  }
  if("dash" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\-.*","",featureDF1$matchName)
  }
  if("numeric" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("[0-9]+","",featureDF1$matchName)
  }
  if("dot" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\..*","",featureDF1$matchName)
  }

  featureDF2$matchName <- toupper(featureDF2$name)
  if("underscore" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\_.*","",featureDF2$matchName)
  }
  if("dash" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\-.*","",featureDF2$matchName)
  }
  if("numeric" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("[0-9]+","",featureDF2$matchName)
  }
  if("dot" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\..*","",featureDF2$matchName)
  }

  .logThis(featureDF1, name = "featureDF1", logFile = logFile)
  .logThis(featureDF2, name = "featureDF2", logFile = logFile)

  #Now Lets see how many matched pairings
  matchP1 <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
  matchP2 <- sum(featureDF2$matchName %in% featureDF1$matchName) / nrow(featureDF2)
  matchP <- max(matchP1, matchP2)

  .logThis(featureDF1$matchName, "featureDF1$matchName", logFile)
  .logThis(featureDF2$matchName, "featureDF2$matchName", logFile)

  if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
    .logStop("Matching of useMatrix1 and useMatrix2 resulted in no mappings!", logFile = logFile)
  }
  if(matchP < 0.05){
    if(force){
      .logStop("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Set force = TRUE to continue!", logFile = logFile)
    }else{
      .logMessage("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Continuing since force = TRUE.", verbose = TRUE, logFile = logFile)
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
  .logDiffTime(main=paste0("Testing ", nrow(mappingDF), " Mappings!"), t1=tstart, verbose=verbose, logFile=logFile)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(main=paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)
  .logThis(knnObj, name = "knnObj", logFile = logFile)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  .logThis(knnObj, name = "knnObjList", logFile = logFile)

  #Get Group Matrices
  .logDiffTime(main="Getting Group Matrix 1", t1=tstart, verbose=verbose, logFile=logFile)
  groupMat1 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF1m, 
    groupList = knnObj, 
    useMatrix = useMatrix1,
    threads = threads,
    verbose = FALSE
  )
 
  .logDiffTime(main="Getting Group Matrix 2", t1=tstart, verbose=verbose, logFile=logFile)
  groupMat2 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF2m, 
    groupList = knnObj, 
    useMatrix = useMatrix2,
    threads = threads,
    verbose = FALSE
  )

  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)

  #We need to divide by number of cells for the mean
  groupMat1 <- t(t(groupMat1) / k)
  groupMat2 <- t(t(groupMat2) / k)

  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)

  #Now we can normalize
  if(log2Norm1){
    if(any(groupMat1 < 0)){
      .logMessage("Some entries in groupMat1 are less than 0, continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.", logFile=logFile)
    }else{
      groupMat1 <- log2(groupMat1 + 1)
    }
  }
  if(log2Norm2){
    if(any(groupMat2 < 0)){
      .logMessage("Some entries in groupMat2 are less than 0, continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.", logFile=logFile)
    }else{
      groupMat2 <- log2(groupMat2 + 1)
    }
  }
 
  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)

  #Row Correlate
  rowTest <- .rowCorTest(
    X = groupMat1,
    Y = groupMat2,
    idxX = mappingDF[,1],
    idxY = mappingDF[,2],
    verbose = verbose,
    padjMethod = "bonferroni",
    logFile = logFile
  )
  .logThis(rowTest, name = "rowTest", logFile = logFile)

  #Output DF
  colnames(featureDF1m) <- paste0(useMatrix1, "_", colnames(featureDF1m))
  colnames(featureDF2m) <- paste0(useMatrix2, "_", colnames(featureDF2m))

  df <- DataFrame(
    cbind(rowTest, featureDF1m[mappingDF[,1],,drop=FALSE], featureDF2m[mappingDF[,2],,drop=FALSE])
  )

  frontOrder <- c(paste0(useMatrix1, "_name"), paste0(useMatrix2, "_name"), "cor", "padj", "pval")
  df <- df[, c(frontOrder, colnames(df)[colnames(df) %ni% frontOrder])]

  .endLogging(logFile = logFile)

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
  verbose = TRUE,
  threads = 1,
  logFile = NULL
  ){
  .logMessage("Getting Correlations...", verbose = verbose, logFile=logFile)
  corTestList <- .safelapply(seq_along(idxX), function(i){
    if(i %% 250 == 0){
      .logMessage("Computing Correlation (",i," of ", length(idxX), ")", logFile=logFile)
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

.checkSeqnames <- function(
  featureDF = NULL, 
  useMatrix = NULL, 
  useSeqnames = NULL, 
  matrixClass = NULL,
  logFile = NULL
  ){

  seqnames <- unique(as.vector(featureDF$seqnames))
  useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
  if(length(useSeqnames)==0){
    useSeqnames <- NULL
  }

  if(!is.null(useSeqnames)){
    if(length(useSeqnames) == 1){
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }else{
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
        "Continuing with first seqname '", seqnames[1], "'!\n",
        "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!", logFile=logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }else{
    if(matrixClass == "Sparse.Assays.Matrix"){
      if(all(seqnames %in% c("deviations", "z"))){
        seqnames <- c("z", "deviations")
      }
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
        "Continuing with first seqname '", seqnames[1], "'!\n",
        "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!", logFile=logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }
  if(!(nrow(featureDF) > 1)){
    .logStop("Less than 1 feature is remaining in featureDF please check input!", logFile=logFile)
  }

  featureDF

}


#' Correlate Trajectories
#' 
#' This function will correlate 2 trajectory matrices from getTrajectory.
#' 
#' @param seTrajectory1 A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param seTrajectory2 A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param corCutOff A numeric describing the cutoff for determining correlated features.
#' @param varCutOff1 The "Variance Quantile Cutoff" to be used for identifying the top variable features across `seTrajectory1`.
#' Only features with a variance above the provided quantile will be retained.
#' @param varCutOff2 The "Variance Quantile Cutoff" to be used for identifying the top variable features across `seTrajectory2`.
#' Only features with a variance above the provided quantile will be retained.
#' @param removeFromName1 A character vector describing how to filter names in matrix 1. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param removeFromName2 A character vector describing how to filter names in matrix 2. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param useRanges A boolean describing whether to use range overlap matching for correlation analysis.
#' @param fix1 A character describing where to resize the coordinates of `seTrajectory1`. Options include "start", "center", "end".
#' @param fix2 A character describing where to resize the coordinates of `seTrajectory2`. Options include "start", "center", "end".
#' @param maxDist A integer specifying the maximum distance between the coordinates of `seTrajectory1` and `seTrajectory2` for 
#' computing correlations.
#' @param log2Norm1 A boolean describing whether to log2 normalize `seTrajectory1`.
#' @param log2Norm2 A boolean describing whether to log2 normalize `seTrajectory2`.
#' @param force A boolean value that determines whether analysis should continue if resizing coordinates in `seTrajectory1` or 
#' `seTrajectory2` does not align with the strandedness. Only when `useRanges = TRUE`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
correlateTrajectories <- function(
  seTrajectory1 = NULL,
  seTrajectory2 = NULL,
  corCutOff = 0.5,
  varCutOff1 = 0.8,
  varCutOff2 = 0.8,
  removeFromName1 = c("underscore", "dash"),
  removeFromName2 = c("underscore", "dash"),
  useRanges = FALSE,
  fix1 = "center",
  fix2 = "start",
  maxDist = 250000,
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  force = FALSE,
  logFile = createLogFile("correlateTrajectories")
  ){

  .validInput(input = seTrajectory1, name = "seTrajectory1", valid = c("SummarizedExperiment"))
  .validInput(input = seTrajectory2, name = "seTrajectory2", valid = c("SummarizedExperiment"))
 	.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
	.validInput(input = varCutOff1, name = "varCutOff1", valid = c("numeric"))
	.validInput(input = varCutOff2, name = "varCutOff2", valid = c("numeric")) 
  .validInput(input = removeFromName1, name = "removeFromName1", valid = c("character", "null"))
  .validInput(input = removeFromName2, name = "removeFromName2", valid = c("character", "null"))
	.validInput(input = useRanges, name = "useRanges", valid = c("boolean"))
	.validInput(input = fix1, name = "fix1", valid = c("character"))
	.validInput(input = fix2, name = "fix2", valid = c("character"))
	.validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = log2Norm1, name = "log2Norm1", valid = c("boolean"))
  .validInput(input = log2Norm2, name = "log2Norm2", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "correlateTrajectories Input-Parameters", logFile=logFile)

  featureDF1 <- rowData(seTrajectory1)
  featureDF2 <- rowData(seTrajectory2)

  .logThis(featureDF1, "featureDF1", logFile = logFile)
  .logThis(featureDF2, "featureDF2", logFile = logFile)

  if("name" %in% colnames(featureDF1)){
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
  }else{
    if(!useRanges){
      .logStop("seTrajectory1 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
    }
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
  }

  if("name" %in% colnames(featureDF2)){
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
  }else{
    if(!useRanges){
      .logStop("seTrajectory2 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
    }
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
  }

  .logThis(rownames(featureDF1), "rownames(featureDF1)", logFile = logFile)
  .logThis(rownames(featureDF2), "rownames(featureDF2)", logFile = logFile)

  if(useRanges){

    if("start" %ni% colnames(featureDF1)){
      .logStop("start is not in seTrajectory1, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
    }

    if("start" %ni% colnames(featureDF2)){
      .logStop("start is not in seTrajectory2, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
    }

    if("strand" %in% colnames(featureDF1)){
      ranges1 <- GRanges(
        seqnames = featureDF1$seqnames, 
        IRanges(
          ifelse(featureDF1$strand == 2, featureDF1$end, featureDF1$start),
          ifelse(featureDF1$strand == 2, featureDF1$start, featureDF1$end)
        ),
        strand = ifelse(featureDF1$strand == 2, "-", "+")
      )
    }else{
      ranges1 <- GRanges(featureDF1$seqnames, IRanges(featureDF1$start, featureDF1$end))
    }
    #mcols(ranges1) <- featureDF1
    names(ranges1) <- rownames(featureDF1)
    rowRanges(seTrajectory1) <- ranges1
    rm(ranges1)

    if("strand" %in% colnames(featureDF2)){
      ranges2 <- GRanges(
        seqnames = featureDF2$seqnames, 
        IRanges(
          ifelse(featureDF2$strand == 2, featureDF2$end, featureDF2$start),
          ifelse(featureDF2$strand == 2, featureDF2$start, featureDF2$end)
        ),
        strand = ifelse(featureDF2$strand == 2, "-", "+")
      )
    }else{
      ranges2 <- GRanges(featureDF2$seqnames, IRanges(featureDF2$start, featureDF2$end))
    }
    #mcols(ranges2) <- featureDF2
    names(ranges2) <- rownames(featureDF2)
    rowRanges(seTrajectory2) <- ranges2
    rm(ranges2)

    .logThis(ranges1, "ranges1", logFile = logFile)
    .logThis(ranges2, "ranges2", logFile = logFile)

    #Find Associations to test
    isStranded1 <- any(as.integer(strand(seTrajectory1)) == 2)
    isStranded2 <- any(as.integer(strand(seTrajectory2)) == 2)

    if(fix1 == "center" & isStranded1){
      if(!force){
        .logStop("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }

    if(fix1 != "center" & !isStranded1){
      if(!force){
        .logStop("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }

    if(fix2 == "center" & isStranded2){
      if(!force){
        .logStop("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }

    if(fix2 != "center" & !isStranded2){
      if(!force){
        .logStop("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }

    #Overlaps
    mappingDF <- DataFrame(
      findOverlaps( 
        resize(rowRanges(seTrajectory1), 1, fix1), 
        .suppressAll(resize(resize(seTrajectory2, 1, fix2), 2 * maxDist + 1, "center")),
        ignore.strand = TRUE
      )
    )

    #Get Distance 
    mappingDF$distance <- distance(
      x = ranges(rowRanges(seTrajectory1)[mappingDF[,1]]), 
      y = ranges(rowRanges(resize(seTrajectory2, 1, fix2))[mappingDF[,2]])
    )

  }else{

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
    if("dot" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\..*","",featureDF1$matchName)
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
    if("dot" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\..*","",featureDF2$matchName)
    }  

    #Now Lets see how many matched pairings
    matchP1 <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
    matchP2 <- sum(featureDF2$matchName %in% featureDF1$matchName) / nrow(featureDF2)
    matchP <- max(matchP1, matchP2)

    .logThis(featureDF1$matchName, "featureDF1$matchName", logFile)
    .logThis(featureDF2$matchName, "featureDF2$matchName", logFile)

    if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
      .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!", logFile = logFile)
      stop("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!")
    }
    if(matchP < 0.05){
      if(!force){
        .logStop("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Set force = TRUE to continue!", logFile=logFile)
      }else{
        .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Continuing since force = TRUE.", verbose=TRUE, logFile=logFile)
      }
    }

    #Create Mappings
    mappingDF <- lapply(seq_len(nrow(featureDF1)), function(x){
      idx <- which(paste0(featureDF2$matchName) %in% paste0(featureDF1$matchName[x]))
      if(length(idx) > 0){
        expand.grid(x, idx)
      }else{
        NULL
      }
    }) %>% Reduce("rbind", .)

  }

  colnames(mappingDF)[1:2] <- c("idx1", "idx2")
  mappingDF <- DataFrame(mappingDF)

  .logThis(mappingDF, "mappingDF", logFile = logFile)

  if(!useRanges){
    mappingDF$matchname1 <- featureDF1$matchName[mappingDF$idx1]
    mappingDF$matchname2 <- featureDF2$matchName[mappingDF$idx2]
    mappingDF$name1 <- rownames(featureDF1)[mappingDF$idx1]
    mappingDF$name2 <- rownames(featureDF2)[mappingDF$idx2]
  }

  mappingDF$Correlation <- rowCorCpp(
    idxX = as.integer(mappingDF[,1]), 
    idxY = as.integer(mappingDF[,2]), 
    X = assays(seTrajectory1)[["mat"]], 
    Y = assays(seTrajectory2)[["mat"]]
  )
  mappingDF$VarAssay1 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory1)[["mat"]]))[as.integer(mappingDF[,1])]
  mappingDF$VarAssay2 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory2)[["mat"]]))[as.integer(mappingDF[,2])]
  mappingDF$TStat <- (mappingDF$Correlation / sqrt((1-mappingDF$Correlation^2)/(ncol(seTrajectory1)-2))) #T-statistic P-value
  mappingDF$Pval <- 2 * pt(-abs(mappingDF$TStat), ncol(seTrajectory1) - 2)
  mappingDF$FDR <- p.adjust(mappingDF$Pval, method = "fdr")

  idxPF <- which(mappingDF$Correlation > corCutOff & mappingDF$VarAssay1 > varCutOff1 & mappingDF$VarAssay2 > varCutOff2)
  .logMessage("Found ", length(idxPF), " Correlated Pairings!", logFile=logFile, verbose=TRUE)

  .logThis(mappingDF[idxPF,], "mappingDF-PF", logFile = logFile)

  out <- SimpleList(
    correlatedMappings = mappingDF[idxPF,],
    allMappings = mappingDF,
    seTrajectory1 = seTrajectory1,
    seTrajectory2 = seTrajectory2
  )

  out

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
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addCoAccessibility <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 100000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1, 
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("addCoAccessibility")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addCoAccessibility Input-Parameters", logFile = logFile)

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

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
  o$correlation <- -999.999
  o$Variability1 <- 0.000
  o$Variability2 <- 0.000

  #Peak Matrix ColSums
  cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

  for(x in seq_along(chri)){
  
    .logDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), t1=tstart, verbose=verbose, logFile=logFile)

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
    corVals <- rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
    .logThis(head(corVals), paste0("SubsetCorVals-", x), logFile = logFile)

    rowVars <- as.numeric(matrixStats::rowVars(groupMat))

    o[idx,]$correlation <- as.numeric(corVals)
    o[idx,]$Variability1 <- rowVars[o[idx,]$idx1]
    o[idx,]$Variability2 <- rowVars[o[idx,]$idx2]

    .logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
    .logThis(o[idx,], paste0("SubsetCoA-", x), logFile = logFile)

  }
  
  o$idx1 <- NULL
  o$idx2 <- NULL
  o <- o[!is.na(o$correlation),]

  o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), length(knnObj) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$VarQuantile1 <- .getQuantiles(o$Variability1)
  o$VarQuantile2 <- .getQuantiles(o$Variability2)

  mcols(peakSet) <- NULL
  o@metadata$peakSet <- peakSet

  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  .endLogging(logFile = logFile)

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
  resolution = 1, 
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
      metadata(coA) <- list()
      mcols(loops) <- coA[,-c(1:3)]
      mcols(loops)$value <- coA$correlation

      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))

      loops <- SimpleList(CoAccessibility = loops)

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

#' Add Peak2GeneLinks to an ArchRProject
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
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param addEmpiricalPval Add empirical p-values based on randomly correlating peaks and genes not on the same seqname.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  addEmpiricalPval = FALSE,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeak2GeneLinks Input-Parameters", logFile = logFile)

  .logDiffTime(main="Getting Available Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)

  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }

  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }

  ArrowFiles <- getArrowFiles(ArchRProj)

  tstart <- Sys.time()

  dfAll <- .safelapply(seq_along(ArrowFiles), function(x){
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e){
      message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(
      cellNames = cNx,
      predictionScore = pSx
    )
  }, threads = threads) %>% Reduce("rbind", .)

  .logDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
    sum(dfAll[,2] < predictionCutoff), 
    nrow(dfAll), 
    round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), t1=tstart, verbose=verbose, logFile=logFile)

  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  .logThis(peakSet, "peakSet", logFile = logFile)

  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  .logThis(geneStart, "geneStart", logFile = logFile)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

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
  .logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  rawMatRNA <- groupMatRNA
  .logThis(groupMatRNA, "groupMatRNA", logFile = logFile)

  #Group Matrix ATAC
  .logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  rawMatATAC <- groupMatATAC
  .logThis(groupMatATAC, "groupMatATAC", logFile = logFile)

  .logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)

  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }

  names(geneStart) <- NULL

  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)

  names(peakSet) <- NULL

  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj
  .logThis(seATAC, "seATAC", logFile = logFile)

  rm(groupMatRNA, groupMatATAC)
  gc()

  #Overlaps
  .logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
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
  if(addEmpiricalPval){
    .logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
    nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
  }

  .logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")  
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart

  if(addEmpiricalPval){
    out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }
  
  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  .safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  .safeSaveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA

  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out

  .logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  .endLogging(logFile = logFile)

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
#' @param varCutOffATAC A numeric describing the minimum variance quantile of the ATAC peak accessibility when selecting links.
#' @param varCutOffRNA A numeric describing the minimum variance quantile of the RNA gene expression when selecting links.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the peak-to-gene links as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export
getPeak2GeneLinks <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = FDRCutOff, name = "FDRCutOff", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
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

    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }

    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }

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

      loops <- SimpleList(Peak2GeneLinks = loops)

      return(loops)

    }else{

      return(p2g)

    }

  }

}

#' @export
peak2GeneHeatmap <- function(...){
    .Deprecated("plotPeak2GeneHeatmap")
    plotPeak2GeneHeatmap(...)
}

#' Plot Peak2Gene Heatmap from an ArchRProject
#' 
#' This function plots side by side heatmaps of linked ATAC and Gene regions from `addPeak2GeneLinks`.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-gene correlation to return.
#' @param FDRCutOff A numeric describing the maximum numeric peak-to-gene false discovery rate to return.
#' @param varCutOffATAC A numeric describing the minimum variance quantile of the ATAC peak accessibility when selecting links.
#' @param varCutOffRNA A numeric describing the minimum variance quantile of the RNA gene expression when selecting links.
#' @param k An integer describing the number of k-means clusters to group peak-to-gene links prior to plotting heatmaps.
#' @param nPlot An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param limitsATAC An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param limitsRNA An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param groupBy The name of the column in `cellColData` to use for labeling KNN groupings. The maximum group appeared in the KNN groupings is used.
#' @param palGroup A color palette describing the colors in `groupBy`. For example, if groupBy = "Clusters" try paletteDiscrete(ArchRProj$Clusters) for a color palette.
#' @param palATAC A color palette describing the colors to be used for the ATAC heatmap. For example, paletteContinuous("solarExtra").
#' @param palRNA A color palette describing the colors to be used for the RNA heatmap. For example, paletteContinuous("blueYellow").
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param returnMatrices A boolean value that determines whether the matrices should be returned with kmeans id versus plotting.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotPeak2GeneHeatmap <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  k = 25,
  nPlot = 25000,
  limitsATAC = c(-2, 2),
  limitsRNA = c(-2, 2),
  groupBy = "Clusters",
  palGroup = NULL,
  palATAC = paletteContinuous("solarExtra"),
  palRNA = paletteContinuous("blueYellow"),
  verbose = TRUE,
  returnMatrices = FALSE,
  seed = 1,
  logFile = createLogFile("plotPeak2GeneHeatmap")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = FDRCutOff, name = "FDRCutOff", valid = c("numeric"))
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = limitsATAC, name = "limitsATAC", valid = c("numeric"))
  .validInput(input = limitsRNA, name = "limitsRNA", valid = c("numeric"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = palGroup, name = "palGroup", valid = c("palette", "null"))
  .validInput(input = palATAC, name = "palATAC", valid = c("palette", "null"))
  .validInput(input = palRNA, name = "palRNA", valid = c("palette", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = returnMatrices, name = "returnMatrices", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "peak2GeneHeatmap Input-Parameters", logFile = logFile)
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    stop("No Peak2GeneLinks Found! Try addPeak2GeneLinks!")
  }
  
  #########################################
  # Get Inputs
  #########################################
  ccd <- getCellColData(ArchRProj, select = groupBy)
  p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]

  if(!is.null(varCutOffATAC)){
    p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
  }

  if(!is.null(varCutOffRNA)){
    p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
  }

  if(nrow(p2g) == 0){
    stop("No peak2genelinks found with your cutoffs!")
  }

  if(!file.exists(metadata(p2g)$seATAC)){
    stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  if(!file.exists(metadata(p2g)$seRNA)){
    stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
  mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  p2g$gene <- rowData(mRNA)$name
  gc()

  mATAC <- assay(mATAC)
  mRNA <- assay(mRNA)

  #########################################
  # Determine Groups from KNN
  #########################################
  .logDiffTime(main="Determining KNN Groups!", t1=tstart, verbose=verbose, logFile=logFile)
  KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
  KNNGroups <- lapply(seq_along(KNNList), function(x){
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
  }) %>% unlist
  cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups)
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
  colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
  colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
  rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
  rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
  rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))

  .logDiffTime(main="Ordering Peak2Gene Links!", t1=tstart, verbose=verbose, logFile=logFile)
  if(!is.null(seed)){
    set.seed(seed)
  }
  k1 <- kmeans(mATAC, k)
  if(nrow(mATAC) > nPlot){
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
  }else{
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
  }
  bS <- .binarySort(t(.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
  rowOrder <- rownames(bS[[1]])
  colOrder <- colnames(bS[[1]])
  kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))

  if(returnMatrices){

    out <- SimpleList(
        ATAC = SimpleList(
          matrix = mATAC[kDF[,2],colOrder],
          kmeansId = kDF[,3],
          colData = cD[colOrder,,drop=FALSE]
        ),
        RNA = SimpleList(
          matrix = mRNA[kDF[,2],colOrder],
          kmeansId = kDF[,3],
          colData = cD[colOrder,,drop=FALSE]
        ),
        Peak2GeneLinks = p2g[kDF[,2],]
      )

    return(out)

  }

  #Log Info
  .logThis(colorMap, "colorMap", logFile = logFile)
  .logThis(colOrder, "colOrder", logFile = logFile)
  .logThis(kDF, "kDF", logFile = logFile)
  .logThis(mATAC, "mATAC", logFile = logFile)
  .logThis(mRNA, "mRNA", logFile = logFile)
  .logThis(cD[colOrder,,drop=FALSE], "cD", logFile = logFile)
  .logThis(mATAC[kDF[,2],colOrder], "mATAC2", logFile = logFile)
  .logThis(mRNA[kDF[,2],colOrder], "mRNA2", logFile = logFile)

  #########################################
  # Plot Heatmaps
  #########################################
  .logDiffTime(main="Constructing ATAC Heatmap!", t1=tstart, verbose=verbose, logFile=logFile)
  htATAC <- .ArchRHeatmap(
    mat = mATAC[kDF[,2],colOrder],
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

  .logDiffTime(main = "Constructing RNA Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
  htRNA <- .ArchRHeatmap(
    mat = mRNA[kDF[,2],colOrder], 
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

  .endLogging(logFile = logFile)

  htATAC + htRNA

}











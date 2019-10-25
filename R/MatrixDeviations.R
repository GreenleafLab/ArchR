#' Add DeviationsMatrix to Arrow Files in ArchRProject
#' 
#' This function for each sample will independently compute counts for each tile
#' per cell and then infer gene activity scores.
#'
#' @param ArchRProj ArchRProject
#' @param annotations annotaions name stored in ArchRProject
#' @param matrixName matrixName to be stored as in Arrow Files
#' @param out save ouptut matrices deviations and/or z
#' @param binarize binarize peaks prior to computing deviations
#' @param threads number of threads for parallel execution
#' @param parallelParam parallel parameters for batch style execution
#' @param force force overwriting previous TileMatrix in ArrowFile
#' @export
addDeviationsMatrix <- function(
  ArchRProj,
  annotations = NULL,
  bdgPeaks = getBdgPeaks(ArchRProj),
  matrixName = NULL,
  out = c("z", "deviations"),
  binarize = FALSE,
  threads = 1,
  parallelParam = NULL,
  force = FALSE,
  ...
  ){

  .requirePackage("SummarizedExperiment")

  set.seed(1)
  tstart <- Sys.time()
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Error Needs to be ArchR Project for Input!")
  }
  ArrowFiles <- getSampleColData(ArchRProj)$ArrowFiles
  threads <- min(length(ArrowFiles), threads)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  ##############################################################
  #Annotations Matrix!
  ##############################################################
  anno <- getAnnotation(ArchRProj, annotations)
  annotationsMatrix <- SummarizedExperiment::assay(readRDS(anno$Matches))
  if(is.null(matrixName)){
    matrixName <- paste0(anno$Name, "Matrix")
  }
  annotationsMatrix <- as(annotationsMatrix, "dgCMatrix")
  rownames(annotationsMatrix) <- NULL
  gc()

  ##############################################################
  #Get Expectations!
  ##############################################################
  ArrowFiles <- getArrowFiles(ArchRProj)
  useMatrix <- "PeakMatrix"
  availableChr <- .availableSeqnames(ArrowFiles, useMatrix)
  
  rS <- suppressMessages(.getRowSums(
      ArrowFiles = ArrowFiles, 
      seqnames = availableChr,
      useMatrix = useMatrix,
      filter0 = FALSE
    ))
  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC

  #Create args list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())

  #Add args to list
  args$annotations <- NULL
  rm(annotations)
  args$annotationsMatrix <- annotationsMatrix
  args$featureDF <- rS
  args$useMatrix <- useMatrix
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName <- matrixName
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addDeviationsMatrix
  args$registryDir <- file.path(getOutputDirectory(ArchRProj), paste0(matrixName,"DeviationsRegistry"))

  #Run With Parallel or lapply
  outList <- .batchlapply(args)
  .messageDiffTime("Completed Computing Deviations!", tstart, addHeader = TRUE)
  gc()

  return(ArchRProj)

}

.addDeviationsMatrix <- function(
  i,
  ArrowFiles, 
  annotationsMatrix,
  out = c("z", "deviations"),
  cellNames = NULL, 
  allCells = NULL,
  featureDF  = NULL,
  bdgPeaks = NULL,
  binarize = FALSE,
  useMatrix = "PeakMatrix",
  matrixName = "Motif", 
  force = FALSE,
  profileMemory = TRUE,
  debug = FALSE,
  tstart = NULL,
  ...
  ){

  gc()

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  ArrowFile <- ArrowFiles[i]
  cellNames <- .availableCells(ArrowFile, subGroup=useMatrix)

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  #Get Matrix and Run ChromVAR!
  .messageDiffTime(sprintf("Computing chromVAR-based deviations %s of %s (see Schep et. al (2017)!", i, length(ArrowFiles)), tstart, addHeader = TRUE)
  dev <- .getMatFromArrow(
    ArrowFile, 
    featureDF = featureDF, 
    binarize = binarize, 
    useMatrix = useMatrix,
    cellNames = cellNames
    ) %>% {.customDeviations(
      countsMatrix = .,
      annotationsMatrix = annotationsMatrix,
      backgroudPeaks = bdgPeaks,
      expectation = featureDF$value/sum(featureDF$value),
      out = out
    )}
  gc()

  #######################################
  # Initialize Matrix Group
  #######################################
  if(length(out)==1){
    featureDF <- data.frame(seqnames = out, idx = seq_len(nrow(dev)), name = rownames(dev), stringsAsFactors = FALSE)
  }else if(length(out)==2){
    featureDF <- rbind(
      data.frame(seqnames = out[1], idx = seq_len(nrow(dev)), name = rownames(dev), stringsAsFactors = FALSE),
      data.frame(seqnames = out[2], idx = seq_len(nrow(dev)), name = rownames(dev), stringsAsFactors = FALSE)
    )
  }else{
    stop("out can only be up to 2 items deviations,z")
  }

  featureDF <- featureDF[order(featureDF[,1]),]
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = "Double",
    cellNames = colnames(dev),
    params = "chromVAR",
    featureDF = featureDF,
    force=TRUE
  )

  #######################################
  # Write Matrices To Arrow
  #######################################
  if("z" %in% tolower(out)){
    o <- .addMatToArrow(
      mat = as(SummarizedExperiment::assays(dev)[["z"]], "dgCMatrix"),
      ArrowFile = ArrowFile,
      Group = paste0(matrixName,"/z"),
      binarize = FALSE,
      addRowSums = FALSE,
      addColSums = FALSE,
      addRowVars = TRUE,
      addRowMeans = TRUE
    )
  }

  if("deviations" %in% tolower(out)){
    o <- .addMatToArrow(
      mat = as(SummarizedExperiment::assays(dev)[["deviations"]], "dgCMatrix"),
      ArrowFile = ArrowFile,
      Group = paste0(matrixName,"/deviations"),
      binarize = FALSE,
      addRowSums = FALSE,
      addColSums = FALSE,
      addRowVars = TRUE,
      addRowMeans = TRUE
    )    
  }

  .messageDiffTime("Finished Computing Deviations!", tstart)
  return(0)

}

############################################################################
# Adapted from chromVAR
############################################################################
.customDeviations <- function(
  countsMatrix,
  annotationsMatrix,
  backgroudPeaks,
  expectation,
  out = c("deviations", "z")
  ){

  tstart <- Sys.time()
  #lets not do this check because we are running on partial matrix
  #if (min(getFragmentsPerPeak(countsMatrix)) <= 0)
  #  stop("All peaks must have at least one fragment in one sample")
  stopifnot(nrow(countsMatrix) == nrow(backgroudPeaks))
  stopifnot(length(expectation) == nrow(countsMatrix))
  colData <- DataFrame(seq_len(ncol(countsMatrix)), row.names = colnames(countsMatrix))[,FALSE]
  norm_expectation <- expectation / sum(expectation) #Double check this sums to 1!
  countsPerSample <- Matrix::colSums(countsMatrix)

  results <- lapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% floor(ncol(annotationsMatrix)/20) == 0){
      .messageDiffTime(sprintf("Computing Deviations for Annotation %s of %s", x, ncol(annotationsMatrix)), tstart)
    }
    if(x %% max(floor(ncol(annotationsMatrix)/20), 10) == 0){
      gc()
    }
    .customDeviationsSingle(
      annotationsVector = annotationsMatrix[, x, drop=FALSE],
      countsMatrix = countsMatrix,
      backgroudPeaks = backgroudPeaks,
      countsPerSample = countsPerSample,
      expectation = norm_expectation,
      out = out
    )
  })
  cn <- colnames(countsMatrix)
  rm(countsMatrix)
  gc()

  if("z" %in% tolower(out)){
    z <- t(vapply(results, function(x) x[["z"]], rep(0, length(cn))))
  }else{
    z <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  if("deviations" %in% tolower(out)){
    dev <- t(vapply(results, function(x) x[["dev"]], rep(0, length(cn))))
  }else{
    dev <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  colnames(z) <- cn
  colnames(dev) <- cn

  #Check First
  nullOverlap <- is.null(results[[1]]$overlap)
  rowData <- lapply(seq_along(results), function(x){
      resx <- results[[x]]
      if(nullOverlap){
        data.frame(fractionMatches = resx$matches)
      }else{
        data.frame(fractionMatches = resx$matches, fractionBackgroundOverlap = resx$overlap)
      }
    }) %>% Reduce("rbind",.)
  rownames(rowData) <- colnames(annotationsMatrix)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      deviations = dev, 
      z = z
    ),
    colData = colData,
    rowData = rowData
  )
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[tolower(out)]
  
  return(se)

}

.customDeviationsSingle <- function(
  annotationsVector,
  countsMatrix,
  countsPerSample,
  backgroudPeaks,
  out = c("deviations", "z"),
  expectation = NULL,
  intermediate_results = FALSE,
  threshold = 1
  ){

  binarizeMat <- function(mat){
    mat@x[mat@x > 0] <- 1
    mat
  }

  if (length(annotationsVector@x) == 0) {
      out <- list(
        z = rep(NA, ncol(countsMatrix)), 
        dev = rep(NA, ncol(countsMatrix)), 
        expFG = NA,
        expBG = NA,
        matches = 0, 
        overlap = NA
        )
    return(out)
  }

  ################################
  # Fore Ground Deviations
  ################################
  .requirePackage("Matrix")
  observed <- as.vector(Matrix::t(annotationsVector) %*% countsMatrix)
  expected <- as.vector(Matrix::t(annotationsVector) %*% expectation %*% countsPerSample)
  observed_deviation <- (observed - expected)/expected

  #Filter those with no matches at all
  fail_filter <- which(expected == 0)
  
  ################################
  # Back Ground Deviations
  ################################
  if("z" %in% tolower(out)){

    #Compute Background Null Per Iteration
    niterations <- ncol(backgroudPeaks)
    sampleMat <- Matrix::sparseMatrix(
        j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
        i = rep(seq_len(niterations), each = length(annotationsVector@x)),
        x = rep(annotationsVector@x, niterations),
        dims = c(niterations, nrow(countsMatrix))
    )  
    sampled <- as.matrix(sampleMat %*% countsMatrix)
    sampledExpected <- sampleMat %*% expectation %*% countsPerSample
    sampledDeviation <- (sampled - sampledExpected)/sampledExpected
    bgOverlap <- Matrix::mean(binarizeMat(sampleMat) %*% binarizeMat(annotationsVector)) / length(annotationsVector@x)
    
    #Summary
    meanSampledDeviation <- Matrix::colMeans(sampledDeviation)
    sdSampledDeviation <- apply(as.matrix(sampledDeviation), 2, sd)

    #Norm Deviation
    normdev <- (observed_deviation - meanSampledDeviation)
    z <- normdev/sdSampledDeviation
    if (length(fail_filter) > 0) {
      z[fail_filter] <- NA
      normdev[fail_filter] <- NA
    }

  }else{

    #Compute Background Null Per Iteration
    niterations <- ncol(backgroudPeaks)
    sampleMat2 <- Matrix::sparseMatrix(
        j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
        i = rep(1, niterations * length(annotationsVector@x)),
        x = rep(annotationsVector@x, niterations),
        dims = c(1, nrow(countsMatrix))
    )
    sampled2 <- (sampleMat2 %*% countsMatrix)[1,]
    sampledExpected2 <- (sampleMat2 %*% expectation %*% countsPerSample)[1,]
    ######################
    # Equivalent to above
    # colMeans(sampled) - colMeans(sampledExpected))/colMeans(sampledExpected)
    ######################
    sampledDeviation2 <- (sampled2 - sampledExpected2)/sampledExpected2
    bgOverlap <- NA

    #Norm Deviation
    normdev <- (observed_deviation - sampledDeviation2)
    z <- NULL
    if (length(fail_filter) > 0) {
      normdev[fail_filter] <- NA
    }

  }

  outList <- list(
    z = z, 
    dev = normdev, 
    matches = length(annotationsVector@x) / nrow(countsMatrix), 
    overlap = bgOverlap
  )

  return(outList)

}

#' @export
addBdgPeaks <- function(
  ArchRProj, 
  niterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
  force = FALSE,
  ...){

  if(!is.null(metadata(getPeakSet(ArchRProj))$bdgPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bdgPeaks)){
      
      stop("Background Peaks Already Exist! set force = TRUE to addBdgPeaks!")

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bdgPeaks <- .getBdgPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = outFile)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set force = TRUE to addBdgPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bdgPeaks <- .getBdgPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = outFile)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bdgPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  metadata(ArchRProj@peakSet)$bdgPeaks <- outFile
  ArchRProj

}

#' @export
getBdgPeaks <- function(
  ArchRProj, 
  niterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  force = FALSE,
  ...){

  if(!is.null(metadata(getPeakSet(ArchRProj))$bdgPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bdgPeaks)){
      
      message("Using Previous Background Peaks!")
      bdgPeaks <- readRDS(metadata(getPeakSet(ArchRProj))$bdgPeaks)

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bdgPeaks <- .getBdgPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = NULL)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set add = TRUE to addBdgPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bdgPeaks <- .getBdgPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = NULL)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bdgPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  bdgPeaks

}

.getBdgPeaks <- function(
  ArchRProj, 
  niterations = 50,
  w = 0.1, 
  binSize = 50, 
  seed = 1, 
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
  ...){

  set.seed(1)
  .requirePackage("chromVAR")

  #Get Expectations
  ArrowFiles <- getArrowFiles(ArchRProj)
  useMatrix <- "PeakMatrix"
  availableChr <- .availableSeqnames(ArrowFiles, useMatrix)
  
  rS <- suppressMessages(.getRowSums(
      ArrowFiles = ArrowFiles, 
      seqnames = availableChr,
      useMatrix = useMatrix,
      filter0 = FALSE
    ))
  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC

  #minimal chromVAR change
  #chromVAR reuiqres a matrix/se of ncol > 1 and with a log10(values) transform removing peaks with 0 reads
  #to disable this we create a column of 1's forcing chromVAR to perform log10(values + 1)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = SimpleList(counts = as.matrix(data.frame(rS$value, 1))),
    rowData = DataFrame(bias = rS$GC)
  )

  bdgPeaks <- chromVAR::getBackgroundPeaks(
    object = se,
    bias = rowData(se)$bias, 
    niterations = niterations, 
    w = w, 
    bs = binSize
  )

  #Save Background Peaks
  if(!is.null(outFile)){
    saveRDS(bdgPeaks, outFile, compress = FALSE)
  }

  return(bdgPeaks)

}

#' @export
annotationEnrichment <- function(
  seMarker = NULL,
  ArchRProj = NULL,
  annotations = NULL,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.1",
  background = "bdgPeaks",
  return = "SE",
  ...){

  if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
    stop("Only markers identified from PeakMatrix can be used!")
  }

  matches <- getMatches(ArchRProj, annotations)
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  mcols(r1) <- NULL

  r2 <- getPeakSet(ArchRProj)
  mcols(r2) <- NULL

  if(!all.equal(paste0(r1), paste0(r2))){
    stop("Peaks from seMarker do not match peakSet in ArchRProj!")
  }

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(tolower(background) %in% c("backgroundpeaks", "bdgpeaks", "background", "bdg")){
    method <- "bdg"
    bdgPeaks <- getBdgPeaks(ArchRProj)
  }else{
    method <- "all"
  }

  enrichList <- lapply(seq_len(ncol(seMarker)), function(x){

    idx <- which(passMat[, x])

    if(method == "bdg"){
      .computeEnrichment(matches, idx, c(idx, as.vector(bdgPeaks[idx,])))
    }else{
      .computeEnrichment(matches, idx, seq_len(nrow(bdgPeaks)))
    }

  }) %>% SimpleList
  names(enrichList) <- colnames(seMarker)


  if(tolower(return)=="se"){

    assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
      d <- lapply(seq_along(enrichList), function(y){
        enrichList[[y]][colnames(matches),x,drop=FALSE]
      }) %>% Reduce("cbind",.)
      colnames(d) <- names(enrichList)
      d
    }) %>% SimpleList
    names(assays) <- colnames(enrichList[[1]])
    assays <- rev(assays)
    se <- SummarizedExperiment::SummarizedExperiment(assays=assays)

  }else{

    enrichList
    
  }

}

.computeEnrichment <- function(matches, compare, background){

  matches <- .getAssay(matches, "matches")
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
     pOut$BackgroundFrequency[x], # Number of all successes in background
     pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
     pOut$nCompare[x], # Number that were drawn
     lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)

  #Minus Log10 FDR
  pOut$mlog10FDR <- -log10(p.adjust(matrixStats::rowMaxs(cbind(10^-pOut$mlog10p, 4.940656e-324)), method = "fdr"))
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

  pOut

}





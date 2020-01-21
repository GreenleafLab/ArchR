####################################################################
# Transcription Factor Deviation Methods
####################################################################

#' Add a matrix of deviations for a given peakAnnotation to Arrow Files in ArchRProject
#' 
#' This function will compute peakAnnotation deviations for each ArrowFiles independently while controlling for global biases (low-memory requirement).
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param peakAnnotation The name of the `peakAnnotation` stored in the `ArchRProject`.
#' @param matrixName The name to be used for storage of the deviations matrix in the provided `ArchRProject`.
#' @param out A string or character vector that indicates whether to save the ouptut matrices as deviations ("deviations") z-scores ("z"), or both (c("deviations","z")).
#' @param binarize A boolean value indicating whether the input matrix should be binarized before calculating deviations. This is often desired when working with insertion counts.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exists in the ArrowFiles asociated with the given `ArchRProject`.
#' @export
addDeviationsMatrix <- function(
  ArchRProj = NULL,
  peakAnnotation = NULL,
  matches = NULL,
  bgdPeaks = getBgdPeaks(ArchRProj),
  matrixName = NULL,
  out = c("z", "deviations"),
  binarize = FALSE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = FALSE,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = peakAnnotation, name = "peakAnnotation", valid = c("character", "null"))
  .validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", "null"))
  .validInput(input = bgdPeaks, name = "bgdPeaks", valid = c("GRanges"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = out, name = "out", valid = c("character"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))

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
  print(matches)
  if(is.null(matches)){
    anno <- getPeakAnnotation(ArchRProj, peakAnnotation)
    matches <- readRDS(anno$Matches)
    if(is.null(matrixName)){
      matrixName <- paste0(anno$Name, "Matrix")
    }
  }else{
    if(is.null(matrixName)){
      matrixName <- paste0("MotifMatrix")
    }    
  }
  annotationsMatrix <- SummarizedExperiment::assay(matches)
  rownames(annotationsMatrix) <- paste0(seqnames(matches), "_", start(matches), "_", end(matches))
  annotationsMatrix <- as(annotationsMatrix, "dgCMatrix")
  #rm(matches)
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
  rownames(rS) <- paste0(rS$seqnames,"_",rS$idx)
  rS <- rS[paste0(seqnames(ArchRProj@peakSet), "_", mcols(ArchRProj@peakSet)$idx),]
  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC
  rownames(rS) <- paste0(rS$seqnames, "_", rS$start, "_", rS$end)

  annotationsMatrix <- annotationsMatrix[rownames(rS), ]

  #Create args list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())

  #Add args to list
  args$peakAnnotation <- NULL
  rm(peakAnnotation)

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
  bgdPeaks = NULL,
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
  sampleName <- .sampleName(ArrowFile)
  prefix <- sprintf("%s (%s of %s)", sampleName, i, length(ArrowFiles))
  cellNames <- .availableCells(ArrowFile, subGroup=useMatrix)

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  #Get Matrix and Run ChromVAR!
  .messageDiffTime(sprintf("chromVAR deviations %s Schep (2017)", prefix), tstart, addHeader = TRUE)
  dev <- .getMatFromArrow(
    ArrowFile, 
    featureDF = featureDF, 
    binarize = binarize, 
    useMatrix = useMatrix,
    cellNames = cellNames
    ) %>% {.customDeviations(
      countsMatrix = .,
      annotationsMatrix = annotationsMatrix,
      prefix = prefix,
      backgroudPeaks = SummarizedExperiment::assay(bgdPeaks),
      expectation = featureDF$rowSums/sum(featureDF$rowSums),
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
    Class = "Assays",
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
  prefix = "",
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

  d <- max(floor(ncol(annotationsMatrix)/20), 1)
  results <- lapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% d == 0){
      .messageDiffTime(sprintf("%s : Deviations for Annotation %s of %s", prefix, x, ncol(annotationsMatrix)), tstart)
    }
    if(x %% max(c(d, 10)) == 0){
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

#' Get Variable Deviations across cells in ArchRProject.
#' 
#' This function will rank the variability of the deviations computed by ArchR and label the top variable annotations.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the `DeviationsMatrix` object stored in the `ArchRProject`. See `addDeviationsMatrix()`.
#' @param plot QQQ WHAT IS AN "ANNOTATION" HERE? A boolean value indicating whether the ranked variability should be plotted for each QQQ annotation.
#' @param n The number of annotations to label with `ggrepel`.
#' @export
getVarDeviations <- function(ArchRProj, name = "MotifMatrix", plot = TRUE, n = 25){

  rowV <- .getRowVars(getArrowFiles(ArchRProj), useMatrix = name, seqnames = "z")
  rowV <- rowV[order(rowV$combinedVars, decreasing=TRUE), ]
  rowV$rank <- seq_len(nrow(rowV))

  if(plot){
    rowV <- data.frame(rowV)
    ggplot(rowV, aes(rank, combinedVars)) +
      geom_point(size = 1) + 
      ggrepel::geom_label_repel(
        data = rowV[rev(seq_len(n)), ], aes(x = rank, y = combinedVars, label = name), 
        size = 1.5,
        nudge_x = 2
      ) + theme_ArchR() + ylab("Variability") + xlab("Rank Sorted Annotations")
  }else{

    rowV

  }

}

#' Add Background Peaks to an ArchRProject
#' 
#' This function will compute background peaks controlling for total accessibility and GC-content and add this information to an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param niterations QQQ SHOULD THIS BE "nIterations" FOR CONSISTENCY? The number of background peaks to sample. See `chromVAR::getBackgroundPeaks()`.
#' @param w QQQ I FEEL LIKE THESE PARAMETERS FOR CHROMVAR NEED TO BE BETTER DESCRIBED. The parameter controlling similarity of background peaks. See `chromVAR::getBackgroundPeaks()`.
#' @param binSize QQQ I FEEL LIKE THESE PARAMETERS FOR CHROMVAR NEED TO BE BETTER DESCRIBED. The precision with which the similarity is computed. See `chromVAR::getBackgroundPeaks()`.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param outFile QQQ The path to save the backgroundPeaks object as a `.RDS` file for the given `ArchRProject`. The default action is to save this file in the `outDir` of the `ArchRProject`.
#' @param force A boolean value indicating whether to force the file indicated by `outFile` to be overwritten if it already exists.
#' @export
addBgdPeaks <- function(
  ArchRProj, 
  niterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
  force = FALSE,
  ...){

  if(!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bgdPeaks)){
      
      stop("Background Peaks Already Exist! set force = TRUE to addBgdPeaks!")

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = outFile)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set force = TRUE to addBgdPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = outFile)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bgdPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  metadata(ArchRProj@peakSet)$bgdPeaks <- outFile
  ArchRProj

}

#' QQQ NEEDS PARAM DEFINITIONS OR NEEDS TO BE HIDDEN. IF YOU MAKE IT HIDDEN, MAKE SURE TO CHANGE THE FUNCTION CALL ABOVE
#' @export
getBgdPeaks <- function(
  ArchRProj, 
  niterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  force = FALSE,
  ...
  ){

  if(!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bgdPeaks)){
      
      message("Using Previous Background Peaks!")
      bgdPeaks <- readRDS(metadata(getPeakSet(ArchRProj))$bgdPeaks)

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = NULL)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set add = TRUE to addBgdPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, niterations=niterations, w=w, binSize=binSize, seed = seed, outFile = NULL)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bgdPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  bgdPeaks

}

.computeBgdPeaks <- function(
  ArchRProj, 
  niterations = 50,
  w = 0.1, 
  binSize = 50, 
  seed = 1, 
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
  ...
  ){

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
    assays = SimpleList(counts = as.matrix(data.frame(rS$rowSums, 1))),
    rowData = DataFrame(bias = rS$GC)
  )

  bgdPeaks <- chromVAR::getBackgroundPeaks(
    object = se,
    bias = rowData(se)$bias, 
    niterations = niterations, 
    w = w, 
    bs = binSize
  )

  bgdPeaks <- SummarizedExperiment(assays = SimpleList(bgdPeaks = bgdPeaks), 
    rowRanges = GRanges(rS$seqnames,IRanges(rS$start,rS$end),value=rS$rowSums,GC=rS$GC))

  #Save Background Peaks
  if(!is.null(outFile)){
    saveRDS(bgdPeaks, outFile, compress = FALSE)
  }

  return(bgdPeaks)

}


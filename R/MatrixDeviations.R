####################################################################
# Transcription Factor Deviation Methods
####################################################################

#' Add a matrix of deviations for a given peakAnnotation to Arrow Files in ArchRProject
#' 
#' This function will compute peakAnnotation deviations for each ArrowFiles independently while controlling for global biases (low-memory requirement).
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param peakAnnotation The name of the `peakAnnotation` stored in the `ArchRProject`.
#' @param matches A custom `peakAnnotation` matches object used as input for the hypergeometric test. See
#' `motifmatchr::matchmotifs()` for additional information.
#' @param bgdPeaks A `SummarizedExperiment` that contains for each peak a set of background peaks matched by biases such as total accessibility
#' and GC nucleotide content. This can be computed using `addBgdPeaks` and accessed by `getBgdPeaks`.
#' @param matrixName The name to be used for storage of the deviations matrix in the provided `ArchRProject`.
#' @param out A string or character vector that indicates whether to save the ouptut matrices as deviations ("deviations")
#' z-scores ("z"), or both (c("deviations","z")).
#' @param binarize A boolean value indicating whether the input matrix should be binarized before calculating deviations.
#' This is often desired when working with insertion counts.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it
#' already exists in the ArrowFiles associated with the given `ArchRProject`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addDeviationsMatrix <- function(
  ArchRProj = NULL,
  peakAnnotation = NULL,
  matches = NULL,
  bgdPeaks = getBgdPeaks(ArchRProj, method = "chromVAR"),
  matrixName = NULL,
  out = c("z", "deviations"),
  binarize = FALSE,
  threads = getArchRThreads(),
  verbose = TRUE,
  parallelParam = NULL,
  force = FALSE,
  logFile = createLogFile("addDeviationsMatrix")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = peakAnnotation, name = "peakAnnotation", valid = c("character", "null"))
  .validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", "null"))
  .validInput(input = bgdPeaks, name = "bgdPeaks", valid = c("SummarizedExperiment"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character", "null"))
  .validInput(input = out, name = "out", valid = c("character"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)

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
  args <- mget(names(formals()),sys.frame(sys.nframe()))

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
  args$logFile <- logFile
  args$registryDir <- file.path(getOutputDirectory(ArchRProj), paste0(matrixName,"DeviationsRegistry"))

  #Remove Project from Args
  args$ArchRProj <- NULL
  args$matches <- NULL

  .logThis(args, "Deviations-Args", logFile = logFile)

  #Run With Parallel or lapply
  outList <- .batchlapply(args)
  .logDiffTime("Completed Computing Deviations!", tstart, addHeader = TRUE, logFile = logFile)

  .endLogging(logFile = logFile)
  gc()

  return(ArchRProj)

}

.addDeviationsMatrix <- function(
  i = NULL,
  ArrowFiles = NULL, 
  annotationsMatrix = NULL,
  out = c("z", "deviations"),
  cellNames = NULL, 
  allCells = NULL,
  featureDF  = NULL,
  bgdPeaks = NULL,
  binarize = FALSE,
  useMatrix = "PeakMatrix",
  matrixName = "Motif", 
  force = FALSE,
  verbose = TRUE,
  tstart = NULL,
  subThreads = 1,
  logFile = NULL
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

  #Check
  # completed <- tryCatch({
  #     h5read(ArrowFile, paste0(matrixName,"/Info/Completed")) #Check if completed
  #     return(TRUE)
  #   },error = function(y){
  #     return(FALSE)
  # })
  # if(completed){
  #   if(!force){
  #     .logMessage(paste0("Previous Run Completed Successfully, to overwrite set force = TRUE! Skipping ", sampleName, " Deviations."), logFile = logFile)
  #     message("Previous Run Completed Successfully, to overwrite set force = TRUE! Skipping ", sampleName, " Deviations.")
  #     return(0)
  #   }else{
  #     .logMessage("Previous Run Completed Successfully, continuing since force = TRUE!", logFile = logFile)
  #     message("Previous Run Completed Successfully, continuing since force = TRUE!")
  #   }
  # }
  o <- h5closeAll()
  o <- .createArrowGroup(ArrowFile = ArrowFile, group = matrixName, force = force, logFile = logFile)

  #Get Matrix and Run ChromVAR!
  .logDiffTime(sprintf("chromVAR deviations %s Schep (2017)", prefix), tstart, addHeader = FALSE, logFile = logFile)
  dev <- tryCatch({
    
    .getMatFromArrow(
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
        out = out,
        verbose = verbose,
        threads = subThreads,
        logFile = logFile
      )}
  
  }, error = function(e){

    errorList <- list(
      annotationsMatrix = annotationsMatrix,
      prefix = prefix,
      backgroudPeaks = SummarizedExperiment::assay(bgdPeaks),
      expectation = featureDF$rowSums/sum(featureDF$rowSums),
      out = out,
      verbose = verbose,
      threads = subThreads,
      logFile = logFile
    )
    
    .logError(e, fn = ".computeDeviations", info = prefix, errorList = errorList, logFile = logFile)

  })
  gc()
  .logThis(dev, "dev", logFile = logFile)

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
  .logThis(featureDF, "featureDF", logFile = logFile)

  Units <- c()
  if("z" %in% tolower(out)){
    Units <- c(Units, "DeviationScores")
  }
  if("deviations" %in% tolower(out)){
    Units <- c(Units, "Deviations")
  }

  #Initialize
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = "Assays",
    Units = Units,
    cellNames = colnames(dev),
    params = "chromVAR",
    featureDF = featureDF,
    force = TRUE
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

  #Add Completion Mark
  o <- h5write(obj = "Finished", file = ArrowFile, name = paste0(matrixName,"/Info/Completed"))

  .logDiffTime("Finished Computing Deviations!", tstart)

  return(0)

}

############################################################################
# Adapted from chromVAR, Approved by Alicia Schep for Modification.
############################################################################
.customDeviations <- function(
  countsMatrix = NULL,
  annotationsMatrix = NULL,
  backgroudPeaks = NULL,
  expectation = NULL,
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE,
  logFile = NULL
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

  .logThis(countsMatrix, paste0(prefix, " : CountsMatrix"))
  .logThis(annotationsMatrix, paste0(prefix, " : annotationsMatrix"))
  .logThis(backgroudPeaks, paste0(prefix, " : backgroudPeaks"))
  .logThis(expectation, paste0(prefix, " : expectation"))

  d <- max(floor(ncol(annotationsMatrix)/20), 1)
  m <- 0
  results <- .safelapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% d == 0){
      .logDiffTime(sprintf("%s : Deviations for Annotation %s of %s", prefix, x, ncol(annotationsMatrix)), tstart, verbose = verbose, logFile = logFile)
      m <- 1 #Print to console
    }
    if(x %% max(floor(d/5), 2) == 0){
      if(m != 1){
        .logDiffTime(sprintf("%s : Deviations for Annotation %s of %s", prefix, x, ncol(annotationsMatrix)), tstart, verbose = FALSE, logFile = logFile)
      }else{
        m <- 0 #Reset
      }
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
      out = out,
      logFile = logFile,
      prefix = prefix
    )
  }, threads = threads)
  cn <- colnames(countsMatrix)
  rm(countsMatrix)
  gc()

  if("z" %in% tolower(out)){
    z <- t(vapply(results, function(x) x[["z"]], rep(0, length(cn))))
    if(length(cn)==1){
      z <- matrix(z, ncol=length(cn))
    }
  }else{
    z <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  if("deviations" %in% tolower(out)){
    dev <- t(vapply(results, function(x) x[["dev"]], rep(0, length(cn))))
    if(length(cn)==1){
      dev <- matrix(dev, ncol=length(cn))
    }
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
  annotationsVector = NULL,
  countsMatrix = NULL,
  countsPerSample = NULL,
  backgroudPeaks = NULL,
  out = c("deviations", "z"),
  expectation = NULL,
  threshold = 1,
  logFile = NULL,
  prefix = ""
  ){

  .binarizeMat <- function(mat = NULL){
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

  outList <- tryCatch({

    ################################
    # Fore Ground Deviations
    ################################
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
      bgOverlap <- Matrix::mean(.binarizeMat(sampleMat) %*% .binarizeMat(annotationsVector)) / length(annotationsVector@x)
      
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

    outList

  }, error = function(e){

    errorList <- list(
      annotationsVector = annotationsVector,
      observed = if(exists("observed", inherits = FALSE)) observed else "observed",
      expected = if(exists("expected", inherits = FALSE)) expected else "expected",
      sampleMat = if(exists("sampleMat", inherits = FALSE)) sampleMat else "sampleMat",
      sampleMat2 = if(exists("sampleMat", inherits = FALSE)) sampleMat2 else "sampleMat2",
      sampledDeviation = if(exists("sampledDeviation", inherits = FALSE)) sampledDeviation else "sampledDeviation",
      sampledDeviation2 = if(exists("sampledDeviation2", inherits = FALSE)) sampledDeviation2 else "sampledDeviation2",
      normdev = if(exists("normdev", inherits = FALSE)) normdev else "normdev",
      z = if(exists("z", inherits = FALSE)) z else "z"
    )

    .logError(e, fn = ".customDeviationsSingle", info = prefix, errorList = errorList, logFile = logFile)

  })

  return(outList)

}

#' Get Variable Deviations across cells in ArchRProject.
#' 
#' This function will rank the variability of the deviations computed by ArchR and label the top variable annotations.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the `DeviationsMatrix` object stored in the `ArchRProject`. See `addDeviationsMatrix()`.
#' @param plot A boolean value indicating whether the ranked variability should be plotted for each peakAnnotation in `DeviationsMatrix`.
#' @param n The number of annotations to label with `ggrepel`.
#' @export
getVarDeviations <- function(ArchRProj = NULL, name = "MotifMatrix", plot = TRUE, n = 25){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = plot, name = "plot", valid = c("boolean"))
  .validInput(input = n, name = "n", valid = c("integer"))

  rowV <- .getRowVars(getArrowFiles(ArchRProj), useMatrix = name, seqnames = "z")
  rowV <- rowV[order(rowV$combinedVars, decreasing=TRUE), ]
  rowV$rank <- seq_len(nrow(rowV))

  print(head(rowV))

  if(plot){
    rowV <- data.frame(rowV)
    ggplot(rowV, aes(x = rank, y = combinedVars, color = combinedVars)) +
      geom_point(size = 1) + 
      scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
      ggrepel::geom_label_repel(
        data = rowV[rev(seq_len(n)), ], aes(x = rank, y = combinedVars, label = name), 
        size = 1.5,
        color = "black",
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
#' @param nIterations The number of background peaks to sample. See `chromVAR::getBackgroundPeaks()`.
#' @param w The parameter controlling similarity of background peaks. See `chromVAR::getBackgroundPeaks()`.
#' @param binSize The precision with which the similarity is computed. See `chromVAR::getBackgroundPeaks()`.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that
#' you can reproduce results downstream.
#' @param outFile The path to save the `backgroundPeaks` object as a `.RDS` file for the given `ArchRProject`. The default action
#' is to save this file in the `outputDirectory` of the `ArchRProject`.
#' @param method A string indicating whether to use chromVAR or ArchR for background peak identification.
#' @param force A boolean value indicating whether to force the file indicated by `outFile` to be overwritten if it already exists.
#' @export
addBgdPeaks <- function(
  ArchRProj = NULL, 
  nIterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  method  = "chromVAR",
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = nIterations, name = "nIterations", valid = c("integer"))
  .validInput(input = w, name = "w", valid = c("numeric"))
  .validInput(input = binSize, name = "binSize", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = outFile, name = "outFile", valid = c("character"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  if ("PeakMatrix" %ni% getAvailableMatrices(ArchRProj)) {
    .logMessage(paste0("PeakMatrix does not exist in the provided ArchRProject. Add a peak matrix using addPeakMatrix(). See available matrix names from getAvailableMatrices()!"), logFile = logFile)
    stop("PeakMatrix does not exist in the provided ArchRProject. Add a peak matrix using addPeakMatrix(). See available matrix names from getAvailableMatrices()!")
  }
  
  if(!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bgdPeaks)){
      
      stop("Background Peaks Already Exist! set force = TRUE to addBgdPeaks!")

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, nIterations=nIterations, w=w, binSize=binSize, seed = seed, outFile = outFile, method = method)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set force = TRUE to addBgdPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, nIterations=nIterations, w=w, binSize=binSize, seed = seed, outFile = outFile, method = method)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bgdPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  metadata(ArchRProj@peakSet)$bgdPeaks <- outFile
  ArchRProj

}


#' Get Background Peaks from an ArchRProject
#' 
#' This function will get/compute background peaks controlling for total accessibility and GC-content from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param nIterations The number of background peaks to sample. See `chromVAR::getBackgroundPeaks()`.
#' @param w The parameter controlling similarity measure of background peaks. See `chromVAR::getBackgroundPeaks()`.
#' @param binSize The precision with which the similarity is computed. See `chromVAR::getBackgroundPeaks()`.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used
#' so that you can reproduce results downstream.
#' @param method A string indicating whether to use chromVAR or ArchR for background peak identification.
#' @param force A boolean value indicating whether to force the file indicated by `outFile` to be overwritten if it already exists.
#' @export
getBgdPeaks <- function(
  ArchRProj = NULL, 
  nIterations = 50, 
  w = 0.1, 
  binSize = 50, 
  seed = 1,
  method  = "chromVAR",
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = nIterations, name = "nIterations", valid = c("integer"))
  .validInput(input = w, name = "w", valid = c("numeric"))
  .validInput(input = binSize, name = "binSize", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  if(!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks) & !force){
    
    if(file.exists(metadata(getPeakSet(ArchRProj))$bgdPeaks)){
      
      message("Using Previous Background Peaks!")
      bgdPeaks <- readRDS(metadata(getPeakSet(ArchRProj))$bgdPeaks)

    }else{

      if(force){
      
        message("Previous Background Peaks file does not exist! Identifying Background Peaks!")
        bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, nIterations=nIterations, w=w, binSize=binSize, seed = seed, outFile = NULL, method = method)
      
      }else{
      
        stop("Previous Background Peaks file does not exist! set add = TRUE to addBgdPeaks!")
      
      }

    }

  }else{
    
    message("Identifying Background Peaks!")
    bgdPeaks <- .computeBgdPeaks(ArchRProj=ArchRProj, nIterations=nIterations, w=w, binSize=binSize, seed = seed, outFile = NULL, method = method)

  }

  if(length(getPeakSet(ArchRProj)) != nrow(bgdPeaks)){
    stop("Number of rows in Background Peaks does not match peakSet!")
  }

  bgdPeaks

}

.computeBgdPeaks <- function(
  ArchRProj = NULL, 
  nIterations = 50,
  w = 0.1, 
  binSize = 50, 
  seed = 1, 
  method = "chromVAR",
  outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds")
  ){

  set.seed(1)
  .requirePackage("chromVAR",source="bioc")

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

  all1 <- all(
    paste0(rS$seqnames, ":", rS$idx) %in% 
    paste0(seqnames(ArchRProj@peakSet), ":", ArchRProj@peakSet$idx)
  )

  all2 <- all(
  paste0(seqnames(ArchRProj@peakSet), ":", ArchRProj@peakSet$idx) %in% 
    paste0(rS$seqnames, ":", rS$idx)
  )

  if(!(all1 & all2)){
    stop("PeakSet in Arrows does not match PeakSet in ArchRProject!
     To try to solve this, try re-running addPeakMatrix(ArchRProj, force=TRUE)")    
  }

  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC

  uniqueDist <- unique(rS$end - rS$start)
  if(length(uniqueDist) > 1){
    if(tolower(method) != "archr"){
      stop("When using non-fixed width peaks, you need to use method = ArchR!")
    }
  }

  #minimal chromVAR change
  #chromVAR reuiqres a matrix/se of ncol > 1 and with a log10(values) transform removing peaks with 0 reads
  #to disable this we create a column of 1's forcing chromVAR to perform log10(values + 1)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = SimpleList(counts = as.matrix(data.frame(rS$rowSums, 1))),
    rowData = DataFrame(bias = rS$GC, start = rS$start, end = rS$end)
  )

  if(tolower(method) == "chromvar"){

    bgdPeaks <- tryCatch({

      chromVAR::getBackgroundPeaks(
        object = se,
        bias = rowData(se)$bias, 
        niterations = nIterations, 
        w = w, 
        bs = binSize
      )

    }, error = function(e){

      message("Error with chromVAR::getBackgroundPeaks! Handling this with a safer method for getting background peaks with ArchR!")

      .ArchRBdgPeaks(
        object = se,
        bias = rowData(se)$bias, 
        nIterations = nIterations
      )

    })

  }else{

      bgdPeaks <- .ArchRBdgPeaks(
        object = se,
        bias = rowData(se)$bias, 
        nIterations = nIterations
      )

  }

  bgdPeaks <- SummarizedExperiment(assays = SimpleList(bgdPeaks = bgdPeaks), 
    rowRanges = GRanges(rS$seqnames,IRanges(rS$start,rS$end),value=rS$rowSums,GC=rS$GC))

  biasDF <- data.frame(
    rowSums = Matrix::rowSums(assay(se)),
    bias = rowData(se)$bias,
    length = rowData(se)$end - rowData(se)$start
  )

  rowData(bgdPeaks)$bgdSumMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),1], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdSumSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),1], nrow = nrow(bgdPeaks))),3)

  rowData(bgdPeaks)$bgdGCMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),2], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdGCSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),2], nrow = nrow(bgdPeaks))),3)

  rowData(bgdPeaks)$bgdLengthMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),3], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdLengthSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),3], nrow = nrow(bgdPeaks))),3)

  #Save Background Peaks
  if(!is.null(outFile)){
    saveRDS(bgdPeaks, outFile, compress = FALSE)
  }

  return(bgdPeaks)

}

.ArchRBdgPeaks <- function(object = NULL, bias = NULL, nIterations = 50){

  .cleanSelf <- function(x){
    xn <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for(i in seq_len(nrow(x))){
        xi <- x[i, ]
        idx <- which(xi != i)
        xn[i, seq_along(idx)] <- xi[idx]
    }
    idx <- which(colSums(xn == 0) > 0)
    if(length(idx) > 0){
      xn <- xn[,-idx]
    }
    xn
  }

  #Bias Dataframe
  biasDF <- data.frame(
    rowSums = Matrix::rowSums(assay(object)),
    bias = bias,
    length = rowData(object)$end - rowData(object)$start
  )

  #Quantile Normalize
  biasDFN <- apply(biasDF, 2, .getQuantiles)

  #Get KNN
  knnObj <- nabor::knn(
    data =  biasDFN,
    k = nIterations + 1
  )[[1]]

  #Filter Self
  knnObj <- .cleanSelf(knnObj)

  #Shuffle
  idx <- seq_len(ncol(knnObj))
  knnObj2 <- matrix(0, nrow = nrow(knnObj), ncol = ncol(knnObj))
  for(x in seq_len(nrow(knnObj2))){
    knnObj2[x,] <- knnObj[x, sample(idx, length(idx))]
  }

  knnObj2

}






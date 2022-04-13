#' Project Bulk ATAC-seq data into single cell subspace
#' 
#' This function will Project Bulk ATAC-seq data into single cell subspace.
#' 
#' @param ArchRProj An `ArchRProject` object containing the dimensionality reduction matrix passed by `reducedDims`.
#' @param seATAC A `SummarizedExperiment` object containing bulk ATAC-seq data.
#' @param reducedDims A string specifying the name of the `reducedDims` object to be used.
#' @param embedding A string specifying the name of the `embedding` object to be used.
#' @param n An integer specifying the number of subsampled "pseudo single cells" per bulk sample.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param threads The number of threads used for parallel execution
#' @param force A boolean value indicating whether to force the projection of bulk ATAC data even if fewer than 25% of the features are present in the bulk ATAC data set.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
#'
projectBulkATAC <- function(
  ArchRProj = NULL,
  seATAC = NULL,
  reducedDims = "IterativeLSI",
  embedding = "UMAP",
  n = 250,
  verbose = TRUE,
  threads = getArchRThreads(),
  force = FALSE,
  logFile = createLogFile("projectBulkATAC")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = seATAC, name = "seATAC", valid = c("SummarizedExperiment"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character", "null"))
  .validInput(input = n, name = "n", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "projectBulkATAC Input-Parameters", logFile = logFile)

  ##################################################
  # Reduced Dimensions
  ##################################################
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, returnMatrix = FALSE)
  .logThis(names(rD), "reducedDimsNames", logFile = logFile)
  .logThis(rD[[1]], "reducedDimsMat", logFile = logFile)
  rDFeatures <- rD[[grep("Features", names(rD))]]
  if("end" %in% colnames(rDFeatures)){
    rDGR <- GRanges(seqnames=rDFeatures$seqnames,IRanges(start=rDFeatures$start, end=rDFeatures$end))
  }else{
    rDGR <- GRanges(seqnames=rDFeatures$seqnames,IRanges(start=rDFeatures$start, width = (rDFeatures$start) / (rDFeatures$idx - 1)))
  }
  .logThis(rDGR, "reducedDimsGRanges", logFile = logFile)
  subATAC <- subsetByOverlaps(seATAC, rDGR, ignore.strand = TRUE)
  subATAC <- subATAC[order(rowSums(as.matrix(.getAssay(subATAC, "counts"))), decreasing = TRUE), ]
  o <- DataFrame(findOverlaps(subATAC, rDGR, ignore.strand = TRUE))
  sumOverlap <- length(unique(o[,2]))
  .logThis(o, "overlapATAC", logFile = logFile)

  if(sumOverlap == 0){
    .logStop(paste0("No overlaps between bulk ATAC data and reduce dimensions feature found.",
      "\nEither recreate counts matrix or most likely these data sets are incompatible!"), logFile = logFile)
  }
  if( (sumOverlap / length(rDGR)) < 0.25 ){
    if(force){
      .logMessage("Less than 25% of the features are present in this bulk ATAC data set! Continuing since force = TRUE!", verbose = TRUE, logFile = logFile)
    }else{
      .logStop("Less than 25% of the features are present in this bulk ATAC data set! Set force = TRUE to continue!", logFile = logFile)
    }
  }
  .logMessage("Overlap Ratio of Reduced Dims Features = ", (sumOverlap / length(rDGR)), verbose = TRUE, logFile = logFile)

  o <- o[!duplicated(o$subjectHits),]
  subATAC <- subATAC[o$queryHits, ]
  rownames(subATAC) <- paste0("f", o$subjectHits)
  .logThis(subATAC, "subsettedATAC", logFile = logFile)

  ##################################################
  # Create Bulk Matrix
  ##################################################
  bulkMat <- .safeSubset(
    mat = .getAssay(subATAC, "counts"), 
    subsetRows = paste0("f", seq_along(rDGR))
  )
  .logThis(bulkMat, "bulkATACMat", logFile = logFile)

  ##################################################
  # Simulate and Project
  ##################################################
  depthN <- round(sum(rD$rowSm / rD$nCol))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments

  simRD <- .safelapply(seq_len(ncol(bulkMat)), function(x){
    .logDiffTime(sprintf("Projecting Sample (%s of %s)", x, ncol(bulkMat)), t1 = tstart, verbose = verbose, logFile = logFile)
    counts <- bulkMat[, x]
    counts <- rep(seq_along(counts), counts)
    simMat <- lapply(seq_len(nRep), function(y){
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[,-1,drop=FALSE]
      simMat[,1] <- simMat[,1] + (y - 1) * n2
      simMat
    }) %>%  Reduce("rbind", .)
    simMat <- Matrix::sparseMatrix(i = simMat[,2], j = simMat[,1], x = rep(1, nrow(simMat)), dims = c(length(rDGR), n2 * nRep))
    simRD <- as.matrix(.projectLSI(simMat, LSI = rD, verbose = FALSE))
    rownames(simRD) <- paste0(colnames(bulkMat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, threads = threads) %>% Reduce("rbind", .)

  if(is.null(embedding)){
    if(rD$scaleDims){
      simRD <- .scaleDims(simRD)
    }
    out <- SimpleList(
      simulatedReducedDims = simRD
    )
    return(out)
  }
  .logThis(simRD, "simulatedReducedDims", logFile = logFile)

  ##################################################
  # Prep Reduced Dims
  ##################################################
  embedding <- getEmbedding(ArchRProj = ArchRProj, embedding = embedding, returnDF = FALSE)
  corCutOff <- embedding$params$corCutOff
  dimsToUse <- embedding$params$dimsToUse
  scaleDims <- embedding$params$scaleDims

  if(is.null(scaleDims)){
    scaleDims <- rD$scaleDims
  }

  simRD <- .scaleDims(simRD)

  if(embedding$params$nc != ncol(simRD)){
    
    if(is.null(dimsToUse)){
      dimsToUse <- seq_len(ncol(rD[[1]]))
    }

    if(!is.null(corCutOff)){
      if(scaleDims){
        corToDepth <- rD$corToDepth$scaled
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      }else{
        corToDepth <- rD$corToDepth$none
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      }
    }

    if(embedding$params$nc != ncol(simRD)){
      .logMessage("Error! Inconsistency found with matching LSI dimensions to those used in addUMAP or addTSNE",
        "\nReturning with simulated reduced dimension coordinates...", verbose = TRUE, logFile = logFile)
      out <- SimpleList(
        simulatedReducedDims = simRD
      )
      return(out)
    }

    simRD <- simRD[, dimsToUse, drop = FALSE]

  }

  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umapModel <- .loadUWOT(embedding$params$uwotModel, embedding$params$nc)

  idx <- sort(sample(seq_len(nrow(rD[[1]])), min(nrow(rD[[1]]), 5000))) #Try to use 5000 or total cells to check validity
  rD2 <- getReducedDims(
    ArchRProj = ArchRProj, 
    reducedDims = reducedDims, 
    dimsToUse = embedding$params$dimsToUse,
    scaleDims = embedding$params$scaleDims,
    corCutOff = embedding$params$corCutOff
  )[idx,,drop=FALSE]

  ##################################################
  # Project UMAP
  ##################################################
  set.seed(1)
  threads2 <- max(floor(threads/2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD2, simRD), 
    model = umapModel, 
    verbose = TRUE, 
    n_threads = threads2
  )
  rownames(simUMAP) <- c(rownames(rD2), rownames(simRD))
  .logThis(simUMAP, "simulatedUMAP", logFile = logFile)

  #Check if the projection matches using previous data
  c1 <- cor(simUMAP[rownames(rD2), 1], embedding[[1]][rownames(rD2),1])
  c2 <- cor(simUMAP[rownames(rD2), 2], embedding[[1]][rownames(rD2),2])
  if(min(c1, c2) < 0.8){
    .logMessage("Warning projection correlation is less than 0.8 (R = ", round(min(c1,c2), 4),").\nThese results may not be accurate because of the lack of heterogeneity in the single cell data.", verbose = TRUE, logFile = logFile)
  }

  dfUMAP <- embedding[[1]]
  colnames(dfUMAP) <- c("UMAP1", "UMAP2")
  colnames(simUMAP) <- c("UMAP1", "UMAP2")
  dfUMAP <- DataFrame(dfUMAP)
  dfUMAP$Type <- Rle("scATAC", lengths = nrow(dfUMAP))

  simUMAP <- DataFrame(simUMAP[rownames(simRD),,drop=FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[,1])

  out <- SimpleList(
    simulatedBulkUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    simulatedReducedDims = simRD
  )
  .endLogging(logFile = logFile)

  return(out)

}



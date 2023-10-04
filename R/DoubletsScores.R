##########################################################################################
# Doublet Identification Methods
##########################################################################################

#' Add Doublet Scores to a collection of ArrowFiles or an ArchRProject
#' 
#' For each sample in the ArrowFiles or ArchRProject provided, this function will independently assign inferred doublet information
#' to each cell. This allows for removing strong heterotypic doublet-based clusters downstream. A doublet results from a droplet that
#' contained two cells, causing the ATAC-seq data to be a mixture of the signal from each cell. 
#'
#' @param input An `ArchRProject` object or a character vector containing the paths to the ArrowFiles to be used.
#' @param useMatrix The name of the matrix to be used for performing doublet identification analyses. Options include "TileMatrix" and "PeakMatrix".
#' @param k The number of cells neighboring a simulated doublet to be considered as putative doublets.
#' @param nTrials The number of times to simulate nCell (number of cells in the sample) doublets to use for doublet simulation when calculating doublet scores.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param LSIMethod A number or string indicating the order of operations in the TF-IDF normalization.
#' Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param scaleDims A boolean that indicates whether to z-score the reduced dimensions for each cell during the LSI
#' method performed for doublet determination. This is useful for minimizing the contribution of strong biases (dominating early PCs)
#' and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation
#' to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param knnMethod The name of the dimensionality reduction method to be used for k-nearest neighbors calculation. Possible values are "UMAP" or "LSI".
#' @param UMAPParams The list of parameters to pass to the UMAP function if "UMAP" is designated to `knnMethod`. See the function `umap` in the uwot package.
#' @param LSIParams The list of parameters to pass to the `IterativeLSI()` function. See `IterativeLSI()`.
#' @param outDir The relative path to the output directory for relevant plots/results from doublet identification.
#' @param threads The number of threads to be used for parallel computing.
#' @param force If the UMAP projection is not accurate (when R < 0.8 for the reprojection of the training data - this occurs when you 
#' have a very homogenous population of cells), setting `force=FALSE` will return -1 for all doubletScores and doubletEnrichments. If you would like to
#' override this (not recommended!), you can bypass this warning by setting `force=TRUE`.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param verbose A boolean value that determines whether standard output is printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addDoubletScores <- function(
  input = NULL,
  useMatrix = "TileMatrix",
  k = 10,
  nTrials = 5,
  dimsToUse = 1:30,
  LSIMethod = 1,
  scaleDims = FALSE,
  corCutOff = 0.75,
  knnMethod = "UMAP",
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "euclidean", verbose = FALSE),
  LSIParams = list(outlierQuantiles = NULL, filterBias = FALSE),
  outDir = getOutputDirectory(input),  
  threads = getArchRThreads(),
  force = FALSE,
  parallelParam = NULL,
  verbose = TRUE,
  logFile = createLogFile("addDoubletScores")
  ){

  .validInput(input = input, name = "input", valid = c("character", "ArchRProject"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = nTrials, name = "nTrials", valid = c("integer"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean"))
  .validInput(input = knnMethod, name = "knnMethod", valid = c("character"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  .validInput(input = LSIParams, name = "LSIParams", valid = c("list"))
  .validInput(input = outDir, name = "outDir", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))


  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addDoubletScores Input-Parameters", logFile = logFile)

  if(tolower(useMatrix) %ni% c("peakmatrix","tilematrix")){
    stop(sprintf("Supported Matrix Names at the moment are PeakMatrix and TileMatrix : ", useMatrix))
  }

  if(inherits(input, "ArchRProject")){
    
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
  
  }else if(inherits(input, "character")){

    ArrowFiles <- input
    allCells <- NULL

  }else{

    stop("Error Unrecognized Input!")

  }

  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addDoubScores
  args$registryDir <- file.path(outDir, "AddDoubletsRegistry")

  #Make Sure these Args are NULL
  args$input <- NULL

  #Run With Parallel or lapply
  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
  outList <- .batchlapply(args, sequential = TRUE)
  names(outList) <- names(ArrowFiles)

  .endLogging(logFile = logFile)

  #Return Output
  if(inherits(input, "ArchRProject")){

    input@cellColData[,"DoubletScore"] <- NA
    input@cellColData[,"DoubletEnrichment"] <- NA

    for(i in seq_along(outList)){
      input@cellColData[names(outList[[i]]$doubletScore), "DoubletScore"] <- outList[[i]]$doubletScore
      input@cellColData[names(outList[[i]]$doubletEnrich), "DoubletEnrichment"] <- outList[[i]]$doubletEnrich
    }

    return(input)

  }else{

    return(outList)

  }

}

.addDoubScores <- function(
  i = NULL,
  ArrowFiles = NULL,
  useMatrix = "TileMatrix",
  allCells = NULL,
  UMAPParams = list(),
  LSIParams = list(),
  nTrials = 5,
  dimsToUse = 1:30,
  corCutOff = 0.75,
  LSIMethod = 1,
  sampleCells = NULL,
  scaleDims = FALSE,
  k = 10,
  nSample = 1000,
  knnMethod = "UMAP",
  outDir = "QualityControl",
  force = FALSE,
  subThreads = 1,
  verbose = TRUE,
  tstart = NULL,
  logFile = NULL
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)
  outDir <- file.path(outDir, sampleName)
  dir.create(outDir, showWarnings = FALSE)
  prefix <- sprintf("%s (%s of %s) : ", sampleName, i, length(ArrowFiles))

  .logDiffTime(sprintf("%s Computing Doublet Statistics", prefix), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)

  #################################################
  # 1. Create ArchRProject For Iterative LSI
  #################################################
  tmpDir <- .tempfile()
  dir.create(tmpDir)
  proj <- suppressMessages(ArchRProject(
    ArrowFiles = ArrowFile,
    outputDirectory = tmpDir,
    copyArrows = FALSE,
    showLogo = FALSE,
    geneAnnotation = .nullGeneAnnotation(), #this doesnt matter just needs to be valid
    genomeAnnotation = .nullGenomeAnnotation() #this doesnt matter just needs to be valid
  ))
  if(is.null(allCells)){
    proj@cellColData <- proj@cellColData[.availableCells(ArrowFile, useMatrix),]
  }else{
    proj@cellColData <- proj@cellColData[which(rownames(proj@cellColData) %in% allCells),]
  }

  #################################################
  # 2. Compute Iterative LSI
  #################################################
  .logDiffTime("Running IterativeLSI", tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
  LSIParams$ArchRProj <- proj
  LSIParams$saveIterations <- FALSE
  LSIParams$useMatrix <- useMatrix
  LSIParams$LSIMethod <- LSIMethod
  LSIParams$dimsToUse <- dimsToUse
  LSIParams$scaleDims <- scaleDims
  LSIParams$corCutOff <- corCutOff
  LSIParams$threads <- subThreads
  LSIParams$verbose <- FALSE
  LSIParams$force  <- TRUE
  LSIParams$logFile <- logFile
  proj <- tryCatch({
    do.call(addIterativeLSI, LSIParams)
  }, error = function(e){
    .logError(e, fn = "addIterativeLSI", info = prefix, errorList = list(ArrowFile = ArrowFile), logFile = logFile)
  })

  #################################################
  # 3. Get LSI Partial Matrix For Simulation
  #################################################
  .logDiffTime("Constructing Partial Matrix for Projection", tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
  LSI <- getReducedDims(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    corCutOff = corCutOff, 
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    returnMatrix = FALSE
  )
  .logThis(LSI, name = paste0(prefix, "LSI Result"), logFile = logFile)

  LSIDims <- seq_len(ncol(LSI[[1]]))
  if(length(LSIDims) < 2){
    .logMessage("Reduced LSI Dims below 2 dimensions, please increase dimsToUse or increase corCutOff!")
    stop("Reduced LSI Dims below 2 dimensions, please increase dimsToUse or increase corCutOff!")
  }
  featureDF <- LSI$LSIFeatures

  mat <- tryCatch({
    .getPartialMatrix(
      ArrowFiles = getArrowFiles(proj),
      featureDF = featureDF,
      threads = subThreads,
      cellNames = rownames(getCellColData(proj)),
      doSampleCells = FALSE,
      verbose = FALSE
    )
  }, error = function(e){
    errorList <- list(
      ArrowFiles = getArrowFiles(proj),
      featureDF = featureDF,
      threads = subThreads,
      cellNames = rownames(getCellColData(proj)),
      doSampleCells = FALSE,
      verbose = FALSE
    )
    .logError(e, fn = "getPartialMatrix", info = prefix, errorList = errorList, logFile = logFile)
  })
  cellNames <- rownames(getCellColData(proj))

  #################################################
  # 4. Run UMAP for LSI-Projection
  #################################################
  .logDiffTime("Running LSI UMAP", tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
  set.seed(1) # Always do this prior to UMAP
  UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="euclidean", verbose=FALSE))
  UMAPParams$X <- LSI$matSVD
  UMAPParams$ret_nn <- TRUE
  UMAPParams$ret_model <- TRUE
  UMAPParams$n_threads <- subThreads
  .logThis(UMAPParams, name = paste0(prefix, "UMAP Params"), logFile = logFile)
  uwotUmap <- tryCatch({
    do.call(uwot::umap, UMAPParams)
  }, error = function(e){
    errorList <- UMAPParams
    .logError(e, fn = "uwot::umap", info = prefix, errorList = errorList, logFile = logFile)
  })

  #################################################
  # 5. Simulate and Project Doublets
  #################################################
  .logDiffTime("Simulating and Projecting Doublets", tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
  simDoubletsSave <- tryCatch({
    .simulateProjectDoublets(
      mat = mat, 
      LSI = LSI, 
      sampleRatio1 = c(1/2), 
      sampleRatio2 = c(1/2), 
      nTrials = nTrials * max( floor(nCells(proj) / nSample), 1 ), 
      nSample = nSample, 
      k = k, 
      uwotUmap = uwotUmap,
      seed = 1, 
      force = force,
      threads = subThreads,
      logFile = logFile,
      prefix = prefix
    )
  }, error = function(e){
    errorList <- list(
      mat = mat, 
      LSI = LSI, 
      sampleRatio1 = c(1/2), 
      sampleRatio2 = c(1/2), 
      nTrials = nTrials * max( floor(nCells(proj) / nSample), 1 ), 
      nSample = nSample, 
      k = k, 
      uwotUmap = uwotUmap,
      seed = 1, 
      force = force,
      threads = subThreads,
      logFile = logFile,
      prefix = prefix
    )
    .logError(e, fn = ".simulateProjectDoublets", info = prefix, errorList = errorList, logFile = logFile)    
  })
  if(tolower(knnMethod)=="lsi"){
    simDoublets <- SimpleList(
      doubletUMAP=simDoubletsSave$doubletUMAP,
      doubletScore=simDoubletsSave$doubletScoreLSI,
      doubletEnrich=simDoubletsSave$doubletEnrichLSI
    )
  }else{
    simDoublets <- SimpleList(
      doubletUMAP=simDoubletsSave$doubletUMAP,
      doubletScore=simDoubletsSave$doubletScoreUMAP,
      doubletEnrich=simDoubletsSave$doubletEnrichUMAP
    )    
  }
  .logThis(simDoublets, name = paste0(prefix, "SimulationResults"), logFile = logFile)

  #################################################
  # 6. Plot / Save Results
  #################################################

  pal <- c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #grey_magma

  #Create Plot DF
  df <- data.frame(row.names = rownames(LSI$matSVD), uwotUmap[[1]], type = "experiment")
  df[,"score"] <- 0
  df[,"enrichment"] <- 0
  df[names(simDoublets$doubletScore),"score"] <- simDoublets$doubletScore
  df[names(simDoublets$doubletScore),"enrichment"] <- simDoublets$doubletEnrich
  
  doubUMAP <- simDoublets$doubletUMAP
  dfDoub <- data.frame(
    row.names = paste0("doublet_", seq_len(nrow(doubUMAP))), 
    .getDensity(doubUMAP[,1], doubUMAP[,2]), 
    type = "simulated_doublet"
  )
  dfDoub <- dfDoub[order(dfDoub$density), , drop = FALSE]
  dfDoub$color <- dfDoub$density
 
  ##################################
  #Save Results
  .logThis(df, name = paste0(prefix, "Sample UMAP"), logFile = logFile)
  .logThis(dfDoub, name = paste0(prefix, "Simulated Doublet UMAP"), logFile = logFile)
  summaryList <- SimpleList(
    originalDataUMAP = df,
    simulatedDoubletUMAP = dfDoub,
    doubletResults = simDoubletsSave
  )

  .safeSaveRDS(summaryList, file.path(outDir, paste0(.sampleName(ArrowFile), "-Doublet-Summary.rds")))
  rm(simDoubletsSave)
  ##################################

  tmpFile <- .tempfile()

  o <- tryCatch({

    #Plot Doublet Summary
    pdf(file.path(outDir, paste0(.sampleName(ArrowFile), "-Doublet-Summary.pdf")), width = 6, height = 6)

    #Plot Doublet Density
    xlim <- range(df$X1) %>% extendrange(f = 0.05)
    ylim <- range(df$X2) %>% extendrange(f = 0.05)
    
    pdensity <- ggplot() + 
      .geom_point_rast2(data = df, aes(x=X1,y=X2),color="lightgrey", size = 0.5) + 
      .geom_point_rast2(data = dfDoub, aes(x=x,y=y,colour=color), size = 0.5) + 
        scale_colour_gradientn(colors = pal) + 
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
        labs(color = "Simulated Doublet Density") +
        guides(fill = "none") + theme_ArchR(baseSize = 10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim, expand = FALSE) +
        ggtitle("Simulated and LSI-Projected Density Overlayed") + theme(legend.direction = "horizontal", 
        legend.box.background = element_rect(color = NA))
    
    # if(!requireNamespace("ggrastr", quietly = TRUE)){
      
    #   message("ggrastr is not available for rastr of points, continuing without rastr!")
    #   message("To install ggrastr try : devtools::install_github('VPetukhov/ggrastr')")

    #   pdensity <- ggplot() + 
    #     geom_point(data = df, aes(x=X1,y=X2),color="lightgrey", size = 0.5) + 
    #     geom_point(data = dfDoub, aes(x=x,y=y,colour=color), size = 0.5) + 
    #       scale_colour_gradientn(colors = pal) + 
    #       xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
    #       guides(fill = "none") + theme_ArchR(baseSize = 10) +
    #       labs(color = "Simulated Doublet Density") +
    #       theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    #             axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    #       coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim, expand = FALSE) +
    #       ggtitle("Simulated and LSI-Projected Doublet Density Overlayed") + theme(legend.direction = "horizontal", 
    #       legend.box.background = element_rect(color = NA))

    # }else{

    #   #.requirePackage("ggrastr", installInfo = "devtools::install_github('VPetukhov/ggrastr')")

    #   pdensity <- ggplot() + 
    #     .geom_point_rast2(data = df, aes(x=X1,y=X2),color="lightgrey", size = 0.5) + 
    #     .geom_point_rast2(data = dfDoub, aes(x=x,y=y,colour=color), size = 0.5) + 
    #       scale_colour_gradientn(colors = pal) + 
    #       xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
    #       labs(color = "Simulated Doublet Density") +
    #       guides(fill = "none") + theme_ArchR(baseSize = 10) +
    #       theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    #             axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    #       coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim, expand = FALSE) +
    #       ggtitle("Simulated and LSI-Projected Density Overlayed") + theme(legend.direction = "horizontal", 
    #       legend.box.background = element_rect(color = NA))

    # }

    #Plot Doublet Score
    pscore <- ggPoint(
      x = df[,1],
      y = df[,2],
      color = .quantileCut(df$score, 0, 0.95),
      xlim = xlim,
      ylim = ylim,
      discrete = FALSE,
      size = 0.5,
      xlab = "UMAP Dimension 1",
      ylab = "UMAP Dimension 2",
      pal = pal,
      title = "Doublet Scores -log10(P-adj.)",
      colorTitle = "Doublet Scores -log10(P-adj.)",
      rastr = TRUE,
      baseSize = 10
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank())

    #Plot Enrichment Summary
    penrich <- ggPoint(
      x = df[,1],
      y = df[,2],
      color = .quantileCut(df$enrichment, 0, 0.95),
      xlim = xlim,
      ylim = ylim,
      discrete = FALSE,
      size = 0.5,
      xlab = "UMAP Dimension 1",
      ylab = "UMAP Dimension 2",
      pal = pal,
      title = "Simulated Doublet Enrichment over Expectation",
      colorTitle = "Doublet Enrichment",
      rastr = TRUE,
      baseSize = 10
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank())


    #1. Doublet Enrichment
    .fixPlotSize(penrich, plotWidth = 6, plotHeight = 6)
    grid::grid.newpage()

    #2. Doublet Scores
    .fixPlotSize(pscore, plotWidth = 6, plotHeight = 6)
    grid::grid.newpage()

    #3. Doublet Density
    .fixPlotSize(pdensity, plotWidth = 6, plotHeight = 6)

    dev.off()

  }, error = function(e){
    errorList <- list(df=df,dfDoub=dfDoub)
    .logError(e, fn = "ggplot", info = prefix, errorList = errorList, logFile = logFile, throwError = FALSE)
  })

  #################################################
  # 7. Add Info To Arrow!
  #################################################
  allCells <- .availableCells(ArrowFile, passQC = FALSE)
  
  allDoubletScores <- rep(-1, length(allCells))
  names(allDoubletScores) <- allCells
  allDoubletScores[names(simDoublets$doubletScore)] <- simDoublets$doubletScore
  
  allDoubletEnrichment <- rep(-1, length(allCells))
  names(allDoubletEnrichment) <- allCells
  allDoubletEnrichment[names(simDoublets$doubletEnrich)] <- simDoublets$doubletEnrich

  o <- h5closeAll()
  h5write(allDoubletScores, file = ArrowFile, "Metadata/DoubletScore")
  h5write(allDoubletEnrichment, file = ArrowFile, "Metadata/DoubletEnrichment")
  o <- h5closeAll()

  out <- SimpleList(doubletScore = simDoublets$doubletScore, doubletEnrich = simDoublets$doubletEnrich)

  return(out)

}

.simulateProjectDoublets <- function(
  mat = NULL, 
  LSI = NULL,
  uwotUmap = NULL,
  sampleRatio1 = c(0.5), 
  sampleRatio2 = c(0.5), 
  nTrials = 100, 
  nSample = 1000, 
  k = 200, 
  knnMethod = "UMAP",
  seed = 1, 
  threads = 1,
  force = FALSE,
  logFile = NULL,
  prefix = ""  
  ){

  .sampleSparseMat <- function(mat = NULL, sampleRatio = 0.5){
    total <- length(mat@x)
    sampleTo <- floor(total * (1-sampleRatio))
    mat@x[sample(seq_len(total), sampleTo)] <- 0
    mat <- drop0(mat)
    mat
  }

  set.seed(seed)
  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
  simLSI <- tryCatch({
    .safelapply(seq_len(nTrials), function(y){
      if(y %% 5 == 0){
        gc()
      }
      lapply(seq_along(sampleRatio1), function(x){
        idx1 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
        idx2 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
        #Simulated Doublet
        simulatedMat <- .sampleSparseMat(mat = mat[,idx1], sampleRatio = sampleRatio1[x]) + 
                        .sampleSparseMat(mat = mat[,idx2], sampleRatio = sampleRatio2[x])
        #Project LSI
        lsiProject <- suppressMessages(.projectLSI(simulatedMat, LSI))
        rownames(lsiProject) <- NULL
        lsiProject
      }) %>% Reduce("rbind", .)
    }, threads = threads) %>% Reduce("rbind", .)

  }, error = function(e){
    .logError(e, fn = "Simulate LSI Project Doublets", info = prefix, errorList = errorList, logFile = logFile)
  })

  .logThis(simLSI, name = paste0(prefix, "SimulatedLSI"), logFile = logFile)

  #Compute original
  ogLSI <- suppressMessages(.projectLSI(mat, LSI))

  .logThis(ogLSI, name = paste0(prefix, "OriginalLSI"), logFile = logFile)
  .logThis(LSI[[1]], name = paste0(prefix, "OriginalReProjectedLSI"), logFile = logFile)

  #Merge
  allLSI <- rbind(simLSI[, LSI$dimsKept, drop = FALSE], ogLSI[, LSI$dimsKept, drop = FALSE])
  nSimLSI <- nrow(simLSI)
  if(nSimLSI==0){
    stop("Simulations must be greater than 0! Please adjust nTrials!")
  }
  rm(simLSI)
  rm(ogLSI)
  gc()

  if(LSI$scaleDims){
    allLSI <- .scaleDims(allLSI)
  }

  #Project UMAP
  set.seed(1) # Always do this prior to UMAP
  umapProject <- tryCatch({
    uwot::umap_transform(
      X = as.matrix(allLSI), 
      model = uwotUmap, 
      verbose = FALSE, 
      n_threads = threads
    )
  }, error = function(e){
    errorList <- list(X = allLSI, model = uwotUmap)
    .logError(e, fn = "uwot::umap_transform", info = prefix, errorList = errorList, logFile = logFile)    
  })

  corProjection <- list(
    LSI = unlist(lapply(seq_len(ncol(allLSI)), function(x) cor(allLSI[-seq_len(nSimLSI), x], LSI$matSVD[, x]) )),
    UMAP =  c(
        dim1 = cor(uwotUmap[[1]][,1], umapProject[-seq_len(nSimLSI), 1]),
        dim2 = cor(uwotUmap[[1]][,2], umapProject[-seq_len(nSimLSI), 2])
      )
  )
  names(corProjection[[1]]) <- paste0("SVD", LSI$dimsKept)

  .logThis(uwotUmap[[1]], name = paste0(prefix, "OriginalUMAP"), logFile = logFile)
  .logThis(umapProject[-seq_len(nSimLSI), ], name = paste0(prefix, "OriginalReProjectedUMAP"), logFile = logFile)

  msg <- paste0(prefix, "UMAP Projection R^2 = ", round(mean(corProjection[[2]])^2, 5))
  .logMessage(msg, logFile = logFile)
  message(msg)

  out <- SimpleList(
    doubletUMAP = umapProject[seq_len(nSimLSI), ],
    projectionCorrelation = corProjection
  )

  if(mean(corProjection[[2]]) < 0.9){
    if(!force){
      msg <- paste0(prefix, "Correlation of UMAP Projection is below 0.9 (normally this is ~0.99)\nThis means there is little heterogeneity in your sample and thus doubletCalling is inaccurate.\nforce = FALSE, thus returning -1 doubletScores and doubletEnrichments!\nSet force = TRUE if you want to continue (not recommended).")
      .logMessage(msg, logFile = logFile)
      message(msg)
      out$doubletEnrichLSI <- rep(-1, nrow(LSI$matSVD))
      out$doubletScoreLSI <- rep(-1, nrow(LSI$matSVD))
      out$doubletEnrichUMAP <- rep(-1, nrow(LSI$matSVD))
      out$doubletScoreUMAP <- rep(-1, nrow(LSI$matSVD))
      return(out)
    }
  }

  ##############################################################################
  # Compute Doublet Scores from LSI (TF-IDF + SVD)
  ##############################################################################
  out2 <- tryCatch({
    
    knnDoub <- .computeKNN(allLSI[-seq_len(nSimLSI),], allLSI[seq_len(nSimLSI),], k)
    
    #Compile KNN Sums
    countKnn <- rep(0, nrow(LSI$matSVD))
    names(countKnn) <- rownames(LSI$matSVD)

    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <-  countKnn[as.integer(names(tabDoub))] + tabDoub
    .logThis(countKnn, paste0(prefix, "Number of Nearby Simulated Doublets"), logFile = logFile)

    nSim <- nrow(LSI$matSVD)
    scaleTo <- 10000
    scaleBy <- scaleTo / nSim

    #P-Values
    pvalBinomDoub <- lapply(seq_along(countKnn), function(x){
      #Round Prediction
      countKnnx <- round(countKnn[x] * scaleBy)
      sumKnnx <- round(sum(countKnn) * scaleBy)
      pbinom(countKnnx - 1, sumKnnx, 1 / scaleTo, lower.tail = FALSE)
    }) %>% unlist

    #Adjust
    padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")

    #Convert To Scores
    doubletScore <- -log10(pmax(padjBinomDoub, 4.940656e-324))
    doubletEnrich <- (countKnn / sum(countKnn)) / (1 / nrow(LSI$matSVD))
    doubletEnrich <- 10000 * doubletEnrich / length(countKnn) #Enrichment Per 10000 Cells in Data Set  
    .logThis(doubletScore, paste0(prefix, "DoubletScoresLSI"), logFile = logFile)
    .logThis(doubletEnrich, paste0(prefix, "DoubletEnrichLSI"), logFile = logFile)

    #Store Results
    out$doubletEnrichLSI <- doubletEnrich
    out$doubletScoreLSI <- doubletScore

    ##############################################################################
    # Compute Doublet Scores from LSI (TF-IDF + SVD) + UMAP Embedding
    ##############################################################################
    knnDoub <- .computeKNN(umapProject[-seq_len(nSimLSI),], umapProject[seq_len(nSimLSI),], k)
    
    #Compile KNN Sums
    countKnn <- rep(0, nrow(LSI$matSVD))
    names(countKnn) <- rownames(LSI$matSVD)

    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <-  countKnn[as.integer(names(tabDoub))] + tabDoub

    nSim <- nrow(LSI$matSVD)
    scaleTo <- 10000
    scaleBy <- scaleTo / nSim

    #P-Values
    pvalBinomDoub <- lapply(seq_along(countKnn), function(x){
      #Round Prediction
      countKnnx <- round(countKnn[x] * scaleBy)
      sumKnnx <- round(sum(countKnn) * scaleBy)
      pbinom(countKnnx - 1, sumKnnx, 1 / scaleTo, lower.tail = FALSE)
    }) %>% unlist

    #Adjust
    padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")

    #Convert To Scores
    doubletScore <- -log10(pmax(padjBinomDoub, 4.940656e-324))
    doubletEnrich <- (countKnn / sum(countKnn)) / (1 / nrow(LSI$matSVD))
    doubletEnrich <- 10000 * doubletEnrich / length(countKnn) #Enrichment Per 10000 Cells in Data Set  
    .logThis(doubletScore, paste0(prefix, "DoubletScoresUMAP"), logFile = logFile)
    .logThis(doubletEnrich, paste0(prefix, "DoubletEnrichUMAP"), logFile = logFile)

    #Store Results
    out$doubletEnrichUMAP <- doubletEnrich
    out$doubletScoreUMAP <- doubletScore

    out

  }, error = function(e){
    errorList <- list(
      prefix = prefix,
      knnDoub = if(exists("knnDoub", inherits = FALSE)) knnDoub else "knnDoublets",
      pvalBinomDoub = if(exists("pvalBinomDoub", inherits = FALSE)) pvalBinomDoub else "Error with Pval!",
      padjBinomDoub = if(exists("padjBinomDoub", inherits = FALSE)) padjBinomDoub else "Error with Padj!",
      out = if(exists("out", inherits = FALSE)) out else "Error with outlist of results!"
    )
    .logError(e, fn = "Simulate LSI Project Doublets", info = prefix, errorList = errorList, logFile = logFile)
  })

  #Save Output
  out2

}

#' Add Demuxlet Results to ArchR Project
#' 
#' This function will read in the .best file output from demuxlet and add the doublet
#' classifications into the cellColData for the ArchR Project
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param bestFiles The file path to the .best files created by Demuxlet. There should be one .best file for each sample in the `ArchRProject`.
#' @param sampleNames The sample names corresponding to the .best files. These must match the sample names present in the `ArchRProject`.
#' @export
addDemuxletResults <- function(ArchRProj = NULL, bestFiles = NULL, sampleNames = NULL){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = bestFiles, name = "bestFiles", valid = c("character"))
  .validInput(input = sampleNames, name = "sampleNames", valid = c("character"))

  .requirePackage("readr", source = "cran")

  if(!all(sampleNames %in% rownames(getSampleColData(ArchRProj)))){
    samples <- sampleNames[sampleNames %ni% rownames(getSampleColData(ArchRProj))]
    warning(sprintf("Sample %s not in sampleNames of ArchRProj!", samples))
  }

  ccd <- getCellColData(ArchRProj)
  ccd[ , "DemuxletClassify"] <- "NotClassified"
  ccd[ , "DemuxletBest"] <- "NotClassified"

  for(x in seq_along(bestFiles)){
    best <- .suppressAll(data.frame(readr::read_tsv(bestFiles[x])))
    classification <- stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[,1]
    cellNames <- paste0(sampleNames[x], "#", best$BARCODE)
    idx <- which(cellNames %in% rownames(ccd))
    ccd[ cellNames[idx], "DemuxletClassify"] <- ifelse(
      stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,1] %in% c("AMB","DBL"),
      stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,1],
      stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,2]
    )
    ccd[ cellNames[idx], "DemuxletBest"] <- gsub("AMB-","",gsub("DBL-","",gsub("SNG-","",best$BEST)))[idx]
  }

  ArchRProj@cellColData <- ccd
  ArchRProj
  
}




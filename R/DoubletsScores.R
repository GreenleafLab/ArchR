#' Add Doublet Scores to a collection of Arrow files or an ArchRProject
#' 
#' For each sample in the Arrow files or ArchRProject provided, this function will independently assign inferred doublet information
#' to each cell. This allows for removing strong heterotypic doublet-based clusters downstream. A doublet results from a droplet that
#' contained two cells, causing the ATAC-seq data to be a mixture of the signal from each cell. 
#'
#' @param input An `ArchRProject` object or a character vector containing the names of ArrowFiles to be used.
#' @param useMatrix QQQ The name of the matrix to be used for performing doublet identification analyses. Options include "TileMatrix", QQQ.
#' @param k The number of cells neighboring a simulated doublet to be considered as putative doublets.
#' @param nTrials QQQ The number of trials (in terms of the number of input cells) to simulate doublets when calculating doublet scores. A value of 5 would utilize 5 trials.
#' @param dimsToUse QQQ A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param corCutOff QQQ A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is QQQ greater than the corCutOff, it will be excluded from analysis.
#' @param knnMethod The name of the dimensionality reduction method to be used for k-nearest neighbors calculation. Possible values are "UMAP" or "SVD".
#' @param UMAPParams The list of parameters to pass to the UMAP function if "UMAP" is designated to `knnMethod`. See the function umap in the uwot package.
#' @param LSIParams QQQ The list of parameters to pass to the IterativeLSI function if QQQ. See IterativeLSI.
#' @param outDir The name or path for the output directory for writing information on doublet identification,
#' @param threads The number threads to be used for parallel computing.
#' @param parallelParam QQQ A list of parameters to be used for batch-style parallel computing.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... additional args
#' @export
addDoubletScores <- function(
  input,
  useMatrix = "TileMatrix",
  k = 10,
  nTrials = 5,
  dimsToUse = 1:25,
  corCutOff = 0.75,
  knnMethod = "UMAP",
  UMAPParams = list(),
  LSIParams = list(sampleCells = NULL),
  outDir = "QualityControl",  
  threads = 1,
  parallelParam = NULL,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  if(tolower(useMatrix) %ni% c("peakmatrix","tilematrix")){
    stop(sprintf("Supported Matrix Names at the moment are PeakMatrix and TileMatrix : ", useMatrix))
  }

  if(inherits(input, "ArchRProject")){
    
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
  
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

  #Run With Parallel or lapply
  outList <- .batchlapply(args, sequential = TRUE)
  names(outList) <- names(ArrowFiles)

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
  i,
  ArrowFiles,
  useMatrix = "TileMatrix",
  allCells = NULL,
  UMAPParams = list(),
  LSIParams = list(),
  nTrials = 5,
  dimsToUse = 1:25,
  corCutOff = 0.75,
  k = 10,
  nSample = 1000,
  knnMethod = "UMAP",
  outDir = "QualityControl",
  #useClusters = FALSE,
  subThreads = 1,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  tstart <- Sys.time()
  useClusters <- FALSE
  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)
  outDir <- file.path(outDir, sampleName)
  dir.create(outDir, showWarnings = FALSE)

  .messageDiffTime(sprintf("Computing Doublet Scores %s (%s of %s)!", sampleName, i, length(ArrowFiles)), tstart, addHeader = TRUE)

  #################################################
  # 1. Create ArchRProject For Iterative LSI
  #################################################
  tmpDir <- .tempfile()
  dir.create(tmpDir)
  proj <- suppressMessages(ArchRProject(
    ArrowFiles = ArrowFile,
    sampleNames = .sampleName(ArrowFile),
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
  LSIParams$ArchRProj <- proj
  LSIParams$saveIterations <- FALSE
  LSIParams$useMatrix <- useMatrix
  LSIParams$dimsToUse <- dimsToUse
  LSIParams$corCutOff <- corCutOff
  LSIParams$threads <- subThreads
  LSIParams$verboseHeader <- verboseHeader
  LSIParams$verboseAll <- verboseAll
  proj <- do.call(addIterativeLSI, LSIParams)
  if(useClusters){
    proj <- addClusters(proj)
  }

  #################################################
  # 3. Get LSI Partial Matrix For Simulation
  #################################################
  .messageDiffTime("Constructing Partial Matrix for Projection", tstart, addHeader = verboseHeader)
  LSI <- getReducedDims(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    corCutOff = 999, 
    dimsToUse = dimsToUse,
    returnMatrix = FALSE
  )
  LSIDims <- intersect(seq_len(ncol(LSI[[1]])), which(LSI$corToDepth < corCutOff))
  if(length(LSIDims) < 2){
    stop("Reduced LSI Dims below 2 dimensions, please increase dimsToUse or increase corCutOff!")
  }
  featureDF <- LSI$LSIFeatures
  mat <- .getPartialMatrix(
      ArrowFiles = getArrowFiles(proj),
      featureDF = featureDF,
      threads = subThreads,
      cellNames = rownames(getCellColData(proj)),
      doSampleCells = FALSE,
      verbose = verboseAll
    )
  cellNames <- rownames(getCellColData(proj))

  #################################################
  # 2. Run UMAP for LSI-Projection
  #################################################
  .messageDiffTime("Running LSI UMAP", tstart, addHeader = verboseHeader)
  set.seed(1) # Always do this prior to UMAP
  UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors=40, min_dist=0.4, metric="euclidean", verbose=FALSE))
  UMAPParams$X <- LSI$matSVD[, LSIDims, drop = FALSE]
  UMAPParams$ret_nn <- TRUE
  UMAPParams$ret_model <- TRUE
  UMAPParams$n_threads <- subThreads
  uwotUmap <- do.call(uwot::umap, UMAPParams)

  #################################################
  # 4. Simulate and Project Doublets
  #################################################
  .messageDiffTime("Simulating and Projecting Doublets", tstart, addHeader = verboseHeader)
  simDoublets <- .simulateProjectDoublets(
    mat = mat, 
    LSI = LSI, 
    LSIDims = LSIDims,
    clusters = if(useClusters) getCellColData(proj, "Clusters", drop = TRUE) else NULL,
    sampleRatio1 = c(1/2), 
    sampleRatio2 = c(1/2), 
    nTrials = floor(nTrials * nCells(proj) / 1000), 
    nSample = nSample, 
    k = k, 
    uwotUmap = uwotUmap,
    knnMethod = knnMethod,
    seed = 1, 
    threads = threads
  )

  #################################################
  # 5. Plot / Save Results
  #################################################

  pal <- c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #grey magma

  #Create Plot DF
  df <- data.frame(row.names = rownames(LSI$matSVD), uwotUmap[[1]], type = "experiment")
  df[,"score"] <- 0
  df[,"enrichment"] <- 0
  df[names(simDoublets$doubletScore),"score"] <- simDoublets$doubletScore
  df[names(simDoublets$doubletScore),"enrichment"] <- simDoublets$doubletEnrich
  
  doubUMAP <- simDoublets$doubletUMAP
  dfDoub <- data.frame(
    row.names = paste0("doublet_", seq_len(nrow(doubUMAP))), 
    ArchR:::.getDensity(doubUMAP[,1], doubUMAP[,2]), 
    type = "simulated_doublet"
  )
  dfDoub <- dfDoub[order(dfDoub$density), , drop = FALSE]
  dfDoub$color <- dfDoub$density
  
  tmpFile <- .tempfile()
  sink(tmpFile)

  #Plot Doublet Summary
  pdf(file.path(outDir, paste0(.sampleName(ArrowFile), "-Doublet-Summary.pdf")), width = 6, height = 6)

  #Plot Doublet Density
  xlim <- range(df$X1) %>% extendrange(f = 0.05)
  ylim <- range(df$X2) %>% extendrange(f = 0.05)
  
  if(!requireNamespace("ggrastr", quietly = TRUE)){
    
    message("ggrastr is not available for rastr of points, continuing without rastr!")

    pdensity <- ggplot() + 
      geom_point(data = df, aes(x=X1,y=X2),color="lightgrey", size = 0.5) + 
      geom_point(data = dfDoub, aes(x=x,y=y,colour=color), size = 0.5) + 
        scale_colour_gradientn(colors = pal) + 
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
        guides(fill = FALSE) + theme_ArchR(baseSize = 6) +
        labs(color = "Simulated Doublet Density") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim, expand = FALSE) +
        ggtitle("Doublet Density Overlayed") + theme(legend.direction = "horizontal", 
        legend.box.background = element_rect(color = NA))

  }else{

    .requirePackage("ggrastr")

    pdensity <- ggplot() + 
      geom_point_rast(data = df, aes(x=X1,y=X2),color="lightgrey", size = 0.5) + 
      geom_point_rast(data = dfDoub, aes(x=x,y=y,colour=color), size = 0.5) + 
        scale_colour_gradientn(colors = pal) + 
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
        labs(color = "Simulated Doublet Density") +
        guides(fill = FALSE) + theme_ArchR(baseSize = 6) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim, expand = FALSE) +
        ggtitle("Doublet Density Overlayed") + theme(legend.direction = "horizontal", 
        legend.box.background = element_rect(color = NA))

  }
  
  print(.fixPlotSize(pdensity, plotWidth = 6, plotHeight = 6))

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
    title = "Doublet Scores -log10(FDR)",
    colorTitle = "Doublet Scores -log10(FDR)",
    rastr = TRUE,
    baseSize = 6
    ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  grid::grid.newpage()
  print(.fixPlotSize(pscore, plotWidth = 6, plotHeight = 6))
  
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
    title = "Doublet Enrichment",
    colorTitle = "Doublet Enrichment",
    rastr = TRUE
    ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  grid::grid.newpage()
  print(.fixPlotSize(penrich, plotWidth = 6, plotHeight = 6))

  dev.off()
  sink()
  file.remove(tmpFile)

  summaryList <- SimpleList(
    originalDataUMAP = df,
    simulatedDoubletUMAP = dfDoub
  )
  saveRDS(summaryList, file.path(outDir, paste0(.sampleName(ArrowFile), "-Doublet-Summary.rds")))

  #################################################
  # 6. Add Info To Arrow!
  #################################################
  allCells <- .availableCells(ArrowFile, passQC = FALSE)
  
  allDoubletScores <- rep(-1, length(allCells))
  names(allDoubletScores) <- allCells
  allDoubletScores[names(simDoublets$doubletScore)] <- simDoublets$doubletScore
  
  allDoubletEnrichment <- rep(-1, length(allCells))
  names(allDoubletEnrichment) <- allCells
  allDoubletEnrichment[names(simDoublets$doubletEnrich)] <- simDoublets$doubletEnrich

  allDoubUMAP1 <- rep(NA, length(allCells))
  names(allDoubUMAP1) <- allCells
  allDoubUMAP1[rownames(df)] <- df[,1]

  allDoubUMAP2 <- rep(NA, length(allCells))
  names(allDoubUMAP2) <- allCells
  allDoubUMAP2[rownames(df)] <- df[,2]

  o <- h5closeAll()
  h5write(allDoubUMAP1, file = ArrowFile, "Metadata/Doublet_UMAP1")
  h5write(allDoubUMAP2, file = ArrowFile, "Metadata/Doublet_UMAP2")
  h5write(allDoubletScores, file = ArrowFile, "Metadata/DoubletScore")
  h5write(allDoubletEnrichment, file = ArrowFile, "Metadata/DoubletEnrichment")
  o <- h5closeAll()

  out <- SimpleList(doubletScore = simDoublets$doubletScore, doubletEnrich = simDoublets$doubletEnrich)

  return(out)

}

.simulateProjectDoublets <- function(
  mat, 
  LSI,
  LSIDims,
  uwotUmap,
  clusters = NULL,
  sampleRatio1 = c(0.5), 
  sampleRatio2 = c(0.5), 
  nTrials = 100, 
  nSample = 1000, 
  k = 200, 
  knnMethod = "UMAP",
  seed = 1, 
  threads = 16
  ){

  .sampleSparseMat <- function(mat, sampleRatio = 0.5){
    total <- length(mat@x)
    sampleTo <- floor(total * (1-sampleRatio))
    mat@x[sample(seq_len(total), sampleTo)] <- 0
    mat <- drop0(mat)
    mat
  }

  set.seed(seed)

  if(is.null(clusters)){

    simLSI <- .safelapply(seq_len(nTrials), function(y){

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
        lsiProject <- suppressMessages(ArchR:::.projectLSI(simulatedMat, LSI))

        lsiProject

      }) %>% Reduce("rbind", .)


    }, threads = threads) %>% Reduce("rbind", .)

  }else{

    comClust <- combn(unique(clusters), 2)

    simLSI <- .safelapply(seq_len(ncol(comClust)), function(y){

      if(y %% 5 == 0){
        gc()
      }

      lapply(seq_along(sampleRatio1), function(x){

        idx1 <- sample(which(clusters==comClust[1,y]), nSample, replace = TRUE)
        idx2 <- sample(which(clusters==comClust[2,y]), nSample, replace = TRUE)

        #Simulated Doublet
        simulatedMat <- .sampleSparseMat(mat = mat[,idx1], sampleRatio = sampleRatio1[x]) + 
                        .sampleSparseMat(mat = mat[,idx2], sampleRatio = sampleRatio2[x])

        #Project LSI
        lsiProject <- suppressMessages(ArchR:::.projectLSI(simulatedMat, LSI))

        lsiProject

      }) %>% Reduce("rbind", .)


    }, threads = threads) %>% Reduce("rbind", .)

  }


  #Project UMAP
  set.seed(1) # Always do this prior to UMAP
  umapProject <- data.frame(uwot::umap_transform(as.matrix(simLSI[, LSIDims, drop = FALSE]), uwotUmap, verbose = FALSE, n_threads = threads))

  #Compute KNN 
  if(toupper(knnMethod) == "SVD"){

    knnDoub <- computeKNN(LSI$matSVD, simLSI, k)

  }else if(toupper(knnMethod) == "UMAP"){

    knnDoub <- computeKNN(uwotUmap[[1]], umapProject, k)

  }else{

    stop("Error KNN Method Not Recognized!")

  }

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
  doubletScore <- -log10(pmax(pvalBinomDoub, 4.940656e-324))
  doubletEnrich <- (countKnn / sum(countKnn)) / (1 / nrow(LSI$matSVD))
  doubletEnrich <- 10000 * doubletEnrich / length(countKnn) #Enrichment Per 10000 Cells in Data Set

  out <- SimpleList(doubletUMAP = umapProject, doubletScore = doubletScore, doubletEnrich = doubletEnrich)

  out

}

#' Add Demuxlet Results to ArchR Project
#' 
#' This function will read in the .best file output from demuxlet and add the doublet
#' classifications into the cellColData for the ArchR Project
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param bestFiles The file path to the .best files created by Demuxlet. There should be one .best file for each sample in the `ArchRProject`.
#' @param sampleNames The sample names corresponding to the .best files. These must match the sample names present in the `ArchRProject`.
#' @param ... additional args
#' @export
addDemuxletResults <- function(ArchRProj, bestFiles, sampleNames, ...){
  
  .requirePackage("readr")

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





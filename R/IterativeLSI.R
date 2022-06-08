##########################################################################################
# LSI Dimensionality Reduction Methods
##########################################################################################

#' Add an Iterative LSI-based dimensionality reduction to an ArchRProject
#' 
#' This function will compute an iterative LSI dimensionality reduction on an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix to retrieve from the ArrowFiles associated with the `ArchRProject`. Valid options are
#' "TileMatrix" or "PeakMatrix".
#' @param name The name to use for storage of the IterativeLSI dimensionality reduction in the `ArchRProject` as a `reducedDims` object.
#' @param iterations The number of LSI iterations to perform.
#' @param clusterParams A list of additional parameters to be passed to `addClusters()` for clustering within each iteration. 
#' These params can be constant across each iteration, or specified for each iteration individually. Thus each param must be of
#' length == 1 or the total number of `iterations` - 1. If you want to use `scran` for clustering, you would pass this as `method="scran"`.
#` PLEASE NOTE - We have updated these params to `resolution=2` and `maxClusters=6`! To use previous settings use `resolution=0.2` and `maxClusters=NULL`.
#' @param firstSelection First iteration selection method for features to use for LSI. Either "Top" for the top accessible/average or "Var" for the top variable features. 
#' "Top" should be used for all scATAC-seq data (binary) while "Var" should be used for all scRNA/other-seq data types (non-binary).
#' @param depthCol A column in the `ArchRProject` that represents the coverage (scATAC = unique fragments, scRNA = unique molecular identifiers) per cell.
#' These values are used to minimize the related biases in the reduction related. For scATAC we recommend "nFrags" and for scRNA we recommend "Gex_nUMI".
#' @param varFeatures The number of N variable features to use for LSI. The top N features will be used based on the `selectionMethod`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param LSIMethod A number or string indicating the order of operations in the TF-IDF normalization.
#' Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param scaleDims A boolean that indicates whether to z-score the reduced dimensions for each cell. This is useful forminimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param binarize A boolean value indicating whether the matrix should be binarized before running LSI. This is often desired when working with insertion counts.
#' @param outlierQuantiles Two numerical values (between 0 and 1) that describe the lower and upper quantiles of bias (number of acessible regions per cell, determined 
#' by `nFrags` or `colSums`) to filter cells prior to LSI. For example a value of c(0.02, 0.98) results in the cells in the bottom 2 percent and upper 98 percent to be 
#' filtered prior to LSI. These cells are then projected back in the LSI subspace. This prevents spurious 'islands' that are identified based on being extremely biased.
#' These quantiles are also used for sub-sampled LSI when determining which cells are used.
#' @param filterBias A boolean indicating whether to drop bias clusters when computing clusters during iterativeLSI.
#' @param sampleCellsPre An integer specifying the number of cells to sample in iterations prior to the last in order to perform a sub-sampled LSI and 
#' sub-sampled clustering. This greatly reduced memory usage and increases speed for early iterations.
#' @param projectCellsPre A boolean indicating whether to reproject all cells into the sub-sampled LSI (see `sampleCellsPre`). Setting this to `FALSE`
#' allows for using the sub-sampled LSI for clustering and variance identification. If `TRUE` the cells are all projected into the sub-sampled LSI
#' and used for cluster and variance identification.
#' @param sampleCellsFinal An integer specifying the number of cells to sample in order to perform a sub-sampled LSI in final iteration.
#' @param selectionMethod The selection method to be used for identifying the top variable features. Valid options are "var" for
#' log-variability or "vmr" for variance-to-mean ratio.
#' @param scaleTo Each column in the matrix designated by `useMatrix` will be normalized to a column sum designated by `scaleTo` prior to
#' variance calculation and TF-IDF normalization.
#' @param totalFeatures The number of features to consider for use in LSI after ranking the features by the total number of insertions.
#' These features are the only ones used throught the variance identification and LSI. These are an equivalent when using a `TileMatrix` to a defined peakSet.
#' @param filterQuantile A number [0,1] that indicates the quantile above which features should be removed based on insertion counts prior
#' to the first iteration of the iterative LSI paradigm. For example, if `filterQuantile = 0.99`, any features above the 99th percentile in
#' insertion counts will be ignored for the first LSI iteration.
#' @param excludeChr A string of chromosomes to exclude for iterativeLSI procedure.
#' @param saveIterations A boolean value indicating whether the results of each LSI iterations should be saved as compressed `.rds` files in
#' the designated `outDir`.
#' @param UMAPParams The list of parameters to pass to the UMAP function if "UMAP" if `saveIterations=TRUE`. See the function `uwot::umap()`.
#' @param nPlot If `saveIterations=TRUE`, how many cells to sample make a UMAP and plot for each iteration.
#' @param outDir The output directory for saving LSI iterations if desired. Default is in the `outputDirectory` of the `ArchRProject`.
#' @param threads The number of threads to be used for parallel computing.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param force A boolean value that indicates whether or not to overwrite relevant data in the `ArchRProject` object.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addIterativeLSI <- function(
  ArchRProj = NULL, 
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
      resolution = c(2), 
      sampleCells = 10000,
      maxClusters = 6,
      n.start = 10
  ),
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:30,
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 10000,
  projectCellsPre = FALSE,
  sampleCellsFinal = NULL,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 500000,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = TRUE,
  UMAPParams = list(
    n_neighbors = 40, 
    min_dist = 0.4, 
    metric = "cosine", 
    verbose = FALSE, 
    fast_sgd = TRUE
  ),
  nPlot = 10000,
  outDir = getOutputDirectory(ArchRProj),
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  force = FALSE,
  logFile = createLogFile("addIterativeLSI")
  ){
  
  if(verbose) message("Checking Inputs...")
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = iterations, name = "iterations", valid = c("integer"))
  .validInput(input = clusterParams, name = "clusterParams", valid = c("list"))
  .validInput(input = varFeatures, name = "varFeatures", valid = c("integer"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  .validInput(input = LSIMethod, name = "LSIMethod", valid = c("integer", "character"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = outlierQuantiles, name = "outlierQuantiles", valid = c("numeric", "null"))
  .validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  .validInput(input = sampleCellsPre, name = "sampleCellsPre", valid = c("integer", "null"))
  .validInput(input = sampleCellsFinal, name = "sampleCellsFinal", valid = c("integer", "null"))
  .validInput(input = selectionMethod, name = "selectionMethod", valid = c("character"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = totalFeatures, name = "totalFeatures", valid = c("integer"))
  .validInput(input = filterQuantile, name = "filterQuantile", valid = c("numeric"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = saveIterations, name = "saveIterations", valid = c("boolean"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = outDir, name = "outDir", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))


  if(varFeatures < 1000){
    stop("Please provide more than 1000 varFeatures!")
  }

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "IterativeLSI Input-Parameters", logFile=logFile)

  .requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()

  if(!is.null(ArchRProj@reducedDims[[name]])){
    if(!force){
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }

  #Set Seed
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  #All the Cell Names
  cellNames <- rownames(getCellColData(ArchRProj))
  if(!is.null(sampleCellsPre)){
    if(length(cellNames) < sampleCellsPre){
      sampleCellsPre <- NULL
    }
  }
  if(!is.null(sampleCellsFinal)){
    if(length(cellNames) < sampleCellsFinal){
      sampleCellsFinal <- NULL
    }
  }

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]

  #Check if Matrix is supported and check type
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }else if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }else{
    tileSize <- NA
  }

  units <- unique(unlist(lapply(ArrowFiles, function(x) h5read(x, paste0(useMatrix, "/Info/Units")))))
  if(length(units) != 1){
    stop("Units of matrices are not identical!")
  }
  if(grepl("log",units,ignore.case=TRUE)){
    stop("Cannot use log transformed values for iterativeLSI!")
  }

  tstart <- Sys.time()  
  ############################################################################################################################
  # Organize Information for LSI
  ############################################################################################################################
  chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)

  if(tolower(firstSelection) == "top"){
    
    if(!binarize){
      matClass <- h5read(ArrowFiles[1], paste0(useMatrix,"/Info/Class"))
      if(matClass != "Sparse.Binary.Matrix"){
        stop("Input matrix is not binarized and binarize != TRUE. Please use binarized data if using top selection for first iteration! Set binarize = TRUE!")
      }
    }

    #Compute Row Sums Across All Samples
    .logDiffTime("Computing Total Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    .logDiffTime("Computing Top Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    rmTop <- floor((1-filterQuantile) * totalFeatures)
    topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else if(tolower(firstSelection) %in% c("var", "variable")){

    if(binarize){
      stop("Please do not binarize data if using variable selection for first iteration! Set binarize = FALSE!")
    }

    if(units %in% "BinarizedCounts"){
      stop("Cannot do variable selection with BinarizedCounts. Set firstSelection = Top!")
    }

    #Compute Row Sums Across All Samples
    .logDiffTime("Computing Variability Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- .getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- .getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    .logDiffTime("Computing Variable Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    if(nFeature > 0.5 * nrow(totalAcc)){
      stop("nFeature for variable selection must be less than 1/2 the total features!")
    }
    topIdx <- head(order(totalAcc$combinedVars, decreasing=TRUE), nFeature)
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else{

    stop("firstSelect method must be Top or Var/Variable!")

  }

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = depthCol)
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(e){
      tryCatch({
        .getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)[ArchRProj$cellNames]
      }, error = function(y){
        stop("Could not determine depth from depthCol or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  ############################################################################################################################
  # LSI Iteration 1
  ############################################################################################################################
  .logDiffTime(paste0("Running LSI (1 of ",iterations,") on Top Features"), tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
  j <- 1

  if(!is.null(clusterParams$sampleCells)){
    if(!is.na(clusterParams$sampleCells[j])){
      sampleJ <- clusterParams$sampleCells[j]
    }else if(!is.na(clusterParams$sampleCells[1])){
      sampleJ <- clusterParams$sampleCells[1]
    }else{
      sampleJ <- sampleCellsPre
    }
  }else{
    sampleJ <- sampleCellsPre 
  }

  outLSI <- .LSIPartialMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures,
    cellNames = cellNames, 
    cellDepth = cellDepth,
    useMatrix = useMatrix,
    sampleNames = getCellColData(ArchRProj)$Sample, 
    LSIMethod = LSIMethod, 
    scaleTo = scaleTo,
    dimsToUse = dimsToUse, 
    binarize = binarize, 
    outlierQuantiles = outlierQuantiles,
    sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
    projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
    threads = threads,
    useIndex = FALSE,
    seed = seed,
    tstart = tstart,
    verbose = verbose,
    logFile = logFile
  )
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()
  .logThis(outLSI, paste0("outLSI-",j), logFile = logFile)

  if(iterations == 1){
    .logDiffTime("Finished Running IterativeLSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    ArchRProj@reducedDims[[name]] <- outLSI
    #.endLogging(logFile = logFile)
    return(ArchRProj)
  }

  #########################
  # Identify LSI Clusters
  #########################
  clusterDF <- .LSICluster(
    outLSI = outLSI,
    filterBias = filterBias,
    cellNames = cellNames,
    cellDepth = cellDepth,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    corCutOff = corCutOff,
    clusterParams = clusterParams,
    j = j,
    verbose = verbose,
    tstart = tstart,
    logFile = logFile
  )
  clusters <- clusterDF$clusters
  nClust <- length(unique(clusters))
  .logDiffTime(sprintf("Identified %s Clusters", nClust), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  .logThis(clusterDF, paste0("clusterDF-",j), logFile = logFile)

  #########################
  # Save LSI Iteration
  #########################
  if(saveIterations){
    .logDiffTime("Saving LSI Iteration", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    .saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims, 
      dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
      nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
  }

  ############################################################################################################################
  # LSI Iteration 2+
  ############################################################################################################################
  variableFeatures <- topFeatures

  while(j < iterations){

    #Jth iteration
    j <- j + 1

    #########################
    # Identify Features for LSI Iteration
    #########################
    variableFeatures <- .identifyVarFeatures(
      outLSI = outLSI,
      clusterDF = clusterDF,
      ArrowFiles = ArrowFiles,
      useMatrix = useMatrix,
      prevFeatures = variableFeatures,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      totalFeatures = totalFeatures,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures,
      tstart = tstart,
      threads = threads,
      verbose = verbose,
      logFile = logFile
    )

    #########################
    # LSI
    #########################
    .logDiffTime(sprintf("Running LSI (%s of %s) on Variable Features", j, iterations), tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
    if(!is.null(clusterParams$sampleCells)){
      if(!is.na(clusterParams$sampleCells[j])){
        sampleJ <- clusterParams$sampleCells[j]
      }else if(!is.na(clusterParams$sampleCells[1])){
        sampleJ <- clusterParams$sampleCells[1]
      }else{
        sampleJ <- sampleCellsPre
      }
    }else{
      sampleJ <- sampleCellsPre 
    }

    #Compute Partial Matrix LSI
    outLSI <- .LSIPartialMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = variableFeatures,
      useMatrix = useMatrix,
      cellNames = cellNames, 
      cellDepth = cellDepth,
      sampleNames = getCellColData(ArchRProj)$Sample, 
      LSIMethod = LSIMethod, 
      scaleTo = scaleTo, 
      dimsToUse = dimsToUse,
      binarize = binarize,
      outlierQuantiles = outlierQuantiles, 
      sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
      projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
      threads = threads,
      useIndex = FALSE,
      seed = seed,
      tstart = tstart,
      verbose = verbose,
      logFile = logFile
    )
    outLSI$scaleDims <- scaleDims
    outLSI$useMatrix <- useMatrix
    outLSI$tileSize <- tileSize
    .logThis(outLSI, paste0("outLSI-",j), logFile = logFile)

    if(j != iterations){

      #########################
      # Identify LSI Clusters
      #########################
      clusterDF <- .LSICluster(
        outLSI = outLSI,
        dimsToUse = dimsToUse,
        scaleDims = scaleDims,
        corCutOff = corCutOff,
        filterBias = filterBias,
        cellNames = cellNames,
        cellDepth = cellDepth,
        j = j,
        clusterParams = clusterParams,
        verbose = verbose,
        tstart = tstart,
        logFile = logFile
      )
      clusters <- clusterDF$clusters
      nClust <- length(unique(clusters))
      .logDiffTime(sprintf("Identified %s Clusters", nClust), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      .logThis(clusterDF, paste0("clusterDF-",j), logFile = logFile)

      #########################
      # Save LSI Iteration
      #########################
      if(saveIterations){
        .logDiffTime("Saving LSI Iteration", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        .saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims, 
          dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
          nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
      }

    }

    gc()

  }

  #Organize Output
  .logDiffTime("Finished Running IterativeLSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  ArchRProj@reducedDims[[name]] <- outLSI

  return(ArchRProj)

}

#########################################################################################
# LSI On Partial Matrix
#########################################################################################
.LSIPartialMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  useMatrix = NULL,
  cellNames = NULL, 
  cellDepth = NULL,
  sampleNames = NULL, 
  dimsToUse = NULL, 
  binarize = TRUE, 
  outlierQuantiles = c(0.02, 0.98),
  LSIMethod = FALSE,
  scaleTo = 10^4,
  sampleCells = 5000, 
  projectAll = TRUE,
  threads = 1, 
  seed = 1,
  useIndex = FALSE, 
  tstart = NULL, 
  verbose = TRUE,
  logFile = NULL
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))

  outLSI2 <- tryCatch({

    if(is.null(sampleCells)){
      
      #Construct Matrix
      .logDiffTime("Creating Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      
      mat <- .getPartialMatrix(
        ArrowFiles = ArrowFiles,
        featureDF = featureDF,
        useMatrix = useMatrix,
        cellNames = cellNames,
        doSampleCells = FALSE,
        threads = threads,
        verbose = FALSE,
        logFile = logFile
      )

      #Compute LSI
      .logDiffTime("Computing LSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      outLSI <- .computeLSI(
       mat = mat, 
       LSIMethod = LSIMethod, 
       scaleTo = scaleTo,
       nDimensions = max(dimsToUse),
       binarize = binarize, 
       outlierQuantiles = outlierQuantiles,
       verbose = FALSE, 
       seed = seed,
       tstart = tstart,
       logFile = logFile
      )
      outLSI$LSIFeatures <- featureDF
      outLSI$corToDepth <- list(
        scaled = abs(cor(.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
        none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
      )

    }else{
     
      sampledCellNames <- .sampleBySample(
        cellNames = cellNames,
        sampleNames = sampleNames,
        cellDepth = cellDepth,
        sampleCells = sampleCells,
        outlierQuantiles = outlierQuantiles,
        factor = 2      
      )
      .logDiffTime(sprintf("Sampling Cells (N = %s) for Estimated LSI", length(sampledCellNames)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)

      #Construct Sampled Matrix
      .logDiffTime("Creating Sampled Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      o <- h5closeAll()
      if(!projectAll){

        mat <- .getPartialMatrix(
          ArrowFiles = ArrowFiles,
          featureDF = featureDF,
          useMatrix = useMatrix,
          cellNames = sampledCellNames,
          doSampleCells = FALSE,
          threads = threads,
          verbose = FALSE,
          logFile = logFile
        )

        #Compute LSI
        .logDiffTime("Computing Estimated LSI (projectAll = FALSE)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        outLSI <- .computeLSI(
         mat = mat, 
         LSIMethod = LSIMethod, 
         scaleTo = scaleTo,
         nDimensions = max(dimsToUse),
         binarize = binarize, 
         outlierQuantiles = outlierQuantiles,
         seed = seed,
         tstart = tstart,
         logFile = logFile
        )
        outLSI$LSIFeatures <- featureDF
        outLSI$corToDepth <- list(
          scaled = abs(cor(.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
          none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
        )

      }else{

        tmpPath <- .tempfile(pattern = "tmp-LSI-PM")

        .logDiffTime(sprintf("Sampling Cells (N = %s) for Estimated LSI", length(sampledCellNames)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        out <- .getPartialMatrix(
            ArrowFiles = ArrowFiles,
            featureDF = featureDF,
            useMatrix = useMatrix,
            cellNames = cellNames,
            doSampleCells = TRUE,
            sampledCellNames = sampledCellNames,
            tmpPath = tmpPath,
            useIndex = useIndex,
            threads = threads,
            verbose = FALSE,
            logFile = logFile
          )
        gc()

        #Perform LSI on Partial Sampled Matrix
        .logDiffTime("Computing Estimated LSI (projectAll = TRUE)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        outLSI <- .computeLSI(
           mat = out$mat, 
           LSIMethod = LSIMethod, 
           scaleTo = scaleTo,
           nDimensions = max(dimsToUse),
           binarize = binarize, 
           outlierQuantiles = outlierQuantiles,
           seed = seed,
           tstart = tstart,
           logFile = logFile
          )

        tmpMatFiles <- out[[2]]
        rm(out)
        gc()

        #Read In Matrices and Project into Manifold
        #Do Threads / 3 just in case of memory here (Needs testing QQQ)
        threads2 <- 1 #max(ceiling(threads / 3), 1)
        .logDiffTime("Projecting Matrices with LSI-Projection (Granja* et al 2019)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        pLSI <- .safelapply(seq_along(tmpMatFiles), function(x){
          .logDiffTime(sprintf("Projecting Matrix (%s of %s) with LSI-Projection", x, length(tmpMatFiles)), tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
          .projectLSI(mat = readRDS(tmpMatFiles[x]), LSI = outLSI, verbose = FALSE, tstart = tstart, logFile = logFile)
        }, threads = threads2) %>% Reduce("rbind", .)

        #Remove Temporary Matrices
        rmf <- file.remove(tmpMatFiles)

        #Set To LSI the SVD Matrices
        outLSI$exlcude <- cellNames[which(cellNames %ni% rownames(pLSI))]
        outLSI$matSVD <- as.matrix(pLSI[cellNames[which(cellNames %in% rownames(pLSI))],])

      }

      outLSI$LSIFeatures <- featureDF
      outLSI$corToDepth <- list(
        scaled = abs(cor(.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
        none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
      )

    }

    outLSI

  }, error = function(e){

    errorList$outLSI <- if(exists("outLSI", inherits = FALSE)) outLSI else "Error with outLSI!"
    errorList$matSVD <- if(exists("outLSI", inherits = FALSE)) outLSI[[1]] else "Error with matSVD!"

    .logError(e, fn = ".LSIPartialMatrix", info = "", errorList = errorList, logFile = logFile)

  })

  return(outLSI2)

}


#########################################################################################
# Sampling Helpers
#########################################################################################
.filterSample <- function(
  x = NULL, 
  n = NULL, 
  vals = x, 
  outlierQuantiles = NULL, 
  factor = 2, 
  ...
  ){
  if(!is.null(outlierQuantiles)){
    quant <- quantile(vals, probs = c(min(outlierQuantiles) / factor, 1 - ((1-max(outlierQuantiles)) / factor)))
    idx <- which(vals >= quant[1] & vals <= quant[2])
  }else{
    idx <- seq_along(x)
  }
  if(length(idx) >= n){
    sample(x = x[idx], size = n)
  }else{
    sample(x = x, size = n)
  }
}

.sampleBySample <- function(
  cellNames = NULL, 
  cellDepth = NULL, 
  sampleNames = NULL, 
  sampleCells = NULL, 
  outlierQuantiles = NULL, 
  factor = 2
  ){

  if(sampleCells < length(cellNames)){

    sampleN <- ceiling(sampleCells * table(sampleNames) / length(sampleNames))
    splitCells <- split(cellNames, sampleNames)
    splitDepth <- split(cellDepth, sampleNames)

    sampledCellNames <- lapply(seq_along(splitCells), function(x){
      .filterSample(
        x = splitCells[[x]], 
        n = sampleN[names(splitCells)[x]], 
        vals = splitDepth[[x]], 
        outlierQuantiles = outlierQuantiles, 
        factor = factor
      )
    }) %>% unlist %>% sort

    return(sampledCellNames)
  
  }else{

    return(cellNames)

  }

}


#########################################################################################
# Save LSI Iteration
#########################################################################################

.saveIteration <- function(
  outLSI = NULL, 
  clusters = NULL, 
  nPlot = NULL, 
  UMAPParams = NULL, 
  ArchRProj = NULL, 
  scaleDims = NULL,
  corCutOff = NULL,
  dimsToUse = NULL,
  j = NULL,
  threads = NULL,
  outDir = NULL,
  logFile = NULL
  ){

    errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))

    o <- tryCatch({

      .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Save iteration", logFile=logFile)

      .requirePackage("uwot", source = "cran")

      if(scaleDims){
        dimsPF <- dimsToUse[which(outLSI$corToDepth$scaled[dimsToUse] <= corCutOff)]
      }else{
        dimsPF <- dimsToUse[which(outLSI$corToDepth$none[dimsToUse] <= corCutOff)]
      }

      if(nrow(outLSI[[1]]) > nPlot){
        saveIdx <- sample(seq_len(nrow(outLSI[[1]])), nPlot)
      }else{
        saveIdx <- seq_len(nrow(outLSI[[1]]))
      }
      
      #Plot Quick UMAP
      UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="cosine", verbose=FALSE, fast_sgd = TRUE))
      if(scaleDims){
        UMAPParams$X <- .scaleDims((outLSI[[1]][saveIdx,,drop=FALSE])[, dimsPF, drop = FALSE])
      }else{
        UMAPParams$X <- (outLSI[[1]][saveIdx,,drop=FALSE])[, dimsPF, drop = FALSE]
      }
      UMAPParams$ret_nn <- FALSE
      UMAPParams$ret_model <- FALSE
      UMAPParams$n_threads <- floor(threads / 2)
      uwotUmap <- do.call(uwot::umap, UMAPParams)
      
      #Plot      
      p1 <- ggPoint(
        uwotUmap[,1], 
        uwotUmap[,2], 
        getCellColData(ArchRProj, select = "Sample")[rownames(outLSI[[1]])[saveIdx],], 
        size = 0.5, 
        title = paste0("SampleName (nCells Plot = ",nrow(UMAPParams$X),")"),
        rastr=TRUE
      )
      
      p2 <- ggPoint(
        uwotUmap[,1], 
        uwotUmap[,2], 
        clusters[saveIdx], 
        size = 0.5, 
        title = paste0("Clusters (nCells Plot = ",nrow(UMAPParams$X),")"),
        rastr=TRUE
      )
      
      p1 <- p1 + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
      p2 <- p2 + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank())

      pdf(file.path(outDir, paste0("Save-LSI-Iteration-",j,".pdf")), width = 6, height = 6)
      .fixPlotSize(p1, plotWidth = 6, plotHeight = 6)
      grid::grid.newpage()
      .fixPlotSize(p2, plotWidth = 6, plotHeight = 6)
      dev.off()

      #Save results
      outj <- SimpleList(LSI = outLSI, clusters = clusters, uwotUmap = uwotUmap)
      .safeSaveRDS(outj, file.path(outDir, paste0("Save-LSI-Iteration-",j,".rds")))
      rm(UMAPParams, uwotUmap)
      gc()

    }, error = function(e){

      .logError(e, fn = ".saveIteration", info = "", errorList = errorList, logFile = logFile, throwError = FALSE)

    })

    return(0)

}

#########################################################################################
# Identify Quick Clusters
#########################################################################################
.LSICluster <- function(
  outLSI = NULL, 
  corCutOff = NULL,
  dimsToUse = NULL, 
  scaleDims = NULL, 
  clusterParams = NULL,
  verbose = NULL,
  j = NULL,
  filterBias = NULL,
  cellNames = NULL,
  cellDepth = NULL,
  tstart = NULL,
  logFile = NULL
  ){

  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Cluster Params", logFile=logFile)
  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))

  df2 <- tryCatch({

    if(scaleDims){
      dimsPF <- dimsToUse[which(outLSI$corToDepth$scaled[dimsToUse] <= corCutOff)]
    }else{
      dimsPF <- dimsToUse[which(outLSI$corToDepth$none[dimsToUse] <= corCutOff)]
    }
    if(length(dimsPF)!=length(dimsToUse)){
      message("Filtering ", length(dimsToUse) - length(dimsPF), " dims correlated > ", corCutOff, " to log10(depth + 1)")
    }
    if(length(dimsPF) < 2){
      stop("Dimensions to use after filtering for correlation to depth lower than 2!")
    }

    #Time to compute clusters
    .logDiffTime("Identifying Clusters", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    parClust <- lapply(clusterParams, function(x){
      if(length(x) > 1){
        return(x[[j]])
      }else{
        return(x[[1]])
      }
    })
    parClust$verbose <- FALSE
    if(scaleDims){
      parClust$input <- .scaleDims(outLSI$matSVD)[, dimsPF, drop = FALSE]
    }else{
      parClust$input <- outLSI$matSVD[, dimsPF, drop = FALSE]
    }

    parClust$input <- as.matrix(parClust$input)

    if(filterBias){
      parClust$testBias <- TRUE
      parClust$filterBias <- TRUE
    }
    parClust$biasVals <- data.frame(row.names = cellNames, x = cellDepth)[rownames(outLSI$matSVD), 1]
    parClust$logFile <- logFile

    clusters <- do.call(addClusters, parClust)
    
    parClust$input <- NULL
    nClust <- length(unique(clusters))  
    
    df <- DataFrame(cellNames = rownames(outLSI$matSVD), clusters = clusters)
    metadata(df)$parClust <- parClust
    df

  }, error = function(e){

    .logError(e, fn = ".LSICluster", info = "", errorList = errorList, logFile = logFile, throwError = FALSE)

  })

  df2

}

#########################################################################################
# Identify Variable Features
#########################################################################################
.identifyVarFeatures <- function(
  outLSI = NULL,
  clusterDF = NULL,
  prevFeatures = NULL,
  ArrowFiles = NULL,
  useMatrix = NULL,
  totalAcc = NULL,
  scaleTo = NULL,
  firstSelection = NULL,
  totalFeatures = NULL,
  selectionMethod = NULL,
  varFeatures = NULL,
  tstart = NULL,
  threads = NULL,
  verbose = NULL,
  logFile = NULL
  ){

  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))

  variableFeatures2 <- tryCatch({

    nClust <- length(unique(clusterDF$clusters))

    if(nClust == 1){

      .logDiffTime("Identified 1 Cluster, Continuing with Previous Features!", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      return(prevFeatures)
    
    }else{

      #Random Sampling for Quick Estimation of Variance
      parClust <- metadata(clusterDF)$parClust
      if(!is.null(parClust$sampleCells)){
        if(is.numeric(parClust$sampleCells)){
          if(floor(parClust$sampleCells) < nrow(outLSI$matSVD)){
            idxSub <- sort(sample(seq_len(nrow(outLSI$matSVD)), floor(parClust$sampleCells)))
          }else{
            idxSub <- seq_len(nrow(outLSI$matSVD))
          }
        }else{
          idxSub <- seq_len(nrow(outLSI$matSVD))
        }
      }else{
        idxSub <- seq_len(nrow(outLSI$matSVD))
      }

      .logDiffTime("Creating Cluster Matrix on the total Group Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      groupList <- SimpleList(split(clusterDF$cellNames, clusterDF$clusters))

      if(tolower(firstSelection) == "top"){
        groupFeatures <- totalAcc[sort(head(order(totalAcc$rowSums, decreasing = TRUE), totalFeatures)),]
      }else if(tolower(firstSelection) %in% c("var", "variable")){
        groupFeatures <- totalAcc
      }

      groupMat <- .getGroupMatrix(
        ArrowFiles = ArrowFiles, 
        featureDF = groupFeatures,
        useMatrix = useMatrix, 
        threads = threads,
        groupList = groupList,
        useIndex = FALSE,
        verbose = FALSE
      )

      .logDiffTime("Computing Variable Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      if(length(varFeatures) > 1){
        nFeature <- varFeatures[j]
      }else{
        nFeature <- varFeatures
      }

      if(tolower(selectionMethod) == "var"){

        #Log-Normalize
        groupMat <- log2(t(t(groupMat) / colSums(groupMat)) * scaleTo + 1)
        var <- matrixStats::rowVars(groupMat)
        idx <- sort(head(order(var,decreasing=TRUE), nFeature))
        variableFeatures <- groupFeatures[idx,]

      }else if(tolower(selectionMethod) == "vmr"){

        #Variance-to-Mean Ratio
        vmr <- matrixStats::rowVars(groupMat) / rowMeans(groupMat)
        idx <- sort(head(order(vmr, decreasing=TRUE), nFeature))
        variableFeatures <- groupFeatures[idx,]

      }else{

        stop("Error selectionMethod is not valid requires var or vmr")

      }

    }

    variableFeatures

  }, error = function(e){

    .logError(e, fn = ".identifyVarFeatures", info = "", errorList = errorList, logFile = logFile)

  })

  return(variableFeatures2)

}

#########################################################################################
# LSI Methods
#########################################################################################
.computeLSI <- function(
  mat = NULL, 
  LSIMethod = 1,
  scaleTo = 10^4,
  nDimensions = 50, 
  binarize = TRUE, 
  outlierQuantiles = c(0.02, 0.98),
  seed = 1, 
  verbose = FALSE, 
  tstart = NULL,
  logFile = NULL
  ){

  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "LSI Parameters", logFile=logFile)

  out2 <- tryCatch({

    set.seed(seed)

    if(is.null(tstart)){
      tstart <- Sys.time()
    }
  
    .logDiffTime(sprintf("Running LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose, logFile = logFile)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        .logDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        mat@x[mat@x > 0] <- 1 
    }

    #Compute Col Sums
    .logDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude, drop = FALSE]
      colSm <- colSm[-exclude]
    }else{
      exclude <- c()
    }

    cn <- colnames(mat)
    filterOutliers <- 0
    if(!is.null(outlierQuantiles)){
      qCS <- quantile(colSm, probs = c(min(outlierQuantiles), max(outlierQuantiles)))
      idxOutlier <- which(colSm <= qCS[1] | colSm >= qCS[2])
      if(length(idxOutlier) > 0){
        .logDiffTime("Filtering Outliers Based On Counts", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        #.safeSaveRDS(mat, "temp.rds", compress = FALSE)
        matO <- mat[, idxOutlier, drop = FALSE]
        mat <- mat[, -idxOutlier, drop = FALSE]
        mat2 <- mat[, head(seq_len(ncol(mat)), 50), drop = FALSE] # A 2nd Matrix to Check Projection is Working
        colSm <- colSm[-idxOutlier]
        filterOutliers <- 1       
      }
    }

    #Clean up zero rows
    .logDiffTime("Removing 0 Sum Rows", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    rowSm <- Matrix::rowSums(mat)
    idx <- which(rowSm > 0)
    mat <- mat[idx, ]
    rowSm <- rowSm[idx]

    #TF - Normalize
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSIMethod == 1 | tolower(LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-LogIDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSIMethod == 2 | tolower(LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(ncol(mat) / rowSm, "sparseVector")

      #TF-IDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * scaleTo + 1)  

    }else if(LSIMethod == 3 | tolower(LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-IDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Calc SVD then LSI
    .logDiffTime("Computing SVD using irlba", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    svd <- irlba::irlba(mat, nDimensions, nDimensions)
    svdDiag <- matrix(0, nrow=nDimensions, ncol=nDimensions)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))

    #Return Object
    out <- SimpleList(
        matSVD = matSVD, 
        rowSm = rowSm, 
        nCol = length(colSm),
        exclude = exclude, 
        idx = idx, 
        svd = svd, 
        binarize = binarize, 
        scaleTo = scaleTo,
        nDimensions = nDimensions,
        LSIMethod = LSIMethod,
        outliers = NA,
        date = Sys.Date(),
        seed = seed
      )

    if(filterOutliers == 1){
      .logDiffTime("Projecting Outliers with LSI-Projection (Granja* et al 2019)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      #Quick Check LSI-Projection Works
      pCheck <- .projectLSI(mat = mat2, LSI = out, verbose = verbose, logFile = logFile)
      pCheck2 <- out[[1]][rownames(pCheck), ]
      pCheck3 <- lapply(seq_len(ncol(pCheck)), function(x){
        cor(pCheck[,x], pCheck2[,x])
      }) %>% unlist
      if(min(pCheck3) < 0.95){
        .logThis(pCheck, "pCheck", logFile=logFile)
        .logThis(pCheck2, "pCheck2", logFile=logFile)
        .logThis(pCheck3, "pCheck3", logFile=logFile)
        warning("Warning with LSI-projection! Cor less than 0.95 of re-projection. Please report this to github with logFile!")
      }
      #Project LSI Outliers
      out$outliers <- colnames(matO)
      outlierLSI <- .projectLSI(mat = matO, LSI = out, verbose = verbose, logFile = logFile)
      allLSI <- rbind(out[[1]], outlierLSI)
      allLSI <- allLSI[cn, , drop = FALSE] #Re-Order Correctly to original
      out[[1]] <- allLSI
    }
    .logDiffTime("Finished LSI (TF-IDF SVD) using irlba", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)

    rm(mat)
    gc()

    out

  }, error = function(e){
    
    errorList <- list(
      mat = mat,
      colSm = if(exists("colSm", inherits = FALSE)) colSm else "Error with colSm!",
      idf = if(exists("idf", inherits = FALSE)) idf else "Error with idf!",
      V = if(exists("V", inherits = FALSE)) V else "Error with V!",
      matSVD = if(exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!",
      out = if(exists("out", inherits = FALSE)) out else "Error with out!"
    )

    .logError(e, fn = ".computeLSI", info = "", errorList = errorList, logFile = logFile)

  })

  out2

}

.projectLSI <- function(
  mat = NULL, 
  LSI = NULL, 
  returnModel = FALSE, 
  verbose = FALSE, 
  tstart = NULL,
  logFile = NULL
  ){   

  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "LSI-Projection Parameters", logFile=logFile)
    
  out2 <- tryCatch({
   
    require(Matrix)
    set.seed(LSI$seed)
    
    if(is.null(tstart)){
      tstart <- Sys.time()
    }

    .logDiffTime(sprintf("Projecting LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose, logFile = logFile)

    #Get Same Features
    .logDiffTime("Subsetting by Non-Zero features in inital Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    mat <- mat[LSI$idx,]

    #Binarize Matrix
    if(LSI$binarize){
        .logDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        mat@x[mat@x > 0] <- 1       
    }
    
    #TF
    .logDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSI$LSIMethod == 1 | tolower(LSI$LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-LogIDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")

      #TF-IDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)  

    }else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      .logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-IDF
      .logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
        .logDiffTime(sprintf("Zeroing %s NA elements", length(idxNA)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        mat[idxNA] <- 0
    }

    #Calc V
    .logDiffTime("Calculating V Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)

    #LSI Diagonal
    .logDiffTime("Computing Projected Coordinates", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))

    if(returnModel){
        .logDiffTime("Calculating Re-Projected Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
        out <- list(matSVD = matSVD, V = V, X = X)
    }else{
        out <- matSVD
    }

    out


  }, error = function(e){

    errorList <- list(
      mat = mat,
      colSm = if(exists("colSm", inherits = FALSE)) colSm else "Error with colSm!",
      idf = if(exists("idf", inherits = FALSE)) idf else "Error with idf!",
      V = if(exists("V", inherits = FALSE)) V else "Error with V!",
      matSVD = if(exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!"
    )

    .logError(e, fn = ".projectLSI", info = "", errorList = errorList, logFile = logFile)

  })

  out2

}



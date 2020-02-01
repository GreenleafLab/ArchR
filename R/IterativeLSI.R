##########################################################################################
# LSI Dimensionality Reduction Methods
##########################################################################################

#' Add an LSI-based dimensionality reduction to an ArchRProject JJJ
#' 
#' This function will compute an LSI dimensionality reduction on an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix to retrieve from the ArrowFiles associated with the `ArchRProject`. Valid options are "TileMatrix" or "PeakMatrix".
#' @param reducedDimsOut The name to use for storage of the LSI dimensionality reduction in the `ArchRProject` as a `reducedDims` object.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean describing whether to rescale the total variance for each principal component. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param LSIMethod A number or string indicating the order of operations in the TF-IDF normalization.
#' Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param binarize A boolean value indicating whether the matrix should be binarized before running LSI. This is often desired when working with insertion counts.
#' @param sampleCells An integer specifying the number of cells to sample in order to perform a sub-sampled LSI and sub-sampled clustering.
#' @param topFeatures The number of N top accessible features to use for LSI.
#' @param scaleTo Each column in the matrix designated by `useMatrix` will be normalized to a column sum designated by `scaleTo` prior to TF-IDF normalization.
#' @param totalFeatures The number of features to consider for use in LSI after ranking the features by the total insertion counts. These are an equivalent when using a `TileMatrix` to a defined peakSet.
#' @param filterQuantile A number [0,1] that indicates the quantile above which features should be removed based on insertion counts prior to the LSI reduction. For example, if `filterQuantile = 0.99`, any features above the 99th percentile in insertion counts will be ignored for LSI reduction.
#' @param runHarmony A boolean value indicating whether harmony-based batch correction should be run on the computed LSI object.
#' @param harmonyParams Additional parameters to be passed to `harmony::HarmonyMatrix()`.
#' @param threads The number of threads to be used for parallel computing.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param force A boolean value that indicates whether or not to overwrite relevant data in the `ArchRProject` object.
#' @export
addLSI <- function(
  ArchRProj = NULL, 
  useMatrix = "TileMatrix",
  reducedDimsOut = "LSI",
  dimsToUse = 1:30,
  scaleDims = TRUE,
  corCutOff = 0.75,
  LSIMethod = 2,
  binarize = TRUE,
  sampleCells = NULL,
  topFeatures = 50000,
  totalFeatures = 500000,
  filterQuantile = 0.995,
  runHarmony = FALSE,
  harmonyParams = list(),
  threads = getArchRThreads(),
  seed = 1,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  force = FALSE
  ){

  args <- mget(names(formals()),sys.frame(sys.nframe()))
  args$iterations <- 1
  args$varFeatures <- args$topFeatures
  args$saveIterations <- FALSE
  do.call(addIterativeLSI, args)

}


#' Add an Iterative LSI-based dimensionality reduction to an ArchRProject
#' 
#' This function will compute an iterative LSI dimensionality reduction on an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix to retrieve from the ArrowFiles associated with the `ArchRProject`. Valid options are "TileMatrix" or "PeakMatrix".
#' @param reducedDimsOut The name to use for storage of the IterativeLSI dimensionality reduction in the `ArchRProject` as a `reducedDims` object.
#' @param iterations The number of LSI iterations to perform.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean describing whether to rescale the total variance for each principal component. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param LSIMethod A number or string indicating the order of operations in the TF-IDF normalization.
#' Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param binarize A boolean value indicating whether the matrix should be binarized before running LSI. This is often desired when working with insertion counts.
#' @param sampleCells An integer specifying the number of cells to sample in order to perform a sub-sampled LSI and sub-sampled clustering.
#' @param varFeatures The number of N variable features to use for LSI. The top N features will be used based on the `selectionMethod`.
#' @param selectionMethod The selection method to be used for identifying the top variable features. Valid options are "var" for log-variability or "vmr" for variance-to-mean ratio.
#' @param scaleTo Each column in the matrix designated by `useMatrix` will be normalized to a column sum designated by `scaleTo` prior to variance calculation and TF-IDF normalization.
#' @param totalFeatures The number of features to consider for use in LSI after ranking the features by the total number of insertions. These features are the only ones used throught the variance identification and LSI. These are an equivalent when using a `TileMatrix` to a defined peakSet.
#' @param filterQuantile A number [0,1] that indicates the quantile above which features should be removed based on insertion counts prior to the first iteration of the iterative LSI paradigm. For example, if `filterQuantile = 0.99`, any features above the 99th percentile in insertion counts will be ignored for the first LSI iteration.
#' @param saveIterations A boolean value indicating whether the results of each LSI iterations should be saved as compressed `.rds` files in the designated `outDir`.
#' @param outDir The output directory for saving LSI iterations if desired. Default is in the `outputDirectory` of the `ArchRProject`.
#' @param clusterParams Additional parameters to be passed to `addClusters()` for clustering within each iteration. These must be either length 1 or the total number of `iterations` - 1.
#' @param runHarmony A boolean value indicating whether harmony-based batch correction should be run during the LSI iterations.
#' @param harmonyParams Additional parameters to be passed to `harmony::HarmonyMatrix()`.
#' @param threads The number of threads to be used for parallel computing.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param force A boolean value that indicates whether or not to overwrite relevant data in the `ArchRProject` object.
#' @export
addIterativeLSI <- function(
  ArchRProj = NULL, 
  useMatrix = "TileMatrix",
  reducedDimsOut = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 25),
  dimsToUse = 1:30,
  scaleDims = TRUE,
  corCutOff = 0.75,
  LSIMethod = 2,
  binarize = TRUE,
  sampleCells = NULL,
  varFeatures = 50000,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 500000,
  filterQuantile = 0.995,
  saveIterations = TRUE,
  outDir = getOutputDirectory(ArchRProj),
  runHarmony = FALSE,
  harmonyParams = list(),
  threads = getArchRThreads(),
  seed = 1,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  force = FALSE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = iterations, name = "iterations", valid = c("integer"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = LSIMethod, name = "LSIMethod", valid = c("integer", "character"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("integer","null"))
  .validInput(input = varFeatures, name = "varFeatures", valid = c("integer"))
  .validInput(input = selectionMethod, name = "selectionMethod", valid = c("character"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = totalFeatures, name = "totalFeatures", valid = c("integer"))
  .validInput(input = filterQuantile, name = "filterQuantile", valid = c("numeric"))
  .validInput(input = saveIterations, name = "saveIterations", valid = c("boolean"))
  .validInput(input = outDir, name = "outDir", valid = c("character"))
  .validInput(input = clusterParams, name = "clusterParams", valid = c("list"))
  .validInput(input = runHarmony, name = "runHarmony", valid = c("boolean"))
  .validInput(input = harmonyParams, name = "harmonyParams", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
  .validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))


  .requirePackage("Matrix")
  tstart <- Sys.time()

  if(!is.null(ArchRProj@reducedDims[[reducedDimsOut]])){
    if(!force){
      stop("Error ReducedDimsOut Already Exists! Set force = TRUE or pick a different name!")
    }
  }

  #Set Seed
  set.seed(seed)
  outDir <- file.path(outDir, reducedDimsOut)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  #All the Cell Names
  cellNames <- rownames(getCellColData(ArchRProj))
  if(!is.null(sampleCells)){
    if(length(cellNames) < sampleCells){
      sampleCells <- NULL
    }
  }

  #Check if Matrix is supported
  stopifnot(any(tolower(useMatrix) %in% c("tilematrix","peakmatrix")))
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
  }
  if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
  }

  tstart <- Sys.time()
  .messageDiffTime(paste0("Running LSI (1 of ",iterations,") on Top Features"), tstart, addHeader = TRUE, verbose = verboseHeader)

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]
  chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)

  #Compute Row Sums Across All Samples
  .messageDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
  totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
  gc()

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = "nFrags")
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(x){
      tryCatch({
        .getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
      }, error = function(y){
        stop("Could not determine depth from nFrags or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  #Identify the top features to be used here
  .messageDiffTime("Computing Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
  nFeature <- varFeatures[1]
  rmTop <- floor((1-filterQuantile) * totalFeatures)
  topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
  topFeatures <- totalAcc[sort(topIdx),]

  #Compute Partial Matrix LSI
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
    sampleCells = sampleCells,
    threads = threads,
    useIndex = FALSE,
    tstart = tstart
    )
  outLSI$scaleDims <- scaleDims
  gc()

  if(runHarmony){
    .messageDiffTime("Harmonizing LSI output on the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    .requirePackage("harmony")
    harmonyParams$data_mat <- outLSI$matSVD
    harmonyParams$meta_data <- data.frame(row.names = rownames(outLSI$matSVD), Group = stringr::str_split(rownames(outLSI$matSVD), pattern = "#", simplify=TRUE)[,1])
    harmonyParams$do_pca <- FALSE
    harmonyParams$vars_use <- "Group"
    harmonyParams$plot_convergence <- FALSE
    harmonyParams$verbose <- verboseAll
    #Harmonize the LSI Results
    outLSI$matSVD <- do.call(HarmonyMatrix, harmonyParams)
  }

  #Time to compute clusters
  .messageDiffTime("Identifying Clusters", tstart, addHeader = verboseAll, verbose = verboseHeader)
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
  parClust <- lapply(clusterParams, function(x) x[[1]])
  if(scaleDims){
      parClust$input <- .scaleDims(outLSI$matSVD)
  }else{
    parClust$input <- outLSI$matSVD[, dimsPF, drop = FALSE]
  }
  parClust$verbose <- verboseAll
  clusters <- do.call(addClusters, parClust)
  
  #Save Output
  if(saveIterations){
    .messageDiffTime("Saving LSI Iteration", tstart, addHeader = verboseAll, verbose = verboseHeader)
    outj <- SimpleList(LSI = outLSI, clusters = clusters, params = parClust[-length(parClust)])
    saveRDS(outj, file.path(outDir, paste0("Save-LSI-Iteration-1.rds")))
  }

  j <- 1
  while(j < iterations){

    #Jth iteration
    j <- j + 1
    
    .messageDiffTime(sprintf("Running LSI (%s of %s) on Variable Features", j, iterations), tstart, addHeader = TRUE, verbose = verboseHeader)
    
    #Create Group Matrix
    .messageDiffTime("Creating Cluster Matrix on the total Group Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    groupList <- SimpleList(split(rownames(outLSI$matSVD), clusters))
    groupFeatures <- totalAcc[sort(head(order(totalAcc$rowSums, decreasing = TRUE), totalFeatures)),]
    groupMat <- .getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = groupFeatures,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = groupList,
      useIndex = FALSE,
      verbose = verboseAll
    )

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

      stop("Error Selection Method is not Valid requires var or vmr")

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
      sampleCells = sampleCells,
      threads = threads,
      useIndex = FALSE,
      tstart = tstart
      )
    outLSI$scaleDims <- scaleDims

    if(runHarmony){
      .messageDiffTime("Harmonizing LSI output on the Variable Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
      harmonyParams$data_mat <- outLSI$matSVD
      harmonyParams$meta_data <- data.frame(row.names = rownames(outLSI$matSVD), Group = stringr::str_split(rownames(outLSI$matSVD), pattern = "#", simplify=TRUE)[,1])
      harmonyParams$do_pca <- FALSE
      harmonyParams$vars_use <- "Group"
      harmonyParams$plot_convergence <- FALSE
      harmonyParams$verbose <- verboseAll
      #Harmonize the LSI Results
      outLSI$matSVD <- do.call(HarmonyMatrix, harmonyParams)
    }

    if(j != iterations){

      #Time to compute clusters
      .messageDiffTime("Identifying Clusters", tstart, addHeader = verboseAll, verbose = verboseHeader)
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
      parClust <- lapply(clusterParams, function(x){
        if(length(x) > 1){
          return(x[[j]])
        }else{
          return(x[[1]])
        }
      })

      if(scaleDims){
          parClust$input <- .scaleDims(outLSI$matSVD)
      }else{
        parClust$input <- outLSI$matSVD[, dimsPF, drop = FALSE]
      }
      parClust$verbose <- verboseAll
      clusters <- do.call(addClusters, parClust)

      #Save Output
      if(saveIterations){
        .messageDiffTime("Saving LSI Iteration", tstart, addHeader = verboseAll, verbose = verboseHeader)
        outj <- SimpleList(LSI = outLSI, clusters = clusters, params = parClust[-length(parClust)])
        saveRDS(outj, file.path(outDir, paste0("Save-LSI-Iteration-",j,".rds")))
      }

    }

    gc()

  }

  #Organize Output
  .messageDiffTime("Finished Running IterativeLSI", tstart, addHeader = verboseAll, verbose = verboseHeader)
  ArchRProj@reducedDims[[reducedDimsOut]] <- outLSI

  return(ArchRProj)

}

.LSIPartialMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  useMatrix = NULL,
  cellNames = NULL, 
  cellDepth = NULL,
  sampleNames = NULL, 
  dimsToUse = NULL, 
  binarize = TRUE, 
  LSIMethod = FALSE,
  scaleTo = 10^4,
  sampleCells = 5000, 
  threads = 1, 
  useIndex = FALSE, 
  tstart = NULL, 
  verboseHeader = TRUE,
  verboseAll = FALSE
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  if(is.null(sampleCells)){
    
    #Construct Matrix
    .messageDiffTime("Creating Partial Matrix of Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    
    mat <- .getPartialMatrix(
      ArrowFiles = ArrowFiles,
      featureDF = featureDF,
      useMatrix = useMatrix,
      cellNames = cellNames,
      doSampleCells = FALSE,
      threads = threads,
      verbose = verboseAll
    )

    #Compute LSI
    .messageDiffTime("Running LSI on the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    outLSI <- .computeLSI(
     mat = mat, 
     LSIMethod = LSIMethod, 
     scaleTo = scaleTo,
     nDimensions = max(dimsToUse),
     binarize = binarize, 
     verbose = verboseAll, 
     tstart = tstart
     )
  
  }else{
   
    set.seed(1)
    .messageDiffTime("Sampling Cells for Estimated LSI", tstart, addHeader = verboseAll, verbose = verboseHeader)
    sampleN <- floor(sampleCells * table(sampleNames) / length(sampleNames))
    splitCells <- split(cellNames, sampleNames)
    sampledCellNames <- lapply(seq_along(splitCells), function(x){
      sample(splitCells[[x]], sampleN[names(splitCells)[x]])
    }) %>% unlist %>% sort

    #Construct Sampled Matrix
    .messageDiffTime("Creating Sampled Partial Matrix of Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    tmpPath <- .tempfile(pattern = "tmp-LSI-PM")
    o <- h5closeAll()
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
        verbose = verboseAll
      )
    gc()

    #Perform LSI on Partial Sampled Matrix
    .messageDiffTime("Running Sampled LSI on the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    outLSI <- .computeLSI(
       mat = out$mat, 
       LSIMethod = LSIMethod, 
       scaleTo = scaleTo,
       nDimensions = max(dimsToUse),
       binarize = binarize, 
       verbose = verboseAll, 
       tstart = tstart
      )
    tmpMatFiles <- out[[2]]
    rm(out)
    gc()

    #Read In Matrices and Project into Manifold
    .messageDiffTime("Projecting Matrices with the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    pLSI <- lapply(seq_along(tmpMatFiles), function(x){
      .projectLSI(mat = readRDS(tmpMatFiles[x]), LSI = outLSI, verbose = FALSE, tstart = tstart)
    }) %>% Reduce("rbind", .)

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

  return(outLSI)

}

.computeLSI <- function(
  mat = NULL, 
  LSIMethod = 1,
  scaleTo = 10^4,
  nDimensions = 50, 
  binarize = TRUE, 
  seed = 1, 
  verbose = TRUE, 
  tstart = NULL
  ){

    set.seed(seed)

    if(is.null(tstart)){
      tstart <- Sys.time()
    }
  
    .messageDiffTime(sprintf("Running LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        .messageDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose)
        mat@x[mat@x > 0] <- 1 
    }

    #Clean up zero rows
    .messageDiffTime("Removing 0 Sum Rows", tstart, addHeader = FALSE, verbose = verbose)
    rowSm <- Matrix::rowSums(mat)
    idx <- which(rowSm > 0)
    mat <- mat[idx,]
    rowSm <- rowSm[idx]

    #TF
    .messageDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose)
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }else{
      exclude <- c()
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSIMethod == 1 | tolower(LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-LogIDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSIMethod == 2 | tolower(LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(ncol(mat) / rowSm, "sparseVector")

      #TF-IDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * scaleTo + 1)  

    }else if(LSIMethod == 3 | tolower(LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-IDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Calc SVD then LSI
    .messageDiffTime("Computing SVD using irlba", tstart, addHeader = FALSE, verbose = verbose)
    svd <- irlba::irlba(mat, nDimensions, nDimensions)
    svdDiag <- matrix(0, nrow=nDimensions, ncol=nDimensions)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    #Return Object
    .messageDiffTime("Finished LSI (TF-IDF SVD) using irlba", tstart, addHeader = FALSE, verbose = verbose)
    out <- SimpleList(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm,
        exclude = exclude, 
        idx = idx, 
        svd = svd, 
        binarize = binarize, 
        scaleTo = scaleTo,
        nDimensions = nDimensions,
        LSIMethod = LSIMethod,
        date = Sys.Date(),
        seed = seed
      )

    rm(mat)
    gc()

    out
}

.projectLSI <- function(
  mat = NULL, 
  LSI = NULL, 
  returnModel = FALSE, 
  verbose = TRUE, 
  tstart = NULL
  ){   
    
    require(Matrix)
    set.seed(LSI$seed)
    
    if(is.null(tstart)){
      tstart <- Sys.time()
    }

    .messageDiffTime(sprintf("Projecting LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose)

    #Get Same Features
    .messageDiffTime("Subsetting by Non-Zero features in inital Matrix", tstart, addHeader = FALSE, verbose = verbose)
    mat <- mat[LSI$idx,]

    #Binarize Matrix
    if(LSI$binarize){
        .messageDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose)
        mat@x[mat@x > 0] <- 1       
    }
    
    #TF
    .messageDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose)
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
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(log(1 + length(LSI$colSm) / LSI$rowSm), "sparseVector")

      #TF-LogIDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(length(LSI$colSm) / LSI$rowSm, "sparseVector")

      #TF-IDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)  

    }else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
      idf   <- as(log(1 + length(LSI$colSm) / LSI$rowSm), "sparseVector")

      #TF-IDF
      .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
        .messageDiffTime(sprintf("Zeroing %s NA elements", length(idxNA)), tstart, addHeader = FALSE, verbose = verbose)
        mat[idxNA] <- 0
    }

    #Calc V
    .messageDiffTime("Calculating V Matrix", tstart, addHeader = FALSE, verbose = verbose)
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)

    #LSI Diagonal
    .messageDiffTime("Computing Projected Coordinates", tstart, addHeader = FALSE, verbose = verbose)
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    if(returnModel){
        .messageDiffTime("Calculating Re-Projected Matrix", tstart, addHeader = FALSE, verbose = verbose)
        X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
        out <- list(matSVD = matSVD, V = V, X = X)
    }else{
        out <- matSVD
    }

    return(out)
}











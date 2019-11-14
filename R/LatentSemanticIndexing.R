#' Compute Iterative LSI
#' 
#' This function will compute an iterative LSI dimensionality reduction
#' on an ArchRProject.
#'
#' @param ArchRProj ArchRProject
#' @param useMatrix use matrix for LSI clustering from Arrow
#' @param reducedDimsOut name of dimensionality reduction to be stored as
#' @param iterations number of LSI iterations to perform
#' @param dimsToUse number of dimensions to compute and use from LSI (TFIDF-SVD) for clustering
#' @param binarize binarize matrix prior to LSI
#' @param sampleCells number of cells to sample for LSI estimation
#' @param varFeatures number of variable features to use for LSI
#' @param selectionMethod selection method for variable features (var or vmr)
#' @param scaleTo scaleTo for Cluster Averages for variance calculation
#' @param totalFeatures number of features to consider (ranked by total number of counts) use for LSI
#' @param filterQuantile filter features for initial LSI that are above this quantile
#' @param saveIterations save LSI iterations as rds in the outDir
#' @param outDir output directory for saving LSI iterations
#' @param clusterParams additional params to pass to IdentifyClusters
#' @param runHarmony run harmony batch correction through the iterations
#' @param harmonyParams additional params to pass to harmony
#' @param threads number of threads for parallel execution
#' @param seed seed for analysis
#' @param verboseHeader verbose sections
#' @param verboseAll verbose sections and subsections
#' @param force verbose sections and subsections
#' @param ... additional args
#' @export
IterativeLSI <- function(
  ArchRProj = NULL, 
  useMatrix = "TileMatrix",
  reducedDimsOut = "IterativeLSI",
  iterations = 3,
  dimsToUse = 1:25,
  binarize = TRUE,
  sampleCells = 5000,
  varFeatures = 50000,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 500000,
  filterQuantile = 0.99,
  saveIterations = TRUE,
  outDir = getOutputDirectory(ArchRProj),
  clusterParams = list(),
  runHarmony = FALSE,
  harmonyParams = list(),
  threads = 1,
  seed = 1,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  force = FALSE,
  ...){
  
  .requirePackage("Matrix")
  tstart <- Sys.time()

  if(!is.null(ArchRProj@reducedDims[[reducedDimsOut]])){
    if(!force){
      stop("Error ReducedDimsOut Already Exists! Set force = TRUE or pick a different name!")
    }
  }

  #What Parameters To Pass
  defaultClustParams <- list(
      method = "Seurat",
      resolution = c(0.4, 0.6),
      n.start = c(10, 10),
      verbose = TRUE
  )
  clusterParams <- .mergeParams(clusterParams, defaultClustParams)

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
  .messageDiffTime(paste0("Computing IterativeLSI on ", useMatrix), tstart, addHeader = TRUE, verbose = verboseHeader)

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]
  chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)

  #Compute Row Sums Across All Samples
  .messageDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
  totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
  gc()

  #Identify the top features to be used here
  .messageDiffTime("Computing Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
  nFeature <- varFeatures[1]
  rmTop <- floor((1-filterQuantile) * totalFeatures)
  topIdx <- head(order(totalAcc$value, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
  topFeatures <- totalAcc[sort(topIdx),]

  #Compute Partial Matrix LSI
  outLSI <- .LSIPartialMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures,
    cellNames = cellNames, 
    sampleNames = getCellColData(ArchRProj)$Sample, 
    dimsToUse = dimsToUse, 
    binarize = binarize, 
    sampleCells = sampleCells,
    threads = threads,
    useIndex = FALSE,
    tstart = tstart
    )
  outLSI$LSIFeatures <- topFeatures
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
  parClust <- lapply(clusterParams, function(x) x[[1]])
  parClust$input <- outLSI$matSVD
  parClust$sampleCells <- sampleCells
  parClust$verbose <- verboseAll
  clusters <- do.call(IdentifyClusters, parClust)
  
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
    
    .messageDiffTime(sprintf("Running LSI %s of %s on Variable Features", j, iterations), tstart, addHeader = TRUE, verbose = verboseHeader)
    
    #Create Group Matrix
    .messageDiffTime("Creating Cluster Matrix on the total Group Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    groupList <- SimpleList(split(rownames(outLSI$matSVD), clusters))
    groupFeatures <- totalAcc[sort(head(order(totalAcc$value, decreasing = TRUE), totalFeatures)),]
    groupMat <- .getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = groupFeatures, 
      threads = threads,
      groupList = SimpleList(split(rownames(outLSI$matSVD), clusters)),
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
      cellNames = cellNames, 
      sampleNames = getCellColData(ArchRProj)$Sample,  
      dimsToUse = dimsToUse, 
      binarize = binarize, 
      sampleCells = sampleCells,
      threads = threads,
      useIndex = FALSE,
      tstart = tstart
      )
    outLSI$LSIFeatures <- variableFeatures

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
      parClust <- lapply(clusterParams, function(x){
        if(length(x) > 1){
          return(x[[j]])
        }else{
          return(x[[1]])
        }
      })
      parClust$input <- outLSI$matSVD
      parClust$sampleCells <- sampleCells
      parClust$verbose <- verboseAll
      clusters <- do.call(IdentifyClusters, parClust)

      #Save Output
      if(saveIterations){
        .messageDiffTime("Saving LSI Iteration", tstart, addHeader = verboseAll, verbose = verboseHeader)
        outj <- SimpleList(LSI = outLSI, clusters = clusters, params = parClust[-length(parClust)])
        saveRDS(outj, file.path(outDir, paste0("Save-LSI-Iteration-",j,".rds")))
      }

    }

  }

  #Organize Output
  .messageDiffTime("Finished Running IterativeLSI", tstart, addHeader = verboseAll, verbose = verboseHeader)
  ArchRProj@reducedDims[[reducedDimsOut]] <- outLSI

  return(ArchRProj)

}

.LSIPartialMatrix <- function(
  ArrowFiles, 
  featureDF, 
  cellNames, 
  sampleNames, 
  dimsToUse, 
  binarize = TRUE, 
  sampleCells = 5000, 
  threads = 1, 
  useIndex = FALSE, 
  tstart = NULL, 
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
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
      cellNames = cellNames,
      doSampleCells = FALSE,
      threads = threads,
      verbose = verboseAll
    )

    #Compute LSI
    .messageDiffTime("Running LSI on the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    outLSI <- computeLSI(mat, nDimensions = max(dimsToUse), binarize = binarize, verbose = verboseAll, tstart = tstart)
  
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
    outLSI <- computeLSI(out$mat, nDimensions = max(dimsToUse), binarize = binarize, verbose = verboseAll, tstart = tstart)
    tmpMatFiles <- out[[2]]
    rm(out)
    gc()

    #Read In Matrices and Project into Manifold
    .messageDiffTime("Projecting Matrices with the Top Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
    pLSI <- lapply(seq_along(tmpMatFiles), function(x){
      projectLSI(mat = readRDS(tmpMatFiles[x]), LSI = outLSI, verbose = FALSE, tstart = tstart)
    }) %>% Reduce("rbind", .)

    #Remove Temporary Matrices
    rmf <- file.remove(tmpMatFiles)

    #Set To LSI the SVD Matrices
    outLSI$exlcude <- cellNames[which(cellNames %ni% rownames(pLSI))]
    outLSI$matSVD <- as.matrix(pLSI[cellNames[which(cellNames %in% rownames(pLSI))],])

  }

  return(outLSI)

}

#' Compute LSI
#' 
#' This function will compute a LSI transform (TF-IDF followed by SVD)
#'
#' @param mat sparseMatrix (dgcMatrix) for LSI
#' @param nDimensions number of LSI dimensions to compute
#' @param binarize binarize matrix prior to LSI
#' @param seed seed for analysis
#' @param verbose verbose
#' @param tstart time stamp to pass
#' @param ... additional args
#' @export
computeLSI <- function(mat, nDimensions = 50, binarize = TRUE, seed = 1, verbose = TRUE, tstart = NULL, ...){

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

    #IDF
    .messageDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose)
    idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

    #TF-IDF
    .messageDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
    mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
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
        nDimensions = nDimensions,
        date = Sys.Date(),
        seed = seed
        )

    rm(mat)
    gc()

    out
}

#' Project LSI
#' 
#' This function will compute a LSI Projection (TF-IDF followed by SVD projection)
#'
#' @param mat sparseMatrix (dgcMatrix) for LSI
#' @param LSI previous LSI transform to project into
#' @param returnModel return projection information
#' @param verbose verbose
#' @param tstart time stamp to pass
#' @param ... additional args
#' @export
projectLSI <- function(mat, LSI, returnModel = FALSE, verbose = TRUE, tstart = NULL, ...){   
    
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

    #IDF
    .messageDiffTime("Computing Inverse Document Frequency of initial Matrix", tstart, addHeader = FALSE, verbose = verbose)
    idf   <- as(log(1 + length(LSI$colSm) / LSI$rowSm), "sparseVector")

    #TF-IDF
    .messageDiffTime("Computing TF-IDF Transform", tstart, addHeader = FALSE, verbose = verbose)
    mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

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




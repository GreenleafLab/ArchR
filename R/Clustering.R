##########################################################################################
# Clustering Methods
##########################################################################################

#' Add cluster information to an ArchRProject
#' 
#' This function will identify clusters from a reduced dimensions object in an ArchRProject or from a supplied reduced dimensions matrix.
#' 
#' @param input Either (i) an `ArchRProject` object containing the dimensionality reduction matrix passed by `reducedDims`
#' or (ii) a dimensionality reduction matrix. This object will be used for cluster identification.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' Not required if input is a matrix.
#' @param name The column name of the cluster label column to be added to `cellColData` if `input` is an `ArchRProject` object.
#' @param sampleCells An integer specifying the number of cells to subsample and perform clustering on. The remaining cells
#' that were not subsampled will be assigned to the cluster of the nearest subsampled cell. This enables a decrease in run time
#' but can sacrifice granularity of clusters.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param method A string indicating the clustering method to be used. Supported methods are "Seurat" and "Scran".
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param knnAssign The number of nearest neighbors to be used during clustering for assignment of outliers (clusters with less than nOutlier cells).
#' @param nOutlier The minimum number of cells required for a group of cells to be called as a cluster. If a group of cells does not reach
#' this threshold, then the cells will be considered outliers and assigned to nearby clusters.
#' @param maxClusters The maximum number of clusters to be called. If the number exceeds this the clusters are merged unbiasedly using hclust and cutree. 
#' This is useful for contraining the cluster calls to be reasonable if they are converging on large numbers. Useful in iterativeLSI as well for initial iteration. Default is set to 25.
#' @param testBias A boolean value that indicates whether or not to test clusters for bias.
#' @param filterBias A boolean value indicates whether or not to filter clusters that are identified as biased.
#' @param biasClusters A numeric value between 0 and 1 indicating that clusters that are smaller than the specified proportion of total cells are
#' to be checked for bias. This should be set close to 0. We recommend a default of 0.01 which specifies clusters below 1 percent of the total cells.
#' @param biasCol The name of a column in `cellColData` that contains the numeric values used for testing bias enrichment.
#' @param biasVals A set of numeric values used for testing bias enrichment if `input` is not an `ArchRProject`.
#' @param biasQuantiles A vector of two numeric values, each between 0 and 1, that describes the lower and upper quantiles of the bias values to use
#' for computing bias enrichment statistics.
#' @param biasEnrich A numeric value that specifies the minimum enrichment of biased cells over the median of the permuted background sets.
#' @param biasProportion A numeric value between 0 and 1 that specifies the minimum proportion of biased cells in a cluster required to determine that the
#' cluster is biased during testing for bias-enriched clusters. 
#' @param biasPval A numeric value between 0 and 1 that specifies the p-value to use when testing for bias-enriched clusters.
#' @param nPerm An integer specifying the number of permutations to perform for testing bias-enriched clusters.
#' @param prefix A character string to be added before each cluster identity. For example, if "Cluster" then cluster results will be "Cluster1", "Cluster2" etc.
#' @param ArchRProj An `ArchRProject` object containing the dimensionality reduction matrix passed by `reducedDims`. This argument can also be supplied as `input`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param tstart A timestamp that is typically passed internally from another function (for ex. "IterativeLSI") to measure how long the clustering analysis
#' has been running relative to the start time when this process was initiated in another function. This argument is rarely manually specified.
#' @param force A boolean value that indicates whether or not to overwrite data in a given column when the value passed to `name` already
#' exists as a column name in `cellColData`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional arguments to be provided to Seurat::FindClusters or scran::buildSNNGraph (for example, knn = 50, jaccard = TRUE)
#' @export
#'
addClusters <- function(
  input = NULL, 
  reducedDims = "IterativeLSI",
  name = "Clusters",
  sampleCells = NULL,
  seed = 1, 
  method = "Seurat",
  dimsToUse = NULL,
  scaleDims = NULL, 
  corCutOff = 0.75,
  knnAssign = 10, 
  nOutlier = 5, 
  maxClusters = 25,
  testBias = TRUE,
  filterBias = FALSE,
  biasClusters = 0.01,
  biasCol = "nFrags",
  biasVals = NULL,
  biasQuantiles = c(0.05, 0.95),
  biasEnrich = 10,
  biasProportion = 0.5,
  biasPval = 0.05,
  nPerm = 500,
  prefix = "C",
  ArchRProj = NULL,
  verbose = TRUE,
  tstart = NULL,
  force = FALSE,
  logFile = createLogFile("addClusters"),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj", "null"))
  if(is(ArchRProj, "ArchRProject")){
    message("When running addClusters 'input' param should be used for 'ArchRProj'. Replacing 'input' param with user 'ArchRPRoj'...")
    input <- ArchRProj
    rm(ArchRProj)
    gc()
  }
  .validInput(input = input, name = "input", valid = c("ArchRProj", "matrix"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = method, name = "method", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = knnAssign, name = "knnAssign", valid = c("integer"))
  .validInput(input = nOutlier, name = "nOutlier", valid = c("integer"))
  .validInput(input = testBias, name = "testBias", valid = c("boolean"))
  .validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  .validInput(input = biasClusters, name = "biasClusters", valid = c("numeric"))
  .validInput(input = biasCol, name = "biasCol", valid = c("character"))
  .validInput(input = biasQuantiles, name = "biasQuantiles", valid = c("numeric"))
  .validInput(input = biasEnrich, name = "biasEnrich", valid = c("numeric"))
  .validInput(input = biasProportion, name = "biasProportion", valid = c("numeric"))
  .validInput(input = biasPval, name = "biasPval", valid = c("numeric"))
  .validInput(input = nPerm, name = "nPerm", valid = c("integer"))
  .validInput(input = prefix, name = "prefix", valid = c("character"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = tstart, name = "tstart", valid = c("timestamp","null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "addClusters Input-Parameters", logFile=logFile)

  if(is.null(tstart)){
      tstart <- Sys.time()
  }

  if(inherits(input, "ArchRProject")){
      #Check
      input <- addCellColData(
          ArchRProj = input, 
          data = rep(NA, nCells(input)), 
          name = name, 
          cells = getCellNames(input), 
          force = force
      )

      if(reducedDims %ni% names(input@reducedDims)){
          stop("Error reducedDims not available!")
      }

      matDR <- getReducedDims(
          ArchRProj = input, 
          reducedDims = reducedDims, 
          dimsToUse = dimsToUse, 
          corCutOff = corCutOff, 
          scaleDims = scaleDims
      )
  
  }else if(inherits(input, "matrix")){
      matDR <- input
  }else{
      stop("Input an ArchRProject or Cell by Reduced Dims Matrix!")
  }

  #Subset Matrix
  set.seed(seed)
  nr <- nrow(matDR)

  if(!is.null(sampleCells)){
    if(sampleCells < nrow(matDR)){
      .logDiffTime("Estimating Clusters by Sampling", tstart, verbose = verbose, logFile = logFile)
      estimatingClusters <- 1
      idx <- sample(seq_len(nrow(matDR)), sampleCells)
      matDRAll <- matDR
      matDR <- matDR[idx,,drop=FALSE]
    }else{
      estimatingClusters <- 0
    }
  }else{
    estimatingClusters <- 0
  }

  #################################################################################
  # Decide on which clustering setup to use
  #################################################################################
  if(grepl("seurat",tolower(method))){
  }else if(grepl("scran",tolower(method))){
  }else{
    stop("Clustering Method Not Recognized!")
  }

  clust <- tryCatch({

    if(grepl("seurat",tolower(method))){

      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clust <- .clustSeurat(mat = matDR, clustParams = clustParams, logFile = logFile)

    }else if(grepl("scran",tolower(method))){

      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clustParams$x <- t(matDR)
      clustParams$d <- ncol(matDR)
      clustParams$k <- ifelse(exists("...$k"), ...$k, 25)
      clust <- .clustScran(clustParams = clustParams, logFile = logFile)

    }

  }, error = function(e){

    errorList <- clustParams
    .logError(e, fn = "runClusters", info = "", errorList = errorList, logFile = logFile)

  })

  #################################################################################
  # If estimating clsuters we will assign to nearest neighbor cluster
  #################################################################################
  if(estimatingClusters == 1){
      
      .logDiffTime("Finding Nearest Clusters", tstart, verbose = verbose, logFile = logFile)
      knnAssigni <- as.matrix(.computeKNN(matDR, matDRAll[-idx,,drop=FALSE], knnAssign))
      clustUnique <- unique(clust)
      clustMatch <- match(clust, clustUnique)
      knnAssigni <- matrix(apply(knnAssigni, 2, function(x) clustMatch[x]), ncol = knnAssign)

      .logDiffTime("Assigning Nearest Clusters", tstart, verbose = verbose, logFile = logFile)
      clustAssign <- lapply(seq_along(clustUnique), function(x){
          rowSums(knnAssigni == x)
      }) %>% Reduce("cbind", .) %>% apply(., 1, which.max)
      clustOld <- clust
      clust <- rep(NA, nr)
      clust[idx] <- clustOld
      clust[-idx] <- clustUnique[clustAssign]
      matDR <- matDRAll
      remove(matDRAll)
      gc()

  }

  #################################################################################
  # Testing Bias
  #################################################################################
  if(testBias){
    if(inherits(input, "ArchRProject")){
      if(is.null(biasVals)){
        biasDF <- getCellColData(input, select = biasCol)
      }else{
        biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
      }
    }else{
      if(!is.null(biasVals)){
        biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
      }else{
        message("No biasVals for testing bias continuing without bias detection")
        testBias <- FALSE
      }
    }
  }

  if(testBias){
    clust <- tryCatch({
      biasDF$Q <- .getQuantiles(biasDF[,1])
      tabClust <- table(clust)
      tabClustP <- tabClust / sum(tabClust)
      idxTest <- which(tabClustP < biasClusters)
      names(clust) <- rownames(matDR)
      if(length(idxTest) > 0){
        .logDiffTime("Testing Biased Clusters", tstart, verbose = verbose, logFile = logFile)
        testDF <- lapply(seq_along(idxTest), function(i){
          clustTesti <- names(tabClustP)[idxTest[i]]
          biasQ <- biasDF[names(clust)[which(clust == clustTesti)], 2]
          biasBgd <- matrix(
            sample(
              x = biasDF[names(clust)[which(clust != clustTesti)], 2],
              size = nPerm * length(biasQ),
              replace = if(nPerm * length(biasQ) > nrow(biasDF[names(clust)[which(clust != clustTesti)], ])) TRUE else FALSE
            ), 
            nrow = length(biasQ), 
            ncol = nPerm
          )
          n1 <- colSums(biasBgd >= max(biasQuantiles))
          n2 <- colSums(biasBgd <= min(biasQuantiles))
          pval1 <- max(sum(sum(biasQ >= max(biasQuantiles)) < n1) * 2, 1) / length(n1)
          pval2 <- max(sum(sum(biasQ <= min(biasQuantiles)) < n2) * 2, 1) / length(n2)
          enrich1 <- sum(biasQ >= max(biasQuantiles)) / max(median(n1), 1)
          enrich2 <- sum(biasQ <= min(biasQuantiles)) / max(median(n2), 1)
          per1 <- sum(biasQ >= max(biasQuantiles)) / length(biasQ)
          per2 <- sum(biasQ <= min(biasQuantiles)) / length(biasQ)
          if(enrich1 > enrich2){
            enrichClust <- enrich1
            enrichPval <- min(pval1, 1)
            enrichPer <- per1
          }else{
            enrichClust <- enrich2
            enrichPval <- min(pval2, 1)
            enrichPer <- per2
          }
          DataFrame(Cluster = clustTesti, enrichClust = enrichClust, enrichPval = enrichPval, enrichProportion = enrichPer)
        }) %>% Reduce("rbind", .)

        clustAssign <- testDF[which(testDF$enrichClust > biasEnrich & testDF$enrichProportion > biasProportion & testDF$enrichPval <= biasPval),1]
        if(length(clustAssign) > 0){
          if(filterBias){
            .logDiffTime(sprintf("Assigning Biased Clusters (n = %s) to Neighbors", length(clustAssign)), tstart, verbose = verbose, logFile = logFile)
            for(i in seq_along(clustAssign)){
              clusti <- clustAssign[i]
              idxi <- which(clust==clusti)
              knni <- .computeKNN(matDR[-idxi,,drop=FALSE], matDR[idxi,,drop=FALSE], knnAssign)
              clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
              clust[idxi] <- clustf
            }
          }else{
            .logDiffTime(sprintf("Identified Biased Clusters (n = %s), set filterBias = TRUE to re-assign these cells: ", length(clustAssign)), tstart, verbose = verbose, logFile = logFile)
            message("Biased Clusters : ", appendLF = FALSE)
            for(i in seq_along(clustAssign)){
              message(clustAssign[i], " ", appendLF = FALSE)
            }
            message("")
          }
        }
      }
      clust
    }, error = function(e){

      errorList <- list(
        idxTest = if(exists("testDF", inherits = FALSE)) fragx else "Error with idxTest!",
        biasDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with biasDF!",
        testDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with testDF!",
        clustAssign = if(exists("idf", inherits = FALSE)) fragx else "Error with clustAssign!"
      )

      .logError(e, fn = "testBias", info = "", errorList = errorList, logFile = logFile)

    })

  }
  
  #################################################################################
  # Test if clusters are outliers identified as cells with fewer than nOutlier
  #################################################################################
  .logDiffTime("Testing Outlier Clusters", tstart, verbose = verbose, logFile = logFile)
  tabClust <- table(clust)
  clustAssign <- which(tabClust < nOutlier)
  if(length(clustAssign) > 0){
      .logDiffTime(sprintf("Assigning Outlier Clusters (n = %s, nOutlier < %s cells) to Neighbors", length(clustAssign), nOutlier), tstart, verbose = verbose, logFile = logFile)
      for(i in seq_along(clustAssign)){
          clusti <- names(clustAssign[i])
          idxi <- which(clust==clusti)
          knni <- .computeKNN(matDR[-idxi,], matDR[idxi,], knnAssign)
          clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
          clust[idxi] <- clustf
      }
  }

  #################################################################################
  # Merging if more than maxClusters
  #################################################################################
  if(!is.null(maxClusters)){
    if(length(unique(clust)) > maxClusters){
      .logDiffTime(sprintf("Identified more clusters than maxClusters allowed (n = %s). Merging clusters to maxClusters (n = %s).\nIf this is not desired set maxClusters = NULL!", length(clustAssign), maxClusters), tstart, verbose = verbose, logFile = logFile)
      meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      ct <- cutree(hc, maxClusters)
      clust <- mapLabels(
        labels = clust, 
        oldLabels = names(ct), 
        newLabels = paste0(prefix, ct)
      )
    }
  }
  
  #################################################################################
  # Renaming Clusters based on Proximity in Reduced Dimensions
  #################################################################################
  .logDiffTime(sprintf("Assigning Cluster Names to %s Clusters", length(unique(clust))), tstart, verbose = verbose, logFile = logFile)
  
  if(length(unique(clust)) > 1){

      meanDR <- t(.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      out <- mapLabels(
        labels = clust, 
        oldLabels = hc$labels[hc$order], 
        newLabels = paste0(prefix, seq_along(hc$labels))
      )

  }else{

      out <- rep(paste0(prefix, "1"), length(clust))

  }

  if(inherits(input, "ArchRProject")){
    input <- .suppressAll(addCellColData(
            input, 
            data = out, 
            name = name, 
            cells = rownames(matDR),
            force = TRUE
        ))
    .logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, logFile = logFile)
    return(input)
  }else if(!inherits(input, "ArchRProject")){
    .logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, logFile = logFile)
    return(out)
  }

}

#Simply a wrapper on Seurats FindClusters
.clustSeurat <- function(mat = NULL, clustParams = NULL, logFile = NULL){

  .requirePackage("Seurat", source = "cran")
  .logDiffTime("Running Seurats FindClusters (Stuart et al. Cell 2019)", clustParams$tstart, verbose=clustParams$verbose, logFile = logFile)

  tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
  colnames(tmp) <- rownames(mat)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))

  clustParams <- tryCatch({

    obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
    obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
    clustParams$object <- obj
    clustParams$reduction <- "pca"
    clustParams$dims <- seq_len(ncol(mat))

    obj <- suppressWarnings(do.call(Seurat::FindNeighbors, clustParams))
    clustParams$object <- obj
    clustParams

  }, error = function(e){

    errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
    .logError(e, fn = "FindNeighbors", info = "", errorList = errorList, logFile = logFile)

  })

  clust <- tryCatch({

    cS <- Matrix::colSums(obj@graphs$RNA_snn)

    if(cS[length(cS)] == 1){

        #Error Handling with Singletons
        idxSingles <- which(cS == 1)
        idxNonSingles <- which(cS != 1)

        rn <- rownames(mat) #original order
        mat <- mat[c(idxSingles, idxNonSingles), ,drop = FALSE]

        tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
        colnames(tmp) <- rownames(mat)
        rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))

        obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
        obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
        clustParams$object <- obj
        clustParams$reduction <- "pca"
        clustParams$dims <- seq_len(ncol(mat))

        obj <- .suppressAll(do.call(Seurat::FindNeighbors, clustParams))
        clustParams$object <- obj

        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))

        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)
        clust <- clust[rn]

    }else{

        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))

        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)

    }

    clust

  }, error = function(e){

    errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
    .logError(e, fn = "FindClusters", info = "", errorList = errorList, logFile = logFile)

  })

  clust

}

.clustScran <- function(clustParams = NULL, logFile = NULL){
  .requirePackage("scran", installInfo='BiocManager::install("scran")')
  .requirePackage("igraph", installInfo='install.packages("igraph")')
  #See Scran Vignette!
  tstart <- clustParams$tstart
  verbose <- clustParams$verbose
  clustParams$tstart <- NULL
  clustParams$verbose <- NULL
  .logDiffTime("Running Scran SNN Graph (Lun et al. F1000Res. 2016)", tstart, verbose=verbose, logFile = logFile)
  snn <- do.call(scran::buildSNNGraph, clustParams)
  .logDiffTime("Identifying Clusters (Lun et al. F1000Res. 2016)", tstart, verbose=verbose, logFile = logFile)
  cluster <- igraph::cluster_walktrap(snn)$membership
  paste0("Cluster", cluster)
}








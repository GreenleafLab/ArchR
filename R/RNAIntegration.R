####################################################################
# Gene Integration Matrix Methods
####################################################################

#' Add a GeneIntegrationMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, for each sample, will independently integrate with a scRNA experiment, compute matched scRNA profiles and
#' then store this in each samples ArrowFile.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param seRNA A scRNA-seq `SummarizedExperiment` or `SeuratObject` to integrated with the scATAC-seq data.
#' @param groupBy QQQ UNCLEAR WHAT PREDICTION SCORES ARE The column name in either `colData` or `metadata` of seRNA to use to compute prediction scores.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`. 
#' @param sampleList QQQ THIS COULD BE BETTER DESCRIBED A named list of sampleNames from `ArchRProj` containing cell names in `seRNA` that are used for constraining data integration.
#' For example when integrating hematopoiesis it is useful to perform integration separately for cells derived from different biological contexts
#' such as "CD34", "BMMC" and "PBMC" single cells. Performing this integration separately not only provides better linking but
#' also increases speed of processing.
#' @param sampleRNA The number of cells from `seRNA` to be used for integration. This is useful for lowering the memory requirements and
#' processing time of `Seurat::FindTransferAnchors()` when integrating datasets.
#' @param nGenes The number of variable genes determined by `Seurat::FindVariableGenes()` to use for integration.
#' @param useImputation A boolean value indicating whether to use imputation for creating the Gene Score Matrix prior to integration.
#' @param reduction The Seurat reduction method to use for integrating modalities. See `Seurat::FindTransferAnchors()` for possible reduction methods.
#' @param addToArrow A boolean value indicating whether to add the log2-normalized transcript counts from the integrated matched RNA to the ArrowFiles.
#' @param threads The number of threads to be used for parallel computing.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exists in the given `input`.
#' @export
addGeneIntegrationMatrix <- function(
  ArchRProj = NULL,
  seRNA = NULL,
  groupBy = NULL,
  reducedDims = "IterativeLSI",
  sampleList = SimpleList(),
  sampleRNA = 5000,
  nGenes = 2000,
  useImputation = TRUE,
  reduction = "cca",
  addToArrow = TRUE,
  threads = getArchRThreads(),
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  force = FALSE,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  tstart <- Sys.time()
  .messageDiffTime("Running Seurat's Integration Stuart et al 2019", tstart, verbose = verboseHeader)

  .requirePackage("Seurat", source = "cran")

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = sampleList, name = "sampleList", valid = c("list"))
  .validInput(input = sampleRNA, name = "sampleRNA", valid = c("integer"))
  .validInput(input = nGenes, name = "nGenes", valid = c("integer"))
  .validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  .validInput(input = reduction, name = "reduction", valid = c("character"))
  .validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
  .validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))

  ArrowFiles <- getArrowFiles(ArchRProj)

  if(!is.null(sampleRNA) & length(sampleList) > 0){
    sampleList2 <- sampleList
    sampleList <- lapply(sampleList, function(x){
      if(length(x) > sampleRNA){
        sample(x, sampleRNA)
      }else{
        x
      }
    }) %>% SimpleList
    names(sampleList) <- names(sampleList2)
    rm(sampleList2)
  }

  idx <- which(names(ArrowFiles) %ni% names(sampleList))
  if(length(idx) > 1){
    for(x in seq_along(idx)){
      if(ncol(seRNA) > sampleRNA){
        sampleList[[names(ArrowFiles)[x]]] <- sample(colnames(seRNA), sampleRNA)
      }else{
        sampleList[[names(ArrowFiles)[x]]] <- colnames(seRNA)
      }
    }
  }
  seRNA <- seRNA[,unique(unlist(sampleList))]

  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupBy %ni% colnames(colData(seRNA))){
      stop("groupBy not in colData of seRNA")
    }
    seuratRNA$Group <- colData(seRNA)[, groupBy, drop = TRUE]
    seuratRNA <- NormalizeData(object = seuratRNA, verbose = verboseAll)
    rm(seRNA)
  }else{
    if(groupBy %ni% colnames(seRNA@meta.data)){
      stop("groupBy not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- seRNA@meta.data[,groupBy]
    seuratRNA <- NormalizeData(object = seuratRNA, verbose = verboseAll)
    rm(seRNA)
  }
  gc()

  #Clean Project For Parallel
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()

  dfAll <- .safelapply(seq_along(ArrowFiles), function(i){

    ArrowFiles <- getArrowFiles(subProj)
    ArrowFile <- ArrowFiles[i]
    sampleName <- .sampleName(ArrowFile)
    matrixName <- "GeneIntegrationMatrix"
    prefix <- sprintf("%s (%s of %s)", sampleName, i, length(ArrowFiles))
    .messageDiffTime(paste0(prefix, " Running Seurat's Integration Stuart et al 2019"), tstart, verbose = verboseHeader)

    #Subset
    subProj@cellColData <- subProj@cellColData[which(subProj$Sample %in% sampleName), ]

    #Determine Number of Cells
    cellsUse <- sampleList[[sampleName]]
    seuratRNA <- seuratRNA[,cellsUse]

    #Determine Features
    geneDF <- .getFeatureDF(ArrowFile, "GeneScoreMatrix")
    geneDF <- geneDF[geneDF$name %in% rownames(seuratRNA), , drop = FALSE]

    #Subet RNA
    seuratRNA <- seuratRNA[rownames(seuratRNA) %in% geneDF$name, ]

    #1. Create Seurat RNA and Normalize
    .messageDiffTime(paste0(prefix, " Identifying Variable Genes for Integration!"), tstart, verbose = verboseHeader)
    seuratRNA <- FindVariableFeatures(object = seuratRNA, nfeatures = nGenes, verbose = verboseAll)
    seuratRNA <- ScaleData(object = seuratRNA, verbose = verboseAll)
    genesUse <- VariableFeatures(object = seuratRNA)

    #2. Get Gene Score Matrix and Create Seurat ATAC
    .messageDiffTime(paste0(prefix, " Creating Seurat ATAC Gene Score Object"), tstart, verbose = verboseHeader)
    mat <- .getPartialMatrix(
      getArrowFiles(subProj)[sampleName], 
      featureDF = geneDF[geneDF$name %in% genesUse,], 
      threads = threads,
      cellNames = subProj$cellNames,
      useMatrix = "GeneScoreMatrix",
      verbose = verboseAll
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]

    # Log-Normalize (its already scaled internally in ArrowFiles)
    mat@x <- log(mat@x + 1) #use natural log
    if(useImputation){
      .messageDiffTime(paste0(prefix, " Imputing Gene Scores"), tstart, verbose = verboseHeader)
      subProj <- .suppressAll(addImputeWeights(subProj))
      mat <- getImputeWeights(subProj)[[1]][colnames(mat), colnames(mat)] %*% Matrix::t(mat)
      mat <- Matrix::t(mat)
    }
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    #Clean Memory
    rm(mat)
    gc()
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = verboseAll)

    #3. Transfer Anchors   
    .messageDiffTime(paste0(prefix, " Seurat Finding Transfer Anchors"), tstart, verbose = verboseHeader)
    transferAnchors <- suppressWarnings(Seurat::FindTransferAnchors(
      reference = seuratRNA, 
      query = seuratATAC, 
      reduction = reduction, 
      features = genesUse,
      ...
    ))

    #Clean Memory
    gc()

    #4. Transfer Data
    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = getReducedDims(
        ArchRProj = subProj, 
        reducedDims = reducedDims
      )[colnames(seuratATAC),,drop=FALSE], 
      key = "LSI_", 
      assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- verboseAll
      
    #Get RNA Labels
    .messageDiffTime(paste0(prefix, " Seurat Transfer Labels"), tstart, verbose = verboseHeader)
    transferParams$refdata <- seuratRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)

    if(addToArrow){
      #Get Matched RNA
      .messageDiffTime(paste0(prefix, " Seurat Transfer Matrix"), tstart, verbose = verboseHeader)
      transferParams$refdata <- GetAssayData(seuratRNA, assay = "RNA", slot = "data")
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }

    #Match results
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC), 
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id
    )
    rownames(matchDF) <- matchDF$cellNames

    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()

    if(addToArrow){

      #5. Organize info for ArchR Arrow

      #Re-Index RNA
      geneDF2 <- geneDF[geneDF$name %in% rownames(matchedRNA), , drop = FALSE]
      splitGeneDF <- S4Vectors::split(geneDF2, geneDF2$seqnames)
      featureDF <- lapply(splitGeneDF, function(x){
        x$idx <- seq_len(nrow(x))
        return(x)
      }) %>% Reduce("rbind", .)
      matchedRNA <- matchedRNA[paste0(featureDF$name), ,drop=FALSE]
      stopifnot(identical(paste0(featureDF$name), paste0(rownames(matchedRNA))))

      #Params
      dfParams <- data.frame(
          reduction = reduction
      )

      ######################################
      # Initialize SP Mat Group
      ######################################
      o <- .initializeMat(
        ArrowFile = ArrowFile,
        Group = matrixName,
        Class = "double",
        Units = "Log2NormCounts",
        cellNames = colnames(matchedRNA),
        params = dfParams,
        featureDF = featureDF,
        force = force
      )

      o <- h5write(
        obj = matchDF[colnames(matchedRNA), "predictionScore"], 
        file = ArrowFile, 
        name = paste0(matrixName, "/Info/predictionScore")
      )

      o <- h5write(
        obj = matchDF[colnames(matchedRNA), "predictedGroup"], 
        file = ArrowFile, 
        name = paste0(matrixName, "/Info/predictedGroup")
      )

      #Clean Memory
      gc()

      #Time to write into Arrow File!
      allChr <- unique(featureDF$seqnames)

      .messageDiffTime(sprintf("%s Adding GeneIntegrationMatrix to ArrowFile!", prefix), tstart, verbose = verboseHeader)

      for(z in seq_along(allChr)){

        chrz <- allChr[z]
        .messageDiffTime(sprintf("Adding GeneIntegrationMatrix to %s for Chr (%s of %s)!", sampleName, z, length(allChr)), tstart, verbose = verboseAll)

        idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
        matz <- matchedRNA[idz, ,drop=FALSE]
        matz@x <- round(matz@x, 3)
        matz <- Matrix::drop0(matz)
        stopifnot(identical(paste0(featureDF$name[idz]), paste0(rownames(matz))))

        #Write sparseMatrix to Arrow File!
        o <- .addMatToArrow(
          mat = matz, 
          ArrowFile = ArrowFile, 
          Group = paste0(matrixName, "/", chrz), 
          binarize = FALSE,
          addColSums = TRUE,
          addRowSums = TRUE,
          addRowVars = TRUE
        )

        #Clean Memory
        rm(matz)

        if(z %% 3 == 0 | z == length(allChr)){
          gc()
        }

      }

    }

    matchDF

  }, threads = 1) %>% Reduce("rbind", .)

  .messageDiffTime("Completed Integration with RNA Matrix", tstart, verbose = verboseHeader)

  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedGroup,
    name = nameGroup,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictionScore,
    name = nameScore,
    force = TRUE
  )

  return(ArchRProj)

}


#' Add Peak2GeneLinks to an ArchRProject JJJ
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#' QQQ IM NOT CLEAR ON HOW THIS IS DIFFERENT FROM addCoAccessibility? THE FUNCTION DESCRIPTION IS IDENTICAL
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims QQQ DOUBLE CHECK A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param knnMethod The method to be used for k-nearest neighbor computations. Options are "nabor", "RANN", and "FNN" and the corresponding package is required.
#' @param threads The number of threads to be used for parallel computing.
#' @export
addPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 150000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  seed = 1, 
  knnMethod = NULL,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = knnMethod, name = "knnMethod", valid = c("character", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  ArchRProj = projHeme5
  reducedDims = "IterativeLSI"
  dimsToUse = 1:30
  scaleDims = NULL
  corCutOff = 0.75
  k = 100
  knnIteration = 500
  overlapCutoff = 0.8
  maxDist = 150000
  scaleTo = 10^4
  log2Norm = TRUE
  predictionCutoff = 0.4
  seed = 1
  knnMethod = NULL
  threads = getArchRThreads()

  ArrowFiles <- getArrowFiles(ArchRProj)

  dfAll <- lapply(seq_along(ArrowFiles), function(x){
    DataFrame(
      cellNames = paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], "GeneIntegrationMatrix/Info/CellNames")),
      predictionScore = h5read(ArrowFiles[x], "GeneIntegrationMatrix/Info/predictionScore")
    )
  }) %>% Reduce("rbind", .)

  keep <- sum(dfAll[,2] > predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  tstart <- Sys.time()

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .messageDiffTime("Computing KNN", tstart)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k, method = knnMethod)

  #Determin Overlap
  .messageDiffTime("Identifying Non-Overlapping KNN pairs", tstart)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .messageDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), tstart)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
  stopifnot(identical(chri,chrj))

  #Create Ranges
  peakSummits <- resize(peakSet, 1, "center")
  geneWindows <- resize(geneStart, maxDist, "center")

  #Create Pairwise Things to Test
  o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
  o <- o[o[,1] != o[,2],]
  o$seqnames <- seqnames(peakSet)[o[,1]]
  o$idx1 <- peakSet$idx[o[,1]]
  o$idx2 <- peakSet$idx[o[,2]]
  o$correlation <- NA

  #Peak Matrix ColSums
  cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

  for(x in seq_along(chri)){
  
    .messageDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), tstart)

    #Features
    featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]),]
    featureDF$seqnames <- chri[x]

    #Group Matrix
    groupMat <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = featureDF, 
      groupList = knnObj, 
      useMatrix = "PeakMatrix",
      threads = threads,
      verbose = FALSE
    )
    
    #Scale
    groupMat <- t(t(groupMat) / gS) * scaleTo

    if(log2Norm){
      groupMat <- log2(groupMat + 1)
    }

    #Correlations
    idx <- BiocGenerics::which(o$seqnames==chri[x])
    o[idx,]$correlation <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = groupMat, Y = groupMat)

  }
  
  o$idx1 <- NULL
  o$idx2 <- NULL
  o <- o[!is.na(o$correlation),]
  mcols(peakSet) <- NULL
  o@metadata$peakSet <- peakSet

  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  ArchRProj

}




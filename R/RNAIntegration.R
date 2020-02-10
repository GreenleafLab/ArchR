####################################################################
# Gene Integration Matrix Methods
####################################################################

#' Add GeneIntegrationMatrix to ArrowFiles or an ArchRProject JJJ
#' 
#' This function, for each sample, will independently compute counts for each tile
#' per cell and then infer gene activity scores.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param seRNA A scRNA `SummarizedExperiment` or `SeuratObject` to be integrated with.
#' @param groupBy A column name in either `colData` or `metadata` of seRNA to compute prediction scores.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`. 
#' @param sampleList A named list of sampleNames from `ArchRProj` containing cell names in `seRNA` that are used for constraining data integration.
#' For example when integrating hematopoiesis it is useful to partition integration for "CD34", "BMMC" and "PBMC" single cells.
#' @param sampleRNA Number of cells from seRNA to be used for integration. This is useful for lowering memory and speeding
#' up results of `Seurat::FindTransferAnchors` when integrating datasets.
#' @param nGenes Number of variable genes determined by `Seurat::FindVariableGenes` to use for integration.
#' @param useImputation Use imputation for Gene Score Matrix prior to integration.
#' @param reduction Seurat reduction method for integrating modalities. See `Seurat::FindTransferAnchors` for reduction methods.
#' @param addToArrow Add integrated matched RNA to ArrowFiles afterwards.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param subThreading A boolean determining whether possible use threads within each multi-threaded subprocess if greater than the number of input samples.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exist in the given `input`.
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
  transferParams = list(),
  force = FALSE,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  tstart <- Sys.time()
  .messageDiffTime("Running Seurat's Integration Stuart et al 2019", tstart, verbose = verboseHeader)

  .requirePackage("Seurat")

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

  if(!is.null(sampleRNA)){
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
    for(x in seq_along(sampleList)){
      if(ncol(seRNA) > sampleRNA){
        sampleList[[names(ArrowFiles)]] <- sample(colnames(seRNA), sampleRNA)
      }else{
        sampleList[[names(ArrowFiles)]] <- colnames(seRNA)
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
    subProj@cellColData <- subProj@cellColData[BiocGenerics::which(subProj$Sample == sampleName), ]

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
    name = "predictedGroup",
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictionScore,
    name = "predictionScore",
    force = TRUE
  )

  return(ArchRProj)

}










































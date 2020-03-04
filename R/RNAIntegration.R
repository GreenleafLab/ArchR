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
  reducedDims = "IterativeLSI",
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  seRNA = NULL,
  groupATAC = NULL,
  groupRNA = NULL,
  groupList = NULL,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  nGenes = 2000,
  useImputation = TRUE,
  reduction = "cca",
  addToArrow = TRUE,
  scaleTo = 10000,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  force = FALSE,
  threads = getArchRThreads(),
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  tstart <- Sys.time()
  .messageDiffTime("Running Seurat's Integration Stuart* et al 2019", tstart, verbose = verboseHeader)

  .requirePackage("Seurat", source = "cran")

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  #JJJ
  .validInput(input = nGenes, name = "nGenes", valid = c("integer"))
  .validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  .validInput(input = reduction, name = "reduction", valid = c("character"))
  .validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
  .validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))

  if(is.null(groupList)){ #If null use all cells (blocking will still occur)
    groupList <- SimpleList()
    groupList[[1]] <- SimpleList(
        ATAC = ArchRProj$cellNames,
        RNA = colnames(seRNA)
    )
  }

  #########################################################################################
  # 1. Check All ATAC is Accounted For!
  #########################################################################################
  .messageDiffTime("Checking ATAC Input", tstart, verbose = verboseHeader)

  if(!is.null(groupATAC)){
    dfATAC <- getCellColData(ArchRProj = ArchRProj, select = groupATAC, drop = FALSE)
  }
  nCell <- rep(0, length(ArchRProj$cellNames))
  names(nCell) <- ArchRProj$cellNames

  groupList <- lapply(seq_along(groupList), function(x){

    ATAC <- groupList[[x]]$ATAC

    if(!is.null(groupATAC)){

      if(any(ATAC %in% dfATAC[,1])){
        idx <- which(ATAC %in% dfATAC[,1])
        ATAC2 <- rownames(dfATAC)[which(dfATAC[,1] %in% ATAC[idx])]
        if(length(idx) == length(ATAC)){
          ATAC <- ATAC2
        }else{
          ATAC <- c(ATAC[-idx], ATAC2)
        }
      }
      
    }

    SimpleList(ATAC = ATAC, RNA = groupList[[x]]$RNA)

  }) %>% SimpleList

  for(i in seq_along(groupList)){
    nCell[groupList[[i]]$ATAC] <- nCell[groupList[[i]]$ATAC] + 1
  }

  if(!all(nCell == 1)){
    stop("Missing ", length(which(nCell == 0)), " Overlapping ", length(which(nCell > 1))," cells from ArchRProj in groupList!")
  }

  #########################################################################################
  # 2. Check All RNA is a Cell Name 
  #########################################################################################
  .messageDiffTime("Checking RNA Input", tstart, verbose = verboseHeader)

  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupBy %ni% colnames(colData(seRNA))){
      stop("groupBy not in colData of seRNA")
    }
    seuratRNA$Group <- colData(seRNA)[, groupBy, drop = TRUE]
    rm(seRNA)
  }else{
    if(groupBy %ni% colnames(seRNA@meta.data)){
      stop("groupBy not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- seRNA@meta.data[,groupBy]
    rm(seRNA)
  }
  gc()

  if(!is.null(groupRNA)){
    dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  }

  groupList <- lapply(seq_along(groupList), function(x){

    RNA <- groupList[[x]]$RNA

    if(!is.null(groupRNA)){

      if(any(RNA %in% dfRNA[,1])){
        idx <- which(RNA %in% dfRNA[,1])
        RNA2 <- rownames(dfRNA)[which(dfRNA[,1] %in% RNA[idx])]
        if(length(idx) == length(RNA)){
          RNA <- RNA2
        }else{
          RNA <- c(RNA[-idx], RNA2)
        }
      }
      
    }

    SimpleList(ATAC = groupList[[x]]$ATAC, RNA = RNA)

  }) %>% SimpleList

  cellRNA <- unlist(lapply(groupList, function(x) x$RNA))
  if(!all(cellRNA %in% colnames(seuratRNA))){
    stop("Found cells for RNA not in colnames(seRNA)! Please retry your input!")
  }

  seuratRNA <- seuratRNA[, unique(cellRNA)]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = verboseAll)

  #########################################################################################
  # 3. Create Integration Blocks
  #########################################################################################
  .messageDiffTime("Creating Integration Blocks", tstart, verbose = verboseHeader)

  blockList <- SimpleList()

  for(i in seq_along(groupList)){

    gLi <- groupList[[i]]

    #######################################
    # ATAC
    #######################################

    if(length(gLi$ATAC) > sampleCellsATAC){
      
      if(!is.null(embeddingATAC)){
        probATAC <- .getDensity(embeddingATAC[gLi$ATAC,1], embeddingATAC[gLi$ATAC,2])$density
        probATAC <- probATAC / max(probATAC)
        cellsATAC <- gLi$ATAC[order(probATAC, decreasing = TRUE)]
      }else{
        cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      }

      cutoffs <- lapply(seq_len(1000), function(x) length(gLi$ATAC) / x) %>% unlist
      blockSize <- ceiling(min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1, length(gLi$ATAC)))
      
      #Density Based Blocking
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)

      blocks <- lapply(seq_len(nBlocks), function(x){
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }) %>% SimpleList
    
    }else{
    
      blocks <- list(gLi$ATAC)
    }

    #######################################
    # RNA
    #######################################

    if(!is.null(embeddingRNA)){
      probRNA <- .getDensity(embeddingRNA[gLi$RNA,1], embeddingRNA[gLi$RNA,2])$density
      probRNA <- probRNA / max(probRNA)
    }else{
      probRNA <- rep(1, length(gLi$RNA))
    }

    blockListi <- lapply(seq_along(blocks), function(x){

      SimpleList(
        ATAC = blocks[[x]],
        RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)) , prob = probRNA)
      )

    }) %>% SimpleList

    blockList <- c(blockList, blockListi)

  }
  rm(groupList)

  #########################################################################################
  # 4. Begin Integration
  #########################################################################################
  .messageDiffTime("Prepping Interation Data", tstart, verbose = verboseHeader)

  #Clean Project For Parallel
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()

  #Gene Score Info
  geneDF <- .getFeatureDF(getArrowFiles(subProj), useMatrix)
  geneDF <- geneDF[geneDF$name %in% rownames(seuratRNA), , drop = FALSE]

  #Re-Index RNA
  splitGeneDF <- S4Vectors::split(geneDF, geneDF$seqnames)
  featureDF <- lapply(splitGeneDF, function(x){
    x$idx <- seq_len(nrow(x))
    return(x)
  }) %>% Reduce("rbind", .)
  dfParams <- data.frame(
      reduction = reduction
  )
  allChr <- unique(featureDF$seqnames)

  #Temp File Prefix
  tmpFile <- ArchR:::.tempfile()
  o <- suppressWarnings(file.remove(paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")))

  if(threads > 1){
    h5disableFileLocking()
  }

  addToArrow <- TRUE

  tstart <- Sys.time()

  threads2 <- max(ceiling(threads/2), 1)

  .messageDiffTime(paste0("Computing Integration in ", length(blockList), " Integration Blocks!"), tstart, verbose = verboseHeader)

  #Integration
  dfAll <- .safelapply(seq_along(blockList), function(i){

    .messageDiffTime(sprintf("Computing Integration : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
    blocki <- blockList[[i]]

    #Subset ATAC
    subProj@cellColData <- subProj@cellColData[blocki$ATAC, ]
    subProj@sampleColData <- subProj@sampleColData[unique(subProj$Sample),,drop=FALSE]

    #Subset RNA
    subRNA <- seuratRNA[, blocki$RNA]

    #Subet RNA
    subRNA <- subRNA[rownames(subRNA) %in% geneDF$name, ]

    ##############################################################################################
    #1. Create Seurat RNA and Normalize
    ##############################################################################################
    .messageDiffTime(sprintf("Identifying Variable Genes for Integration : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
    subRNA <- FindVariableFeatures(object = subRNA, nfeatures = nGenes, verbose = verboseAll)
    subRNA <- ScaleData(object = subRNA, verbose = verboseAll)
    genesUse <- VariableFeatures(object = subRNA)

    ##############################################################################################
    #2. Get Gene Score Matrix and Create Seurat ATAC
    ##############################################################################################
    .messageDiffTime(sprintf("Creating ATAC Gene Score Matrix : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
    mat <- .getPartialMatrix(
      getArrowFiles(subProj), 
      featureDF = geneDF[geneDF$name %in% genesUse,], 
      threads = 1,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = verboseAll
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]

    #Impute Matrix (its already scaled internally in ArrowFiles)
    if(useImputation){
      .messageDiffTime(sprintf("Imputing ATAC Gene Score Matrix : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$threads <- 1
      subProj <- do.call(addImputeWeights, imputeParams)
      mat <- imputeMatrix(mat = mat, imputeWeights = getImputeWeights(subProj), verbose = verboseAll)
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]]))) #Clean Up Space
    }

    #Log-Normalize 
    mat <- log(mat + 1) #use natural log
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    
    #Clean Memory
    rm(mat)
    
    #Set Default Assay
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = verboseAll)

    ##############################################################################################
    #3. Transfer Anchors  
    ############################################################################################## 
    .messageDiffTime(sprintf("Seurat FindTransferAnchors : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
    gc()
    verbose <- verboseAll #Force This to Be a Parameter Input
    transferAnchors <- suppressWarnings(Seurat::FindTransferAnchors(
      reference = subRNA, 
      query = seuratATAC, 
      reduction = reduction, 
      features = genesUse,
      verbose = verbose
      ...
    ))

    ##############################################################################################
    #4. Transfer Data
    ##############################################################################################
    .messageDiffTime(sprintf("Seurat TransferData Labels : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
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
    
    #Group
    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)

    #RNA Names
    transferParams$refdata <- colnames(subRNA)
    rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[,1]

    if(addToArrow){
      .messageDiffTime(sprintf("Seurat TransferData Matrix : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)
      transferParams$refdata <- GetAssayData(subRNA, assay = "RNA", slot = "data")
      gc()
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }

    #Match results
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC), 
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id,
      predictedCell = rnaLabels2
    )
    rownames(matchDF) <- matchDF$cellNames

    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()

    ##############################################################################################
    #5. Add To Temp Hdf5
    ##############################################################################################

    if(addToArrow){

      .messageDiffTime(sprintf("Transferring Paired RNA to Temp File : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)

      #Quickly Write to A Temp Hdf5 File Split By Sample to Then Enable Writing to Each Arrow File

      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- ArchR:::.safeSubset( #If Rownames disappeared this will catch that!
        mat = matchedRNA, 
        subsetRows = paste0(featureDF$name), 
        subsetCols = matchDF$cellNames
      )

      for(z in seq_along(uniqueSamples)){

        mat <- matchedRNA[, which(sampleNames == uniqueSamples[z]), drop = FALSE]
        Group <- uniqueSamples[z]

        o <- tryCatch({h5delete(tmpFilei, paste0(Group))}, error = function(x){})
        o <- h5createGroup(tmpFilei, paste0(Group))

        #Convert Columns to Rle
        j <- Rle(findInterval(seq(mat@x)-1, mat@p[-1]) + 1)

        #Info
        lengthRle <- length(j@lengths)
        lengthI <- length(mat@i)

        #Create Data Set
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/i"), storage.mode = "integer", 
          dims = c(lengthI, 1), level = 0))

        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jLengths"), storage.mode = "integer", 
          dims = c(lengthRle, 1), level = 0))

        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jValues"), storage.mode = "integer", 
          dims = c(lengthRle, 1), level = 0))

        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group, "/x"), storage.mode = "double", 
          dims = c(lengthI, 1), level = 0))

        #Write Data Set
        o <- .suppressAll(h5write(obj = mat@i + 1, file = tmpFilei, name = paste0(Group,"/i")))
        o <- .suppressAll(h5write(obj = j@lengths, file = tmpFilei, name = paste0(Group,"/jLengths")))
        o <- .suppressAll(h5write(obj = j@values, file = tmpFilei, name = paste0(Group,"/jValues")))
        o <- .suppressAll(h5write(obj = mat@x, file = tmpFilei, name = paste0(Group, "/x")))
        o <- .suppressAll(h5write(obj = colnames(mat), file = tmpFilei, name = paste0(Group, "/cellNames")))
        #Row Names is always the same

      }

      rm(matchedRNA, mat, j)

    }

    .messageDiffTime(sprintf("Completed Integration : Block (%s of %s)", i, length(blockList)), tstart, verbose = verboseHeader)

    gc()

    matchDF

  }, threads = threads2) %>% Reduce("rbind", .)

  ##############################################################################################
  #5. Read sub-matrices and store in ArrowFiles
  ##############################################################################################

  if(addToArrow){

    .messageDiffTime("Transferring Data to ArrowFiles", tstart, verbose = verboseHeader)

    integrationFiles <- paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")

    if(!all(file.exists(integrationFiles))){
      stop("Something went wrong with integration as not all temporary files containing integrated RNA exist!")
    }

    h5list <- .safelapply(seq_along(integrationFiles), function(x){
      h5ls(integrationFiles[x])
    }, threads = threads)

    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)

    o <- ArchR:::.safelapply(seq_along(allSamples), function(y){

      sample <- allSamples[y]

      prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))

      .messageDiffTime(sprintf("%s Getting GeneIntegrationMatrix From TempFiles!", prefix), tstart, verbose = verboseHeader)

      sampleIF <- lapply(seq_along(h5list), function(x){
        if(any(h5list[[x]]$group==paste0("/",sample))){
          integrationFiles[x]
        }else{
          NULL
        }
      }) %>% unlist

      sampleMat <- lapply(seq_along(sampleIF), function(x){

        cellNames <- .h5read(sampleIF[x], paste0(sample, "/cellNames"))

        mat <- sparseMatrix(
          i = .h5read(sampleIF[x], paste0(sample, "/i"))[,1], 
          j = as.vector(
            Rle(
              .h5read(sampleIF[x], paste0(sample, "/jValues"))[,1], 
              .h5read(sampleIF[x], paste0(sample, "/jLengths"))[,1]
            )
          ), 
          x = .h5read(sampleIF[x], paste0(sample, "/x"))[,1],
          dims = c(nrow(featureDF), length(cellNames))
        )
        colnames(mat) <- cellNames

        mat
        
      }) %>% Reduce("cbind", .)

      sampleMat@x <- exp(sampleMat@x) - 1 #Back To Counts
      sampleMat <- .normalizeCols(sampleMat, scaleTo = scaleTo) #Scale to 10,000
      sampleMat <- drop0(sampleMat) # Drop 0's
      rownames(sampleMat) <- paste0(featureDF$name)
      sampleMat <- sampleMat[,ArchRProj$cellNames[BiocGenerics::which(ArchRProj$Sample == sample)], drop = FALSE]

      ######################################
      # Initialize SP Mat Group
      ######################################
      o <- .initializeMat(
        ArrowFile = ArrowFiles[sample],
        Group = matrixName,
        Class = "double",
        Units = "Log2NormCounts",
        cellNames = colnames(sampleMat),
        params = dfParams,
        featureDF = featureDF,
        force = force
      )

      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictionScore"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictionScore")
      )

      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedGroup"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedGroup")
      )

      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedCell"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedCell")
      )

      .messageDiffTime(sprintf("%s Adding GeneIntegrationMatrix to ArrowFile!", prefix), tstart, verbose = verboseHeader)

      for(z in seq_along(allChr)){

        chrz <- allChr[z]

        .messageDiffTime(sprintf("Adding GeneIntegrationMatrix to %s for Chr (%s of %s)!", sample, z, length(allChr)), tstart, verbose = verboseAll)

        idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
        matz <- sampleMat[idz, ,drop=FALSE]
        stopifnot(identical(paste0(featureDF$name[idz]), paste0(rownames(matz))))

        #Write sparseMatrix to Arrow File!
        o <- .addMatToArrow(
          mat = matz, 
          ArrowFile = ArrowFiles[sample], 
          Group = paste0(matrixName, "/", chrz), 
          binarize = FALSE,
          addColSums = TRUE,
          addRowSums = TRUE,
          addRowVarsLog2 = TRUE
        )

        #Clean Memory
        rm(matz)

        if(z %% 3 == 0 | z == length(allChr)){
          gc()
        }

      }

      0

    }, threads = threads)

    o <- suppressWarnings(file.remove(integrationFiles))

  }

  .messageDiffTime("Completed Integration with RNA Matrix", tstart, verbose = verboseHeader)

  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedCell,
    name = nameCell,
    force = TRUE
  )

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




####################################################################
# Gene Integration Matrix Methods
####################################################################

#' Add a GeneIntegrationMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, will integrate multiple subsets of scATAC cells with a scRNA experiment, compute matched scRNA profiles and
#' then store this in each samples ArrowFile.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of a matrix in the `ArchRProject` containing gene scores to be used for RNA integration.
#' @param matrixName The name to use for the output matrix containing scRNA-seq integration to be stored in the `ArchRProject`.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' This `reducedDims` will be used in weighting the transfer of data to scRNA to scATAC. See `Seurat::TransferData` for more info.
#' @param seRNA A `SeuratObject` or a scRNA-seq `SummarizedExperiment` (cell x gene) to be integrated with the scATAC-seq data.
#' @param groupATAC A column name in `cellColData` of the `ArchRProj` that will be used to determine the subgroupings specified in `groupList`.
#' This is used to constrain the integration to occur across biologically relevant groups.
#' @param groupRNA A column name in either `colData` (if `SummarizedExperiment`) or `metadata` (if `SeuratObject`) of `seRNA` that 
#' will be used to determine the subgroupings specified in `groupList`. This is used to constrain the integration to occur across biologically relevant groups.
#' Additionally this groupRNA is used for the `nameGroup` output of this function.
#' @param groupList A list of cell groupings for both ATAC-seq and RNA-seq cells to be used for RNA-ATAC integration.
#' This is used to constrain the integration to occur across biologically relevant groups. The format of this should be a list of groups 
#' with subgroups of ATAC and RNA specifying cells to integrate from both platforms. 
#' For example `groupList` <- list(groupA = list(ATAC = cellsATAC_A, RNA = cellsRNA_A), groupB = list(ATAC = cellsATAC_B, RNA = cellsRNA_B))
#' @param sampleCellsATAC An integer describing the number of scATAC-seq cells to be used for integration. 
#' This number will be evenly sampled across the total number of cells in the ArchRProject.
#' @param sampleCellsRNA An integer describing the number of scRNA-seq cells to be used for integration.
#' @param embeddingATAC A `data.frame` of cell embeddings such as a UMAP for scATAC-seq cells to be used for density sampling. The `data.frame` object
#' should have a row for each single cell described in `row.names` and 2 columns, one for each dimension of the embedding.
#' @param embeddingRNA A `data.frame` of cell embeddings such as a UMAP for scRNA-seq cells to be used for density sampling. The `data.frame` object
#' should have a row for each single cell described in `row.names` and 2 columns, one for each dimension of the embedding.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param plotUMAP A boolean determining whether to plot a UMAP for each integration block.
#' @param UMAPParams The list of parameters to pass to the UMAP function if "plotUMAP = TRUE". See the function `umap` in the uwot package.
#' @param nGenes The number of variable genes determined by `Seurat::FindVariableGenes()` to use for integration.
#' @param useImputation A boolean value indicating whether to use imputation for creating the Gene Score Matrix prior to integration.
#' @param reduction The Seurat reduction method to use for integrating modalities. See `Seurat::FindTransferAnchors()` for possible reduction methods.
#' @param addToArrow A boolean value indicating whether to add the log2-normalized transcript counts from the integrated matched RNA to the Arrow files.
#' @param scaleTo Each column in the integrated RNA matrix will be normalized to a column sum designated by `scaleTo` prior to adding to Arrow files.
#' @param genesUse If desired a character vector of gene names to use for integration instead of determined ones from Seurat::variableGenes.
#' @param nameCell A column name to add to `cellColData` for the predicted scRNA-seq cell in the specified `ArchRProject`. This is useful for identifying which cell was closest to the scATAC-seq cell.
#' @param nameGroup A column name to add to `cellColData` for the predicted scRNA-seq group in the specified `ArchRProject`. See `groupRNA` for more details.
#' @param nameScore A column name to add to `cellColData` for the predicted scRNA-seq score in the specified `ArchRProject`. These scores represent
#' the assignment accuracy of the group in the RNA cells. Lower scores represent ambiguous predictions and higher scores represent precise predictions. 
#' @param transferParams Additional params to be passed to `Seurat::TransferData`.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exists in the given `input`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional params to be added to `Seurat::FindTransferAnchors`
#' @export
addGeneIntegrationMatrix <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = NULL,
  groupATAC = NULL,
  groupRNA = NULL,
  groupList = NULL,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  plotUMAP = TRUE,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
  nGenes = 2000,
  useImputation = TRUE,
  reduction = "cca",
  addToArrow = TRUE,
  scaleTo = 10000,
  genesUse = NULL,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  threads = getArchRThreads(),
  verbose = TRUE,
  force = FALSE,
  logFile = createLogFile("addGeneIntegrationMatrix"),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  .validInput(input = groupATAC, name = "groupATAC", valid = c("character", "null"))
  .validInput(input = groupRNA, name = "groupRNA", valid = c("character"))
  .validInput(input = groupList, name = "groupList", valid = c("list", "null"))
  .validInput(input = sampleCellsATAC, name = "sampleCellsATAC", valid = c("integer", "null"))
  .validInput(input = sampleCellsRNA, name = "sampleCellsRNA", valid = c("integer", "null"))
  .validInput(input = embeddingATAC, name = "embeddingATAC", valid = c("data.frame", "null"))
  .validInput(input = embeddingRNA, name = "embeddingRNA", valid = c("data.frame", "null"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = plotUMAP, name = "plotUMAP", valid = c("boolean"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))  
  .validInput(input = nGenes, name = "nGenes", valid = c("integer"))
  .validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  .validInput(input = reduction, name = "reduction", valid = c("character"))
  .validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = genesUse, name = "genesUse", valid = c("character", "null"))
  .validInput(input = nameCell, name = "nameCell", valid = c("character"))
  .validInput(input = nameGroup, name = "nameGroup", valid = c("character"))
  .validInput(input = nameScore, name = "nameScore", valid = c("character"))
  .validInput(input = transferParams, name = "transferParams", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logDiffTime("Running Seurat's Integration Stuart* et al 2019", tstart, verbose = verbose, logFile = logFile)

  .requirePackage("Seurat", source = "cran")

  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Input-Parameters", logFile=logFile)

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
  .logDiffTime("Checking ATAC Input", tstart, verbose = verbose, logFile = logFile)

  if (useMatrix %ni% getAvailableMatrices(ArchRProj)) {
    .logMessage(paste0("Matrix ", useMatrix, " does not exist in the provided ArchRProject. See available matrix names from getAvailableMatrices()!"), logFile = logFile)
    stop("Matrix name provided to useMatrix does not exist in ArchRProject!")
  }
  
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
    .logMessage(paste0("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!"), logFile = logFile)
    stop("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!")
  }

  #########################################################################################
  # 2. Check All RNA is a Cell Name 
  #########################################################################################
  .logDiffTime("Checking RNA Input", tstart, verbose = verbose, logFile = logFile)

  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupRNA %ni% colnames(colData(seRNA))){
      .logMessage("groupRNA not in colData of seRNA", logFile = logFile)
      stop("groupRNA not in colData of seRNA")
    }
    seuratRNA$Group <- paste0(colData(seRNA)[, groupRNA, drop = TRUE])
    rm(seRNA)
  }else{
    if(groupRNA %ni% colnames(seRNA@meta.data)){
      .logMessage("groupRNA not in meta.data of Seurat Object", logFile = logFile)
      stop("groupRNA not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- paste0(seRNA@meta.data[,groupRNA])
    rm(seRNA)
  }

  if("RNA" %in% names(seuratRNA@assays)){
    DefaultAssay(seuratRNA) <- "RNA"
  }else{
    stop("'RNA' is not present in Seurat Object's Assays! Please make sure that this assay is present!")
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
    .logMessage("Found cells for RNA not in colnames(seRNA)! Please retry your input!", logFile = logFile)
    stop("Found cells for RNA not in colnames(seRNA)! Please retry your input!")
  }

  seuratRNA <- seuratRNA[, unique(cellRNA)]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = FALSE)

  #########################################################################################
  # 3. Create Integration Blocks
  #########################################################################################

  #Check Gene Names And Seurat RowNames
  geneDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  sumOverlap <- sum(unique(geneDF$name) %in% unique(rownames(seuratRNA)))
  if(sumOverlap < 5){
    stop("Error not enough overlaps (",sumOverlap,") between gene names from gene scores (ArchR) and rna matrix (seRNA)!")
  }
  .logDiffTime(paste0("Found ", sumOverlap, " overlapping gene names from gene scores and rna matrix!"), tstart, verbose = TRUE, logFile = logFile)

  .logDiffTime("Creating Integration Blocks", tstart, verbose = verbose, logFile = logFile)

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
  .logDiffTime("Prepping Interation Data", tstart, verbose = verbose, logFile = logFile)

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
  tmpFile <- .tempfile()
  o <- suppressWarnings(file.remove(paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")))

  if(threads > 1){
    h5disableFileLocking()
  }

  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)

  #Create Output Directory
  outDir1 <- getOutputDirectory(ArchRProj)
  outDir2 <- file.path(outDir1, "RNAIntegration")
  outDir3 <- file.path(outDir2, matrixName)
  dir.create(outDir1, showWarnings = FALSE)
  dir.create(outDir2, showWarnings = FALSE)
  dir.create(outDir3, showWarnings = FALSE)
  prevFiles <- list.files(outDir3, full.names = TRUE)
  prevFiles <- .suppressAll(file.remove(prevFiles))

  tstart <- Sys.time()

  threads2 <- max(ceiling(threads * 0.75), 1) #A Little Less here for now

  .logDiffTime(paste0("Computing Integration in ", length(blockList), " Integration Blocks!"), tstart, verbose = verbose, logFile = logFile)

  #Integration
  dfAll <- .safelapply(seq_along(blockList), function(i){

    prefix <- sprintf("Block (%s of %s) :", i , length(blockList))

    .logDiffTime(sprintf("%s Computing Integration", prefix), tstart, verbose = verbose, logFile = logFile)
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
    .logDiffTime(sprintf("%s Identifying Variable Genes", prefix), tstart, verbose = verbose, logFile = logFile)
    subRNA <- FindVariableFeatures(object = subRNA, nfeatures = nGenes, verbose = FALSE)
    subRNA <- ScaleData(object = subRNA, verbose = FALSE)
    if(is.null(genesUse)){
      genesUse <- VariableFeatures(object = subRNA)
    }

    ##############################################################################################
    #2. Get Gene Score Matrix and Create Seurat ATAC
    ##############################################################################################
    .logDiffTime(sprintf("%s Getting GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
    mat <- .getPartialMatrix(
      getArrowFiles(subProj), 
      featureDF = geneDF[geneDF$name %in% genesUse,], 
      threads = 1,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = FALSE
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]
    .logThis(mat, paste0("GeneScoreMat-Block-",i), logFile=logFile)

    #Impute Matrix (its already scaled internally in ArrowFiles)
    if(useImputation){
      .logDiffTime(sprintf("%s Imputing GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$reducedDims <- reducedDims
      imputeParams$dimsToUse <- dimsToUse
      imputeParams$scaleDims <- scaleDims
      imputeParams$corCutOff <- corCutOff
      imputeParams$threads <- 1
      imputeParams$logFile <- logFile
      subProj <- suppressMessages(do.call(addImputeWeights, imputeParams))
      mat <- suppressMessages(imputeMatrix(mat = mat, imputeWeights = getImputeWeights(subProj), verbose = FALSE, logFile = logFile))
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]]))) #Clean Up Space
      .logThis(mat, paste0("GeneScoreMat-Block-Impute-",i), logFile=logFile)
    }

    #Log-Normalize 
    mat <- log(mat + 1) #use natural log
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    
    #Clean Memory
    rm(mat)
    
    #Set Default Assay
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = FALSE)

    ##############################################################################################
    #3. Transfer Anchors  
    ############################################################################################## 
    .logDiffTime(sprintf("%s Seurat FindTransferAnchors", prefix), tstart, verbose = verbose, logFile = logFile)
    transferAnchors <- .retryCatch({ #This sometimes can crash in mclapply so we can just add a re-run parameter
      gc()
      Seurat::FindTransferAnchors(
        reference = subRNA, 
        query = seuratATAC, 
        reduction = reduction, 
        features = genesUse,
        verbose = FALSE,
        ...
      )
    }, maxAttempts = 2, logFile = logFile)
    .logThis(paste0(utils::capture.output(transferAnchors),collapse="\n"), paste0("transferAnchors-",i), logFile=logFile)

    ##############################################################################################
    #4. Transfer Data
    ##############################################################################################
    rDSub <- rD[colnames(seuratATAC),,drop=FALSE]
    .logThis(rDSub, paste0("rDSub-", i), logFile = logFile)
    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = rDSub, 
      key = "LSI_", 
      assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- FALSE
    transferParams$dims <- seq_len(ncol(rDSub))
    
    #Group
    .logDiffTime(sprintf("%s Seurat TransferData Cell Group Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)

    #RNA Names
    .logDiffTime(sprintf("%s Seurat TransferData Cell Names Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- colnames(subRNA)
    rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[,1]

    if(addToArrow){
      .logDiffTime(sprintf("%s Seurat TransferData GeneMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
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

    .logDiffTime(sprintf("%s Saving TransferAnchors Joint CCA", prefix), tstart, verbose = verbose, logFile = logFile)
    jointCCA <- DataFrame(transferAnchors@object.list[[1]]@reductions$cca@cell.embeddings)
    jointCCA$Assay <- ifelse(endsWith(rownames(jointCCA), "_reference"), "RNA", "ATAC")
    jointCCA$Group <- NA
    jointCCA$Score <- NA
    jointCCA[paste0(colnames(subRNA), "_reference"), "Group"] <- subRNA$Group
    jointCCA[paste0(matchDF$cellNames, "_query"), "Group"] <- matchDF$predictedGroup
    jointCCA[paste0(matchDF$cellNames, "_query"), "Score"] <- matchDF$predictionScore
    .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))

    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()

    ##############################################################################################
    #5. Add To Temp Hdf5
    ##############################################################################################

    if(addToArrow){

      .logDiffTime(sprintf("%s Transferring Paired RNA to Temp File", prefix), tstart, verbose = verbose, logFile = logFile)

      #Quickly Write to A Temp Hdf5 File Split By Sample to Then Enable Writing to Each Arrow File

      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- .safeSubset( #If Rownames disappeared this will catch that!
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

    .logDiffTime(sprintf("%s Completed Integration", prefix), tstart, verbose = verbose, logFile = logFile)

    gc()

    matchDF$Block <- Rle(i)
    matchDF

  }, threads = threads2) %>% Reduce("rbind", .)

  ##############################################################################################
  #5. Plot UMAPs for Co-Embeddings from CCA
  ##############################################################################################
  if(plotUMAP){
    
    for(i in seq_along(blockList)){

      o <- tryCatch({

        prefix <- sprintf("Block (%s of %s) :", i , length(blockList))

        .logDiffTime(sprintf("%s Plotting Joint UMAP", prefix), tstart, verbose = verbose, logFile = logFile)

        jointCCA <- readRDS(file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))

        set.seed(1) # Always do this prior to UMAP
        UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="cosine", verbose=FALSE))
        UMAPParams$X <- as.data.frame(jointCCA[, grep("CC_", colnames(jointCCA))])
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- tryCatch({
          do.call(uwot::umap, UMAPParams)
        }, error = function(e){
          errorList <- UMAPParams
          .logError(e, fn = "uwot::umap", info = prefix, errorList = errorList, logFile = logFile)
        })

        #Add UMAP and Save Again
        jointCCA$UMAP1 <- uwotUmap[,1]
        jointCCA$UMAP2 <- uwotUmap[,2]
        .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))

        p1 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Assay,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())

        p2 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Group, 
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by scRNA Group"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())

        pdf(file.path(outDir3, paste0("Save-Block", i,"-JointCCA-UMAP.pdf")), width = 12, height = 6, useDingbats = FALSE)
        ggAlignPlots(p1,p2,type="h")
        dev.off()

      }, error = function(e){

      })

    }

  }

  ##############################################################################################
  #6. Read sub-matrices and store in ArrowFiles
  ##############################################################################################

  if(addToArrow){

    .logDiffTime("Transferring Data to ArrowFiles", tstart, verbose = verbose, logFile = logFile)

    matrixName <- .isProtectedArray(matrixName)

    integrationFiles <- paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")

    if(!all(file.exists(integrationFiles))){
      .logMessage("Something went wrong with integration as not all temporary files containing integrated RNA exist!", logFile = logFile)
      stop("Something went wrong with integration as not all temporary files containing integrated RNA exist!")
    }

    h5list <- .safelapply(seq_along(integrationFiles), function(x){
      h5ls(integrationFiles[x])
    }, threads = threads)

    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)

    o <- .safelapply(seq_along(allSamples), function(y){

      sample <- allSamples[y]

      prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))

      .logDiffTime(sprintf("%s Getting GeneIntegrationMatrix From TempFiles!", prefix), tstart, verbose = verbose, logFile = logFile)

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
      o <- .createArrowGroup(ArrowFile = ArrowFiles[sample], group = matrixName, force = force)

      o <- .initializeMat(
        ArrowFile = ArrowFiles[sample],
        Group = matrixName,
        Class = "double",
        Units = "NormCounts",
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

      .logDiffTime(sprintf("%s Adding GeneIntegrationMatrix to ArrowFile!", prefix), tstart, verbose = verbose, logFile = logFile)

      for(z in seq_along(allChr)){

        chrz <- allChr[z]

        .logDiffTime(sprintf("Adding GeneIntegrationMatrix to %s for Chr (%s of %s)!", sample, z, length(allChr)), tstart, verbose = FALSE, logFile = logFile)

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
          addRowVarsLog2 = TRUE,
          logFile = logFile
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

  .logDiffTime("Completed Integration with RNA Matrix", tstart, verbose = verbose, logFile = logFile)

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

  .endLogging(logFile = logFile)

  return(ArchRProj)

}

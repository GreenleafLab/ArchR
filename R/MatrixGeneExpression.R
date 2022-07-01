####################################################################
# Gene Activity Score Methods
####################################################################

#' Add Gene Expression Matrix to ArrowFiles or an ArchRProject
#' 
#' This function, for each sample, will add gene expression values from a paired scATAC-seq + scRNA-seq
#' multi modal assay to the ArrowFiles or ArchRProject.
#'
#' @param input An `ArchRProject` object or character vector of ArrowFiles.
#' @param seRNA A a scRNA-seq `SummarizedExperiment` (cell x gene) to be integrated with the scATAC-seq data. 
#' Cell names from this object much match those of the cell names in the ArrowFiles/ArchRProject. We will add support shortly
#' for Seurat Objects (see `Seurat::as.SingleCellExperiment`). The provided values `MUST` be in counts (integer), not log transformed. 
#' @param chromSizes A GRanges object of the chromosome lengths. See `getChromSizes` for more info.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from this analysis.
#' @param scaleTo Each column in the calculated gene score matrix will be normalized to a column sum designated by `scaleTo`.
#' @param verbose A boolean describing whether to print to console messages of progress.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param strictMatch A boolean value indicating whether every cell in `input` must be represented in `seRNA`. If set to `FALSE`,
#' and this `GeneExpressionMatrix` is used for certain downstream analyses such as `addIterativeLSI()`, then errors may occur
#' because not all cells will have relevant information.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exist in the given `input`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addGeneExpressionMatrix <- function(
  input = NULL,
  seRNA = NULL,
  chromSizes = getChromSizes(input),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  strictMatch = FALSE,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
  ){

  .validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment"))
  .validInput(input = chromSizes, name = "chromSizes", valid = c("granges"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = strictMatch, name = "strictMatch", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))


  if(inherits(input, "ArchRProject")){
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
    if(is.null(chromSizes)){
      chromSizes <- getChromSizes(input)
    }
  }else if(inherits(input, "character")){
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
    if(is.null(chromSizes)){
      chromSizes <- getChromSizes()
    }
  }else{
    stop("Error Unrecognized Input!")
  }
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addGeneExpressionMatrix Input-Parameters", logFile = logFile)

  cellsInArrows <- unlist(lapply(ArrowFiles, .availableCells), use.names=FALSE)
  if(!is.null(allCells)){
    cellsInArrows <- allCells
  }

  overlap <- sum(cellsInArrows %in% colnames(seRNA)) / length(cellsInArrows)
  .logMessage("Overlap w/ scATAC = ", round(overlap,3), logFile = logFile, verbose = TRUE)

  if(overlap == 0){
    stop("No overlapping cell names found between ArrowFiles and seRNA object! Cell names in ArrowFiles must match colnames in seRNA!")
  } else if(overlap != 1) {
    if(strictMatch){
      stop("Error! 'strictMatch = TRUE' and not all cells in input are represented in the provided gene expression seRNA. To proceed, please subset your ArchRProject using the subsetArchRProject() function to contain only cells present in seRNA or set 'strictMatch = FALSE'.")
    } else {
      .logMessage("Warning! Not all cells in input exist in seRNA! This may cause downstream issues with functions that require information from all cells. For example, addIterativeLSI() will not work on this GeneExpressionMatrix! To remove these mis-matched cells, subset your ArchRProject using the subsetArchRProject() function to contain only cells present in seRNA and set 'strictMatch = TRUE'", logFile = logFile, verbose = TRUE)
    }
  }

  splitCells <- split(cellsInArrows, stringr::str_split(cellsInArrows, pattern = "#", simplify=TRUE)[,1])
  overlapPerSample <- unlist(lapply(splitCells, function(x) sum(x %in% colnames(seRNA))))
  .logMessage("Overlap Per Sample w/ scATAC : ", paste(paste(names(overlapPerSample), round(overlapPerSample,3), sep = "="), collapse=","), logFile = logFile, verbose = TRUE)

  #Get QC Info
  assay(seRNA) <- Matrix::Matrix(assay(seRNA), sparse=TRUE)
  nUMI <- Matrix::colSums(assay(seRNA))
  mb <- assay(seRNA)
  mb@x[mb@x > 0] <- 1
  nGenes <- Matrix::colSums(mb)
  rm(mb)
  MitoRatio <- Matrix::colSums(assay(seRNA)[grep("^MT", rownames(assay(seRNA))),]) / nUMI
  RiboRatio <- Matrix::colSums(assay(seRNA)[grep("^RP", rownames(assay(seRNA))),]) / nUMI
  qcInfo <- DataFrame(nUMI = nUMI, nGenes = nGenes, MitoRatio = MitoRatio, RiboRatio = RiboRatio)
  colnames(qcInfo) <- paste0("Gex_", colnames(qcInfo))
  
  #Filter seRNA
  seRNA <- seRNA[BiocGenerics::which(seqnames(seRNA) %bcin% seqnames(chromSizes))]
  seRNA <- seRNA[BiocGenerics::which(seqnames(seRNA) %bcni% excludeChr)]

  #Dedup
  idxDup <- which(rownames(seRNA) %in% rownames(seRNA[duplicated(rownames(seRNA))]))
  names(idxDup) <- rownames(seRNA)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(seRNA[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    seRNA <- seRNA[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  #Add Index To RNA Ranges
  features <- rowRanges(seRNA)
  features <- sort(sortSeqlevels(features), ignore.strand = TRUE)
  features <- split(features, seqnames(features))
  features <- lapply(features, function(x){
  mcols(x)$idx <- seq_along(x)
    return(x)
  })
  features <- Reduce("c",features)
  rowData(seRNA)$idx <- features[rownames(seRNA)]$idx

  .logThis(qcInfo, "qcInfo", logFile = logFile)

   #Add args to list
  args <- mget(names(formals()), sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addGeneExpressionMat
  args$registryDir <- file.path(outDir, "addGeneExpressionMatRegistry")
  args$qcInfo <- qcInfo
  args$seRNA <- seRNA

  #Remove Input from args
  args$input <- NULL
  args$chromSizes <- NULL
  args$strictMatch <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

    .endLogging(logFile = logFile)

    #Return Output
    if(inherits(input, "ArchRProject")){

      qcInfo <- qcInfo[rownames(qcInfo) %in% input$cellNames, ]

      for(i in seq_len(ncol(qcInfo))){
    input <- addCellColData(
      ArchRProj = input, 
      data = as.vector(qcInfo[,i]), 
      name = paste0(colnames(qcInfo)[i]), 
      cells = paste0(rownames(qcInfo)), 
      force = force
    )
      }
    
      return(input)

    }else{

      return(unlist(outList))

    }


} 

.addGeneExpressionMat <- function(
  i = NULL,
  ArrowFiles = NULL, 
  seRNA = NULL,
  qcInfo = NULL,
  excludeChr = NULL,
  scaleTo = NULL,
  cellNames = NULL, 
  allCells = NULL,
  tstart = NULL,
  subThreads = 1,
  force = FALSE,
  verbose = TRUE,
  logFile = NULL
  ){

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)
  
  #Check
  matrixName <- "GeneExpressionMatrix"
  o <- h5closeAll()
  o <- .createArrowGroup(ArrowFile = ArrowFile, group = matrixName, force = force, logFile = logFile)
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
 
  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  #Identify Overlapping Cells
  cellNames <- cellNames[cellNames %in% colnames(seRNA)]
  seRNA <- seRNA[, cellNames]

  dfParams <- data.frame(
    scaleTo = scaleTo, 
    exclude = excludeChr,
    stringsAsFactors = FALSE
  )

  featureDF <- data.frame(
  seqnames = paste0(seqnames(seRNA)), 
  idx = mcols(seRNA)$idx, 
  start = start(seRNA), 
  end = end(seRNA), 
  name = rownames(seRNA),
    strand = as.integer(strand(seRNA)),
  stringsAsFactors = FALSE
  )

  .logThis(featureDF, "featureDF", logFile = logFile)

  ######################################
  # Initialize SP Mat Group
  ######################################
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = "double",
    Units = "NormCounts",
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = TRUE
  )

  ######################################
  # Normalize and Insert Log2 Normalized Counts
  ######################################

  assay(seRNA) <- .normalizeCols(assay(seRNA), scaleTo = scaleTo)

  uniqueChr <- unique(featureDF$seqnames)

  for(z in seq_along(uniqueChr)){

  o <- tryCatch({

    o <- h5closeAll()
    chr <- uniqueChr[z]
    matz <- assay(seRNA[BiocGenerics::which(seqnames(seRNA)==chr), ])
    .logDiffTime(sprintf("Adding %s to %s for Chr (%s of %s)!", sampleName, matrixName, z, length(uniqueChr)), tstart, verbose = verbose, logFile = logFile)

      #Write sparseMatrix to Arrow File!
      o <- .addMatToArrow(
        mat = matz, 
        ArrowFile = ArrowFile, 
        Group = paste0(matrixName, "/", chr), 
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = TRUE #add for integration analyses
      )
      gc()

      if(z %% 3 == 0 | z == length(uniqueChr)){
        gc()
      }

    }, error = function(e){

      errorList <- list(
        ArrowFile = ArrowFile,
        chr = chr,
        mat = if(exists("matz", inherits = FALSE)) matz else "matz"
      )

      .logError(e, fn = ".addGeneExpressionMat AddToArrow", info = sampleName, errorList = errorList, logFile = logFile)

    })


  }

  #Add Info To Arrow Files
  allCells <- .availableCells(ArrowFile, passQC = FALSE)

  qcInfoi <- qcInfo[rownames(qcInfo) %in% colnames(seRNA), ]

  for(i in seq_len(ncol(qcInfo))){

    infoi <- rep(-1, length(allCells))
    names(infoi) <- allCells
    infoi[rownames(qcInfoi)] <- qcInfoi[,i]

    o <- h5closeAll()
    h5write(infoi, file = ArrowFile, paste0("Metadata/", colnames(qcInfoi)[i]))
    o <- h5closeAll()

  }

  ArrowFile

}



#' @include HelperUtils.R
NULL

#' @export
exportGroupSE <- function(...){
    .Deprecated("getGroupSE")
    getGroupSE(...)
}

#' Export Group Summarized Experiment 
#' 
#' This function will group, summarize and export a summarized experiment for a assay in a ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the matrix in the ArrowFiles. See getAvailableMatrices to see options
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param divideN A boolean describing whether to divide by the number of cells.
#' @param scaleTo Depth normalize to this value if not NULL.
#' @param threads An integer specifying the number of threads for parallel.
#' @param verbose A boolean specifying to print messages during computation.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Get Group SE
#' se <- getGroupSE(proj, useMatrix = "PeakMatrix", groupBy = "Clusters")
#'
#' @export
getGroupSE <- function(
  ArchRProj = NULL,
  useMatrix = NULL,
  groupBy = "Sample",
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = divideN, name = "divideN", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getGroupSE Input-Parameters", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)
  featureDF <- .getFeatureDF(ArrowFiles, subGroup = useMatrix)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)
  if(!.isDiscrete(Groups)){
    .logStop("groupBy must be a discrete variable!", logFile = logFile)
  }
  Cells <- ArchRProj$cellNames

  .logMessage("Getting Group Matrix", logFile=logFile)
  groupMat <- tryCatch({
    .getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = split(Cells, Groups),
      useIndex = FALSE,
      verbose = verbose
    )
  }, error = function(e){
    errorList <- list(
      ArrowFiles = ArrowFiles, 
      featureDF = featureDF,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = split(Cells, Groups),
      useIndex = FALSE,
      verbose = verbose
    )
    .logError(e, fn = ".getGroupMatrix", info = "", errorList = errorList, logFile = logFile)
  })

  if(divideN){
    .logMessage("Normalizing by number of Cells", logFile=logFile)
    nCells <- table(Groups)[colnames(groupMat)]
    groupMat <- t(t(groupMat) / as.vector(nCells))
  }

  #Normalize
  if(!is.null(scaleTo)){
    .logMessage("Depth Normalizing", logFile=logFile)
    groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
  }

  assayList <- SimpleList(
    matrix = groupMat
  )
  names(assayList) <- useMatrix
  rm(groupMat)
  .logThis(assayList, "assayList", logFile = logFile)

  ccd <- getCellColData(ArchRProj)
  .logThis(ccd, "ccd", logFile = logFile)

  idx <- lapply(seq_len(ncol(ccd)), function(x){
    is.numeric(ccd[,x])
  }) %>% unlist %>% which

  if(length(idx) > 0){
    cD <- lapply(seq_along(idx), function(x){
      unlist(lapply(split(ccd[Cells, idx[x]], Groups), median))[colnames(assayList[[1]])]
    }) %>% Reduce("cbind", .) %>% DataFrame
    colnames(cD) <- colnames(ccd)[idx]
    rownames(cD) <- colnames(assayList[[1]])
    cD$nCells <- as.vector(table(Groups)[colnames(assayList[[1]])])
  }else{
    cD <- DataFrame(
      rownames = colnames(assayList[[1]]), 
      nCells = as.vector(table(Groups)[colnames(assayList[[1]])])
    )
  }

  .logThis(cD, "cD", logFile = logFile)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assayList,
    colData = cD,
    rowData = featureDF
  )

  .logThis(se, "se", logFile = logFile)

  .endLogging(logFile = logFile)

  se

}

#' Export PseudoBulk Group Summarized Experiment 
#' 
#' This function will determine cell groups for pseudobulk, summarize and export a summarized experiment for a assay in a ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the matrix in the ArrowFiles. See getAvailableMatrices to see options
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param divideN A boolean describing whether to divide by the number of cells.
#' @param scaleTo Depth normalize to this value if not NULL.
#' @param useLabels A boolean value indicating whether to use sample labels to create sample-aware subgroupings during as pseudo-bulk replicate generation.
#' @param sampleLabels The name of a column in `cellColData` to use to identify samples. In most cases, this parameter should be left as `NULL` and you
#' should only use this parameter if you do not want to use the default sample labels stored in `cellColData$Sample`. However, if your individual Arrow
#' files do not map to individual samples, then you should set this parameter to accurately identify your samples. This is the case in (for example)
#' multiplexing applications where cells from different biological samples are mixed into the same reaction and demultiplexed based on a lipid barcode or genotype.
#' @param minCells The minimum number of cells required in a given cell group to permit insertion coverage file generation.
#' @param maxCells The maximum number of cells to use during insertion coverage file generation.
#' @param minReplicates The minimum number of pseudo-bulk replicates to be generated.
#' @param maxReplicates The maximum number of pseudo-bulk replicates to be generated.
#' @param sampleRatio The fraction of the total cells that can be sampled to generate any given pseudo-bulk replicate.
#' @param verbose A boolean specifying to print messages during computation.
#' @param threads An integer specifying the number of threads for parallel.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Get Group SE
#' se <- getPBGroupSE(proj, useMatrix = "PeakMatrix", groupBy = "Clusters")
#'
#' @export
getPBGroupSE <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  divideN = TRUE,
  scaleTo = 10000,
  useLabels = TRUE,
  sampleLabels = "Sample",
  minCells = 40,
  maxCells = 500,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getPBGroupSE")
  ){

  #Get PB
  cellGroups <- suppressMessages(addGroupCoverages(
    ArchRProj = ArchRProj,
    groupBy = groupBy,
    useLabels = useLabels,
    sampleLabels = sampleLabels,
    minCells = minCells,
    maxCells = maxCells,
    maxFragments = 10^9,
    minReplicates = minReplicates,
    maxReplicates = maxReplicates,
    sampleRatio = sampleRatio,
    returnGroups = TRUE,
    force = TRUE,
    logFile = logFile
  ))
  labeledCells <- unlist(unlist(cellGroups,use.names=TRUE), use.names=TRUE)
  labeledGroups <- names(labeledCells)
  names(labeledGroups) <- labeledCells
  ArchRProj@cellColData$TMP_PB_12312412312312312 <- paste0(as.vector(labeledGroups[rownames(ArchRProj@cellColData)]))

  #Matrix
  sePB <- getGroupSE(
    ArchRProj = ArchRProj,
    useMatrix = useMatrix,
    groupBy = "TMP_PB_12312412312312312",
    divideN = divideN,
    scaleTo = scaleTo,
    threads = threads,
    verbose = verbose,
    logFile = logFile
  )
  sePB <- sePB[,colnames(sePB) != "NA"]
  colData(sePB)$Group <- stringr::str_split(colnames(sePB), pattern="\\.Rep", simplify=TRUE)[,1]

  sePB

}

#' Export Group BigWigs
#' 
#' This function will group, summarize and export a bigwig for each group in an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality. Accepted values are
#' "None", "ReadsInTSS", "nCells", "ReadsInPromoter", or "nFrags".
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param maxCells Maximum number of cells used for each bigwig.
#' @param ceiling Maximum contribution of accessibility per cell in each tile.
#' @param verbose A boolean specifying to print messages during computation.
#' @param threads An integer specifying the number of threads for parallel.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Get Group BW
#' bw <- getGroupBW(proj, groupBy = "Clusters")
#'
#' @export
getGroupBW <- function(
  ArchRProj = NULL,
  groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "useMatrix", valid = c("character"))
  .validInput(input = normMethod, name = "groupBy", valid = c("character"))
  .validInput(input = tileSize, name = "divideN", valid = c("integer"))
  .validInput(input = maxCells, name = "scaleTo", valid = c("integer", "null"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()

  normMethods <- c("None", "ReadsInTSS", "nCells", "ReadsInPromoter", "nFrags")
  
  if(tolower(normMethod) %ni% tolower(normMethods)){
    stop(paste0("normMethod (",normMethod,") not in supported normMethods : ", paste0(normMethods, collapse=", ")))
  }

  .startLogging(logFile = logFile)
  .logThis(normMethod, "normMethod", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)

  if(tolower(normMethod) %in% tolower(c("ReadsInTSS", "ReadsInPromoter", "nFrags"))){
    normBy <- getCellColData(ArchRProj = ArchRProj, select = normMethod)
  }else{
    normBy <- NULL
  }

  if(!.isDiscrete(Groups)){
    stop("groupBy must be a discrete variable!")
  }

  Cells <- ArchRProj$cellNames

  cellGroups <- split(Cells, Groups)

  if(!is.null(maxCells)){
  
    gnames <- names(cellGroups)
  
    cellGroups <- lapply(seq_along(cellGroups), function(x){
      if(length(cellGroups[[x]]) > maxCells){
        sample(cellGroups[[x]], maxCells)
      }else{
        cellGroups[[x]]
      }
    })
  
    names(cellGroups) <- gnames
  
  }

  bwDir1 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs")
  bwDir2 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs", groupBy)

  dir.create(bwDir1, showWarnings = FALSE)
  dir.create(bwDir2, showWarnings = FALSE)

  o <- suppressWarnings(file.remove(list.files(bwDir2, full.names = TRUE)))

  #cellsInArrow <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  cellsInArrow <- split(
    rownames(getCellColData(ArchRProj)), 
    stringr::str_split(rownames(getCellColData(ArchRProj)), pattern="\\#", simplify=TRUE)[,1]
  )
  availableChr <- .availableSeqnames(head(getArrowFiles(ArchRProj)))
  chromLengths <- getChromLengths(ArchRProj)
  chromSizes <- getChromSizes(ArchRProj)
  tiles <- unlist(slidingWindows(chromSizes, width = tileSize, step = tileSize))

  #H5 File Lock Check
  h5lock <- setArchRLocking()
  if(h5lock){
    threads <- 1
  }else{
    if(threads > 1){
      message("subThreading Enabled since ArchRLocking is FALSE see `addArchRLocking`")
    }
    threads <- threads
  }

  covFiles <- c()

  for(i in seq_along(cellGroups)){

    o <- .createGroupBW(
      i = i, 
      cellGroups = cellGroups,
      ArrowFiles = getArrowFiles(ArchRProj), 
      cellsInArrow = cellsInArrow, 
      availableChr = availableChr,
      chromLengths = chromLengths,
      normMethod = normMethod,
      normBy = normBy,
      ceiling = ceiling,
      tiles = tiles, 
      tileSize = tileSize,
      bwDir = bwDir2,
      tstart = tstart, 
      verbose = verbose,
      logFile = logFile,
      threads = threads
    )

    covFiles <- c(covFiles, o)

  }

  .endLogging(logFile = logFile)

  covFiles

}

.createGroupBW <- function(
  i = NULL, 
  cellGroups = NULL,
  ArrowFiles = NULL, 
  cellsInArrow = NULL, 
  availableChr = NULL,
  chromLengths = NULL, 
  tiles = NULL,
  ceiling = NULL,
  tileSize = 100,
  normMethod = NULL,
  normBy = NULL,
  bwDir = "bigwigs",
  tstart = NULL, 
  verbose = TRUE,
  logFile = NULL,
  threads = 1
  ){

  .logDiffTime(sprintf("%s (%s of %s) : Creating BigWig for Group", names(cellGroups)[i], i, length(cellGroups)), tstart, logFile = logFile, verbose = verbose)

  #Cells
  cellGroupi <- cellGroups[[i]]
  #print(sum(normBy[cellGroupi, 1]))

  #Bigwig File!
  covFile <- file.path(bwDir, paste0(make.names(names(cellGroups)[i]), "-TileSize-",tileSize,"-normMethod-",normMethod,"-ArchR.bw"))
  rmf <- .suppressAll(file.remove(covFile))

  covList <- .safelapply(seq_along(availableChr), function(k){

    it <- 0

    for(j in seq_along(ArrowFiles)){
      cellsInI <- sum(cellsInArrow[[names(ArrowFiles)[j]]] %in% cellGroupi)
      if(cellsInI > 0){
        it <- it + 1
        if(it == 1){
          fragik <- .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi)
        }else{
          fragik <- c(fragik, .getFragsFromArrow(ArrowFiles[j], chr = availableChr[k], out = "GRanges", cellNames = cellGroupi))
        }
      }
    }

    tilesk <- tiles[BiocGenerics::which(seqnames(tiles) %bcin% availableChr[k])]

    if(length(fragik) == 0){

      tilesk$reads <- 0

    }else{

      #N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      if (nTiles%%1 != 0) {
          nTiles <- trunc(nTiles) + 1
      }

      #Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$RG, cellGroupi)
      
      mat <- Matrix::sparseMatrix(
          i = c(trunc(start(fragik) / tileSize), trunc(end(fragik) / tileSize)) + 1,
          j = as.vector(c(matchID, matchID)),
          x = rep(1,  2*length(fragik)),
          dims = c(nTiles, length(cellGroupi))
        )

      if(!is.null(ceiling)){
        mat@x[mat@x > ceiling] <- ceiling
      }
      
      mat <- Matrix::rowSums(mat)
      
      rm(fragik, matchID)
         
      tilesk$reads <- mat

      if(tolower(normMethod) %in% c("readsintss", "readsinpromoter", "nfrags")){
        tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
      }else if(tolower(normMethod) %in% c("ncells")){
        tilesk$reads <- tilesk$reads / length(cellGroupi)
      }else if(tolower(normMethod) %in% c("none")){
      }else{
        stop("NormMethod not recognized!")
      }

    }

    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]

    tilesk

  }, threads = threads)

  names(covList) <- availableChr

  covList <- as(covList, "RleList")

  rtracklayer::export.bw(object = covList, con = covFile)

  return(covFile)

}

#' Export Group Fragment Files
#' 
#' This function will group export fragment files for each group in an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and their fragments exported to `outputDirectory`/GroupFragments.
#' @param threads An integer specifying the number of threads for parallel.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Get Group BW
#' frags <- getGroupFragments(proj, groupBy = "Clusters")
#'
#' @export
getGroupFragments <- function(
  ArchRProj = NULL,
  groupBy = "Clusters",
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupFragments")
  ){

  #Cell Col Data
  ccd <- getCellColData(ArchRProj = ArchRProj)

  #Get Groups
  Groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)

  #Cell Split
  cellGroups <- split(getCellNames(ArchRProj), Groups)

  #Outdir
  outDir <- file.path(getOutputDirectory(ArchRProj), "GroupFragments")
  dir.create(outDir, showWarnings = FALSE)

  #Read Fragments From Each Sample Export
  outList <- .safelapply(seq_along(cellGroups), function(x){

    message("Export Fragments : ", x, " of ", length(cellGroups))

    #Get Fragments
    frags <- suppressMessages(getFragmentsFromProject(
      ArchRProj = ArchRProj,
      cellNames = cellGroups[[x]],
      logFile = logFile,
    ) %>% Reduce("c", .))

    #Export Fragments
    dt <- data.frame(
      V1 = seqnames(frags),
      V2 = start(frags) - 1,
      V3 = end(frags),
      V4 = mcols(frags)$RG
    ) %>% data.table

    #Write Bgzip Etc
    groupFile <- file.path(outDir, paste0(groupBy, ".", names(cellGroups)[x], ".tsv"))
    data.table::fwrite(dt, groupFile, sep = "\t", col.names = FALSE)
    Rsamtools::bgzip(groupFile)
    file.remove(groupFile)
    .fileRename(paste0(groupFile, ".bgz"), paste0(groupFile, ".gz"))
    paste0(groupFile, ".gz")

  }, threads = threads)

  #Return
  unlist(outList)

}

#' Export Group Fragment Files from a Project
#' 
#' This function will group export fragment files for each user-specified
#' group in an ArchRProject and output them under a directory. 
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and their fragments exported to `outputDirectory`/GroupFragments.
#' @param outDir the directory to output the group fragment files.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Create directory for fragments
#' getGroupFragmentsFromProj(proj, groupBy = "Clusters", outDir = "./Shiny/Fragments")
#'
.getGroupFragsFromProj <- function(ArchRProj = NULL,
                                   groupBy = NULL,
                                   outDir = file.path("Shiny", "fragments")) {
  dir.create(outDir, showWarnings = FALSE)
  
  # find barcodes of cells in that groupBy.
  cellGroups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, cellGroups)
  
  # outputs unique cell groups (e.g. cluster).
  groupIDs <- names(cellGroups)
  
  
   .safelapply(seq_along(groupIDs), function(x)){
    cat("Making fragment file for cluster:", groupIDs[x], "\n")
    # get GRanges with all fragments for that cluster
    cellNames = cellGroups[[groupIDs[x]]]
    fragments <-
      getFragmentsFromProject(ArchRProj = ArchRProj, cellNames = cellNames)
    fragments <- unlist(fragments, use.names = FALSE)
    # filter Fragments
    fragments <-
      GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
    saveRDS(fragments, file.path(outDir, paste0(groupIDs[x], "_cvg.rds")))
  }
}

#' Export Cluster Coverage from an ArchRProject
#' 
#' This function will group export fragment files for each user-specified
#' group in an ArchRProject and output them under a directory. 
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. 
#'        All insertions in a single bin will be summed.
#' @param scaleFactor A numeric scaling factor to weight genes based on the inverse of there length
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param outDir the directory to output the group fragment files.
#'
.getClusterCoverage <- function(ArchRProj = NULL,
                                tileSize = 100,
                                scaleFactor = 1,
                                groupBy = "Clusters",
                                outDir = file.path(getOutputDirectory(ArchRProj), "Shiny", "coverage")) {
  fragFiles = list.files(path = file.path(getOutputDirectory(ArchRProj), "Shiny", "fragments"), full.names = TRUE)
  dir.create(outDir, showWarnings = FALSE)
  
  # find barcodes of cells in that groupBy.
  cellGroups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters.
  groupIDs <- names(cellGroups)
  
  chrRegions <- getChromSizes(ArchRProj)
  genome <- getGenome(ArchRProj)
  
  for (file in fragFiles) {
    fragments <- readRDS(file)
    left <- GRanges(seqnames = seqnames(fragments),
                    ranges = IRanges(start(fragments), width = 1))
    right <- GRanges(seqnames = seqnames(fragments),
                     ranges = IRanges(end(fragments), width = 1))
    # call sort() after sortSeqlevels() to sort also the ranges in addition to the chromosomes.
    insertions <- c(left, right) %>% sortSeqlevels() %>% sort()
    
    groupID <- file %>% basename() %>% gsub(".{4}$", "", .)
    # binnedCoverage
    message("Creating bins for group ", groupID, "...")
    bins <-
      unlist(slidingWindows(chrRegions, width = tileSize, step = tileSize))
    
    message("Counting overlaps for group ", groupID, "...")
    bins$reads <-
      countOverlaps(
        bins,
        insertions,
        maxgap = -1L,
        minoverlap = 0L,
        type = "any"
      )
    addSeqLengths(bins, genome)
    
    groupReadsInTSS <-
      ArchRProj@cellColData$ReadsInTSS[cells %in% cellGroups$groupID]

    binnedCoverage <- coverage(bins, weight = bins$reads * scaleFactor )
    saveRDS(binnedCoverage, file.path(outDir, paste0(groupID, "_cvg.rds")))
  }  
}

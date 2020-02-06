####################################################################
# Gene Activity Score Methods
####################################################################

#' Add GeneScoreMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, for each sample, will independently compute counts for each tile
#' per cell and then infer gene activity scores.
#'
#' @param input An `ArchRProject` object or character vector of ArrowFiles.
#' @param genes A stranded `GRanges` object containing the ranges associated with all gene start and end coordinates. 
#' @param geneModel A string giving a "gene model function" used for weighting peaks for gene score calculation. This string should be a function of `x`, where `x` is the stranded distance from the transcription start site of the gene. 
#' @param matrixName The name to be used for storage of the gene activity score matrix in the provided `ArchRProject` or ArrowFiles.
#' @param upstream The number of basepairs upstream of the transcription start site to consider for gene activity score calculation.
#' @param downstream The number of basepairs downstream of the transcription start site to consider for gene activity score calculation.
#' @param tileSize The size of the tiles used for binning counts prior to gene activity score calculation.
#' @param ceiling The maximum counts per tile allowed. This is used to prevent large biases in tile counts.
#' @param useGeneBoundaries A boolean value indicating whether gene boundaries should be employed during gene activity score calculation. Gene boundaries refers to the process of preventing tiles from contributing to the gene score of a given gene if there is a second gene's transcription start site between the tile and the gene of interest.
#' @param scaleTo Each column in the calculated gene score matrix will be normalized to a column sum designated by `scaleTo`.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from this analysis.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist that may be extremeley over-represented and thus biasing the geneScores for genes nearby that locus.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param subThreading A boolean determining whether possible use threads within each multi-threaded subprocess if greater than the number of input samples.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exist in the given `input`.
#' @export
addGeneScoreMatrix <- function(
  input = NULL,
  genes = if(inherits(input, "ArchRProject")) getGenes(input) else NULL,
  geneModel = "max(exp(-abs(x)/5000), exp(-1))",
  matrixName = "GeneScoreMatrix",
  upstream = c(5000, 100000),
  downstream = c(5000, 100000),
  tileSize = 500,
  ceiling = 4,
  useGeneBoundaries = TRUE,
  scaleTo = 10000,
  excludeChr = c("chrY","chrM"),
  blacklist = if(inherits(input, "ArchRProject")) getBlacklist(input) else NULL,
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = FALSE
  ){

  .validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  .validInput(input = genes, name = "genes", valid = c("GRanges"))
  .validInput(input = geneModel, name = "geneModel", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer"))
  .validInput(input = useGeneBoundaries, name = "useGeneBoundaries", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  matrixName <- .isProtectedArray(matrixName, exclude = "GeneScoreMatrix")

  if(inherits(input, "ArchRProject")){
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
  }else if(inherits(input, "character")){
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
  }else{
    stop("Error Unrecognized Input!")
  }
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  #Valid GRanges
  genes <- .validGRanges(genes)

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addGeneScoreMat
  args$registryDir <- file.path(outDir, "GeneScoresRegistry")
  
  if(subThreading){
    h5disableFileLocking()
  }else{
    args$threads <- length(inputFiles)
  }

  #Remove Input from args
  args$input <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(subThreading){
    h5enableFileLocking()
  }

  if(inherits(input, "ArchRProject")){

    return(input)

  }else{

    return(unlist(outList))

  }

}

.addGeneScoreMat <- function(
  i = NULL,
  ArrowFiles = NULL,
  genes = NULL,
  geneModel = "max(exp(-abs(x)/5000), exp(-1))",
  matrixName = "GeneScoreMatrix",
  upstream = c(5000, 100000),
  downstream = c(5000, 100000),
  tileSize = 500,
  ceiling = 4,
  useGeneBoundaries = TRUE,
  scaleTo = 10000,
  excludeChr = c("chrY","chrM"),
  blacklist = NULL,
  cellNames = NULL,
  allCells = NULL,
  force = FALSE,
  tmpFile = NULL,
  subThreads = 1
  ){

  .validInput(input = i, name = "i", valid = c("integer"))
  .validInput(input = ArrowFiles, name = "ArrowFiles", valid = c("character"))
  .validInput(input = genes, name = "genes", valid = c("GRanges"))
  .validInput(input = geneModel, name = "geneModel", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer"))
  .validInput(input = useGeneBoundaries, name = "useGeneBoundaries", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character", "null"))
  .validInput(input = allCells, name = "allCells", valid = c("character", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = tmpFile, name = "tmpFile", valid = c("character", "null"))

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)

  if(is.null(tmpFile)){
    tmpFile <- .tempfile(pattern = paste0("tmp-", .sampleName(ArrowFile)))
  }

  #Check
  if(!suppressMessages(h5createGroup(file = ArrowFile, matrixName))){
    if(force){
      o <- h5delete(file = ArrowFile, name = matrixName)
      o <- h5createGroup(ArrowFile, matrixName)
    }else{
      stop(matrixName, " Already Exists!, set force = TRUE to override!")
    }
  }

  o <- h5closeAll()

  #Add Gene Index
  geneStart <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  seqlevels(geneStart) <- as.character(unique(seqnames(geneStart)))
  geneStart <- geneStart[!is.na(mcols(geneStart)$symbol)]
  geneStart <- resize(geneStart, 1, "start")
  strand(geneStart) <- "*"
  geneStart <- sort(sortSeqlevels(geneStart))
  geneStart <- split(geneStart, seqnames(geneStart))
  geneStart <- lapply(geneStart, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })

  #Blacklist Split
  if(!is.null(blacklist)){
    blacklist <- split(blacklist, seqnames(blacklist))
  }

  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }
  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  tstart <- Sys.time()


  #########################################################################################################
  #First we will write gene scores to a temporary path! rhdf5 delete doesnt actually delete the memory!
  #########################################################################################################
  totalGS <- .safelapply(seq_along(geneStart), function(z){

    #Get Gene Starts
    geneStarti <- geneStart[[z]]
    geneStarti <- geneStarti[order(geneStarti$idx)]
    chri <- paste0(unique(seqnames(geneStarti)))
    .messageDiffTime(sprintf("Creating Temp GeneScoreMatrix for %s, Chr (%s of %s)!", sampleName, z, length(geneStart)), tstart)

    #Read in Fragments
    frag <- .getFragsFromArrow(ArrowFile, chr = chri, out = "IRanges", cellNames = cellNames)
    fragSt <- trunc(start(frag)/tileSize) * tileSize
    fragEd <- trunc(end(frag)/tileSize) * tileSize
    fragBC <- rep(S4Vectors::match(mcols(frag)$RG, cellNames), 2)
    rm(frag)
    gc()

    #Unique Inserts
    uniqIns <- sort(unique(c(fragSt,fragEd)))

    #Construct tile by cell mat!
    matGS <- Matrix::sparseMatrix(
        i = match(c(fragSt, fragEd), uniqIns),
        j = as.vector(fragBC),
        x = rep(1,  2*length(fragSt)),
        dims = c(length(uniqIns), length(cellNames))
      )  
    
    if(!is.null(ceiling)){
      matGS@x[matGS@x > ceiling] <- ceiling
    }

    #Unique Tiles
    uniqueTiles <- IRanges(start = uniqIns, width = tileSize)
    
    #Clean Memory
    rm(uniqIns, fragSt, fragEd, fragBC)
    gc() 

    #Time to Overlap Gene Windows
    if(useGeneBoundaries){

      s <- pmax(
          c(1, start(geneStarti)[-length(geneStarti)] + tileSize), 
          start(geneStarti) - max(downstream)
        )

      e <- pmin(
          c(start(geneStarti)[-1] - tileSize, start(geneStarti)[length(geneStarti)] + max(upstream)), 
          start(geneStarti) + max(upstream)
        )

      extenedGeneStart <- IRanges(start = s, end = pmax(s,e)) #handle negative widths!

      idx <- which(width(extenedGeneStart) < (min(upstream) + min(downstream)))
      
      extenedGeneStart[idx] <- ranges(suppressWarnings(extendGR(geneStarti[idx], upstream = min(upstream), downstream = min(downstream))))
      
    }else{

      extenedGeneStart <- ranges(suppressWarnings(extendGR(geneStarti, upstream = max(upstream), downstream = max(downstream))))

    }

    tmp <- suppressWarnings(findOverlaps(extenedGeneStart, uniqueTiles))
    x <- distance(ranges(geneStarti)[queryHits(tmp)], uniqueTiles[subjectHits(tmp)])

    #Determine Sign for Distance relative to strand
    isMinus <- BiocGenerics::which(strand(geneStarti) == "-")
    signDist <- sign(start(uniqueTiles)[subjectHits(tmp)] - start(ranges(geneStarti))[queryHits(tmp)])
    signDist[isMinus] <- signDist[isMinus] * -1

    #Correct the orientation for the distance!
    x <- x * signDist

    #Evaluate Input Model
    x <- eval(parse(text=geneModel))

    #Remove Blacklisted Tiles!
    if(!is.null(blacklist)){
      blacklisti <- blacklist[[chri]]
      if(is.null(blacklisti) | length(blacklisti) > 0){
        tilesBlacklist <- 1 * (!overlapsAny(uniqueTiles, ranges(blacklisti)))
        if(length(tilesBlacklist) > 0){
          x <- x * tilesBlacklist[subjectHits(tmp)] #Multiply Such That All Blacklisted Tiles weight is now 0!
        }
      }
    }

    #Clean Memory
    rm(isMinus, signDist, extenedGeneStart, uniqueTiles)
    gc()

    #Creating Sparse Matrix
    tmp <- Matrix::sparseMatrix(
      i = queryHits(tmp), 
      j = subjectHits(tmp), 
      x = x, 
      dims = c(length(geneStarti), nrow(matGS)))

    #Calculate Gene Scores
    matGS <- tmp %*% matGS
    colnames(matGS) <- cellNames
    totalGSz <- Matrix::colSums(matGS)

    #Save tmp file
    saveRDS(matGS, file = paste0(tmpFile, "-", chri, ".rds"), compress = FALSE)

    #Clean Memory
    rm(matGS, tmp)
    gc()

    totalGSz

  }, threads = subThreads) %>% Reduce("+", .)


  #########################################################################################################
  #Organize info for ArchR Arrow
  #########################################################################################################
  featureDF <- Reduce("c",geneStart) %>% 
    {data.frame(
      row.names=NULL,
      seqnames=as.character(seqnames(.)),
      start=start(.),
      name=mcols(.)$symbol,
      idx=mcols(.)$idx,
      stringsAsFactors=FALSE)}

  dfParams <- data.frame(
      upstream = upstream,
      downstream = downstream,
      scaleTo = scaleTo,
      tileSize = tileSize,
      ceiling = ceiling,
      geneModel = geneModel,
      stringsAsFactors=FALSE
    )

  ######################################
  # Initialize SP Mat Group
  ######################################
  o <- .initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = "double",
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = force
  )

  #Clean Memory
  rm(dfParams, featureDF, genes)
  gc()

  #Normalize and add to Arrow File!
  for(z in seq_along(geneStart)){

    #Get Chromosome
    chri <- paste0(unique(seqnames(geneStart[[z]])))

    .messageDiffTime(sprintf("Adding GeneScoreMatrix to %s for Chr (%s of %s)!", sampleName, z, length(geneStart)), tstart)

    #Re-Create Matrix for that chromosome!
    matGS <- readRDS(paste0(tmpFile, "-", chri, ".rds"))
    file.remove(paste0(tmpFile, "-", chri, ".rds"))

    #Normalize
    matGS@x <- as.numeric(scaleTo * matGS@x/rep.int(totalGS, Matrix::diff(matGS@p)))

    #Round to Reduce Digits After Final Normalization
    matGS@x <- round(matGS@x, 2)
    matGS <- Matrix::drop0(matGS)

    #Write sparseMatrix to Arrow File!
    o <- .addMatToArrow(
      mat = matGS, 
      ArrowFile = ArrowFile, 
      Group = paste0(matrixName, "/", chri), 
      binarize = FALSE,
      addColSums = TRUE,
      addRowSums = TRUE
      )
    gc()

    #Clean Memory
    rm(matGS)
    gc()

  }

  return(ArrowFile)

}


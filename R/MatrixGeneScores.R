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
#' @param geneModel A string giving a "gene model function" used for weighting peaks for gene score calculation. This string
#' should be a function of `x`, where `x` is the stranded distance from the transcription start site of the gene. 
#' @param matrixName The name to be used for storage of the gene activity score matrix in the provided `ArchRProject` or ArrowFiles.
#' @param extendUpstream The minimum and maximum number of basepairs upstream of the transcription start site to consider for gene
#' activity score calculation.
#' @param extendDownstream The minimum and maximum number of basepairs downstream of the transcription start site or transcription termination site 
#' (based on 'useTSS') to consider for gene activity score calculation.
#' @param useGeneBoundaries A boolean value indicating whether gene boundaries should be employed during gene activity score
#' calculation. Gene boundaries refers to the process of preventing tiles from contributing to the gene score of a given gene
#' if there is a second gene's transcription start site between the tile and the gene of interest.
#' @param geneUpstream An integer describing the number of bp upstream the gene to extend the gene body. This effectively makes the gene body larger as there
#' are proximal peaks that should be weighted equally to the gene body. This parameter is used if 'useTSS=FALSE'.
#' @param geneDownstream An integer describing the number of bp downstream the gene to extend the gene body.This effectively makes the gene body larger as there
#' are proximal peaks that should be weighted equally to the gene body. This parameter is used if 'useTSS=FALSE'.
#' @param useTSS A boolean describing whether to build gene model based on gene TSS or the gene body.
#' @param extendTSS A boolean describing whether to extend the gene TSS. By default useTSS uses the 1bp TSS while this parameter enables the extension of this
#' region with 'geneUpstream' and 'geneDownstream' respectively.
#' @param tileSize The size of the tiles used for binning counts prior to gene activity score calculation.
#' @param ceiling The maximum counts per tile allowed. This is used to prevent large biases in tile counts.
#' @param geneScaleFactor A numeric scaling factor to weight genes based on the inverse of there length i.e. [(Scale Factor)/(Gene Length)]. This
#' is scaled from 1 to the scale factor. Small genes will be the scale factor while extremely large genes will be closer to 1. This scaling helps with
#' the relative gene score value.
#' @param scaleTo Each column in the calculated gene score matrix will be normalized to a column sum designated by `scaleTo`.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from this analysis.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist that may be extremeley over-represented and thus
#' biasing the geneScores for genes nearby that locus.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param subThreading A boolean determining whether possible use threads within each multi-threaded subprocess if greater than the number of input samples.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exist in the given `input`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addGeneScoreMatrix <- function(
  input = NULL,
  genes = getGenes(input),
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "GeneScoreMatrix",
  extendUpstream = c(1000, 100000),
  extendDownstream = c(1000, 100000),
  geneUpstream = 5000, #New Param
  geneDownstream = 0, #New Param
  useGeneBoundaries = TRUE,
  useTSS = FALSE, #New Param
  extendTSS = FALSE,
  tileSize = 500,
  ceiling = 4,
  geneScaleFactor = 5, #New Param
  scaleTo = 10000,
  excludeChr = c("chrY", "chrM"),
  blacklist = getBlacklist(input),
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = FALSE,
  logFile = createLogFile("addGeneScoreMatrix")
  ){

  .validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  .validInput(input = genes, name = "genes", valid = c("GRanges"))
  .validInput(input = geneModel, name = "geneModel", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = extendUpstream, name = "extendUpstream", valid = c("integer"))
  .validInput(input = extendDownstream, name = "extendDownstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = ceiling, name = "ceiling", valid = c("integer"))
  .validInput(input = useGeneBoundaries, name = "useGeneBoundaries", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

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

  if(inherits(mcols(genes)$symbol, "list") | inherits(mcols(genes)$symbol, "SimpleList")){
    stop("Found a list in genes symbol! This is an incorrect format. Please correct your genes!")
  }
  if(!any(colnames(mcols(genes)) == "symbol")) {
    stop("No symbol column in genes! A column named symbol is exected in the GRanges object passed to the genes parameter!")
  }

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addGeneScoreMatrix Input-Parameters", logFile = logFile)

  #Valid GRanges
  genes <- .validGRanges(genes)

  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addGeneScoreMat
  args$registryDir <- file.path(outDir, "GeneScoresRegistry")
  args$logFile <- logFile

  if(subThreading){
    h5disableFileLocking()
  }else{
    args$threads <- min(length(ArrowFiles), threads)
  }

  #Remove Input from args
  args$input <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(subThreading){
    h5enableFileLocking()
  }

  .endLogging(logFile = logFile)

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
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "GeneScoreMatrix",
  extendUpstream = c(1000, 100000),
  extendDownstream = c(1000, 100000),
  geneUpstream = 5000, #New Param
  geneDownstream = 0, #New Param
  useGeneBoundaries = TRUE,
  useTSS = FALSE, #New Param
  extendTSS = FALSE,
  tileSize = 500,
  ceiling = 4,
  geneScaleFactor = 5, #New Param
  scaleTo = 10000,
  excludeChr = c("chrY","chrM"),
  blacklist = NULL,
  cellNames = NULL,
  allCells = NULL,
  force = FALSE,
  tmpFile = NULL,
  subThreads = 1,
  tstart = NULL,
  logFile = NULL
  ){

  .validInput(input = i, name = "i", valid = c("integer"))
  .validInput(input = ArrowFiles, name = "ArrowFiles", valid = c("character"))
  .validInput(input = genes, name = "genes", valid = c("GRanges"))
  .validInput(input = geneModel, name = "geneModel", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = extendUpstream, name = "extendUpstream", valid = c("integer"))
  .validInput(input = extendDownstream, name = "extendDownstream", valid = c("integer"))
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

  if(inherits(mcols(genes)$symbol, "list") | inherits(mcols(genes)$symbol, "SimpleList")){
    stop("Found a list in genes symbol! This is an incorrect format. Please correct your genes!")
  }

  ArrowFile <- ArrowFiles[i]
  sampleName <- .sampleName(ArrowFile)

  if(is.null(tmpFile)){
    tmpFile <- .tempfile(pattern = paste0("tmp-", .sampleName(ArrowFile)))
  }

  #Check
  o <- h5closeAll()
  o <- .createArrowGroup(ArrowFile = ArrowFile, group = matrixName, force = force, logFile = logFile)

  geneRegions <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]

  #Create Gene Regions Then Remove Strand Column
  if(useTSS){
    .logMessage(paste0(sampleName, " .addGeneScoreMat useTSS = TRUE"))
    distMethod <- "GenePromoter"
    geneRegions$geneStart <- start(GenomicRanges::resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(GenomicRanges::resize(geneRegions, 1, "end"))
    geneRegions <- GenomicRanges::resize(geneRegions, 1, "start")
    if(extendTSS){
      geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)
    }
    geneRegions$geneWeight <- geneScaleFactor
  }else{
    .logMessage(paste0(sampleName, " .addGeneScoreMat useTSS = FALSE"))
    distMethod <- "GeneBody"
    geneRegions$geneStart <- start(GenomicRanges::resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(GenomicRanges::resize(geneRegions, 1, "end"))
    geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)
    m <- 1 / width(geneRegions)
    geneRegions$geneWeight <- 1 + m * (geneScaleFactor - 1) / (max(m) - min(m))
  }

  .logDiffTime(sprintf("Computing Gene Scores using distance relative to %s! ", distMethod), tstart, logFile = logFile)

  #Add Gene Index For ArrowFile
  geneRegions <- sort(sortSeqlevels(geneRegions), ignore.strand = TRUE)
  .logThis(geneRegions, paste0(sampleName, " .addGeneScoreMat geneRegions"), logFile = logFile)
  
  geneRegions <- split(geneRegions, seqnames(geneRegions))
  geneRegions <- lapply(geneRegions, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })

  #Blacklist Split
  if(!is.null(blacklist)){
    if(length(blacklist) > 0){
      blacklist <- split(blacklist, seqnames(blacklist))
    }
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
  totalGS <- .safelapply(seq_along(geneRegions), function(z){

    totalGSz <- tryCatch({

      .logDiffTime(sprintf("Creating Temp GeneScoreMatrix for %s, Chr (%s of %s)!", sampleName, z, length(geneRegions)), 
        tstart, verbose = FALSE, logFile = logFile)

      #Get Gene Starts
      geneRegionz <- geneRegions[[z]]
      geneRegionz <- geneRegionz[order(geneRegionz$idx)]
      chrz <- paste0(unique(seqnames(geneRegionz)))

      #Read in Fragments
      frag <- .getFragsFromArrow(ArrowFile, chr = chrz, out = "IRanges", cellNames = cellNames)
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

        geneStartz <- start(GenomicRanges::resize(geneRegionz, 1, "start"))
        geneEndz <- start(GenomicRanges::resize(geneRegionz, 1, "end"))

        pminGene <- pmin(geneStartz, geneEndz)
        pmaxGene <- pmax(geneStartz, geneEndz)

        idxMinus <- BiocGenerics::which(strand(geneRegionz) != "-")
    
        pReverse <- rep(max(extendDownstream), length(pminGene))
        pReverse[idxMinus] <- rep(max(extendUpstream), length(idxMinus))

        pReverseMin <- rep(min(extendDownstream), length(pminGene))
        pReverseMin[idxMinus] <- rep(min(extendUpstream), length(idxMinus))

        pForward <- rep(max(extendUpstream), length(pminGene))
        pForward[idxMinus] <- rep(max(extendDownstream), length(idxMinus))      

        pForwardMin <- rep(min(extendUpstream), length(pminGene))
        pForwardMin[idxMinus] <- rep(min(extendDownstream), length(idxMinus))      

        ################################################################
        #We will test when genes pass by another gene promoter
        ################################################################

        #Start of Range is based on the max observed gene ranged <- direction
        s <- pmax(
          c(1, pmaxGene[-length(pmaxGene)] + tileSize), 
          pminGene - pReverse
        )
        s <- pmin(pminGene - pReverseMin, s)

        #End of Range is based on the max observed gene ranged -> direction
        e <- pmin(
            c(pminGene[-1] - tileSize, pmaxGene[length(pmaxGene)] + pForward[length(pmaxGene)]), 
            pmaxGene + pForward
          )
        e <- pmax(pmaxGene + pForwardMin, e)

        extendedGeneRegion <- IRanges(start = s, end = e)

        idx1 <- which(pminGene - pReverseMin < start(extendedGeneRegion))
        if(length(idx1) > 0){
          stop("Error in gene boundaries minError")
        }

        idx2 <- which(pmaxGene + pForwardMin > end(extendedGeneRegion))
        if(length(idx2) > 0){
          stop("Error in gene boundaries maxError")
        }
       
       rm(s, e, pReverse, pReverseMin, pForward, pForwardMin, geneStartz, geneEndz, pminGene, pmaxGene)

      }else{

        extendedGeneRegion <- ranges(suppressWarnings(extendGR(geneRegionz, upstream = max(extendUpstream), downstream = max(extendDownstream))))

      }

      tmp <- suppressWarnings(findOverlaps(extendedGeneRegion, uniqueTiles))
      x <- distance(ranges(geneRegionz)[queryHits(tmp)], uniqueTiles[subjectHits(tmp)])

      #Determine Sign for Distance relative to strand (Directionality determined based on dist from gene start)
      isMinus <- BiocGenerics::which(strand(geneRegionz) == "-")
      signDist <- sign(start(uniqueTiles)[subjectHits(tmp)] - start(GenomicRanges::resize(geneRegionz,1,"start"))[queryHits(tmp)])
      signDist[isMinus] <- signDist[isMinus] * -1

      #Correct the orientation for the distance!
      x <- x * signDist

      #Evaluate Input Model
      x <- eval(parse(text=geneModel))

      #Get Gene Weights Related to Gene Width
      x <- x * mcols(geneRegionz)$geneWeight[queryHits(tmp)]

      #Remove Blacklisted Tiles!
      if(!is.null(blacklist)){
        if(length(blacklist) > 0){
          blacklistz <- blacklist[[chrz]]
          if(!is.null(blacklistz) | length(blacklistz) > 0){
            tilesBlacklist <- 1 * (!overlapsAny(uniqueTiles, ranges(blacklistz)))
            if(sum(tilesBlacklist == 0) > 0){
              x <- x * tilesBlacklist[subjectHits(tmp)] #Multiply Such That All Blacklisted Tiles weight is now 0!
            }
          }
        }
      }

      #Creating Sparse Matrix
      tmp <- Matrix::sparseMatrix(
        i = queryHits(tmp), 
        j = subjectHits(tmp), 
        x = x, 
        dims = c(length(geneRegionz), nrow(matGS))
      )

      #Calculate Gene Scores
      matGS <- tmp %*% matGS
      colnames(matGS) <- cellNames

      totalGSz <- Matrix::colSums(matGS)

      #Save tmp file
      .safeSaveRDS(matGS, file = paste0(tmpFile, "-", chrz, ".rds"), compress = FALSE)

      #Clean Memory
      rm(isMinus, signDist, extendedGeneRegion, uniqueTiles)
      rm(matGS, tmp)
      gc()

      totalGSz
   
    }, error = function(e){

      errorList <- list(
        ArrowFile = ArrowFile,
        geneRegions = geneRegions,
        blacklist = blacklist,
        chr = chrz,
        totalGSz = if(exists("totalGSz", inherits = FALSE)) totalGSz else "totalGSz",
        matGS = if(exists("matGS", inherits = FALSE)) matGS else "matGS"
      )

      .logError(e, fn = ".addGeneScoreMat TmpGS", info = sampleName, errorList = errorList, logFile = logFile)

    })

    totalGSz

  }, threads = subThreads) %>% Reduce("+", .)


  #########################################################################################################
  #Organize info for ArchR Arrow
  #########################################################################################################
  featureDF <- Reduce("c",geneRegions) %>% 
    {data.frame(
      row.names=NULL,
      seqnames=as.character(seqnames(.)),
      start=mcols(.)$geneStart,
      end=mcols(.)$geneEnd,
      strand=as.integer(strand(.)),
      name=mcols(.)$symbol,
      idx=mcols(.)$idx,
      stringsAsFactors=FALSE)}
  .logThis(featureDF, paste0(sampleName, " .addGeneScoreMat FeatureDF"), logFile = logFile)

  dfParams <- data.frame(
      extendUpstream = extendUpstream,
      extendDownstream = extendDownstream,
      geneUpstream = extendUpstream,
      geneDownstream = extendDownstream,
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
    Units = "NormCounts",
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = TRUE
  )

  #Clean Memory
  rm(dfParams, featureDF, genes)
  gc()

  #Normalize and add to Arrow File!
  for(z in seq_along(geneRegions)){

    o <- tryCatch({

      #Get Chromosome
      chrz <- paste0(unique(seqnames(geneRegions[[z]])))

      .logDiffTime(sprintf("Adding GeneScoreMatrix to %s for Chr (%s of %s)!", sampleName, z, length(geneRegions)), 
        tstart, verbose = FALSE, logFile = logFile)

      #Re-Create Matrix for that chromosome!
      matGS <- readRDS(paste0(tmpFile, "-", chrz, ".rds"))
      file.remove(paste0(tmpFile, "-", chrz, ".rds"))

      #Normalize
      matGS@x <- as.numeric(scaleTo * matGS@x/rep.int(totalGS, Matrix::diff(matGS@p)))

      #Round to Reduce Digits After Final Normalization
      matGS@x <- round(matGS@x, 3)
      matGS <- Matrix::drop0(matGS)

      #Write sparseMatrix to Arrow File!
      o <- .addMatToArrow(
        mat = matGS, 
        ArrowFile = ArrowFile, 
        Group = paste0(matrixName, "/", chrz), 
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = TRUE #add for integration analyses
      )

      #Clean Memory
      rm(matGS)

      if(z %% 3 == 0 | z == length(geneRegions)){
        gc()
      }

    }, error = function(e){

      errorList <- list(
        ArrowFile = ArrowFile,
        geneRegions = geneRegions,
        blacklist = blacklist,
        chr = chrz,
        mat = if(exists("mat", inherits = FALSE)) mat else "mat"
      )

      .logError(e, fn = ".addGeneScoreMat AddToArrow", info = sampleName, errorList = errorList, logFile = logFile)

    })

  }

  return(ArrowFile)

}


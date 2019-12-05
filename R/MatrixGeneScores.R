#' Add GeneScoreMatrix to Arrows/ArchRProject
#' 
#' This function for each sample will independently compute counts for each tile
#' per cell and then infer gene activity scores.
#'
#' @param input ArchRProject or ArrowFiles
#' @param genes genes as a GRanges object
#' @param geneModel gene model as a string for weighting peaks for gene score calculation (function of x)
#' @param upstream upstream the Gene Start to consider for calculation
#' @param downstream downstream the Gene Start to consider for calculation
#' @param tileSize tileSize for binning counts prior to gene score calculation
#' @param ceiling ceiling of read counts per tile (prevent huge biases)
#' @param scaleTo scale gene scores to
#' @param excludeChr exclude chromosomes from this analysis
#' @param blacklist blacklist GRanges used to remove tiles prior to calculation
#' @param threads number of threads
#' @param parallelParam parallel parameters for batch style execution
#' @param force force overwriting previous TileMatrix in ArrowFile
#' @export
addGeneScoreMatrix <- function(
  input = NULL,
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
  threads = 1,
  parallelParam = NULL,
  force = FALSE,
  ...
  ){

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

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  if(inherits(input, "ArchRProject")){

    return(input)

  }else{

    return(unlist(outList))

  }

}

.addGeneScoreMat <- function(
  i,
  ArrowFiles,
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
  ...
  ){

  ArrowFile <- ArrowFiles[i]

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

  totalGS <- rep(0, length(cellNames))
  names(totalGS) <- cellNames

  #########################################################################################################
  #First we will write gene scores to a temporary path! rhdf5 delete doesnt actually delete the memory!
  #########################################################################################################
  for(z in seq_along(geneStart)){

    #Get Gene Starts
    geneStarti <- geneStart[[z]]
    geneStarti <- geneStarti[order(geneStarti$idx)]
    chri <- paste0(unique(seqnames(geneStarti)))
    .messageDiffTime(sprintf("Creating Temporary Gene Score Matrix for Chromosome %s of %s!", z, length(geneStart)), tstart)

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
      
      extenedGeneStart[idx] <- ranges(suppressWarnings(extendGRanges(geneStarti[idx], upstream = min(upstream), downstream = min(downstream))))
      
    }else{

      extenedGeneStart <- ranges(suppressWarnings(extendGRanges(geneStarti, upstream = max(upstream), downstream = max(downstream))))

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
    totalGS <- totalGS + Matrix::colSums(matGS)

    #Save tmp file
    saveRDS(matGS, file = paste0(tmpFile, "-", chri, ".rds"), compress = FALSE)

    #Clean Memory
    rm(matGS, tmp)
    gc()

  }


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

    .messageDiffTime(sprintf("Adding Normalized Gene Score Matrix for Chromosome %s of %s to Arrow File!", z, length(geneStart)), tstart)

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




# #' Add GeneScoreMatrix to Arrows/ArchRProject
# #' 
# #' This function for each sample will independently compute counts for each tile
# #' per cell and then infer gene activity scores.
# #'
# #' @param input ArchRProject or ArrowFiles
# #' @param genes genes as a GRanges object
# #' @param geneModel gene model as a string for weighting peaks for gene score calculation (function of x)
# #' @param upstream upstream the Gene Start to consider for calculation
# #' @param downstream downstream the Gene Start to consider for calculation
# #' @param tileSize tileSize for binning counts prior to gene score calculation
# #' @param ceiling ceiling of read counts per tile (prevent huge biases)
# #' @param scaleTo scale gene scores to
# #' @param excludeChr exclude chromosomes from this analysis
# #' @param blacklist blacklist GRanges used to remove tiles prior to calculation
# #' @param threads number of threads
# #' @param parallelParam parallel parameters for batch style execution
# #' @param force force overwriting previous TileMatrix in ArrowFile
# #' @export
# addGeneScoreMatrix <- function(
#   input = NULL,
#   genes = NULL,
#   geneModel = "exp(-abs(x)/10000)",
#   matrixName = "GeneScoreMatrix",
#   upstream = 100000,
#   downstream = 100000,
#   tileSize = 500,
#   ceiling = 4,
#   scaleTo = 10000,
#   excludeChr = c("chrY","chrM"),
#   blacklist = NULL,
#   threads = 1,
#   parallelParam = NULL,
#   force = FALSE,
#   ...
#   ){

#   matrixName <- .isProtectedArray(matrixName, exclude = "GeneScoreMatrix")

#   if(inherits(input, "ArchRProject")){
#     ArrowFiles <- getArrowFiles(input)
#     allCells <- rownames(getCellColData(input))
#     outDir <- getOutputDirectory(input)
#   }else if(inherits(input, "character")){
#     outDir <- ""
#     ArrowFiles <- input
#     allCells <- NULL
#   }else{
#     stop("Error Unrecognized Input!")
#   }
#   if(!all(file.exists(ArrowFiles))){
#     stop("Error Input Arrow Files do not all exist!")
#   }

#   #Valid GRanges
#   genes <- .validGRanges(genes)

#   #Add args to list
#   args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
#   args$ArrowFiles <- ArrowFiles
#   args$allCells <- allCells
#   args$X <- seq_along(ArrowFiles)
#   args$FUN <- .addGeneScoreMat
#   args$registryDir <- file.path(outDir, "GeneScoresRegistry")

#   #Run With Parallel or lapply
#   outList <- .batchlapply(args)

#   if(inherits(input, "ArchRProject")){

#     return(input)

#   }else{

#     return(unlist(outList))

#   }

# }
  
# .addGeneScoreMat <- function(
#   i,
#   ArrowFiles,
#   genes,
#   matrixName = "GeneScoreMatrix",
#   cellNames = NULL,
#   allCells = NULL,
#   upstream = 100000,
#   downstream = 100000,
#   scaleTo = 10000,
#   tileSize = 200,
#   ceiling = 4,
#   blacklist = NULL,
#   geneModel = "exp(-abs(x)/10000)",
#   excludeChr = c("chrY","chrM"),
#   force = FALSE,
#   tmpFile = NULL,
#   ...
#   ){

#   ArrowFile <- ArrowFiles[i]

#   if(is.null(tmpFile)){
#     tmpFile <- .tempfile(pattern = paste0("tmp-", .sampleName(ArrowFile)))
#   }

#   #Check
#   if(!suppressMessages(h5createGroup(file = ArrowFile, matrixName))){
#     if(force){
#       o <- h5delete(file = ArrowFile, name = matrixName)
#       o <- h5createGroup(ArrowFile, matrixName)
#     }else{
#       stop(matrixName, " Already Exists!, set force = TRUE to override!")
#     }
#   }

#   o <- h5closeAll()

#   #Add Gene Index
#   geneStart <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
#   geneStart <- sort(sortSeqlevels(geneStart))
#   seqlevels(geneStart) <- as.character(unique(seqnames(geneStart)))
#   geneStart <- geneStart[!is.na(mcols(geneStart)$symbol)]
#   geneStart <- resize(geneStart, 1, "start")
#   geneStart <- split(geneStart, seqnames(geneStart))
#   geneStart <- lapply(geneStart, function(x){
#     mcols(x)$idx <- seq_along(x)
#     return(x)
#   })

#   #Blacklist Split
#   if(!is.null(blacklist)){
#     blacklist <- split(blacklist, seqnames(blacklist))
#   }

#   #Get all cell ids before constructing matrix
#   if(is.null(cellNames)){
#     cellNames <- .availableCells(ArrowFile)
#   }
#   if(!is.null(allCells)){
#     cellNames <- cellNames[cellNames %in% allCells]
#   }

#   tstart <- Sys.time()

#   totalGS <- rep(0, length(cellNames))
#   names(totalGS) <- cellNames

#   #########################################################################################################
#   #First we will write gene scores to a temporary path! rhdf5 delete doesnt actually delete the memory!
#   #########################################################################################################
#   for(z in seq_along(geneStart)){

#     #Get Gene Starts
#     geneStarti <- geneStart[[z]]
#     mcols(geneStarti)$idx <- seq_along(geneStarti)
#     chri <- paste0(unique(seqnames(geneStarti)))
#     .messageDiffTime(sprintf("Creating Temporary Gene Score Matrix for Chromosome %s of %s!", z, length(geneStart)), tstart)

#     #Read in Fragments
#     frag <- .getFragsFromArrow(ArrowFile, chr = chri, out = "IRanges", cellNames = cellNames)
#     fragSt <- trunc(start(frag)/tileSize) * tileSize
#     fragEd <- trunc(end(frag)/tileSize) * tileSize
#     fragBC <- rep(S4Vectors::match(mcols(frag)$RG, cellNames), 2)
#     rm(frag)
#     gc()

#     #Unique Inserts
#     uniqIns <- sort(unique(c(fragSt,fragEd)))

#     #Construct tile by cell mat!
#     matGS <- Matrix::sparseMatrix(
#         i = match(c(fragSt, fragEd), uniqIns),
#         j = as.vector(fragBC),
#         x = rep(1,  2*length(fragSt)),
#         dims = c(length(uniqIns), length(cellNames))
#       )  
    
#     if(!is.null(ceiling)){
#       matGS@x[matGS@x > ceiling] <- ceiling
#     }

#     #Unique Tiles
#     uniqueTiles <- IRanges(start = uniqIns, width = tileSize)
    
#     #Clean Memory
#     rm(uniqIns, fragSt, fragEd, fragBC)
#     gc() 

#     #Time to Overlap Gene Windows
#     extenedGeneStart <- ranges(suppressWarnings(extendGRanges(geneStarti, upstream = upstream, downstream = downstream))) #Warning if beyond chromosome this doesnt matter for this analysis
#     tmp <- suppressWarnings(findOverlaps(extenedGeneStart, uniqueTiles))
#     x <- distance(ranges(geneStarti)[queryHits(tmp)], uniqueTiles[subjectHits(tmp)])

#     #Determine Sign for Distance relative to strand
#     isMinus <- BiocGenerics::which(strand(geneStarti) == "-")
#     signDist <- sign(start(uniqueTiles)[subjectHits(tmp)] - start(ranges(geneStarti))[queryHits(tmp)])
#     signDist[isMinus] <- signDist[isMinus] * -1

#     #Correct the orientation for the distance!
#     x <- x * signDist

#     #Evaluate Input Model
#     x <- eval(parse(text=geneModel))

#     #Remove Blacklisted Tiles!
#     if(!is.null(blacklist)){
#       blacklisti <- blacklist[[chri]]
#       if(is.null(blacklisti) | length(blacklisti) > 0){
#         tilesBlacklist <- 1 * (!overlapsAny(uniqueTiles, ranges(blacklisti)))
#         if(length(tilesBlacklist) > 0){
#           x <- x * tilesBlacklist[subjectHits(tmp)] #Multiply Such That All Blacklisted Tiles weight is now 0!
#         }
#       }
#     }

#     #Clean Memory
#     rm(isMinus, signDist, extenedGeneStart, uniqueTiles)
#     gc()

#     #Creating Sparse Matrix
#     tmp <- Matrix::sparseMatrix(
#       i = queryHits(tmp), 
#       j = subjectHits(tmp), 
#       x = x, 
#       dims = c(length(geneStarti), nrow(matGS)))

#     #Calculate Gene Scores
#     matGS <- tmp %*% matGS
#     colnames(matGS) <- cellNames
#     totalGS <- totalGS + Matrix::colSums(matGS)

#     #Save tmp file
#     saveRDS(matGS, file = paste0(tmpFile, "-", chri, ".rds"), compress = FALSE)

#     #Clean Memory
#     rm(matGS, tmp)
#     gc()

#   }


#   #########################################################################################################
#   #Organize info for ArchR Arrow
#   #########################################################################################################
#   featureDF <- Reduce("c",geneStart) %>% 
#     {data.frame(
#       row.names=NULL,
#       seqnames=as.character(seqnames(.)),
#       start=start(.),
#       name=mcols(.)$symbol,
#       idx=mcols(.)$idx,
#       stringsAsFactors=FALSE)}

#   dfParams <- data.frame(
#       upstream = upstream,
#       downstream = downstream,
#       scaleTo = scaleTo,
#       tileSize = tileSize,
#       ceiling = ceiling,
#       geneModel = geneModel,
#       stringsAsFactors=FALSE
#     )

#   ######################################
#   # Initialize SP Mat Group
#   ######################################
#   o <- .initializeMat(
#     ArrowFile = ArrowFile,
#     Group = matrixName,
#     Class = "double",
#     cellNames = cellNames,
#     params = dfParams,
#     featureDF = featureDF,
#     force = force
#   )

#   #Clean Memory
#   rm(dfParams, featureDF, genes)
#   gc()

#   #Normalize and add to Arrow File!
#   for(z in seq_along(geneStart)){

#     #Get Chromosome
#     chri <- paste0(unique(seqnames(geneStart[[z]])))

#     .messageDiffTime(sprintf("Adding Normalized Gene Score Matrix for Chromosome %s of %s to Arrow File!", z, length(geneStart)), tstart)

#     #Re-Create Matrix for that chromosome!
#     matGS <- readRDS(paste0(tmpFile, "-", chri, ".rds"))
#     file.remove(paste0(tmpFile, "-", chri, ".rds"))

#     #Normalize
#     matGS@x <- as.numeric(scaleTo * matGS@x/rep.int(totalGS, Matrix::diff(matGS@p)))

#     #Round to Reduce Digits After Final Normalization
#     matGS@x <- round(matGS@x, 2)
#     matGS <- Matrix::drop0(matGS)

#     #Write sparseMatrix to Arrow File!
#     o <- .addMatToArrow(
#       mat = matGS, 
#       ArrowFile = ArrowFile, 
#       Group = paste0(matrixName, "/", chri), 
#       binarize = FALSE,
#       addColSums = TRUE,
#       addRowSums = TRUE
#       )
#     gc()

#     #Clean Memory
#     rm(matGS)
#     gc()

#   }

#   return(ArrowFile)

# }




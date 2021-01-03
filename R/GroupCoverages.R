#' Add Group Coverages to an ArchRProject object
#' 
#' This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates and then
#' merge these replicates into a single insertion coverage file.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping multiple cells together prior to generation of the insertion coverage file.
#' @param useLabels A boolean value indicating whether to use sample labels to create sample-aware subgroupings during as pseudo-bulk replicate generation.
#' @param minCells The minimum number of cells required in a given cell group to permit insertion coverage file generation.
#' @param maxCells The maximum number of cells to use during insertion coverage file generation.
#' @param maxFragments The maximum number of fragments per cell group to use in insertion coverage file generation. This prevents the generation
#' of excessively large files which would negatively impact memory requirements.
#' @param minReplicates The minimum number of pseudo-bulk replicates to be generated.
#' @param maxReplicates The maximum number of pseudo-bulk replicates to be generated.
#' @param sampleRatio The fraction of the total cells that can be sampled to generate any given pseudo-bulk replicate.
#' @param kmerLength The length of the k-mer used for estimating Tn5 bias.
#' @param threads The number of threads to be used for parallel computing.
#' @param returnGroups A boolean value that indicates whether to return sample-guided cell-groupings without creating coverages.
#' This is used mainly in `addReproduciblePeakSet()` when MACS2 is not being used to call peaks but rather peaks are called from a
#' TileMatrix (`peakMethod = "Tiles"`).
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value that indicates whether or not to overwrite the relevant data in the `ArchRProject` object if
#' insertion coverage / pseudo-bulk replicate information already exists.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addGroupCoverages <- function(
  ArchRProj = NULL,
  groupBy = "Clusters",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 500,
  maxFragments = 25*10^6,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = FALSE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = useLabels, name = "useLabels", valid = c("boolean"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = maxCells, name = "maxCells", valid = c("integer"))
  .validInput(input = maxFragments, name = "maxFragments", valid = c("integer"))
  .validInput(input = minReplicates, name = "minReplicates", valid = c("integer"))
  .validInput(input = maxReplicates, name = "maxReplicates", valid = c("integer"))
  .validInput(input = sampleRatio, name = "sampleRatio", valid = c("numeric"))
  .validInput(input = kmerLength, name = "kmerLength", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = returnGroups, name = "returnGroups", valid = c("boolean"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam","null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  if(minReplicates < 2){
    stop("minReplicates must be at least 2!")
  }

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addGroupCoverages Input-Parameters", logFile = logFile)

  Params <- SimpleList(
    groupBy = groupBy,
    minCells = minCells,
    maxCells = maxCells,
    minReplicates = minReplicates,
    sampleRatio = sampleRatio,
    kmerLength = kmerLength
  )

  if(is.null(ArchRProj@projectMetadata$GroupCoverages)){
    ArchRProj@projectMetadata$GroupCoverages <- SimpleList()
  }

  if(!returnGroups){
    if(!is.null(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])){
      if(!force){
        stop("Group Coverages Already Computed, Set force = TRUE to continue!")
      }
    }
  }else{
    if(!is.null(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])){
      if(!force){
        message("Group Coverages Already Computed Returning Groups, Set force = TRUE to Recompute!")
        return(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])
      }
    }      
  }

  #####################################################
  #Groups 
  #####################################################
  cellNames <- rownames(getCellColData(ArchRProj))
  groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(any(is.na(groups))){
      cellNames <- cellNames[!is.na(groups)]
      groups <- groups[!is.na(groups)]
  }
  uniqueGroups <- gtools::mixedsort(unique(groups))
  tableGroups <- table(groups)[uniqueGroups]

  #####################################################
  #Create Cell Groups
  #####################################################
  cellGroups <- lapply(seq_along(uniqueGroups), function(x){
      subColDat <- getCellColData(ArchRProj)[which(groups==uniqueGroups[x]),]
      cellNamesx <- rownames(subColDat)
      #if(length(cellNamesx) < minCells){
      #  outListx <- SimpleList(LowCellGroup = cellNamesx) or NULL
      #}
      if(useLabels){
        sampleLabelsx <- paste0(subColDat$Sample)
      }else{
        sampleLabelsx <- NULL
      }
      outListx <- .identifyGroupsForPseudoBulk(
        cells = cellNamesx, 
        sampleLabels = sampleLabelsx,
        useLabels = useLabels,
        minCells = minCells, 
        maxCells = maxCells,
        minReplicates = minReplicates, 
        maxReplicates = maxReplicates,
        sampleRatio = sampleRatio,
        prefix = sprintf("%s (%s of %s) :", uniqueGroups[x], x, length(uniqueGroups)),
        logFile = logFile
      )
      if(is.null(outListx)){
        return(NULL)
      }
      if(is.null(names(outListx))){
          names(outListx) <- paste0("Rep", seq_along(outListx))
      }else if(any(names(outListx)=="")){
          names(outListx)[which(names(outListx)=="")] <- paste0("Rep", which(names(outListx)==""))
      }
      outListx
  }) %>% SimpleList
  names(cellGroups) <- uniqueGroups
  Params$cellGroups <- cellGroups

  #####################################################
  #Check For Max Fragments!
  #####################################################
  it <- 0
  for(i in seq_along(cellGroups)){
    for(j in seq_along(cellGroups[[i]])){
      if(sum(getCellColData(ArchRProj, "nFrags")[cellGroups[[i]][[j]],]) > maxFragments){
        it <- it + 1
        nFrags <- getCellColData(ArchRProj, "nFrags")[cellGroups[[i]][[j]],]
        cells <- cellGroups[[i]][[j]][order(nFrags)]
        nFrags <- nFrags[order(nFrags)]
        cellGroups[[i]][[j]] <- cells[which(cumsum(nFrags) < maxFragments)]
      }
    }
  }
  if(it > 0){
    .logDiffTime(sprintf("Further Sampled %s Groups above the Max Fragments!", it), tstart)
  }

  if(returnGroups){
    return(cellGroups)
  }

  #####################################################
  # Arguments for Coverages
  #####################################################

  dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages"), showWarnings = FALSE)
  dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages", groupBy), showWarnings = FALSE)

  unlistGroups <- lapply(seq_along(cellGroups), function(x){
    if(is.null(cellGroups[[x]])){
      NULL
    }else{
      names(cellGroups[[x]]) <- paste0(names(cellGroups)[x], "._.", names(cellGroups[[x]]))
      cellGroups[[x]]
    }
  }) %>% SimpleList %>%unlist()

  args <- list()
  args$X <- seq_along(unlistGroups)
  args$FUN <- .createCoverages
  args$cellGroups <- unlistGroups
  args$genome <- getGenome(ArchRProj)
  args$kmerLength <- kmerLength
  args$ArrowFiles <- getArrowFiles(ArchRProj)
  args$availableChr <- .availableSeqnames(getArrowFiles(ArchRProj))
  args$chromLengths <- getChromLengths(ArchRProj)
  args$cellsInArrow <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  args$covDir <- file.path(getOutputDirectory(ArchRProj), "GroupCoverages", groupBy)
  args$parallelParam <- parallelParam
  args$threads <- threads
  args$verbose <- verbose
  args$tstart <- tstart
  args$logFile <- logFile
  args$registryDir <- file.path(getOutputDirectory(ArchRProj), "GroupCoverages", "batchRegistry")

  #####################################################
  # Batch Apply to Create Insertion Coverage Files
  #####################################################

  #Disable Hdf5 File Locking
  h5disableFileLocking()

  #Batch Apply
  .logDiffTime(sprintf("Creating Coverage Files!"), tstart, addHeader = FALSE)
  batchOut <- .batchlapply(args)
  coverageFiles <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$covFile) %>% unlist
  nCells <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$nCells) %>% unlist
  nFragments <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$nFragments) %>% unlist

  #Add To Project
  coverageMetadata <- DataFrame(
    Group = stringr::str_split(names(unlistGroups), pattern = "\\._.", simplify=TRUE)[,1],
    Name = names(unlistGroups), 
    File = coverageFiles,
    nCells = nCells,
    nInsertions = nFragments * 2
  )

  #####################################################
  # Compute Kmer Bias for each coverage file!
  #####################################################
  .logDiffTime(sprintf("Adding Kmer Bias to Coverage Files!"), tstart, addHeader = FALSE)
  o <- .addKmerBiasToCoverage(
    coverageMetadata = coverageMetadata, 
    genome = getGenome(ArchRProj), 
    kmerLength = kmerLength, 
    threads = threads,
    verbose = FALSE,
    logFile = logFile
  )

  ArchRProj@projectMetadata$GroupCoverages[[groupBy]] <- SimpleList(Params = Params, coverageMetadata = coverageMetadata)

  #Enable Hdf5 File Locking
  h5enableFileLocking()

  .logDiffTime(sprintf("Finished Creation of Coverage Files!"), tstart, addHeader = FALSE)
  .endLogging(logFile = logFile)

  ArchRProj

}

#####################################################################################################
# Creating Insertion (1bp) Coverage Hdf5 Files for downstream group analyses
#####################################################################################################

.createCoverages <- function(
  i = NULL, 
  cellGroups,
  kmerBias = NULL, 
  kmerLength = 6, 
  genome = NULL,
  ArrowFiles = NULL, 
  cellsInArrow = NULL, 
  availableChr = NULL,
  chromLengths = NULL, 
  covDir = NULL, 
  tstart = NULL, 
  subThreads = 1,
  verbose = TRUE,
  logFile = NULL
  ){

  prefix <- sprintf("Group %s (%s of %s) :", names(cellGroups)[i], i, length(cellGroups))

  #Cells
  cellGroupi <- cellGroups[[i]]
  
  #Coverage File!
  namei <- make.names(names(cellGroups)[i]) #Maybe naming convention is weird
  covFile <- file.path(covDir, paste0(namei, ".insertions.coverage.h5"))
  rmf <- .suppressAll(file.remove(covFile))

  .logDiffTime(sprintf("%s Creating Group Coverage File : %s", prefix, basename(covFile)), tstart, verbose = verbose, logFile = logFile)

  #Dealing with sampling w/o replacement!
  tableGroupi <- table(cellGroupi)

  .logThis(cellGroupi, paste0(prefix, " cellGroupi"), logFile = logFile)

  .logMessage(paste0("Number of Cells = ", length(cellGroups[[i]])), logFile = logFile)

  #Create Hdf5 File!
  o <- tryCatch({
    o <- h5closeAll()
    o <- h5createFile(covFile)   
  }, error = function(e){
    rmf <- .suppressAll(file.remove(covFile))
    o <- h5closeAll()
    o <- h5createFile(covFile)
  })

  if(file.exists(covFile)){
    .logMessage("Coverage File Exists!", logFile = logFile)
  }else{
    .logMessage("Coverage File Does Not Exist!", logFile = logFile)
  }

  o <- h5closeAll()
  o <- h5createGroup(covFile, paste0("Coverage"))
  .logMessage("Added Coverage Group", logFile = logFile)

  o <- h5closeAll()
  o <- h5createGroup(covFile, paste0("Metadata"))
  .logMessage("Added Metadata Group", logFile = logFile)

  o <- h5closeAll()
  o <- h5write(obj = "ArrowCoverage", file = covFile, name = "Class")
  .logMessage("Added ArrowCoverage Class", logFile = logFile)

  o <- h5closeAll()
  o <- h5createGroup(covFile, paste0("Coverage/Info"))
  .logMessage("Added Coverage/Info", logFile = logFile)

  o <- h5closeAll()
  o <- h5write(as.character(cellGroupi), covFile, "Coverage/Info/CellNames")
  .logMessage("Added Coverage/Info/CellNames", logFile = logFile)

  #We need to dump all the cells into a coverage file
  nFragDump <- 0
  nCells <- c()
  for(k in seq_along(availableChr)){
    
    .logDiffTime(sprintf("%s Processed Fragments Chr (%s of %s)", prefix, k, length(availableChr)), tstart, verbose = FALSE, logFile = logFile)
    
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

    #Dealing with sampling w/o replacement!
    matchRG <- as.vector(S4Vectors::match(mcols(fragik)$RG, names(tableGroupi)))
    fragik <- rep(fragik, tableGroupi[matchRG])
    nCells <- c(nCells, unique(mcols(fragik)$RG))

    #Compute Rle Coverage
    covk <- coverage(IRanges(start = c( start(fragik), end(fragik) ), width = 1), width = chromLengths[availableChr[k]])
        nFragDump <- nFragDump + length(fragik)
    rm(fragik)

    #Write To Hdf5
    chrLengths <- paste0("Coverage/",availableChr[k],"/Lengths")
    chrValues <- paste0("Coverage/",availableChr[k],"/Values")
    lengthRle <- length(covk@lengths)
    o <- h5createGroup(covFile, paste0("Coverage/",availableChr[k]))
    o <- .suppressAll(h5createDataset(covFile, chrLengths, storage.mode = "integer", dims = c(lengthRle, 1), level = 0))
    o <- .suppressAll(h5createDataset(covFile, chrValues, storage.mode = "integer", dims = c(lengthRle, 1), level = 0))
    o <- h5write(obj = covk@lengths, file = covFile, name = chrLengths)
    o <- h5write(obj = covk@values, file = covFile, name = chrValues)

    gc()

  }

  if(length(unique(cellGroupi)) != length(unique(nCells))){
    .logMessage(paste0("Not all cells (", length(unique(cellGroupi)), ") were found for coverage creation (", length(unique(nCells)), ")!"), logFile = logFile)
    stop("Not all cells (", length(unique(cellGroupi)), ") were found for coverage creation (", length(unique(nCells)), ")!")
  }

  out <- list(covFile = covFile, nCells = length(cellGroupi), nFragments = nFragDump)

  return(out)

}

#####################################################################################################
# Creating Groups of Cells For Pseudobulk Coverage Files
#####################################################################################################

.identifyGroupsForPseudoBulk <- function(
  cells = NULL,
  sampleLabels = NULL,
  useLabels = TRUE,
  minCells = 50,
  maxCells = 500,
  filterGroups = FALSE,
  minReplicates = 2,
  maxReplicates = NULL,
  sampleRatio = 0.8,
  prefix = NULL,
  logFile = NULL
  ){

  if(is.null(sampleLabels)){
    sampleLabels <- rep("A", length(cells))
  }else{
    if(length(cells) != length(sampleLabels)){
      .logMessage("Length of cells need to be same length as sample labels!", logFile = logFile)
      stop("Length of cells need to be same length as sample labels!")
    }
  }
  nCells <- length(cells)
  nCellsPerSample <- table(sampleLabels)
  nCellsPerSample <- nCellsPerSample[sample(seq_along(nCellsPerSample), length(nCellsPerSample))]
  #Samples Passing Min Filter
  samplesPassFilter <- sum(nCellsPerSample >= minCells)   
  samplesThatCouldBeMergedToPass <- floor(sum(nCellsPerSample[nCellsPerSample < minCells]) / minCells)

  errorList <- mget(names(formals()),sys.frame(sys.nframe()))
  errorList$nCells <- nCells
  errorList$nCellsPerSample <- nCellsPerSample
  errorList$samplesPassFilter <- samplesPassFilter
  errorList$samplesThatCouldBeMergedToPass <- samplesThatCouldBeMergedToPass

  cellGroupsPass2 <- tryCatch({

    if(nCells >= minCells * minReplicates & useLabels){
        ############################################################
        # Identifying High-Quality peaks when Cells and Fragments are abundant
        ############################################################
        #Samples Passing Min Filter
        samplesPassFilter <- sum(nCellsPerSample >= minCells) 
        samplesThatCouldBeMergedToPass <- floor(sum(nCellsPerSample[nCellsPerSample < minCells]) / minCells)
        #First Group Cells By Sample
        cellGroups <- split(cells, sampleLabels)
        #Identify Samples That Pass Min Cells
        samples <- names(nCellsPerSample)[nCellsPerSample > minCells]
        cellGroupsPass <- cellGroups[samples]
        #Samples That Do Not Pass
        if(!all(names(cellGroups) %in% names(cellGroupsPass))){
          cellGroupsNotPass <- cellGroups[names(cellGroups) %ni% samples]
        }else{
          cellGroupsNotPass <- list()
        }
        if(samplesPassFilter >= minReplicates){
          ############################################################
          # If there are at least minReplicates with > minCells
          ############################################################
          #If we look at the remaining cells and see that there are enough to make an additional replicate
          nCellsRemaining <- length(unlist(cellGroupsNotPass))
          if(nCellsRemaining >= minCells){
            cellGroupsPass$Other <- unlist(cellGroupsNotPass)
          }
        }else if(samplesPassFilter + samplesThatCouldBeMergedToPass >= minReplicates){
          cellsRemaining <- unlist(cellGroupsNotPass, use.names = FALSE)
          nGroupsRemaining <- minReplicates - samplesPassFilter
          cellsRemaining <- sample(cellsRemaining, length(cellsRemaining))
          cellGroupsSample <- split(cellsRemaining, ceiling(seq_along(cellsRemaining)/ ceiling(length(cellsRemaining)/nGroupsRemaining)))
          #Add to Groups Pass QC
          if(samplesPassFilter == 0){
            cellGroupsPass <- cellGroupsSample
            names(cellGroupsPass) <- paste0("Rep", seq_along(cellGroupsPass))
          }else{
            cellGroupsPass <- append(cellGroupsPass, cellGroupsSample)
          }
        }else{
          cellsRemaining <- unlist(cellGroupsNotPass, use.names = FALSE)
          nCellsRemaining <- length(cellsRemaining)
          cellsNeeded <- minCells * (minReplicates - samplesPassFilter) - nCellsRemaining
          cellsFromPass <- sample(unlist(cellGroupsPass, use.names = FALSE), cellsNeeded)
          cellGroupsPass <- lapply(cellGroupsPass, function(x){
            x[x %ni% cellsFromPass]
          })
          cellGroupsPass[[length(cellGroupsPass) + 1]] <- c(cellsRemaining, cellsFromPass)
          names(cellGroupsPass) <- paste0("Rep", seq_along(cellGroupsPass))
        }
    }else{
        ############################################################
        # Identifying High-Quality peaks when Cells are not abundant
        ############################################################
        if(nCells >= minCells / sampleRatio){
            ############################################################
            # When there are more cells than min cells but not enough for robust reproducibility
            ############################################################
            cellGroupsPass <- .leastOverlapCells(x = cells, n = minReplicates, nSample = minCells) #length(cells) * sampleRatio)
        }else{
            ############################################################
            # Sampling With Replacement, Not Super Desirable
            ############################################################
            if(filterGroups){
              return(NULL)
            }else{            
              cellGroupsPass <- .leastOverlapCells(x = cells, n = minReplicates, nSample = minCells, replace = TRUE) #length(cells) * sampleRatio)
            }
        }
    }
    cellGroupsPass <- as(cellGroupsPass, "SimpleList")

    for(i in seq_along(cellGroupsPass)){
        if(length(cellGroupsPass[[i]]) > maxCells){
            cellGroupsPass[[i]] <- sample(cellGroupsPass[[i]], maxCells)
        }
    }

    #Rank Group by Unique # of Cells
    nUnique <- lapply(cellGroupsPass, function(x){
      length(unique(x))
    }) %>% unlist

    cellGroupsPass <- cellGroupsPass[order(nUnique, decreasing = TRUE)]

    if(!is.null(maxReplicates)){
      if(length(cellGroupsPass) > maxReplicates){
        cellGroupsPass <- cellGroupsPass[seq_len(maxReplicates)]
      }    
    }

    cellGroupsPass

  }, error = function(e){

   .logError(e, fn = ".identifyGroupsForPseudoBulk", info = prefix, errorList = errorList, logFile = logFile) 

  })

  .logMessage(paste0(prefix, " CellGroups N = ", length(cellGroupsPass2)), logFile = logFile)
  #.logThis(cellGroupsPass2, paste0(prefix, " CellGroups"), logFile = logFile)

  return(cellGroupsPass2)
  
}

.leastOverlapCells <- function(x = NULL, n = 2, nSample = 0.8 * length(l), iterations = 100, replace = FALSE){   
    set.seed(1)
    maxMat <- matrix(0, nrow = length(x), ncol = n)
    for(i in seq_len(iterations)){
      currentMat <- matrix(0, nrow = length(x), ncol = n)
      for(j in seq_len(n)){
        currentMat[sample(seq_along(x), nSample, replace = replace), j] <- 1
      }
      disti <- max(dist(t(currentMat), method = "euclidean"))
      if(i==1){
        maxMat <- currentMat
        maxDist <- disti
      }else{
        if(disti > maxDist){
          maxMat <- currentMat
          maxDist <- disti        
        }
      }
    }
    out <- lapply(seq_len(ncol(maxMat)), function(i){
        x[which(maxMat[,i]==1)]
    })
    return(out)
}

#####################################################################################################
# Add Kmer Tn5 Bias Values to Each Coverage File!
#####################################################################################################
.addKmerBiasToCoverage <- function(
  coverageMetadata = NULL,
  genome = NULL,
  kmerLength = NULL,
  threads = NULL,
  verbose = TRUE,
  tstart = NULL,
  logFile = NULL
  ){

  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "kmerBias-Parameters", logFile = logFile)
  
  .requirePackage(genome)
  .requirePackage("Biostrings", source = "bioc")
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  coverageFiles <- coverageMetadata$File
  names(coverageFiles) <- coverageMetadata$Name
  availableChr <- .availableSeqnames(coverageFiles, "Coverage")

  biasList <- .safelapply(seq_along(availableChr), function(x){
    .logMessage(sprintf("Kmer Bias %s (%s of %s)", availableChr[x], x, length(availableChr)), logFile = logFile)
    message(availableChr[x]," ", appendLF = FALSE)
    chrBS <- BSgenome[[availableChr[x]]]
    exp <- Biostrings::oligonucleotideFrequency(chrBS, width = kmerLength)
    obsList <- lapply(seq_along(coverageFiles), function(y){
      .logMessage(sprintf("Coverage File %s (%s of %s)", availableChr[x], y, length(coverageFiles)), logFile = logFile)
      tryCatch({
        obsx <- .getCoverageInsertionSites(coverageFiles[y], availableChr[x]) %>%
          {BSgenome::Views(chrBS, IRanges(start = . - floor(kmerLength/2), width = kmerLength))} %>%
          {Biostrings::oligonucleotideFrequency(., width = kmerLength, simplify.as="collapsed")}
        tryCatch({gc()}, error = function(e){})
        obsx
      }, error = function(e){
        errorList <- list(
          y = y, 
          coverageFile = coverageFiles[y], 
          chr = availableChr[x], 
          iS = tryCatch({
            .getCoverageInsertionSites(coverageFiles[y], availableChr[x])
            }, error = function(e){
              "Error .getCoverageInsertionSites"
            })
        )
        .logError(e, fn = ".addKmerBiasToCoverage", info = "", errorList = errorList, logFile = logFile)
      })
    }) %>% SimpleList
    names(obsList) <- names(coverageFiles)
    SimpleList(expected = exp, observed = obsList)
  }, threads = threads) %>% SimpleList
  names(biasList) <- availableChr
  .logMessage("Completed Kmer Bias Calculation", logFile = logFile)

  #Summarize Bias
  for(i in seq_along(biasList)){
    if(i == 1){
      expAll <- biasList[[i]]$expected
      obsAll <- biasList[[i]]$observed
    }else{
      expAll <- expAll + biasList[[i]]$expected 
      for(j in seq_along(obsAll)){
        obsAll[[j]] <- obsAll[[j]] + biasList[[i]]$observed[[names(obsAll)[j]]]
      }
    }
  }

  #Write Bias to Coverage Files
  for(i in seq_along(coverageFiles)){

    .logMessage(sprintf("Adding Kmer Bias (%s of %s)", i, length(coverageFiles)), logFile = logFile)

    obsAlli <- obsAll[[names(coverageFiles)[i]]]
    if(!identical(names(expAll), names(obsAlli))){
      .logMessage("Kmer Names in Exp and Obs not Identical!", logFile = logFile)
      stop("Kmer Names in Exp and Obs not Identical!")
    }
    o <- h5createGroup(coverageFiles[i], "KmerBias")
    o <- h5createGroup(coverageFiles[i], "KmerBias/Info")
    o <- h5write(obj = genome, file = coverageFiles[i], name = "KmerBias/Info/Genome")
    o <- h5write(obj = kmerLength, file = coverageFiles[i], name = "KmerBias/Info/KmerLength")
    o <- h5write(obj = paste0(names(obsAlli)), file = coverageFiles[i], name = "KmerBias/Kmer")
    o <- h5write(obj = obsAlli, file = coverageFiles[i], name = "KmerBias/ObservedKmers")
    o <- h5write(obj = expAll, file = coverageFiles[i], name = "KmerBias/ExpectedKmers")

  }

  return(0)

}

#####################################################################################################
# Get Coverage Metadata and Params from ArchR Project
#####################################################################################################

.getCoverageMetadata <- function(ArchRProj = NULL, groupBy = NULL, useGroups = NULL, minCells = NULL){
  coverageMetadata <- ArchRProj@projectMetadata$GroupCoverages[[groupBy]]$coverageMetadata
  if(is.null(coverageMetadata)){
    stop("No Coverage Metadata found for : ", groupBy, ". Please run addGroupCoverages!")
  }
  if(!is.null(useGroups)){
    if(sum(coverageMetadata[,1] %in% useGroups) == 0){
      stop("No Groups found matching useGroups!")
    }
    coverageMetadata <- coverageMetadata[coverageMetadata[,1] %in% useGroups, , drop = FALSE]
  }
  if(!is.null(minCells)){
    coverageMetadata <- coverageMetadata[coverageMetadata$nCells >= minCells, , drop = FALSE]
  }
  if(nrow(coverageMetadata)==0){
    stop("No coverages remain!")
  }
  coverageMetadata
}

.getCoverageParams <- function(ArchRProj = NULL, groupBy = NULL, useGroups = NULL){
  coverageParams <- ArchRProj@projectMetadata$GroupCoverages[[groupBy]]$Params
  if(is.null(coverageParams)){
    stop("No Coverage Metadata found for : ", groupBy)
  }
  coverageParams
}

#####################################################################################################
# Create Coverage Rle List of all chr
#####################################################################################################

.getCoverageRle <- function(coverageFile = NULL, allChr = NULL){
  cov <- lapply(seq_along(allChr), function(x){
    Rle(
      lengths = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Lengths")), 
      values = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Values"))
    )
  }) %>% {as(., "RleList")}
  names(cov) <- allChr
  cov
}

#####################################################################################################
# Get All Non-Zero Insertion Sites and N
#####################################################################################################

.getCoverageInsertionSites <- function(coverageFile = NULL, chr = NULL){
    cov <- Rle(
      lengths = h5read(coverageFile, paste0("Coverage/", chr, "/Lengths")), 
      values = h5read(coverageFile, paste0("Coverage/", chr, "/Values"))
    )
    rV <- runValue(cov)
    cov <- ranges(cov)
    mcols(cov)$values <- rV
    cov <- cov[mcols(cov)$values > 0]
    cov <- unlist(start(slidingWindows(rep(cov, mcols(cov)$values), width = 1, step = 1)))
    cov
}

#####################################################################################################
# Write Coverage To Bed File for MACS2
#####################################################################################################

.writeCoverageToBed <- function(coverageFile = NULL, out = NULL, excludeChr = NULL, logFile = NULL){
  rmf <- .suppressAll(file.remove(out))
  allChr <- .availableSeqnames(coverageFile, "Coverage")
  if(!is.null(excludeChr)){
    allChr <- allChr[allChr %ni% excludeChr]
  }
  if(length(allChr)==0){
    stop("No Chromosomes in Coverage after Excluding Chr!")
  }
  ##Note that there was a bug with data.table vs data.frame with stacking
  for(x in seq_along(allChr)){
    o <- tryCatch({
      iS <- .getCoverageInsertionSites(coverageFile = coverageFile, chr = allChr[x])
      if(x == 1) .logThis(iS, "InsertionSites", logFile = logFile)
      iS <- data.table(seqnames = allChr[x], start = iS - 1L, end = iS)
      if(x == 1) .logThis(iS, "InsertionSites-DT", logFile = logFile)
      data.table::fwrite(iS, out, sep = "\t", col.names = FALSE, append = TRUE)
    }, error = function(e){
      errorList <- list(
        x = x, 
        coverageFile = coverageFile, 
        allChr = allChr, 
        iS = if(exists("iS", inherits = FALSE)) iS else "Error with insertion sites!"
      )
      .logError(e, fn = ".writeCoverageToBed", info = "", errorList = errorList, logFile = logFile)
    })
  }
  if(!file.exists(out)){
    stop("Error in writing coverage to bedfile")
  }
  out
}


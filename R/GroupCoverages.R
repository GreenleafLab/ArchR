#' Add Group Coverages to an ArchRProject object
#' 
#' This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates and then merge these replicates into a single insertion coverage file.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping multiple cells together prior to generation of the insertion coverage file.
#' @param useLabels A boolean value indicating whether to use sample labels to create sample-aware subgrouping during as pseudo-bulk replicate generation.
#' @param minCells The minimum number of cells required in a given cell group to permit insertion coverage file generation.
#' @param maxCells The maximum number of cells to use during insertion coverage file generation.
#' @param maxFragments The maximum number of fragments per cell group to use in insertion coverage file generation. This prevents the generation of excessively large files which would negatively impact memory requirements.
#' @param minReplicates The minimum number of pseudo-bulk replicates to be generated.
#' @param maxReplicates The maximum number of pseudo-bulk replicates to be generated.
#' @param sampleRatio The fraction of the total cells that can be sampled to generate any given pseudo-bulk replicate.
#' @param kmerLength The length of the kmer used for estimating Tn5 bias.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value that indicates whether or not to overwrite the relevant data in the `ArchRProject` object if insertion coverage / pseudo-bulk replicate information already exists.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... additional args
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
  threads = 1,
  parallelParam = "mclapply",
  force = FALSE,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  if(verboseAll){
    verboseHeader <- TRUE
  }

  tstart <- Sys.time()
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

  if(!is.null(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])){
    if(!force){
      stop("Group Coverages Already Computed, Set force = TRUE to continue!")
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
      if(length(cellNamesx) < minCells){
        return(NULL)
      }
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
        sampleRatio = sampleRatio
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
    .messageDiffTime(sprintf("Further Sampled %s Groups above the Max Fragments!", it), tstart)
  }

  #####################################################
  # Arguments for Coverages
  #####################################################

  dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages"), showWarnings = FALSE)
  dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages", groupBy), showWarnings = FALSE)

  args <- list()
  args$X <- seq_along(unlist(cellGroups))
  args$FUN <- .createCoverages
  args$cellGroups <- unlist(cellGroups)
  args$genome <- getGenome(ArchRProj)
  args$kmerLength <- kmerLength
  args$ArrowFiles <- getArrowFiles(ArchRProj)
  args$availableChr <- .availableSeqnames(getArrowFiles(ArchRProj))
  args$chromLengths <- getChromLengths(ArchRProj)
  args$cellsInArrow <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  args$covDir <- file.path(getOutputDirectory(ArchRProj), "GroupCoverages", groupBy)
  args$parallelParam <- parallelParam
  args$threads <- threads
  args$verbose <- verboseAll
  args$tstart <- tstart
  args$registryDir <- file.path(getOutputDirectory(ArchRProj), "GroupCoverages", "batchRegistry")

  #####################################################
  # Batch Apply to Create Insertion Coverage Files
  #####################################################

  #Disable Hdf5 File Locking
  h5disableFileLocking()

  #Batch Apply
  .messageDiffTime(sprintf("Creating Coverage Files!"), tstart, addHeader = verboseAll)
  batchOut <- .batchlapply(args)
  coverageFiles <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$covFile) %>% unlist
  nCells <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$nCells) %>% unlist
  nFragments <- lapply(seq_along(batchOut),function(x) batchOut[[x]]$nFragments) %>% unlist

  #Enable Hdf5 File Locking
  h5enableFileLocking()

  #Add To Project
  coverageMetadata <- DataFrame(
    Group = stringr::str_split(names(unlist(cellGroups)), pattern = "\\.", simplify=TRUE)[,1],
    Name = names(unlist(cellGroups)), 
    File = coverageFiles,
    nCells = nCells,
    nInsertions = nFragments * 2
  )

  #####################################################
  # Compute Kmer Bias for each coverage file!
  #####################################################
  .messageDiffTime(sprintf("Adding Kmer Bias to Coverage Files!"), tstart, addHeader = verboseAll)
  o <- .addKmerBiasToCoverage(
    coverageMetadata = coverageMetadata, 
    genome = getGenome(ArchRProj), 
    kmerLength = kmerLength, 
    threads = threads,
    verbose = verboseAll
    )

  ArchRProj@projectMetadata$GroupCoverages[[groupBy]] <- SimpleList(Params = Params, coverageMetadata = coverageMetadata)

  ArchRProj

}

#####################################################################################################
# Creating Insertion (1bp) Coverage Hdf5 Files for downstream group analyses
#####################################################################################################

.createCoverages <- function(
  i, 
  cellGroups,
  kmerBias = NULL, 
  kmerLength = 5, 
  genome,
  ArrowFiles, 
  cellsInArrow, 
  availableChr,
  chromLengths, 
  covDir, 
  tstart, 
  verbose = TRUE,
  ...
  ){

  #Cells
  cellGroupi <- cellGroups[[i]]
  
  #Dealing with sampling w/o replacement!
  tableGroupi <- table(cellGroupi)

  #Coverage File!
  covFile <- file.path(covDir, paste0(names(cellGroups)[i], ".insertions.coverage.h5"))
  rmf <- .suppressAll(file.remove(covFile))

  #Create Hdf5 File!
  o <- h5createFile(covFile)
  o <- h5createGroup(covFile, paste0("Coverage"))
  o <- h5createGroup(covFile, paste0("Metadata"))
  o <- h5write(obj = "ArrowCoverage", file = covFile, name = "Class")

  o <- h5createGroup(covFile, paste0("Coverage/Info"))
  o <- h5write(as.character(cellGroupi), covFile, "Coverage/Info/CellNames")

  #We need to dump all the cells into a coverage file
  nFragDump <- 0
  nCells <- c()
  for(k in seq_along(availableChr)){
    
    if(k %% 3 == 0){
      .messageDiffTime(sprintf("Group %s of %s, Read Fragments %s of %s!", i, 
        length(cellGroups), k, length(availableChr)), tstart, verbose = verbose)
    }
    
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
    nCells <- c(nCells, unique(runValue(mcols(fragik)$RG)))

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
    stop("Not all cells (", length(unique(cellGroupi)), ") were found for coverage creation (", length(unique(nCells)), ")!")
  }

  out <- list(covFile = covFile, nCells = length(cellGroupi), nFragments = nFragDump)

  return(out)

}

#####################################################################################################
# Creating Groups of Cells For Pseudobulk Coverage Files
#####################################################################################################

.identifyGroupsForPseudoBulk <- function(
  cells, sampleLabels = NULL, useLabels = TRUE,
  minCells = 50, maxCells = 500, filterGroups = FALSE,
  minReplicates = 2, maxReplicates = NULL, sampleRatio = 0.8){

    .leastOverlapCells <- function(x, n = 2, nSample = 0.8 * length(l), iterations = 100, replace = FALSE){   
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

    if(is.null(sampleLabels)){
      sampleLabels <- rep("A", length(cells))
    }else{
      if(length(cells) != length(sampleLabels)){
        stop("Length of cells need to be same length as sample labels!")
      }
    }
    nCells <- length(cells)
    nCellsPerSample <- table(sampleLabels)
    nCellsPerSample <- nCellsPerSample[sample(seq_along(nCellsPerSample), length(nCellsPerSample))]
    #Samples Passing Min Filter
    samplesPassFilter <- sum(nCellsPerSample >= minCells)   
    samplesThatCouldBeMergedToPass <- floor(sum(nCellsPerSample[nCellsPerSample < minCells]) / minCells)
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
    if(!is.null(maxReplicates)){
      if(length(cellsGroupPass) > maxReplicates){
        cellsGroupPass <- cellsGroupPass[seq_len(maxReplicates)]
      }    
    }

  return(cellGroupsPass)
  
}

#####################################################################################################
# Add Kmer Tn5 Bias Values to Each Coverage File!
#####################################################################################################
.addKmerBiasToCoverage <- function(coverageMetadata, genome, kmerLength, threads, verbose = TRUE, tstart = NULL){
  
  .requirePackage(genome)
  .requirePackage("Biostrings")
  BSgenome <- eval(parse(text = genome))
  BSgenome <- .validBSgenome(BSgenome)

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  coverageFiles <- coverageMetadata$File
  names(coverageFiles) <- coverageMetadata$Name
  availableChr <- .availableSeqnames(coverageFiles, "Coverage")

  biasList <- .safelapply(seq_along(availableChr), function(x){
    .messageDiffTime(sprintf("Computing Kmer Bias Chr %s of %s!", x, length(availableChr)), tstart, verbose=verbose)
    chrBS <- BSgenome[[availableChr[x]]]
    exp <- Biostrings::oligonucleotideFrequency(chrBS, width = kmerLength)
    obsList <- lapply(seq_along(coverageFiles), function(y){
      obsx <- .getCoverageInsertionSites(coverageFiles[y], availableChr[x]) %>%
        {BSgenome::Views(chrBS, IRanges(start = . - floor(kmerLength/2), width = kmerLength))} %>%
        {Biostrings::oligonucleotideFrequency(., width = kmerLength, simplify.as="collapsed")}
        gc()
        obsx
    }) %>% SimpleList
    names(obsList) <- names(coverageFiles)
    SimpleList(expected = exp, observed = obsList)
  }, threads = threads) %>% SimpleList
  names(biasList) <- availableChr

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
    obsAlli <- obsAll[[names(coverageFiles)[i]]]
    if(!identical(names(expAll), names(obsAlli))){
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

.getCoverageMetadata <- function(ArchRProj, groupBy, useGroups = NULL, minCells = NULL){
  coverageMetadata <- ArchRProj@projectMetadata$GroupCoverages[[groupBy]]$coverageMetadata
  if(is.null(coverageMetadata)){
    stop("No Coverage Metadata found for : ", groupBy)
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

.getCoverageParams <- function(ArchRProj, groupBy, useGroups = NULL){
  coverageParams <- ArchRProj@projectMetadata$GroupCoverages[[groupBy]]$Params
  if(is.null(coverageParams)){
    stop("No Coverage Metadata found for : ", groupBy)
  }
  coverageParams
}

#####################################################################################################
# Create Coverage Rle List of all chr
#####################################################################################################

.getCoverageRle <- function(coverageFile, allChr){
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

.getCoverageInsertionSites <- function(coverageFile, chr){
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

.writeCoverageToBed <- function(coverageFile, out, excludeChr = NULL){
  rmf <- .suppressAll(file.remove(out))
  allChr <- .availableSeqnames(coverageFile, "Coverage")
  if(!is.null(excludeChr)){
    allChr <- allChr[allChr %ni% excludeChr]
  }
  if(length(allChr)==0){
    stop("No Chromosomes in Coverage after Excluding Chr!")
  }
  for(x in seq_along(allChr)){
    .getCoverageInsertionSites(coverageFile, allChr[x]) %>% 
      {data.frame(seqnames = allChr[x], start = . - 1L, end = .)} %>% 
      {data.table::fwrite(., out, sep = "\t", col.names = FALSE, append = TRUE)}
  }
  out
}




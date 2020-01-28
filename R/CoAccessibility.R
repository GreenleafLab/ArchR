##########################################################################################
# Co-accessibility Methods
##########################################################################################

#' Add Peak Co-Accessibility to an ArchRProject
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo A numeric value indicating what to depth-normalize the accessibility of a single-cell group prior to computing co-accessibility correlations.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in cluster determination. It is recommended to keep track of the seed used so that you can reproduce results downstream.
#' @param knnMethod The method to be used for k-nearest neighbor computations. Options are "nabor", "RANN", and "FNN" and the corresponding package is required.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... additional args
#' @export
addCoAccessibility <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:25,
  corCutOff = 0.75,
  k = 200, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 100000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1, 
  knnMethod = NULL,
  threads = getArchRThreads(),
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = knnMethod, name = "knnMethod", valid = c("character", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

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
  peakWindows <- resize(peakSummits, maxDist, "center")

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

#' Get the peak coAccessibility from an ArchRProject
#' 
#' This function gets the peak set as a GRanges object from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param corCutOff A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to regurn CoAccessibility as loops designed for `ArchRBrowser` or `ArchRRegionTrack`.
#' @export
getCoAccessibility <- function(ArchRProj = NULL, corCutOff = 0.5, resolution = 1000, returnLoops = TRUE){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = "numeric")
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){

    return(NULL)
  }

  
  if(is.null(metadata(ArchRProj@peakSet)$CoAccessibility)){
  
    return(NULL)
  
  }else{
   
    coA <- metadata(ArchRProj@peakSet)$CoAccessibility
    coA <- coA[coA$correlation >= corCutOff,,drop=FALSE]

    if(returnLoops){
      
      peakSummits <- resize(metadata(coA)$peakSet, 1, "center")
      summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
    
      loops <- constructGR(
        seqnames = seqnames(peakSummits[coA[,1]]),
        start = summitTiles[coA[,1]],
        end = summitTiles[coA[,2]]
      )
      mcols(loops)$value <- coA$correlation

      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) >= resolution]
      loops <- sort(sortSeqlevels(loops))

      loops <- GenomicRangesList(CoAccessibility = loops)

      return(loops)

    }else{

      return(coA)

    }

  }

}


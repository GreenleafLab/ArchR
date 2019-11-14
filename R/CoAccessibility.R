#' Add Peak Co-Accessibility to ArchR Project
#' 
#' This function will randomly group cells and compute correlations of knn groupings
#'
#' @param ArchRProj ArchRProject
#' @param reducedDims reduced dimensions for KNN groupings
#' @param k k-nearest neighbors
#' @param knnIteration number of KNN groupings to test overlapCutoff
#' @param overlapCutoff overlap maximum between group and previous groups to be added to group list
#' @param maxDist maximum distance in bp between peaks for co-accessibility
#' @param scaleTo scale group accessibility to prior to computing correlations
#' @param log2Norm log2 normalize prior to computing correlations
#' @param seed seed for sampling
#' @param knnMethod method for KNN computations
#' @param threads number of threads
#' @param ... additional args
#' @export
addCoAccessibility <- function(
  ArchRProj,
  reducedDims = "IterativeLSI",
  k = 100, 
  knnIteration = 5000, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1, 
  knnMethod = "nabor", 
  ...
  ){

  tstart <- Sys.time()

  set.seed(seed)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims)

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = ifelse(nrow(rD) >= knnIteration, FALSE, TRUE))

  #KNN Matrix
  .messageDiffTime("Computing KNN", tstart)
  knnObj <- computeKNN(data = rD, query = rD[idx,], k = k, method = knnMethod)

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
  peakSet <- getPeakSet(ArchRProj)
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
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]])))

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
  o@metadata$peakSet <- peakSet

  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  ArchRProj

}







#' Extend Filter then Normalize Scores for Summits
#' @export
computeCoAccessibility <- function(
  ArchRProject,
  k = 100, 
  knnIteration = 5000, 
  overlapCutoff = 0.8, 
  seed = 1, 
  knnMethod = "nabor", 
  maxDist = 250000,
  ...){

  set.seed(seed)

  #Get Matrix List
  matrixDF <- validInputMatrix(proj, useMatrix = "PeakMatrix")

  #Check for existence
  stopifnot(all(apply(matrixDF, 2, function(x) all(file.exists(paste0(x))))))

  #Subsample
  idx <- sample(seq_len(nrow(mat)), knnIteration)

  #KNN Matrix
  knnMat <- computeKnn(data = mat, query = mat[idx,], method = knnMethod)

  #Determin Overlap
  keepKnn <- determineOverlapCpp(k1, floor(cutOffOverlap*k))

  #Keep Above Cutoff
  knnMat <- knnMat[keepKnn==0,]
  message("Identified ", nrow(knnMat), " Groupings!")

  # #Time to compute partial group matrix for coAccessibility!
  # featureDF <- DataFrame()
  #   groupMat <- constructGroupMatrix(
  #     inputFilesMatrix = matrixDF, 
  #     featureDF = , 
  #     groupList = groupList
  #   )

  #Work in progress!
  keepKnn

}
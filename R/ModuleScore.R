#' Add Module Scores to an ArchRProject
#' 
#' This function computes imputations weights that describe each cell as a linear combination of many cells based on a MAGIC diffusion matrix.
#'
#' RRR
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addModuleScore <- function(
  ArchRProj = NULL,
  useMatrix = NULL,
  name = "Module",
  features = NULL,
  nBin = 25,
  nBgd = 100,
  seed = 1,
  threads = getArchRThreads(),
  logFile = createLogFile("addModuleScore")
  ){

  if(!is.null(seed)) set.seed(seed)

  #Get Feature DF
  featureDF <- ArchR:::.getFeatureDF(head(getArrowFiles(ArchRProj),2), subGroup=useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, ":", featureDF$idx)
    featureDF$Match <- seq_len(nrow(featureDF))

    if(useMatrix %ni% getAvailableMatrices(ArchRProj)){
      stop("useMatrix not in available matrices! See getAvailableMatrices!")
    }

  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    if(!all(unlist(lapply(unlist(features), function(x) grepl(":",x))))){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!", logFile = logFile)
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }

  if(grepl(":",unlist(features)[1])){

    sname <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,1]
    name <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,2]

    idx <- lapply(seq_along(name), function(x){
      ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }else{

    idx <- lapply(seq_along(unlist(features)), function(x){
      ix <- which(tolower(unlist(features)[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", unlist(features)[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }

  if(is.null(names(features))){
    names(features) <- paste0(name, seq_along(features))
  }else{
    names(features) <- paste0(name, ".", names(features))
  }

  featuresUse <- featureDF[idx,]
  featuresUse$Module <- Rle(stack(features)[,2])

  #Get Averages
  rS <- ArchR:::.getRowSums(ArrowFiles = getArrowFiles(ArchRProj), useMatrix = useMatrix)
  rS <- rS[order(rS[,3]), ]
  rS$Bins <- Rle(ggplot2::cut_number(x = rS[,3] + rnorm(length(rS[,3]))/1e30, n = nBin, labels = FALSE, right = FALSE))
  rS$Match <- match(paste0(rS$seqnames, ":", rS$idx), rownames(featureDF))
  
  if(nBgd > min(rS$Bins@lengths)){
    stop("nBgd must be lower than ", min(rS$Bins@lengths), "!")
  }

  idxMatch <- match(paste0(featuresUse$seqnames, ":", featuresUse$idx), paste0(rS$seqnames, ":", rS$idx))
  featuresUse$Bins <- as.vector(rS$Bins[idxMatch])
  
  #MakeLists
  featureList <- split(featuresUse$Match, featuresUse$Module)
  moduleList <- split(featuresUse$Bins, featuresUse$Module)
  binList <- split(rS$Match, rS$Bins)

  dfM <- lapply(seq_along(featureList), function(x){
    message("Computing Module ",x, " of ", length(featureList))
    binx <- binList[moduleList[[x]]]
    idxFgd <- featureList[[x]]
    idxBgd <- unlist(lapply(binx, function(x) sample(x, nBgd)), use.names=FALSE)
    m <- ArchR:::.getPartialMatrix(
      ArrowFiles = getArrowFiles(ArchRProj),
      featureDF = featureDF[c(idxFgd, idxBgd), ],
      useMatrix = useMatrix,
      cellNames = ArchRProj$cellNames,
      threads = threads,
      verbose = FALSE,
      doSampleCells = FALSE
    )
    Matrix::colMeans(m[seq_along(idxFgd), ]) - Matrix::colMeans(m[-seq_along(idxFgd), ])
  }) %>% Reduce("cbind", .)

  for(x in seq_len(ncol(dfM))){
    ArchRProj <- addCellColData(ArchRProj, data = dfM[,x], name=names(featureList)[x], cells=rownames(dfM), force = TRUE)
  }

  ArchRProj

}



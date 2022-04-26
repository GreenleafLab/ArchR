#' Add Module Scores to an ArchRProject
#' 
#' This function calculates a module score from a set of features across all cells. This allows for
#' grouping of multiple features together into a single quantitative measurement. Currently, this
#' function only works for modules derived from the `GeneScoreMatrix`. Each module is added as a
#' new column in `cellColData`
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the matrix to be used for calculation of the module score. See `getAvailableMatrices()` to view available options.
#' @param name The name to be given to the designated module. If `features` is a list, this name will be prepended to the feature set names given in the list as shown below.
#' @param features A list of feature names to be grouped into modules. For example, `list(BScore = c("MS4A1", "CD79A", "CD74"),	TScore = c("CD3D", "CD8A", "GZMB", "CCR7", "LEF1"))`.
#' Each named element in this list will be stored as a separate module. The examples given in these parameters would yield two modules called `Module.Bscore` and `Module.Tscore`.
#' If the elements of this list are not named, they will be numbered in order, i.e. `Module1`, `Module2`.
#' @param nBin The number of bins to use to divide all features for identification of signal-matched features for background calculation
#' @param nBgd The number of background features to use for signal normalization.
#' @param seed A number to be used as the seed for random number generation required when sampling cells for the background set. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
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

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = features, name = "features", valid = c("list"))
  .validInput(input = nBin, name = "nBin", valid = c("integer"))
  .validInput(input = nBgd, name = "nBgd", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer","null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character", "null"))
  
  if(useMatrix %ni% getAvailableMatrices(ArchRProj)){
      stop("useMatrix not in available matrices! See getAvailableMatrices!")
  }
  
  if(!is.null(seed)) set.seed(seed)

  #Get Feature DF
  featureDF <- ArchR:::.getFeatureDF(head(getArrowFiles(ArchRProj),2), subGroup=useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, ":", featureDF$idx)
    featureDF$Match <- seq_len(nrow(featureDF))

  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    if(!all(unlist(lapply(unlist(features), function(x) grepl(":",x))))){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!", logFile = logFile)
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }

  #Figure out the index numbers of the selected features within the given matrix
  if(grepl(":",unlist(features)[1])){

    sname <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,1]
    name <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,2]

    idx <- lapply(seq_along(name), function(x){
      ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See available features using getFeatures()", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }else{

    idx <- lapply(seq_along(unlist(features)), function(x){
      ix <- which(tolower(unlist(features)[x]) == tolower(featureDF$name))[1]
      if(is.na(ix)){
        .logStop(sprintf("FeatureName (%s) does not exist! See available features using getFeatures()", unlist(features)[x]), logFile = logFile)
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

  #Get average values for all features and then order the features based on their average values
  #so that the features can be binned into nBins
  rS <- ArchR:::.getRowSums(ArrowFiles = getArrowFiles(ArchRProj), useMatrix = useMatrix)
  rS <- rS[order(rS[,3]), ]
  rS$Bins <- Rle(ggplot2::cut_number(x = rS[,3] + rnorm(length(rS[,3]))/1e30, n = nBin, labels = FALSE, right = FALSE))
  rS$Match <- match(paste0(rS$seqnames, ":", rS$idx), rownames(featureDF))
  
  #check that the number of selected background features isnt bigger than the size of each bin
  if(nBgd > min(rS$Bins@lengths)){
    stop("nBgd must be lower than ", min(rS$Bins@lengths), "!")
  }

  #Match the indicies across the different vectors
  idxMatch <- match(paste0(featuresUse$seqnames, ":", featuresUse$idx), paste0(rS$seqnames, ":", rS$idx))
  featuresUse$Bins <- as.vector(rS$Bins[idxMatch])
  
  #Make lists
  featureList <- split(featuresUse$Match, featuresUse$Module) #feature indicies per module
  moduleList <- split(featuresUse$Bins, featuresUse$Module) #bins for each feature per module
  binList <- split(rS$Match, rS$Bins) #list of all indicies for each bin

  #calculate the module score by normalizing to a background set of features
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
  })
  
  if (length(features) > 1) {
    dfM <- Reduce("cbind", dfM)
  } else {
    dfM <- as.data.frame(dfM[[1]], row.names = names(dfM), drop = FALSE)
  }
  
  #add the module scores as new columns in cellColData
  for(x in seq_len(ncol(dfM))){
    ArchRProj <- addCellColData(ArchRProj, data = dfM[,x], name=names(featureList)[x], cells=rownames(dfM), force = TRUE)
  }

  ArchRProj

}



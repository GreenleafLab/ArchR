#' Plot a Heatmap of Features across a Trajectory
#' 
#' This function will plot a heatmap of the results from getTrajectory
#' 
#' @param sePositiveTF A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param varCutOff The "Variance Quantile Cutoff" to be used for identifying the top variable features across the given trajectory.
#' Only features with a variance above the provided quantile will be retained.
#' @param maxFeatures The maximum number of features, ordered by variance, to consider from `useMatrix` when generating a trajectory.
#' This prevents smoothing a large number number of features which can be very time consuming.
#' @param scaleRows A boolean value that indicates whether row-wise z-scores should be computed on the matrix provided by `sePositiveTF`.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies
#' rownames from `sePositiveTF` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param labelMarkers A character vector listing the `rownames` of `sePositiveTF` that should be labeled on the side of the heatmap.
#' @param labelTop A number indicating how many of the top N features, based on variance, in `sePositiveTF` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param labelCols A boolean value that indicates whether all columns should be labeled at the bottom of the heatmap.
#' @param rowOrder If wanting to set the order of rows to be plotted, the indices (integer or character correpsonding 
#' to rownmaes) can be provided here.
#' @param useSeqnames A character vector that indicates which `seqnames` should be plotted in the heatmap. Features from
#' `seqnames` that are not listed will be ignored. In the context of a `Sparse.Assays.Matrix`, such as a matrix containing chromVAR
#' deviations, the `seqnames` do not correspond to chromosomes, rather they correspond to the sub-portions of the matrix, for example
#' raw deviations ("deviations") or deviation z-scores ("z") for a chromVAR deviations matrix.
#' @param returnMatrix A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param force If useSeqnames is longer than 1 if matrixClass is "Sparse.Assays.Matrix" to continue. This is not recommended because these matrices
#' can be in different units.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotPositiveTFHeatmap <- function(
    sePositiveTF = NULL,
    varCutOff = 0.9,
    maxFeatures = 25000,
    scaleRows = TRUE,
    limits = c(-1.5, 1.5),
    grepExclude = NULL,
    pal = NULL,
    labelMarkers = NULL,
    labelTop = 50,
    labelRows = FALSE,
    labelCols = FALSE,
    rowOrder = NULL, 
    useSeqnames = NULL,
    returnMatrix = FALSE,
    force = FALSE,
    logFile = createLogFile("plotPositiveTFHeatmap")
){
  
  ArchR:::.validInput(input = sePositiveTF, name = "sePositiveTF", valid = c("SummarizedExperiment"))
  ArchR:::.validInput(input = varCutOff, name = "varCutOff", valid = c("numeric", "null"))
  ArchR:::.validInput(input = maxFeatures, name = "maxFeatures", valid = c("integer", "null"))
  ArchR:::.validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
  ArchR:::.validInput(input = limits, name = "limits", valid = c("numeric"))
  ArchR:::.validInput(input = grepExclude, name = "grepExclude", valid = c("character", "null"))
  ArchR:::.validInput(input = pal, name = "pal", valid = c("palette", "null"))
  ArchR:::.validInput(input = labelMarkers, name = "labelMarkers", valid = c("character", "null"))
  ArchR:::.validInput(input = labelTop, name = "labelTop", valid = c("integer"))
  ArchR:::.validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
  ArchR:::.validInput(input = labelCols, name = "labelCols", valid = c("boolean"))
  ArchR:::.validInput(input = rowOrder, name = "rowOrder", valid = c("vector", "null"))
  ArchR:::.validInput(input = useSeqnames, name = "useSeqnames", valid = c("character", "null"))
  ArchR:::.validInput(input = returnMatrix, name = "returnMatrix", valid = c("boolean"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotTrajectoryHeatmap Input-Parameters", logFile = logFile)

  
  if(!is.null(useSeqnames)){
    sePositiveTF <- sePositiveTF[paste0(rowData(sePositiveTF)$seqnames) %in% paste0(useSeqnames), ]
  }
  
  if(nrow(sePositiveTF) == 0){
    ArchR:::.logStop("No features left in sePositiveTF, please check input!", logFile = logFile)
  }
  
  mat <- assay(sePositiveTF)
  
  if(!is.null(grepExclude)){
    idxExclude <- grep(grepExclude, rownames(mat))
    if(length(idxExclude) > 0){
      mat <- mat[-grep(grepExclude, rownames(mat)), , drop = FALSE]
    }
  }
  
  #Rows with NA
  rSNA <- rowSums(is.na(mat))
  if(sum(rSNA > 0) > 0){
    ArchR:::.logMessage("Removing rows with NA values...", verbose = TRUE, logFile = logFile)
    mat <- mat[rSNA == 0, ]#Remove NA Rows
  }
  ArchR:::.logThis(mat, "mat-pre", logFile = logFile)
  varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))
  ArchR:::.logThis(varQ, "varQ", logFile = logFile)
  orderedVar <- FALSE
  if(is.null(rowOrder)){
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if(is.null(varCutOff) & is.null(maxFeatures)){
      n <- nrow(mat)
    }else if(is.null(varCutOff)){
      n <- maxFeatures
    }else if(is.null(maxFeatures)){
      n <- (1-varCutOff) * nrow(mat)
    }else{
      n <- min((1-varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n ),]
  }
  ArchR:::.logThis(mat, "mat-post", logFile = logFile)
  
  #rownames(mat) <- rowData(sePositiveTF)$name
  
  if(!is.null(labelTop)){
    if(orderedVar){
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    }else{
      idxLabel <- rownames(mat)[order(varQ,decreasing=TRUE)][seq_len(labelTop)]
    }
  }else{
    idxLabel <- NULL
  }
  ArchR:::.logThis(idxLabel, "idxLabel", logFile = logFile)
  
  if(!is.null(labelMarkers)){
    idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
    idxLabel2 <- idxLabel2[idxLabel2 > 0]
  }else{
    idxLabel2 <- NULL
  }
  ArchR:::.logThis(idxLabel2, "idxLabel2", logFile = logFile)
  
  idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])
  ArchR:::.logThis(idxLabel, "idxLabel", logFile = logFile)
  
  if(scaleRows){
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    ArchR:::.logThis(mat, "mat-zscores", logFile = logFile)
  }
  
  if(nrow(mat) == 0){
    stop("No Features Remaining!")
  }
  
  if(is.null(pal)){
    if(is.null(metadata(sePositiveTF)$Params$useMatrix)){
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }else if(tolower(metadata(sePositiveTF)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blueYellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }
  }
  
  if(!is.null(rowOrder)){
    idx <- rowOrder
  }else{
    idx <- order(apply(mat, 1, which.max))
  }
  ArchR:::.logThis(idx, "idx", logFile = logFile)
  
  ht <- tryCatch({
    
    ArchR:::.ArchRHeatmap(
      mat = mat[idx, ],
      scale = FALSE,
      limits = c(min(mat), max(mat)),
      color = pal, 
      clusterCols = FALSE, 
      clusterRows = FALSE,
      labelRows = labelRows,
      labelCols = labelCols,
      customRowLabel = match(idxLabel, rownames(mat[idx,])),
      showColDendrogram = TRUE,
      name = names(sePositiveTF@assays),
      draw = FALSE
    )
    
  }, error = function(e){
    
    errorList = list(
      mat = mat[idx, ],
      scale = FALSE,
      limits = c(min(mat), max(mat)),
      color = pal, 
      clusterCols = FALSE, 
      clusterRows = FALSE,
      labelRows = labelRows,
      labelCols = labelCols,
      customRowLabel = match(idxLabel, rownames(mat[idx,])),
      showColDendrogram = TRUE,
      name = metadata(sePositiveTF)$Params$useMatrix,
      draw = FALSE
    )
    
    ArchR:::.logError(e, fn = ".ArchRHeatmap", info = "", errorList = errorList, logFile = logFile)      
    
  })
  
  ArchR:::.endLogging(logFile = logFile)
  
  if(returnMatrix){
    return(mat[idx, ])
  }else{
    return(ht)
  }
  
}
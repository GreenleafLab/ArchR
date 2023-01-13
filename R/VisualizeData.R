####################################################################
# Visualization Methods
####################################################################

#' Visualize an Embedding from ArchR Project
#' 
#' This function will plot an embedding stored in an ArchRProject
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the embedding stored in the `ArchRProject` to be plotted. See `computeEmbedding()` for more information.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData") or by
#' a data matrix in the corresponding ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName/rowname of the data matrix to be used for plotting. 
#' For example if colorBy is "cellColData" then `name` refers to a column name in the `cellcoldata` (see `getCellcoldata()`). If `colorBy`
#' is "GeneScoreMatrix" then `name` refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values (if continuous) in plotting.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells values.
#' See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring cells. Typically created using `paletteDiscrete()` or `paletteContinuous()`.
#' To make a custom palette, you must construct this following strict specifications. If the coloring is for discrete data (i.e. "Clusters"),
#' then this palette must be a named vector of colors where each color is named for the corresponding group (e.g. `"C1" = "#F97070"`). If the coloring
#' for continuous data, then it just needs to be a vector of colors. If you are using `pal` in conjuction with `highlightCells`, your palette
#' must be a named vector with two entries, one named for the value of the cells in the `name` column of `cellColData` and the other named
#' "Non.Highlighted". For example, `pal=c("Mono" = "green", "Non.Highlighted" = "lightgrey")` would be used to change the color of cells with the value
#' "Mono" in the `cellColData` column indicated by `name`. Because of this, the cells indicated by `highlightCells` must also match this value in the `name` column.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param sampleCells A numeric describing number of cells to use for plot. If using impute weights, this will occur after imputation.
#' @param highlightCells A character vector of cellNames describing which cells to hightlight if using `plotAs = "points"` (default if discrete). 
#' The remainder of cells will be colored light gray.
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels, just the
#' internal portions of the plot.
#' @param quantCut If this is not `NULL`, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values. 
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the lower threshold and y is 
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the 2.5th percentile and 97.5 percentile of values and set
#' values below/above to the value of the 2.5th and 97.5th percentile values respectively.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a discrete color set is desired.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a continuous color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster being
#' uniformly present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x- and y-axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted. By default
#' if `colorBy` is numeric, then `plotAs` is set to "hex".
#' @param Shiny A boolean value that tells the function is calling for Shiny or not.                              
#' @param threads The number of threads to be used for parallel computing.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' 
#' @examples
#'
#' #Get Test Project
#' proj <- getTestProject()
#' 
#' #Plot UMAP
#' p <- plotEmbedding(proj, name = "Clusters")
#' 
#' #PDF
#' plotPDF(p, name = "UMAP-Clusters", ArchRProj = proj)
#'
#' @export
plotEmbedding <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "Sample",
  log2Norm = NULL,
  imputeWeights = if(!grepl("coldata",tolower(colorBy[1]))) getImputeWeights(ArchRProj),
  pal = NULL,
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  Shiny = FALSE,
  matrices = NULL,
  imputeMatrices = NULL,
  embeddingDF = NULL,
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding"),
  ...
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = embedding, name = "reducedDims", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = pal, name = "pal", valid = c("palette", "null"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = sampleCells, name = "sampleCells", valid = c("numeric", "null"))
  .validInput(input = highlightCells, name = "highlightCells", valid = c("character", "null"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = quantCut, name = "quantCut", valid = c("numeric", "null"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character", "null"))
  .validInput(input = continuousSet, name = "continuousSet", valid = c("character", "null"))
  .validInput(input = randomize, name = "randomize", valid = c("boolean"))
  .validInput(input = keepAxis, name = "keepAxis", valid = c("boolean"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character", "null"))
  .validInput(input = Shiny, name = "Shiny", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .requirePackage("ggplot2", source = "cran")
  
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)
  
  ##############################
  # Get Embedding
  ##############################
  .logMessage("Getting Embedding", logFile = logFile)
  if(Shiny){
     df <- embeddingDF
  } else{
    df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  }

  if(!all(rownames(df) %in% ArchRProj$cellNames)){
    stop("Not all cells in embedding are present in ArchRProject!")
  }
  .logThis(df, name = "Embedding data.frame", logFile = logFile)
  
  if(!is.null(sampleCells)){
    if(sampleCells < nrow(df)){
      if(!is.null(imputeWeights)){
        stop("Cannot sampleCells with imputeWeights not equal to NULL at this time!")
      }
      df <- df[sort(sample(seq_len(nrow(df)), sampleCells)), , drop = FALSE]
    }
  }
  
  #Parameters
  plotParams <- list()
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize
  
  #Additional Params!
  plotParams$xlabel <- gsub("_", " ",stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,2])
  plotParams$ylabel <- gsub("_", " ",stringr::str_split(colnames(df)[2],pattern="#",simplify=TRUE)[,2])
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  
  #Check if Cells To Be Highlighted
  if(!is.null(highlightCells)){
    highlightPoints <- match(highlightCells, rownames(df), nomatch = 0)
    if(any(highlightPoints==0)){
      stop("highlightCells contain cells not in Embedding cellNames! Please make sure that these match!")
    }
  }
  
  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }

  if(!Shiny){
    allColorBy <- .availableArrays(head(getArrowFiles(ArchRProj), 2))
  } else {
    allColorBy <-  matrices$allColorBy
  }
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]
  .logMessage(paste0("ColorBy = ", colorBy), logFile = logFile)
  
  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    colorList <- lapply(seq_along(name), function(x){
      colorParams <- list()
      colorParams$color <- as.vector(getCellColData(ArchRProj, select = name[x], drop = FALSE)[rownames(df), 1])
      colorParams$discrete <- .isDiscrete(colorParams$color)
      colorParams$continuousSet <- "solarExtra"
      colorParams$discreteSet <- "stallion"
      colorParams$title <- paste(plotParams$title, " colored by\ncolData : ", name[x])
      if(!is.null(continuousSet)){
        colorParams$continuousSet <- continuousSet
      }
      if(!is.null(discreteSet)){
        colorParams$discreteSet <- discreteSet
      }
      if(x == 1){
        .logThis(colorParams, name = "ColorParams 1", logFile = logFile)
      }
      
      if(!is.null(imputeWeights)){
        if(getArchRVerbose()) message("Imputing Matrix")
        colorMat <- matrix(colorParams$color, nrow=1)
        colnames(colorMat) <- rownames(df)
        colorMat <- imputeMatrix(mat = colorMat, imputeWeights = imputeWeights, logFile = logFile)
        colorParams$color <- as.vector(colorMat)
      }
      colorParams
    }) 
  }else{# plotting embedding for matrix instead of col in cellcoldata
    suppressMessages(message(logFile))
    
    if(!Shiny){
      units <- tryCatch({
      .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
    },error=function(e){
      "values"
    })
    }else{
        units <- ArchRProj@projectMetadata[["units"]]
    }
    
    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
    }
    
    if(is.null(log2Norm)){
      log2Norm <- FALSE
    }
    
    if(!Shiny){
      colorMat <- .getMatrixValues(
      ArchRProj = ArchRProj, 
      name = name, 
      matrixName = colorBy, 
      log2Norm = FALSE, 
      threads = threads,
      logFile = logFile
    )
    }else{ 
      #get values from pre-saved list
      colorMat = tryCatch({
      t(as.matrix(matrices[[colorBy]][name,]))
      }, warning = function(warning_condition) {
      message(paste("name doesn't exist:", name))
      message(warning_condition)
      return(NULL)    
      }, error = function(error_condition) {
      message(paste("name doesn't exist:", name))
      message(error_condition)
      return(NA)
      }, finally={
      })
     rownames(colorMat)=name
    }
    
    if(!all(rownames(df) %in% colnames(colorMat))){
      .logMessage("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.", logFile = logFile)
      stop("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.")
    }
    
    colorMat <- colorMat[,rownames(df), drop=FALSE]
    
    .logThis(colorMat, "colorMat-Before-Impute", logFile = logFile)
    
    if(!is.null(imputeWeights)){
      if(getArchRVerbose()) message("Imputing Matrix")
        if(!Shiny){
          colorMat <- imputeMatrix(mat = as.matrix(colorMat), imputeWeights = imputeWeights, logFile = logFile)
        }else{
          colorMat <- imputeMatrices[[colorBy]][name,]
        }
        if(!inherits(colorMat, "matrix")){
          colorMat <- matrix(colorMat, ncol = nrow(df))
          colnames(colorMat) <- rownames(df)
        }     
    }
    
    .logThis(colorMat, "colorMat-After-Impute", logFile = logFile)
    
    colorList <- lapply(seq_len(nrow(colorMat)), function(x){
      colorParams <- list()
      colorParams$color <- colorMat[x, ]
      colorParams$discrete <- FALSE
      colorParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name[x])
      if(tolower(colorBy) == "genescorematrix"){
        colorParams$continuousSet <- "horizonExtra"
      }else{
        colorParams$continuousSet <- "solarExtra"
      }
      if(!is.null(continuousSet)){
        colorParams$continuousSet <- continuousSet
      }
      if(!is.null(discreteSet)){
        colorParams$discreteSet <- discreteSet
      }
      if(x == 1){
        .logThis(colorParams, name = "ColorParams 1", logFile = logFile)
      }
      colorParams
    })
  }
  
  if(getArchRVerbose()) message("Plotting Embedding")
  
  ggList <- lapply(seq_along(colorList), function(x){
    
    if(getArchRVerbose()) message(x, " ", appendLF = FALSE)
    
    plotParamsx <- .mergeParams(colorList[[x]], plotParams)
    
    if(plotParamsx$discrete){
      plotParamsx$color <- paste0(plotParamsx$color)
    }
    
    if(!plotParamsx$discrete){
      
      if(!is.null(quantCut)){
        plotParamsx$color <- .quantileCut(plotParamsx$color, min(quantCut), max(quantCut))
      }
      
      plotParamsx$pal <- paletteContinuous(set = plotParamsx$continuousSet)
      
      if(!is.null(pal)){
        
        plotParamsx$pal <- pal
        
      }
      
      if(is.null(plotAs)){
        plotAs <- "hexplot"
      }
      
      if(!is.null(log2Norm)){
        if(log2Norm){
          plotParamsx$color <- log2(plotParamsx$color + 1)
          plotParamsx$colorTitle <- paste0("Log2(",units," + 1)")
        }else{
          plotParamsx$colorTitle <- units
        }
      }
      
      if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
        
        plotParamsx$discrete <- NULL
        plotParamsx$continuousSet <- NULL
        plotParamsx$rastr <- NULL
        plotParamsx$size <- NULL
        plotParamsx$randomize <- NULL
        
        .logThis(plotParamsx, name = paste0("PlotParams-", x), logFile = logFile)
        gg <- do.call(ggHex, plotParamsx)
        
      }else{
        
        if(!is.null(highlightCells)){
          plotParamsx$highlightPoints <- highlightPoints
        }
        
        .logThis(plotParamsx, name = paste0("PlotParams-", x), logFile = logFile)
        gg <- do.call(ggPoint, plotParamsx)
        
      }
      
    }else{
      
      if(!is.null(pal)){
        plotParamsx$pal <- pal
      }
      
      if(!is.null(highlightCells)){
        plotParamsx$highlightPoints <- highlightPoints
      }
      
      .logThis(plotParamsx, name = paste0("PlotParams-", x), logFile = logFile)
      gg <- do.call(ggPoint, plotParamsx)
      
    }
    
    if(!keepAxis){
      gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }
    
    gg
    
  })
  names(ggList) <- name
  if(getArchRVerbose()) message("")
  
  if(length(ggList) == 1){
    ggList <- ggList[[1]]
  }
  
  .endLogging(logFile = logFile)
  
  ggList
  
}



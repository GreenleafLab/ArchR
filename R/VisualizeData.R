####################################################################
# Save Visualization Methods
####################################################################

#' Plot PDF in outputDirectory of an ArchRProject
#' 
#' This function will save a plot or set of plots as a PDF file in the outputDirectory of a given ArchRProject.
#' 
#' @param ... vector of plots to be plotted (if input is a list use plotList instead)
#' @param name The file name to be used for the output PDF file.
#' @param width The width in inches to be used for the output PDF file.
#' @param height The height in inches to be used for the output PDF.
#' @param ArchRProj An `ArchRProject` object to be used for retrieving the desired `outputDirectory` which will be used to store the output
#' plots in a subfolder called "plots".
#' @param addDOC A boolean variable that determines whether to add the date of creation to the end of the PDF file name. This is useful
#' for preventing overwritting of old plots.
#' @param useDingbats A boolean variable that determines wheter to use dingbats characters for plotting points.
#' @param plotList A `list` of plots to be printed to the output PDF file. Each element of `plotList` should be a printable plot formatted
#' object (ggplot2, plot, heatmap, etc).
#' @export
plotPDF <- function(
  ...,
  name = "Plot",
  width = 6,
  height = 6,
  ArchRProj = NULL,
  addDOC = TRUE,
  useDingbats = FALSE,
  plotList = NULL
  ){

  #Validate
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = width, name = "width", valid = "numeric")
  .validInput(input = height, name = "height", valid = "numeric")
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject", "null"))
  .validInput(input = addDOC, name = "addDOC", valid = "boolean")
  .validInput(input = useDingbats, name = "useDingbats", valid = "boolean")
  .validInput(input = plotList, name = "plotList", valid = c("list","null"))
  #########

  if(is.null(plotList)){
    plotList <- list(...)
    plotList2 <- list()
    for(i in seq_along(plotList)){
      if(inherits(plotList[[i]], "list")){
        for(j in seq_along(plotList[[i]])){
          plotList2[[length(plotList2) + 1]] <- plotList[[i]][[j]]
        }
      }else{
        plotList2[[length(plotList2) + 1]] <- plotList[[i]]
      }
    }
    plotList <- plotList2
    rm(plotList2)
    gc()
  }else{
    plotList2 <- list()
    for(i in seq_along(plotList)){
      if(inherits(plotList[[i]], "list")){
        for(j in seq_along(plotList[[i]])){
          plotList2[[length(plotList2) + 1]] <- plotList[[i]][[j]]
        }
      }else{
        plotList2[[length(plotList2) + 1]] <- plotList[[i]]
      }
    }
    plotList <- plotList2
    rm(plotList2)
    gc()
  }
  
  name <- gsub("\\.pdf", "", name)
  if(is.null(ArchRProj)){
    outDir <- "Plots"
  }else{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
    outDir <- file.path(getOutputDirectory(ArchRProj), "Plots")
  }
  
  dir.create(outDir, showWarnings = FALSE)
  if(addDOC){
    doc <- gsub(":","-",stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2])
    filename <- file.path(outDir, paste0(name, "_Date-", Sys.Date(), "_Time-", doc, ".pdf"))
  }else{
    filename <- file.path(outDir, paste0(name, ".pdf"))
  }

  o <- tryCatch({

    pdf(filename, width = width, height = height, useDingbats = useDingbats)
    for(i in seq_along(plotList)){
      
      if(inherits(plotList[[i]], "gg")){

        if(inherits(plotList[[i]], "patchwork")){

          if(getArchRVerbose()) message("Plotting Patchwork!")
          print(plotList[[i]])
        
        }else{

          if(getArchRVerbose()) message("Plotting Ggplot!")

          if(!is.null(attr(plotList[[i]], "ratioYX"))){
            .fixPlotSize(plotList[[i]], plotWidth = width, plotHeight = height, height = attr(plotList[[i]], "ratioYX"), newPage = FALSE)
          }else{
            .fixPlotSize(plotList[[i]], plotWidth = width, plotHeight = height, newPage = FALSE)
          }

        }

        if(i != length(plotList)){
          grid::grid.newpage()
        }
      
      }else if(inherits(plotList[[i]], "gtable")){

        if(getArchRVerbose()) message("Plotting Gtable!")
        
        print(grid::grid.draw(plotList[[i]]))
        if(i != length(plotList)){
          grid::grid.newpage()
        }
      }else if(inherits(plotList[[i]], "HeatmapList") | inherits(plotList[[i]], "Heatmap") ){ 

        if(getArchRVerbose()) message("Plotting ComplexHeatmap!")

        padding <- 15
        draw(plotList[[i]], 
          padding = unit(c(padding, padding, padding, padding), "mm"), 
          heatmap_legend_side = "bot", 
          annotation_legend_side = "bot"
        )

      }else{

        if(getArchRVerbose()) message("Plotting Other")
       
        print(plotList[[i]])

      }

    }
    dev.off()


  }, error = function(x){

    if(getArchRVerbose()) message(x)

  })

  return(invisible(0))

}

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
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override discreteSet/continuousSet for coloring vector.
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
#' @param threads The number of threads to be used for parallel computing.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional parameters to pass to `ggPoint()` or `ggHex()`.
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
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .requirePackage("ggplot2", source = "cran")

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)

  ##############################
  # Get Embedding
  ##############################
  .logMessage("Getting UMAP Embedding", logFile = logFile)
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)

  if(!all(rownames(df) %in% ArchRProj$cellNames)){
    stop("Not all cells in embedding are present in ArchRProject!")
  }

  .logThis(df, name = "Embedding data.frame", logFile = logFile)
  if(!is.null(sampleCells)){
    if(sampleCells < nrow(df)){
      if(!is.null(imputeWeights)){
        stop("Cannot sampleCells with imputeWeights not equalt to NULL at this time!")
      }
      df <- df[sort(sample(seq_len(nrow(df)), sampleCells)), , drop = FALSE]
    }
  }

  #Parameters
  plotParams <- list(...)
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

  #Check if Cells To Be Highlighed
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
  allColorBy <-  c("colData", "cellColData", .availableArrays(head(getArrowFiles(ArchRProj), 2)))
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


  }else{

    suppressMessages(message(logFile))

    units <- tryCatch({
        .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
      },error=function(e){
        "values"
    })
    
    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
    }

    if(is.null(log2Norm)){
      log2Norm <- FALSE
    }

    colorMat <- .getMatrixValues(
      ArchRProj = ArchRProj, 
      name = name, 
      matrixName = colorBy, 
      log2Norm = FALSE, 
      threads = threads,
      logFile = logFile
    )

    if(!all(rownames(df) %in% colnames(colorMat))){
      .logMessage("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.", logFile = logFile)
      stop("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.")
    }

    colorMat <- colorMat[,rownames(df), drop=FALSE]

    .logThis(colorMat, "colorMat-Before-Impute", logFile = logFile)

    if(!is.null(imputeWeights)){
      if(getArchRVerbose()) message("Imputing Matrix")
      colorMat <- imputeMatrix(mat = as.matrix(colorMat), imputeWeights = imputeWeights, logFile = logFile)
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

      plotParamsx$color <- .quantileCut(plotParamsx$color, min(quantCut), max(quantCut))

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


#' Visualize Groups from ArchR Project
#' 
#' This function will group, summarize and then plot data from an ArchRProject for visual comparison.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing and plotting.
#' @param colorBy A string indicating whether the numeric values to be used in the violin plot should be from a column in
#' `cellColData` ("cellColData") or from a data matrix in the ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName/rowname of the data matrix to be used for plotting. 
#' For example if `colorBy` is "cellColData" then `name` refers to a column name in the cellcoldata (see `getCellcoldata()`). If `colorBy`
#' is "GeneScoreMatrix" then `name` refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells values. See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param maxCells The maximum cells to consider when making the plot.
#' @param quantCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution of values.
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(a,b) where `a` is the upper threshold and
#' `b` is the lower threshold. For example, quantCut = c(0.025,0.975) will take the top and bottom 2.5 percent of values and set them
#' to the value of the 97.5th and 2.5th percentile values respectively.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values (if continuous) in plotting.
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override discreteSet/continuousSet for coloring vector.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` if a discrete color set is desired.
#' @param ylim A vector of two numeric values indicating the lower and upper bounds of the y-axis on the plot.
#' @param size The numeric size of the points to be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param ridgeScale The scale factor for the relative heights of each ridge when making a ridgeplot with `ggridges`.
#' @param plotAs A string that indicates whether a rigdge plot ("ridges") should be plotted or a violin plot ("violin") should be plotted.
#' @param threads The number of threads to be used for parallel computing.
#' @param ... Additional parameters to pass to `ggGroup()`.
#' @export
plotGroups <- function(
  ArchRProj = NULL,
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment",
  imputeWeights = if(!grepl("coldata",tolower(colorBy[1]))) getImputeWeights(ArchRProj),
  maxCells = 1000,
  quantCut = c(0.002, 0.998),
  log2Norm = NULL,
  pal = NULL,
  discreteSet = "stallion",
  ylim = NULL, 
  size = 0.5, 
  baseSize = 6, 
  ratioYX = NULL,
  ridgeScale = 2,
  plotAs = "ridges",
  threads = getArchRThreads(),
  ...
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = maxCells, name = "maxCells", valid = c("integer"))
  .validInput(input = quantCut, name = "quantCut", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric", "null"))
  .validInput(input = ridgeScale, name = "ridgeScale", valid = c("numeric"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  .requirePackage("ggplot2", source = "cran")

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }
  allColorBy <-  c("colData", "cellColData", .availableArrays(head(getArrowFiles(ArchRProj), 2)))
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]

  groups <- getCellColData(ArchRProj, groupBy, drop = FALSE)
  groupNames <- groups[,1]
  names(groupNames) <- rownames(groups)
  groupNames2 <- gtools::mixedsort(unique(groupNames))


  plotParams <- list(...)

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
      
    colorList <- lapply(seq_along(name), function(x){
      colorParams <- list()
      colorParams$color <- as.vector(getCellColData(ArchRProj, select = name[x], drop = TRUE))
      if(!is.numeric(colorParams$color)){
        stop(paste0("colorBy = cellColData, name = ", name[x], " : name must correspond to a numeric column!"))
      }
      if(!is.null(discreteSet)){
        colorParams$pal <- paletteDiscrete(values = groupNames2, set = discreteSet)
      }
      if(!is.null(pal)){
        colorParams$pal <- pal
      }
      colorParams
    })

  }else{

    units <- tryCatch({
        .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
      },error=function(e){
        "values"
    })
    
    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
    }

    if(is.null(log2Norm)){
      log2Norm <- FALSE
    }

    colorMat <- .getMatrixValues(
      ArchRProj = ArchRProj, 
      name = name, 
      matrixName = colorBy, 
      log2Norm = FALSE, 
      threads = threads
    )

    if(!is.null(imputeWeights)){
      colorMat <- imputeMatrix(mat = as.matrix(colorMat), imputeWeights = imputeWeights)
      if(!inherits(colorMat, "matrix")){
        colorMat <- matrix(colorMat, ncol = nCells(ArchRProj))
        colnames(colorMat) <- ArchRProj$cellNames
      }
    }

    colorList <- lapply(seq_len(nrow(colorMat)), function(x){
      colorParams <- list()
      colorParams$color <- colorMat[x, ]
      if(!is.null(discreteSet)){
        colorParams$pal <- suppressMessages(paletteDiscrete(values = groupNames2, set = discreteSet))
      }
      if(!is.null(pal)){
        colorParams$pal <- pal
      }
      colorParams
    })

  }

  if(!is.null(maxCells)){
    splitGroup <- split(names(groupNames), groupNames)
    useCells <- lapply(splitGroup, function(x){
      if(length(x) > maxCells){
        sample(x, maxCells)
      }else{
        x
      }
    }) %>% unlist %>% as.vector
    idx <- match(useCells, names(groupNames))
  }else{
    idx <- seq_along(groupNames)
  }

  pl <- lapply(seq_along(colorList), function(x){

    if(getArchRVerbose()) message(paste0(x, " "), appendLF = FALSE)

    if(is.null(ylim)){
      ylim <- range(colorList[[x]]$color,na.rm=TRUE) %>% extendrange(f = 0.05)
    }

    plotParamsx <- plotParams
    plotParamsx$x <- groupNames[idx]
    if(!is.null(quantCut)){
      plotParamsx$y <- .quantileCut(colorList[[x]]$color[idx], min(quantCut), max(quantCut))
    }else{
      plotParamsx$y <- colorList[[x]]$color[idx]
    }
    plotParamsx$xlabel <- groupBy
    plotParamsx$ylabel <- name[x]
    plotParamsx$baseSize <- baseSize
    plotParamsx$ridgeScale <- ridgeScale
    plotParamsx$ratioYX <- ratioYX
    plotParamsx$size <- size
    plotParamsx$plotAs <- plotAs
    plotParamsx$pal <- colorList[[x]]$pal

    p <- do.call(ggGroup, plotParamsx)

    p

  })

  names(pl) <- name
  if(getArchRVerbose()) message("")
  
  if(length(name)==1){
    pl[[1]]
  }else{
    pl
  }

}

.getMatrixValues <- function(
  ArchRProj = NULL, 
  name = NULL, 
  matrixName = NULL, 
  log2Norm = FALSE, 
  threads = getArchRThreads(),
  logFile = NULL
  ){
  
  o <- h5closeAll()

  .logMessage("Getting Matrix Values...", verbose = TRUE, logFile = logFile)

  featureDF <- .getFeatureDF(head(getArrowFiles(ArchRProj), 2), matrixName)
  .logThis(featureDF, "FeatureDF", logFile = logFile)

  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(matrixName, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    if(!all(unlist(lapply(name, function(x) grepl(":",x))))){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!", logFile = logFile)
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }

  if(grepl(":",name[1])){

    sname <- stringr::str_split(name,pattern=":",simplify=TRUE)[,1]
    name <- stringr::str_split(name,pattern=":",simplify=TRUE)[,2]

    idx <- lapply(seq_along(name), function(x){
      ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }else{
    
    idx <- lapply(seq_along(name), function(x){
      ix <- which(tolower(name[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }
  .logThis(idx, "idx", logFile = logFile)

  if(any(is.na(idx))){
    .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", paste0(name[which(is.na(idx))], collapse=",")), logFile = logFile)
  }

  featureDF <- featureDF[idx, ,drop=FALSE]
  .logThis(featureDF, "FeatureDF-Subset", logFile = logFile)

  #Get Values for FeatureName
  cellNamesList <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  
  values <- .safelapply(seq_along(cellNamesList), function(x){
    if(getArchRVerbose()) message(x, " ", appendLF = FALSE)
    valuesx <- tryCatch({
      o <- h5closeAll()
      ArrowFile <- getSampleColData(ArchRProj)[names(cellNamesList)[x],"ArrowFiles"]
      valuesx <- .getMatFromArrow(
          ArrowFile = ArrowFile, 
          featureDF = featureDF,
          binarize = FALSE, 
          useMatrix = matrixName, 
          cellNames = cellNamesList[[x]],
          threads = 1
        )
      colnames(valuesx) <- cellNamesList[[x]]
      valuesx
    }, error = function(e){
      errorList <- list(
        x = x,
        ArrowFile = ArrowFile,
        ArchRProj = ArchRProj, 
        cellNames = ArchRProj$cellNames, 
        cellNamesList = cellNamesList, 
        featureDF = featureDF
      )
      .logError(e, fn = ".getMatFromArrow", info = "", errorList = errorList, logFile = logFile)  
    })
    valuesx
  }, threads = threads) %>% Reduce("cbind", .)
  values <- values[, ArchRProj$cellNames, drop = FALSE]
  if(getArchRVerbose()) message("")
  gc()
  .logThis(values, "Feature-Matrix", logFile = logFile)

  if(!inherits(values, "matrix")){
    values <- matrix(as.matrix(values), ncol = nCells(ArchRProj))
    colnames(values) <- ArchRProj$cellNames
  }

  #Values Summary
  if(!is.null(log2Norm)){
    if(log2Norm){
      if(getArchRVerbose()) message("Log2 Normalizing...")
      values <- log2(values + 1)
    }
  }

  rownames(values) <- name

  return(values)

}

.fixPlotSize <- function(
  p = NULL, 
  plotWidth = unit(6, "in"),
  plotHeight = unit(6, "in"),
  margin = 0.25,
  height = 1,
  it = 0.05,
  newPage = FALSE
  ){

  .requirePackage("grid", source = "cran")
  .requirePackage("gridExtra", source = "cran")

  if(!inherits(plotWidth, "unit")){
    plotWidth <- unit(plotWidth, "in") 
  }

  if(!inherits(plotHeight, "unit")){
    plotHeight <- unit(plotHeight, "in") 
  }

  #adapted from https://github.com/jwdink/egg/blob/master/R/set_panel_size.r
  g <- ggplotGrob(p)
  
  legend <- grep("guide-box", g$layout$name)
  if(length(legend)!=0){
    gl <- g$grobs[[legend]]
    g <- ggplotGrob(p + theme(legend.position = "none"))
  }else{
    gl <- NULL
    g <- ggplotGrob(p)
  }

  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])

  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  pw <- convertWidth(plotWidth, unitTo = "in", valueOnly = TRUE)
  ph <- convertWidth(plotHeight, unitTo = "in", valueOnly = TRUE)

  pw <- pw * 0.95
  ph <- ph * 0.95

  x <- 0
  width <- 1
  sm <- FALSE
  
  while(!sm){
    
    x <- x + it

    w <- unit(x * width, "in")
    h <- unit(x * height / width, "in")
    m <- unit(x * margin / width, "in")

    g$widths[panel_index_w] <-  rep(w, nw)
    g$heights[panel_index_h] <- rep(h, nh)

    sw <- convertWidth(
      x = sum(g$widths)  + m, 
      unitTo = "in", 
      valueOnly = TRUE
    )

    sh <- convertHeight(
      x = sum(g$heights) + m, 
      unitTo = "in", 
      valueOnly = TRUE
    )
    
    sm <- sw > pw | sh > ph

  }

  if(length(legend)!=0){

    sgh <- convertHeight(
      x = sum(g$heights), 
      unitTo = "in", 
      valueOnly = TRUE
    )

    sgw <- convertWidth(
      x = sum(g$widths), 
      unitTo = "in", 
      valueOnly = TRUE
    )

    slh <- convertHeight(
      x = sum(gl$heights), 
      unitTo = "in", 
      valueOnly = TRUE
    )

    slw <- convertWidth(
      x = sum(gl$widths), 
      unitTo = "in", 
      valueOnly = TRUE
    )

    size <- 6
    wh <- 0.1
    it <- 0

    while(slh > 0.2 * ph | slw > pw){

      it <- it + 1

      if(it > 3){
        break
      }

      size <- size * 0.8
      wh <- wh * 0.8

      gl <- ggplotGrob(
        p + theme(
              legend.key.width = unit(wh, "cm"),
              legend.key.height = unit(wh, "cm"),
              legend.spacing.x = unit(0, 'cm'),
              legend.spacing.y = unit(0, 'cm'),
              legend.text = element_text(size = max(size, 2))
              ) + guides(fill = guide_legend(ncol = 4), color =  guide_legend(ncol = 4))
        )$grobs[[legend]]

      slh <- convertHeight(
        x = sum(gl$heights), 
        unitTo = "in", 
        valueOnly = TRUE
      )

      slw <- convertWidth(
        x = sum(gl$widths), 
        unitTo = "in", 
        valueOnly = TRUE
      )

    }

    p <- grid.arrange(g, gl, ncol=1, nrow=2, 
      heights = unit.c(unit(sgh,"in"), unit(min(slh, 0.2 * pw), "in")),
      newpage = newPage
    )

  }else{

    p <- grid.arrange(g, newpage = newPage)

  }


  invisible(p)

}


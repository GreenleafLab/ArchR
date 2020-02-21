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
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = embedding, name = "reducedDims", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = quantCut, name = "quantCut", valid = c("numeric"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character", "null"))
  .validInput(input = continuousSet, name = "continuousSet", valid = c("character", "null"))
  .validInput(input = randomize, name = "randomize", valid = c("boolean"))
  .validInput(input = keepAxis, name = "keepAxis", valid = c("boolean"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character", "null"))

  .requirePackage("ggplot2", source = "cran")

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)

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

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }
  allColorBy <-  c("colData", "cellColData", .availableArrays(getArrowFiles(ArchRProj)))
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
      
    colorList <- lapply(seq_along(name), function(x){
      colorParams <- list()
      colorParams$color <- as.vector(getCellColData(ArchRProj, select = name[x], drop = TRUE))
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
      colorParams
    })

  }else{
    
    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
    }

    colorMat <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)[,rownames(df), drop=FALSE]

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
      colorParams
    })

  }

  message("Plotting Embedding")

  ggList <- lapply(seq_along(colorList), function(x){

    message(x, " ", appendLF = FALSE)

    plotParamsx <- .mergeParams(colorList[[x]], plotParams)

    if(plotParamsx$discrete){
      plotParamsx$color <- paste0(plotParamsx$color)
    }

    if(!plotParamsx$discrete){

      plotParamsx$color <- .quantileCut(plotParamsx$color, min(quantCut), max(quantCut))

      if(!is.null(imputeWeights)){
        imputeWeights <- imputeWeights$Weights[rownames(df), rownames(df)]
        plotParamsx$color <- (imputeWeights %*% as(as.matrix(plotParamsx$color), "dgCMatrix"))[,1] 
      }

      plotParamsx$pal <- paletteContinuous(set = plotParamsx$continuousSet)

      if(!is.null(pal)){

        plotParamsx$pal <- pal
        
      }

      if(is.null(plotAs)){
        plotAs <- "hexplot"
      }

      if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){

        plotParamsx$discrete <- NULL
        plotParamsx$continuousSet <- NULL
        plotParamsx$rastr <- NULL
        plotParamsx$size <- NULL
        plotParamsx$randomize <- NULL

        gg <- do.call(ggHex, plotParamsx)

      }else{

        gg <- do.call(ggPoint, plotParamsx)

      }

    }else{
      
      if(!is.null(pal)){
        plotParamsx$pal <- pal
      }

      gg <- do.call(ggPoint, plotParamsx)

    }

    if(!keepAxis){
      gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }

    gg

  })
  names(ggList) <- name
  message("")

  if(length(ggList) == 1){
    ggList <- ggList[[1]]
  }

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
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values (if continuous) in plotting.
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override discreteSet/continuousSet for coloring vector.
#' @param ylim A vector of two numeric values indicating the lower and upper bounds of the y-axis on the plot.
#' @param size The numeric size of the points to be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param ridgeScale The scale factor for the relative heights of each ridge when making a ridgeplot with `ggridges`.
#' @param plotAs A string that indicates whether a rigdge plot ("ridges") should be plotted or a violin plot ("violin") should be plotted.
#' @param ... Additional parameters to pass to `ggGroup()`.
#' @export
plotGroups <- function(
  ArchRProj = NULL, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment",
  imputeWeights = if(!grepl("coldata",tolower(colorBy[1]))) getImputeWeights(ArchRProj),
  log2Norm = NULL,
  pal = NULL,
  ylim = NULL, 
  size = 0.5, 
  baseSize = 6, 
  ratioYX = NULL,
  ridgeScale = 1,
  plotAs = "ridges",
  ...
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric", "null"))
  .validInput(input = ridgeScale, name = "ridgeScale", valid = c("numeric"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character"))

  .requirePackage("ggplot2", source = "cran")

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }
  allColorBy <-  c("colData", "cellColData", .availableArrays(getArrowFiles(ArchRProj)))
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]

  groups <- getCellColData(ArchRProj, groupBy, drop = FALSE)
  groupNames <- groups[,1]
  names(groupNames) <- rownames(groups)

  plotParams <- list(...)

  pl <- lapply(seq_along(name), function(x){

    message(paste0(x, " "), appendLF = FALSE)

    if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
      values <- getCellColData(ArchRProj, name[x], drop = TRUE)
    }else{
      if(tolower(colorBy) == "genescorematrix"){
        if(is.null(log2Norm)){
          log2Norm <- TRUE
        }
      }
      values <- .getMatrixValues(ArchRProj, name = name[x], matrixName = colorBy, log2Norm = log2Norm)[1, names(groupNames)]
    }

    if(!is.null(imputeWeights)){
      imputeWeights <- imputeWeights$Weights[names(groupNames), names(groupNames)]
      values <- (imputeWeights %*% as(as.matrix(values), "dgCMatrix"))[,1] 
    }

    if(is.null(ylim)){
      ylim <- range(values) %>% extendrange(f = 0.05)
    }

    plotParamsx <- plotParams
    plotParamsx$x <- groupNames
    plotParamsx$y <- values
    plotParamsx$xlabel <- groupBy
    plotParamsx$ylabel <- name[x]
    plotParamsx$baseSize <- baseSize
    plotParamsx$ridgeScale <- ridgeScale
    plotParamsx$ratioYX <- ratioYX
    plotParamsx$size <- size
    plotParamsx$plotAs <- plotAs

    p <- do.call(ggGroup, plotParamsx)

    p

  })
  names(pl) <- name

  message("\n")
  
  if(length(name)==1){
    pl[[1]]
  }else{
    pl
  }

}

.getMatrixValues <- function(ArchRProj = NULL, name = NULL, matrixName = NULL, log2Norm = TRUE){
  
  o <- h5closeAll()

  message("Getting Matrix Values...")

  featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), matrixName)

  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(matrixName, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    if(!all(unlist(lapply(name, function(x) grepl(":",x))))){
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }

  if(grepl(":",name[1])){

    sname <- stringr::str_split(name,pattern=":",simplify=TRUE)[,1]
    name <- stringr::str_split(name,pattern=":",simplify=TRUE)[,2]

    idx <- lapply(seq_along(name), function(x){
      ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        stop(sprintf("FeatureName (%s) does not exist! See availableFeatures", name[x]))
      }
      ix
    }) %>% unlist

  }else{
    
    idx <- lapply(seq_along(name), function(x){
      ix <- which(tolower(name[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        stop(sprintf("FeatureName (%s) does not exist! See availableFeatures", name[x]))
      }
      ix
    }) %>% unlist

  }

  if(any(is.na(idx))){
    stop(sprintf("FeatureName (%s) does not exist! See availableFeatures", name[which(is.na(idx))]))
  }

  featureDF <- featureDF[idx, ,drop=FALSE]

  #Get Values for FeatureName
  cellNamesList <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  
  values <- lapply(seq_along(cellNamesList), function(x){
    message(x, " ", appendLF = FALSE)
    o <- h5closeAll()
    ArrowFile <- getSampleColData(ArchRProj)[names(cellNamesList)[x],"ArrowFiles"]
    valuesx <- .getMatFromArrow(
        ArrowFile = ArrowFile, 
        featureDF = featureDF,
        binarize = FALSE, 
        useMatrix = matrixName, 
        cellNames = cellNamesList[[x]]
      )
    colnames(valuesx) <- cellNamesList[[x]]
    valuesx
  }) %>% Reduce("cbind", .)
  message("")
  gc()

  #Values Summary
  if(!is.null(log2Norm)){
    if(log2Norm){
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

.isDiscrete <- function(x = NULL){
  is.factor(x) || is.character(x) || is.logical(x)
}



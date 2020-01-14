####################################################################
# Visualization Methods
####################################################################

#' Visualize Embedding from ArchR Project
#' 
#' This function will plot an embedding that was created from
#' computeEmbedding
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the embedding to plot. See `ArchR::computeEmbedding()` for more information.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in cellColData ("cellColData") or by a data matrix in the ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrx", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName in the data matrix.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values (if continuous) in plotting.
#' @param imputeWeights imputation weights for imputing numerical values for each cell as a linear combination of other cells values (see add/getImutationWeights).
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized with ggrastr. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param quantCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values. 
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the upper threshold and y is 
#' the lower threshold. For example, quantileCut = c(0.975,0.025) will take the top and bottom 2.5% of values and set them to the value of 
#' the 97.5th and 2.5th percentile values respectively.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing colorBy in the embedding.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing colorBy in the embedding.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted.
#' @param plotParams Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' @param ... additional args
#' @export
plotEmbedding <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  colorBy = "colData",
  name = "Sample",
  log2Norm = NULL,
  imputeWeights = NULL,
  pal = NULL,
  size = 0.5,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  plotAs = NULL,
  plotParams = list(),
  ...
  ){

  .requirePackage("ggplot2")

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)

  #Parameters
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


  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
      
    colorList <- lapply(seq_along(name), function(x){
      colorParams <- list()
      colorParams$color <- as.vector(getCellColData(ArchRProj)[,name[x]])
      colorParams$discrete <- .isDiscrete(colorParams$color)
      colorParams$continuousSet <- "solar_extra"
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
        colorParams$continuousSet <- "horizon_extra"
      }else{
        colorParams$continuousSet <- "solar_extra"
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

  ggList <- lapply(seq_along(colorList), function(x){

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

  if(length(ggList) == 1){
    ggList <- ggList[[1]]
  }

  ggList

}


#' Visualize Groups from ArchR Project
#' 
#' This function will plot an embedding that was created from computeEmbedding
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy use groupings in cellColData for summarizing and plotting
#' @param colorBy A string indicating whether numeric values in violin plot should be from a column in cellColData ("cellColData") or by a data matrix in the ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrx", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName in the data matrix.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param ylim A vector of two numeric values indicating the lower and upper bounds of the y-axis on the plot.
#' @param size The numeric size of the points to be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param addPoints A boolean value that indicates whether points should be added to the plot using `geom_quasirandom()`
#' @param ... additional args
#' @export
plotGroups <- function(
  ArchRProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment",
  imputeWeights = NULL, 
  log2Norm = NULL,
  pal = NULL,
  ylim = NULL, 
  size = 0.5, 
  baseSize = 6, 
  ratioYX = 0.5, 
  addPoints = FALSE, 
  ...
  ){
  
  .requirePackage("ggplot2")

  groups <- getCellColData(ArchRProj, groupBy, drop = FALSE)
  groupNames <- groups[,1]
  names(groupNames) <- rownames(groups)
  
  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    values <- getCellColData(ArchRProj, name, drop = TRUE)
  }else{
    if (tolower(colorBy) == "genescorematrix"){
      if(is.null(log2Norm)){
        log2Norm <- TRUE
      }
    }
    values <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)[1, names(groupNames)]
  }

  if(!is.null(imputeWeights)){
    imputeWeights <- imputeWeights$Weights[names(groupNames), names(groupNames)]
    values <- (imputeWeights %*% as(as.matrix(values), "dgCMatrix"))[,1] 
  }

  if(is.null(ylim)){
    ylim <- range(values) %>% extendrange(f = 0.05)
  }

  p <- ggViolin(
    x = groupNames, 
    y = values, 
    xlabel = groupBy, 
    ylabel = name, 
    baseSize = baseSize, 
    ratioYX = ratioYX,
    size = size,
    addPoints = addPoints
    )

  p

}

#' @export
.getMatrixValues <- function(ArchRProj, name, matrixName, log2Norm = TRUE){
  
  o <- h5closeAll()

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
        stop(sprintf("FeatureName does not exist for %s! See availableFeatures", name))
      }
      ix
    }) %>% unlist

  }else{
    
    idx <- lapply(seq_along(name), function(x){
      ix <- which(tolower(name[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        stop(sprintf("FeatureName does not exist for %s! See availableFeatures", name))
      }
      ix
    }) %>% unlist

  }

  if(is.na(idx)){
    stop("name is not in featureNames!")
  }

  featureDF <- featureDF[idx, ,drop=FALSE]

  #Get Values for FeatureName
  cellNamesList <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  
  values <- lapply(seq_along(cellNamesList), function(x){
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


#' @export
.fixPlotSize <- function(
  p = NULL, 
  plotWidth = unit(6, "in"),
  plotHeight = unit(6, "in"),
  margin = 0.25,
  height = 1,
  it = 0.05,
  newPage = FALSE,
  ...
  ){

  .requirePackage("grid")
  .requirePackage("gridExtra")

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

.isDiscrete <- function(x){
  is.factor(x) || is.character(x) || is.logical(x)
}













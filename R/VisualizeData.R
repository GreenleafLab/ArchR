#' Visualize Embedding from ArchR Project
#' 
#' This function will plot an embedding that was created from
#' computeEmbedding
#'
#' @param ArchRProj ArchRProject
#' @param embedding embedding to visualize (see computeEmbedding)
#' @param colorBy colorBy cellColData or Arrays in Arrows (ie GeneScoreMatrix)
#' @param name name of column in cellColData or Feature in Array in Arrows
#' @param log2Norm log2 Normalize features if they are continuous
#' @param pal custom palette to use for plotting
#' @param size size of points in plot
#' @param rastr rastr points in plot
#' @param quantCut quantile cut of continuous features
#' @param quantHex quantile evaluation for each hex in geom_hex
#' @param discreteSet discrete palette for visualizing embedding
#' @param continuousSet continuous palette for visualizing embedding
#' @param randomize randomize points prior to plotting
#' @param keepAxis keep x and y axis for plot
#' @param baseSize base size for text in plot
#' @param plotContinuous how to plot continuous features (points and hex)
#' @param plotParams additional params to pass to ggPoint/ggHex
#' @param ... additional args
#' @export
VisualizeEmbedding <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  colorBy = "colData",
  name = "Sample",
  log2Norm = NULL,
  pal = NULL,
  size = 0.5,
  rastr = TRUE,
  quantCut = c(0.05, 0.95),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  plotContinuous = NULL,
  plotParams = list(),
  ...
  ){

  .requirePackage("ggplot2")

  ##############################
  # Plot Helpers
  ##############################
  .quantileCut <- function (x, lo = 0, hi = 0.975, rm0 = TRUE){
    q <- quantile(x, probs = c(lo, hi), na.rm = TRUE)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(x)
  }

  .summarizeHex <- function(x){
    quantile(x, quantHex)
  }

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, return = "df")

  #Parameters
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    
    plotParams$color <- as.vector(getCellColData(ArchRProj)[,name])
    plotParams$discrete <- .isDiscrete(plotParams$color)
    plotParams$continuousSet <- "solar_extra"
    plotParams$discreteSet <- "stallion"
    plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
    if(is.null(plotContinuous)){
      plotContinuous <- "hexplot"
    }

  }else{
    if (tolower(colorBy) == "genescorematrix"){
      if(is.null(log2Norm)){
        log2Norm <- TRUE
      }
      plotParams$continuousSet <- "white_blue_purple"
    }else{
      plotParams$continuousSet <- "solar_extra"
    }
    plotParams$color <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)
    plotParams$discrete <- FALSE
    plotParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name)
    if(is.null(plotContinuous)){
      plotContinuous <- "hexplot"
    }
    if(plotContinuous=="hexplot"){
      plotParams$fun <- .summarizeHex
    }

  }

  #Additional Params!
  plotParams$xlabel <- gsub("_", " ",stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,2])
  plotParams$ylabel <- gsub("_", " ",stringr::str_split(colnames(df)[2],pattern="#",simplify=TRUE)[,2])

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  
  if(plotParams$discrete){
    plotParams$color <- paste0(plotParams$color)
  }

  if(!plotParams$discrete){
    plotParams$color <- .quantileCut(plotParams$color, min(quantCut), max(quantCut))
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)
    if(tolower(plotContinuous) == "hex" | tolower(plotContinuous) == "hexplot"){
      out <- do.call(ggHex, plotParams)
    }else{
      out <- do.call(ggPoint, plotParams)
    }
  }else{
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  out

}


#' Visualize Groups from ArchR Project
#' 
#' This function will plot an embedding that was created from
#' computeEmbedding
#'
#' @param ArchRProj ArchRProject
#' @param groupBy use groupings in cellColData for summarizing and plotting
#' @param colorBy colorBy cellColData or Arrays in Arrows (ie GeneScoreMatrix)
#' @param name name of column in cellColData or Feature in Array in Arrows
#' @param pal custom palette to use for plotting
#' @param ylim limits for features in plot
#' @param size size of points in ggplot
#' @param baseSize rastr points in ggplot
#' @param ratioYX ratio of Y axis to X axis
#' @param points add points to plot using quasirandom
#' @param ... additional args
#' @export
VisualizeGroups <- function(
  ArchRProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment", 
  log2Norm = NULL,
  pal = NULL,
  ylim = NULL, 
  size = 0.5, 
  baseSize = 6, 
  ratioYX = NULL, 
  points = FALSE, 
  ...
  ){
  
  .requirePackage("ggplot2")

  groupNames <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  
  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    values <- getCellColData(ArchRProj, name, drop = TRUE)
  }else{
    if (tolower(colorBy) == "genescorematrix"){
      if(is.null(log2Norm)){
        log2Norm <- TRUE
      }
    }
    values <- .getMatrixValues(ArchRProj, name = name, matrixName = colorBy, log2Norm = log2Norm)
  }

  if(is.null(ylim)){
    ylim <- range(values) %>% extendrange(f = 0.05)
  }

  if(is.null(ratioYX)){
    ratioYX <-  sqrt(length(unique(groupNames)) / 2)
  }

  p <- ggViolin(
    x = groupNames, 
    y = values, 
    xlabel = groupBy, 
    ylabel = name, 
    baseSize = baseSize, 
    ratioYX = ratioYX * length(unique(groupNames)) / diff(ylim),
    size = size,
    points = points
    )

  p

}

#' @export
.getMatrixValues <- function(ArchRProj, name, matrixName, log2Norm = TRUE){
  
  o <- h5closeAll()
  featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), matrixName)
  if(grepl(":",name)){
    sname <- stringr::str_split(name,pattern=":",simplify=TRUE)[1,1]
    name <- stringr::str_split(name,pattern=":",simplify=TRUE)[1,2]
    idx <- intersect(which(tolower(name) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname) == tolower(featureDF$seqnames)))
  }else{
    idx <- which(tolower(name) == tolower(featureDF$name))[1]
  }
  if(length(idx)==0){
    stop(sprintf("FeatureName does not exist for %s! See availableFeatures", name))
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
  values <- values[1,]
  if(!is.null(log2Norm)){
    if(log2Norm){
      values <- log2(values + 1)
    }
  }

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
    g <- ggplotGrob(p)
  }

  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])

  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  pw <- convertWidth(plotWidth, unitTo = "in", valueOnly = TRUE)
  ph <- convertWidth(plotHeight, unitTo = "in", valueOnly = TRUE)

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

    slh <- convertHeight(
      x = sum(gl$heights), 
      unitTo = "in", 
      valueOnly = TRUE
    )

    p <- grid.arrange(g, gl, ncol=1, nrow=2, heights = unit.c(unit(sgh,"in"), unit(slh, "in")), newpage = newPage)

  }else{

    p <- grid.arrange(g, newpage = newPage)

  }


  invisible(p)

}

.isDiscrete <- function(x){
  is.factor(x) || is.character(x) || is.logical(x)
}













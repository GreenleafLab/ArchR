####################################################################################################
# Applications of Markers!
####################################################################################################

#' Plot a Heatmap of Identified Marker Features
#' 
#' This function will plot a heatmap of the results from markerFeatures
#' 
#' @param seMarker QQQ A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted in the heatmap. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param log2Norm A boolean value indicating whether a log2 transformation whould be performed on the values in `seMarker` prior to plotting. Should be set to `TRUE` for counts-based assays (but not assays like `DeviationsMatrix`).
#' @param scaleTo QQQ scale to prior to log2 Normalization, if log2Norm is FALSE this does nothing
#' @param scaleRows A boolean value that indicates whether the heatmap should display row-wise z-scores instead of raw values.
#' @param limits A numeric vector of two numbers that represent the lower and upper color limits of the heatmap color scheme.
#' @param grepExclude QQQ A character vector or string that indicates the `rownames` or a specific pattern that identifies rownames from `seMarker` to be excluded from the heatmap.
#' @param pal QQQ The name or numeric index of a custom palette from `ArchRPalettes` to use for the heatmap
#' @param binaryClusterRows A boolean value that indicates whether a binary sorting algorithm should be used for fast clustering of heatmap rows.
#' @param labelMarkers A character vector listing the `rownames` of `seMarker` that should be labeled on the side of the heatmap.
#' @param labelTop A boolean value that indicates whether the top features for each column in `seMarker` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMat QQQ A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param invert QQQ
#' @param ... additional args
#' @export
markerHeatmap <- function(
  seMarker,
  cutOff = "FDR <= 0.001 & Log2FC >= 0.1",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  limits = c(-2,2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  labelMarkers = NULL,
  labelTop = NULL,
  labelRows = FALSE,
  returnMat = FALSE,
  ...
  ){

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  #Now Get Values
  mat <- SummarizedExperiment::assays(seMarker)[["Mean"]]
  idx <- which(rowSums(passMat, na.rm = TRUE) > 0 & matrixStats::rowVars(mat) != 0)
  if(log2Norm){
    mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
  }
  mat <- mat[idx,]
  passMat <- passMat[idx,]

  if(scaleRows){
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }

  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }

  #add rownames
  rd <- SummarizedExperiment::rowData(seMarker)[idx,]
  if(is.null(rd$name)){
    rn <- paste0(rd$seqnames,":",rd$start,"-",rd$end)
  }else{
    if(sum(duplicated(rd$name)) > 0){
      rn <- paste0(rd$seqnames,":",rd$name)
    }else{
      rn <- rd$name
    }
  }
  rownames(mat) <- rn
  rownames(passMat) <- rn

  #identify to remove
  if(!is.null(grepExclude) & !is.null(rownames(mat))){
    idx2 <- which(!grepl(grepExclude, rownames(mat)))
    mat <- mat[idx2,]
  }

  if(nrow(mat)==0){
    stop("No Makers Found!")
  }

  if(!is.null(labelTop)){
    spmat <- passMat / rowSums(passMat)
    idx2 <- lapply(seq_len(ncol(spmat)), function(x){
      head(order(spmat[,x], decreasing = TRUE), labelTop)
    }) %>% unlist %>% unique %>% sort
    mat <- mat[idx2,]
    labelRows <- TRUE
  }

  if(binaryClusterRows){
    bS <- .binarySort(mat, lmat = passMat[rownames(mat), colnames(mat)])
    mat <- bS[[1]][,colnames(mat)]
    clusterRows <- FALSE
    clusterCols <- bS[[2]]
  }else{
    clusterRows <- TRUE
    clusterCols <- TRUE
  }

  if(!is.null(labelMarkers)){
    mn <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
    mn <- mn[mn > 0]
  }else{
    mn <- NULL
  }

  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }

  message(sprintf("Identified %s markers!", nrow(mat)))

  if(is.null(pal)){
    if(is.null(metadata(seMarker)$Params$useMatrix)){
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }else if(tolower(metadata(seMarker)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blue_yellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }
  }

  ht <- .ArchRHeatmap(
    mat = mat,
    scale = FALSE,
    limits = c(min(mat), max(mat)),
    color = pal, 
    clusterCols = clusterCols, 
    clusterRows = clusterRows,
    labelRows = labelRows,
    labelCols = TRUE,
    customRowLabel = mn,
    showColDendrogram = TRUE,
    draw = FALSE,
    ...
  )

  if(returnMat){
    return(mat)
  }else{
    return(ht)
  }

}

########################################################################################################
# Helpers for Nice Heatmap with Bioconductors ComplexHeamtap
########################################################################################################

.ArchRHeatmap <- function(
  mat, 
  scale = FALSE,
  limits = c(min(mat), max(mat)),
  colData = NULL, 
  color = paletteContinuous(set = "solar_extra", n = 100),
  clusterCols = TRUE,
  clusterRows = FALSE,
  labelCols = FALSE,
  labelRows = FALSE,
  colorMap = NULL,
  useRaster = TRUE,
  rasterQuality = 5,
  split = NULL,
  fontsize = 6,
  colAnnoPerRow = 4,
  showRowDendrogram = FALSE,
  showColDendrogram = FALSE,
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.75,
  rasterDevice = "png",
  padding = 45,
  borderColor = NA,
  draw = TRUE,
  name = ""){
  
  #Packages
  .requirePackage("ComplexHeatmap")
  .requirePackage("circlize")
  
  #Z-score
  if (scale) {
    message("Scaling Matrix...")
    mat <- .rowZscores(mat, limit = FALSE)
    name <- paste0(name," Z-Scores")
  }
  
  #Get A Color map if null
  if (is.null(colorMap)) {
    colorMap <- .colorMapAnno(colData)
  }
  
  #Prepare ColorMap format for Complex Heatmap
  if (!is.null(colData)){
    colData = data.frame(colData)
    colorMap <- .colorMapForCH(colorMap, colData) #change
    showLegend <- .checkShowLegend(colorMap[match(names(colorMap), colnames(colData))]) #change
  }else {
    colorMap <- NULL
    showLegend <- NULL
  }
  
  #Prepare Limits if needed
  breaks <- NULL
  if (!is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    breaks <- seq(min(limits), max(limits), length.out = length(color))
    color <- circlize::colorRamp2(breaks, color)
  }
  
  if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
    anno_check_version_rows <- ComplexHeatmap::anno_mark
    anno_check_version_cols <- ComplexHeatmap::anno_mark
  }else{
    anno_check_version_rows <- ComplexHeatmap::row_anno_link
    anno_check_version_cols <- ComplexHeatmap::column_anno_link
  }

  #Annotation Heatmap
  if(!is.null(colData) & !is.null(customColLabel)){
    message("Adding Annotations...")
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customRowLabel]
    }
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        ),
      link = anno_check_version_cols(
        at = customColLabel, labels = customColLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs)

    )
  }else if(!is.null(colData)){
    message("Adding Annotations...")
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        )
    )
  }else if(is.null(colData) & !is.null(customColLabel)){
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customRowLabel]
    }
    message("Adding Annotations...")
    ht1Anno <- HeatmapAnnotation(
      link = anno_check_version_cols(
        at = customColLabel, labels = customColLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs)
    )
  }else{
    ht1Anno <- NULL
  }

  message("Preparing Main Heatmap...")
  ht1 <- Heatmap(
    
    #Main Stuff
    matrix = mat,
    name = name,
    col = color, 
    
    #Heatmap Legend
    heatmap_legend_param = 
      list(color_bar = "continuous", 
           legend_direction = "horizontal",
           legend_width = unit(5, "cm")
      ), 
    rect_gp = gpar(col = borderColor), 
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols, 
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontsize), 
    column_names_max_height = unit(100, "mm"),
    
    #Row Options
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontsize), 
    cluster_rows = clusterRows, 
    show_row_dend = showRowDendrogram, 
    clustering_method_rows = "ward.D2",
    split = split, 
    
    #Annotation
    top_annotation = ht1Anno, 

    #Raster Info
    use_raster = useRaster, 
    raster_device = rasterDevice, 
    raster_quality = rasterQuality
  )

  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    ht1 <- ht1 + rowAnnotation(link = 
        anno_check_version_rows(at = customRowLabel, labels = customRowLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs))
  }

  if(draw){
    draw(ht1, 
      padding = unit(c(padding, padding, padding, padding), "mm"), 
      heatmap_legend_side = "bot", 
      annotation_legend_side = "bot")
  }else{
    ht1
  }

}

.colorMapForCH <- function(colorMap, colData){
  colorMap <- colorMap[which(names(colorMap) %in% colnames(colData))]
  colorMapCH <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete")){
      colorx <- colorMap[[x]]
    }else{
      vals <- colData[[names(colorMap)[x]]][!is.na(colData[[names(colorMap)[x]]])]
      s <-  seq(min(vals), max(vals), length.out = length(colorMap[[x]]))
      colorx <- circlize::colorRamp2(s, colorMap[[x]])
    }
    if(any(is.na(names(colorx)))){
      names(colorx)[is.na(names(colorx))] <- paste0("NA",seq_along(names(colorx)[is.na(names(colorx))]))
    }
    return(colorx)
  })
  names(colorMapCH) <- names(colorMap)
  return(colorMapCH)
}

.checkShowLegend <- function(colorMap, max_discrete = 30){
  show <- lapply(seq_along(colorMap), function(x){
      if(attr(colorMap[[x]],"discrete") && length(unique(colorMap[[x]])) > max_discrete){
        sl <- FALSE
      }else{
        sl <- TRUE
      }
      return(sl)
    }) %>% unlist
  names(show) <- names(colorMap)
  return(show)
}

.colorMapAnno <- function(colData, customAnno = NULL, discreteSet = "stallion", continuousSet = "solar_extra"){
  discreteCols <- sapply(colData,function(x) !is.numeric(x))
  if(!is.null(customAnno)){
    colorMap <- lapply(seq_along(discreteCols),function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      if(length(which(customAnno[,1] %in% names(discreteCols[x]))) > 0){
        if(length(which(customAnno[,2] %in% names(colors))) > 0){
          customAnnox <- customAnno[which(customAnno[,2] %in% names(colors)),]
          colors[which(names(colors) %in% customAnnox[,2])] <- paste0(customAnnox[match(names(colors),customAnnox[,2]),3])
        }
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }else{
    colorMap <- lapply(seq_along(discreteCols), function(x){
      if(discreteCols[x]){
       colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
       names(colors) <- unique(colData[[names(discreteCols[x])]])
       attr(colors, "discrete") <- TRUE
      }else{
       colors <- paletteContinuous(set = continuousSet)
       attr(colors, "discrete") <- FALSE
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }

}

.binarySort <- function(m, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){

  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }

  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  return(list(mat = m, hclust = hc))

}

#' QQQ
#' 
#' This function QQQ
#' 
#' @param seMarker QQQ A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param ArchRProj An `ArchRProject` object.
#' @param annotations QQQ
#' @param matches QQQ
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted in the heatmap. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param background QQQ
#' @param ... additional args
#' @export
markerAnnoEnrich <- function(
  seMarker = NULL,
  ArchRProj = NULL,
  annotations = NULL,
  matches = NULL,
  cutOff = "FDR <= 0.01 & Log2FC >= 0",
  background = "bdgPeaks",
  ...){

  tstart <- Sys.time()
  if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
    stop("Only markers identified from PeakMatrix can be used!")
  }

  if(is.null(matches)){
    matches <- getMatches(ArchRProj, annotations)
  }
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  mcols(r1) <- NULL

  r2 <- getPeakSet(ArchRProj)
  mcols(r2) <- NULL

  if(length(which(paste0(seqnames(r1),start(r1),end(r1), sep = "_") %ni% paste0(seqnames(r2),start(r2),end(r2),sep="_"))) != 0){
    stop("Peaks from matches do not match peakSet in ArchRProj!")
  }

  r3 <- GRanges(rowData(seMarker)$seqnames,IRanges(rowData(seMarker)$start, rowData(seMarker)$end))
  mcols(r3) <- NULL
  rownames(matches) <- paste0(seqnames(matches),start(matches),end(matches),sep="_")
  matches <- matches[paste0(seqnames(r3),start(r3),end(r3), sep = "_"), ]

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(tolower(background) %in% c("backgroundpeaks", "bdgpeaks", "background", "bdg")){
    method <- "bdg"
    bdgPeaks <- SummarizedExperiment::assay(getBdgPeaks(ArchRProj))
  }else{
    method <- "all"
  }

  enrichList <- lapply(seq_len(ncol(seMarker)), function(x){
    .messageDiffTime(sprintf("Computing Enrichments %s of %s",x,ncol(seMarker)),tstart)
    idx <- which(passMat[, x])
    if(method == "bdg"){
      .computeEnrichment(matches, idx, c(idx, as.vector(bdgPeaks[idx,])))
    }else{
      .computeEnrichment(matches, idx, seq_len(nrow(matches)))
    }
  }) %>% SimpleList
  names(enrichList) <- colnames(seMarker)

  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
      enrichList[[y]][colnames(matches),x,drop=FALSE]
    }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- colnames(enrichList[[1]])
  assays <- rev(assays)
  out <- SummarizedExperiment::SummarizedExperiment(assays=assays)

  out

}

.computeEnrichment <- function(matches, compare, background){

  matches <- .getAssay(matches,  grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
     pOut$BackgroundFrequency[x], # Number of all successes in background
     pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
     pOut$nCompare[x], # Number that were drawn
     lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)

  #Minus Log10 FDR
  pOut$mlog10FDR <- -log10(p.adjust(matrixStats::rowMaxs(cbind(10^-pOut$mlog10p, 4.940656e-324)), method = "fdr"))
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

  pOut

}









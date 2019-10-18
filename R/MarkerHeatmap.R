#' Plot a Heatmap of Identified Marker Features
#' 
#' This function will plot a heatmap of the results from markerFeatures
#' 
#' @param seMarker Summarized Experiment result from markerFeatures
#' @param FDR False-Discovery Rate Cutoff to Be called a Marker
#' @param log2FC Log2 Fold Change Cutoff to Be called a Marker
#' @param log2Norm log2 Normalization prior to plotting set true for counting assays (not DeviationsMatrix!)
#' @param scaleTo scale to prior to log2 Normalization, if log2Norm is FALSE this does nothing
#' @param scaleRows compute row z-scores on matrix
#' @param limits heatmap color limits 
#' @param grepExclude remove features by grep
#' @param pal palette for heatmap, default will use solar_extra
#' @param binaryClusterRows fast clustering implementation for row clustering by binary sorting
#' @param labelMarkers label specific markers by name on heatmap (matches rownames of seMarker)
#' @param labelTop label the top features for each column in seMarker
#' @param labelRows label all rows
#' @param returnMat return final matrix that is used for plotting heatmap
#' @param ... additional args
#' @export
markerHeatmap <- function(
  seMarker, 
  FDR = 0.001, 
  log2FC = 0.1, 
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

  passMat <- SummarizedExperiment::assays(seMarker)[["Log2FC"]] >= log2FC & SummarizedExperiment::assays(seMarker)[["FDR"]] <= FDR
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
      pal <- paletteContinuous(set = "viridis", n = 100)
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
    ...
  )

  if(returnMat){
    return(mat)
  }else{
    return(0)
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

.binarySort <- function(m, scale = FALSE, cutOff = 1, lmat = NULL){

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
  hc <- hclust(dist(m))
  colIdx <- hc$order
  m <- t(m[colIdx,])
  lmat <- t(lmat[colIdx,])

  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  return(list(mat = m, hclust = hc))

}










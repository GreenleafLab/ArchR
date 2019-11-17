#' Plot ArchR Region Track
#' 
#' This function will plot the coverage at an input region
#'
#' @param ArchRProj ArchRProject
#' @param region GRanges region that will be plotted in (if more that one first will be selected)
#' @param groupBy use groupings for bulk/scTrack
#' @param useGroups select a subset of groups for plotting
#' @param useCoverages use group coverages for track plotting
#' @param plotSummary summary of region track to be plotted
#' @param sizes sizes corresponding to plotSummary
#' @param features GRanges features to be plotted (ie getPeakSet(ArchRProj))
#' @param geneSymbol if region is null plotting can be centered at gene start site corresponding to the gene symbol
#' @param upstream bp upstream of geneStart to extend
#' @param downstream bp downstream of geneStart to extend
#' @param tileSize with of tiles to plot bulk/scTrack
#' @param normMethod normMethod normalization column in cellColData to normalize bulkTrack
#' @param threads number of threads for parallel execution
#' @param ylim y-limits for bulkTrack
#' @param baseSize size of font in plot
#' @param borderWidth border width in plot
#' @param tickWidth axis tick width in plot
#' @param geneAnno geneAnnotation for geneTrack
#' @param title verbose sections
#' @param ... additional args
#' @export
ArchRRegionTrack <- function(
  ArchRProj, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL,
  useCoverages = FALSE,
  plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
  sizes = c(10, 2, 4),
  features = getPeakSet(ArchRProj),
  geneSymbol = NULL,
  upstream = 50000,
  downstream = 50000,
  tileSize = 100, 
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnno = getGeneAnnotation(ArchRProj),
  title = "",
  ...
  ){
  
  tstart <- Sys.time()

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  .messageDiffTime("Validating Region", tstart)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnno$genes
      region <- region[which(tolower(mcols(region)$symbol) == tolower(geneSymbol))]
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGRanges(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)[1]
  plotList <- list()

  ##########################################################
  # Bulk Tracks
  ##########################################################
  if("bulktrack" %in% tolower(plotSummary)){
    .messageDiffTime("Adding Bulk Tracks", tstart)
    plotList$bulktrack <- .bulkTracks(
      ArchRProj = ArchRProj, 
      region = region, 
      tileSize = tileSize, 
      groupBy = groupBy,
      threads = threads, 
      ylim = ylim,
      baseSize = baseSize,
      borderWidth = borderWidth,
      tickWidth = tickWidth,
      facetbaseSize = facetbaseSize,
      normMethod = normMethod,
      geneAnno = geneAnno,
      title = title,
      useGroups = useGroups,
      useCoverages = useCoverages,
      tstart = tstart) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
  }
  
  ##########################################################
  # Feature Tracks
  ##########################################################
  if("featuretrack" %in% tolower(plotSummary)){
    .messageDiffTime("Adding Feature Tracks", tstart)
    if(!is.null(features)){
      plotList$featuretrack <- .featureTracks(
          features = features, 
          region = region, 
          hideX = TRUE, 
          title = "Peaks") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }
  }

  ##########################################################
  # Gene Tracks
  ##########################################################
  if("genetrack" %in% tolower(plotSummary)){
    .messageDiffTime("Adding Gene Tracks", tstart)
    plotList$genetrack <- .geneTracks(
      geneAnnotation = geneAnno, 
      region = region, 
      title = "Genes") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
  }

  ##########################################################
  # Time to plot
  ##########################################################
  plotSummary <- tolower(plotSummary)
  sizes <- sizes[order(plotSummary)]
  plotSummary <- plotSummary[order(plotSummary)]

  nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
  if(any(nullSummary)){
    sizes <- sizes[-which(nullSummary)]
  }

  .messageDiffTime("Plotting", tstart)
  ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE)

}

#######################################################
# Bulk Aggregated ATAC Track Methods
#######################################################

.bulkTracks <- function(
  ArchRProj, 
  region = NULL, 
  tileSize = 100, 
  groupBy = "Clusters",
  useGroups = NULL,
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnno = getGeneAnnotation(ArchRProj),
  title = "",
  useCoverages = TRUE,
  tstart = NULL,
  ...
  ){

  .requirePackage("ggplot2")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  if(useCoverages){
    df <- .groupRegionSumCoverages(
        ArchRProj = ArchRProj, 
        groupBy = groupBy, 
        normMethod = normMethod,
        region = region, 
        tileSize = tileSize, 
        verbose = verbose
      )
  }else{
    df <- .groupRegionSumArrows(
      ArchRProj = ArchRProj, 
      groupBy = groupBy, 
      normMethod = normMethod,
      region = region, 
      tileSize = tileSize, 
      verbose = verbose
    )
  }

  ######################################################
  # Plot Track
  ######################################################
  if(!is.null(ylim)){
    ylim <- quantile(df$y, ylim)
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }else{
    ylim <- c(0,quantile(df$y, probs=c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }
  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  df$group <- factor(df$group, levels = uniqueGroups)
  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  pal <- suppressWarnings(paletteDiscrete(values = uniqueGroups))
  
  #Plot Track
  p <- ggplot(df, aes_string("x","y", color = "group", fill = "group")) + 
    geom_area(stat = "identity") + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage (Normalized ATAC Insertions Range %s - %s by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) +
    guides(fill = FALSE, colour = FALSE) + ggtitle(title)

  p

}

##############################################################################
# Create Average Tracks from Coverages
##############################################################################
.groupRegionSumCoverages <- function(ArchRProj, groupBy, useGroups = NULL, region, tileSize, normMethod, verbose){

  coverageMetadata <- .getCoverageMetadata(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    useGroups = useGroups
  )

  cellGroups <- .getCoverageParams(
    ArchRProj = ArchRProj, 
    groupBy = groupBy
  )[["cellGroups"]] %>% unlist
  
  groupRegionRle <- .groupRegionCoverages(
    coverageMetadata = coverageMetadata, 
    region = region, 
    tileSize = tileSize, 
    buffer = tileSize * 5, 
    threads = threads
  )
  groupNames <- names(groupRegionRle)
  
  #Normalization 
  g <- names(unlist(cellGroups, use.names = TRUE))
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = FALSE)[unlist(cellGroups),]
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = FALSE)[unlist(cellGroups),]
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }
  
  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors

  #Normalize
  groupRegionRle <- lapply(seq_along(groupRegionRle), function(x){
      groupRegionRle[[x]] * scaleFactors[names(groupRegionRle)[x]]
  })

  #Group And Average
  groupRegionRle <- split(groupRegionRle, coverageMetadata$Group)
  groupRegionList <- lapply(seq_along(groupRegionRle), function(x){
    Reduce("+", groupRegionRle[[x]]) / length(groupRegionRle[[x]])
  })
  names(groupRegionList) <- names(groupRegionRle)

  #Tile Region
  tileSize <- floor(tileSize / 2)
  regionTiles <- seq(trunc(start(region) / tileSize) - 1, trunc(end(region) / tileSize) + 1) * tileSize

  plotDF <- lapply(seq_along(groupRegionList), function(x){
    data.frame(x = regionTiles, y = as.vector(groupRegionList[[x]][regionTiles]), group = names(groupRegionList)[x])
  }) %>% Reduce("rbind", .)

  plotDF

}

.groupRegionCoverages <- function(coverageMetadata, region, tileSize = 100, buffer = 1000, threads = 1){
  
  region <- .validGRanges(region[1])
  coverageFiles <- coverageMetadata$File
  names(coverageFiles) <- coverageMetadata$Name
  
  covList <- .safelapply(seq_along(coverageFiles), function(x){
    .getCoverageFromRegion(coverageFiles[x], region, tileSize, buffer)
  }, threads = threads) %>% {as(.,"RleList")}
  names(covList) <- names(coverageFiles)
  
  covList

}

.getCoverageFromRegion <- function(coverageFile, region, tileSize, buffer){
  chr <- as.character(seqnames(region))
  cov <- Rle(
    lengths = h5read(coverageFile, paste0("Coverage/",chr,"/Lengths")), 
    values = h5read(coverageFile, paste0("Coverage/",chr,"/Values"))
  )
  w <- sum(runLength(cov))
  idx <- cumsum(runLength(cov)) %>% {which(.  >= start(region) - buffer & . <= end(region) + buffer)}
  runValue(cov)[-idx] <- 0
  covRanges <- ranges(cov)
  mcols(covRanges)$values <- runValue(cov)
  covRanges <- covRanges[mcols(covRanges)$values > 0]
  start(covRanges) <- trunc(start(covRanges) / tileSize) * tileSize
  end(covRanges) <- (trunc(end(covRanges) / tileSize) + 1) * tileSize - 1
  o <- coverage(covRanges, weight = mcols(covRanges)$values, width = w)
}


##############################################################################
# Create Average Tracks from Arrows
##############################################################################
.groupRegionSumArrows <- function(ArchRProj, groupBy, region, tileSize, normMethod, verbose){

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  tabGroups <- table(cellGroups)
  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    gmi <- .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
    )
  }, threads = threads) %>% Reduce("+" , .)

  #Plot DF
  df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
  df$y <- groupMat[cbind(df[,1], df[,2])]

  #Minus 1 Tile Size
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0

  #Plus 1 Size
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0

  #Create plot DF
  df <- rbind(df, dfm1, dfp1)
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0,]
  df$x <- regionTiles[df$row]
  df$group <- uniqueGroups[df$col]

  #Add In Ends
  dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
  )

  dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
  )
  
  #Final output
  plotDF <- rbind(df,dfs,dfe)
  plotDF <- df[order(df$group,df$x),]
  plotDF <- df[,c("x", "y", "group")]
  
  #Normalization 
  g <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }

  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors
  matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
  plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])

  return(plotDF)

}

.regionSumArrows <- function(ArrowFile, region, regionTiles, tileSize, cellNames, cellGroups, uniqueGroups){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = FALSE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = c(matchID[ids], matchID[ide]),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 0] <- 1

  #Create Group Matrix
  groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
  colnames(groupMat) <- uniqueGroups
  uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
  for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
  }

  return(groupMat)

}

#######################################################
# Gene Tracks
#######################################################
.geneTracks <- function(
  geneAnnotation, 
  region, 
  baseSize = 9, 
  borderWidth = 0.4, 
  title = "Genes",
  geneWidth = 2, 
  exonWidth = 4, 
  labelSize = 2,
  colorMinus = "dodgerblue2",
  colorPlus = "red",
  ...
  ){

  .requirePackage("ggplot2")
  .requirePackage("ggrepel")

  #only take first region
  region <- ArchR::.validGRanges(region)
  region <- subsetSeqnames(region[1],as.character(seqnames(region[1])))

  genes <- sort(sortSeqlevels(geneAnnotation$genes), ignore.strand = TRUE)
  exons <- sort(sortSeqlevels(geneAnnotation$exons), ignore.strand = TRUE)
  genesO <- data.frame(subsetByOverlaps(genes, region, ignore.strand = TRUE))

  if(nrow(genesO) > 0){

    #Identify Info for Exons and Genes
    exonsO <- data.frame(subsetByOverlaps(exons, region, ignore.strand = TRUE))
    exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol),]
    genesO$facet = title
    genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, start(region)))
    genesO$end <- matrixStats::rowMins(cbind(genesO$end, end(region)))

    #Collapse Iteratively
    #backwards iteration so that the last value chosen is the lowest cluster possible to fit in.
    genesO$cluster <- 0
    for(i in seq_len(nrow(genesO))){
      if(i==1){
        genesO$cluster[i] <- 1
      }else{
        for(j in seq_len(max(genesO$cluster))){
          jEnd <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], rev(genesO$cluster))]
          if(genesO$start[i] > jEnd + median(genesO$width)){
            genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j],rev(genesO$cluster))]
          }
        }
        if(genesO$cluster[i]==0){
          genesO$cluster[i] <- genesO$cluster[i-1] + 1
        }
      }
    }
    exonsO$cluster <- genesO$cluster[match(exonsO$symbol, genesO$symbol)]
    pal <- c("-"=colorMinus,"+"=colorPlus,"*"=colorPlus)
    
    p <- ggplot(data = genesO, aes(color = strand, fill = strand)) +
      facet_grid(facet~.) +
      #################################################
      #Limits
      #################################################
      ylim(c(0.5, max(genesO$cluster) + 0.5)) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
      #################################################
      #Segment for Not Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)!="-"),], 
        aes(x = start, xend = end, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segment for Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)=="-"),], 
        aes(x = end, xend = start, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segement for Exons
      #################################################
      geom_segment(data = exonsO, aes(x = start, xend = end, y = cluster, 
        yend = cluster, color = strand),size=exonWidth) +
      #################################################
      #Colors
      #################################################
      scale_color_manual(values = pal, guide = FALSE) + 
      scale_fill_manual(values = pal) +
      #################################################
      #Theme
      #################################################
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(angle = 0)) +
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = FALSE) + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand!="-"),], 
        aes(x = start, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          size = labelSize, direction = "x")
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand=="-"),], 
        aes(x = end, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          size = labelSize, direction = "x")
    }

    p <- p + theme(legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA), legend.position="none")

  }else{

    #create empty plot
    df <- data.frame(facet = "GeneTrack", start = 0, end = 0, strand = "*", symbol = "none")
    pal <- c("*"=colorPlus)
    p <- ggplot(data = df, aes(start, end, fill = strand)) + geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  return(p)

}

#######################################################
# Feature Tracks
#######################################################
.featureTracks <- function(
  features, 
  region, 
  title = "FeatureTrack", 
  pal = NULL,
  baseSize = 9, 
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  ...
  ){

  .requirePackage("ggplot2")

  #only take first region
  region <- ArchR::.validGRanges(region)
  region <- subsetSeqnames(region[1],as.character(seqnames(region[1])))

  if(!inherits(features,"GRangesList") & !inherits(features,"GenomicRangesList")){
    features <- ArchR::.validGRanges(features)
    featureList <- GenomicRanges::GenomicRangesList(features)
    names(featureList) <- "FeatureTrack"
    hideY <- TRUE
  }else{
    featureList <- features
    hideY <- FALSE
  }
  featureList <- featureList[rev(seq_along(featureList))]

  featureO <- lapply(seq_along(featureList), function(x){
    featurex <- featureList[[x]]
    namex <- names(featureList)[x]
    mcols(featurex) <- NULL
    sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
    if(length(sub) > 0){
      data.frame(sub, name = namex)
    }else{
      empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
      data.frame(empty, name = namex)
    }

  })

  featureO <- Reduce("rbind", featureO)
  featureO$facet <- title

  if(is.null(pal)){
    pal <- paletteDiscrete(set = "stallion", rev(unique(paste0(featureO$name))))
  }

  p <- ggplot(data = featureO, aes(color = name)) +
    facet_grid(facet~.) +
    geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
    ylab("") + xlab("") + 
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_color_manual(values = pal) +
    theme(legend.text = element_text(size = baseSize)) + 
    theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
    guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(angle = 0), strip.background = element_blank())

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  return(p)

}



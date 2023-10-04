#' Plot a TSS Enrichment Plot for Each Sample
#' 
#' This function will plot a TSS enrichment plot for each sample. Cells in `ArchRProject` are the only ones
#' used when making this plot.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param chromSizes A GRanges object of the chromosome lengths. See `getChromSizes` for more info.
#' @param TSS A `GRanges` object containing the locations of stranded TSS regions. The default behavior is to try to retrieve
#' this information from the `geneAnnotation` stored in the `ArchRProject`.
#' @param flank An integer that specifies how far in bp (+/-) to extend the TSS for plotting.
#' @param norm An integer that specifies the number of base pairs from the ends of the flanks to be used for normalization. 
#' For example if `flank=2000` and `norm=100`, the TSS insertions will be normalized by +/- 1900-2000 bp from the TSS.
#' @param smooth An integer that indicates the smoothing window (in basepairs) to be applied to the TSS plot.
#' @param pal A color palette representing the groups from groupBy in TSS plot.
#' @param returnDF A boolean value that indicates whether to return a `data.frame` containing the plot information
#' instead of plotting the TSS enrichment plot.
#' @param threads An integer specifying the number of threads to use for calculation. By default this uses the number of threads set by `addArchRThreads()`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotTSSEnrichment <- function(
  ArchRProj = NULL,
  groupBy = "Sample",
  chromSizes = getChromSizes(ArchRProj),
  TSS = getTSS(ArchRProj),
  flank = 2000, 
  norm = 100, 
  smooth = 11,
  pal = NULL,
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotTSSEnrichment")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = TSS, name = "TSS", valid = c("granges"))
  .validInput(input = flank, name = "flank", valid = c("integer"))
  .validInput(input = norm, name = "norm", valid = c("integer"))
  .validInput(input = smooth, name = "smooth", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotTSSEnrichment Input-Parameters", logFile = logFile)

  chr <- paste0(seqnames(chromSizes))
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  .logThis(chr, paste0("chr"), logFile = logFile)
  TSS <- sort(sortSeqlevels(TSS))
  splitTSS <- split(GenomicRanges::resize(TSS,1,"start"), seqnames(TSS))[chr]
  .logThis(splitTSS, paste0("splitTSS"), logFile = logFile)
  window <- 2 * flank + 1
  groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = FALSE)
  uniqGroups <- gtools::mixedsort(unique(groups[,1]))

  if(threads > 1){
     h5disableFileLocking()
  }

  dfTSS <- .safelapply(seq_along(uniqGroups), function(z){

    .logDiffTime(paste0(uniqGroups[z], " Computing TSS (",z," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    cellx <- rownames(groups)[which(paste0(groups[,1]) == uniqGroups[z])]

    for(k in seq_along(chr)){

      #TSS for Chr
      TSSi <- splitTSS[[chr[k]]]

      #Set TSS To be a dummy chr1
      TSSi <- GRanges(seqnames=rep("chr1",length(TSSi)), ranges = ranges(TSSi), strand = strand(TSSi))
      .logThis(TSSi, paste0(uniqGroups[z], " : TSSi : ", chr[k]), logFile = logFile)

      #Extract Fragments
      covi <- suppressMessages(getFragmentsFromProject(
        ArchRProj = ArchRProj,
        subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[k]],
        cellNames = cellx,
        logFile = logFile
      ) %>% unlist(use.names = FALSE))
      .logThis(covi, paste0(uniqGroups[z], " : Fragments : ", chr[k]), logFile = logFile)
     
      #Get Insertions
      covi <- sort(c(start(covi), end(covi)))
      .logThis(covi, paste0(uniqGroups[z], " : Insertions : ", chr[k]), logFile = logFile)

      #IRanges
      covi <- IRanges(start = covi, width = 1)
      .logThis(covi, paste0(uniqGroups[z], " : Insertions2 : ", chr[k]), logFile = logFile)

      #Coverage
      covi <- IRanges::coverage(covi)
      .logThis(covi, paste0(uniqGroups[z], " : Cov : ", chr[k]), logFile = logFile)

      #Compute Sum
      sumTSSi <- rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
      .logThis(sumTSSi, paste0(uniqGroups[z], " : SumTSS 1 : ", chr[k]), logFile = logFile)

      if(k == 1){
        sumTSS <- sumTSSi
      }else{
        sumTSS <- sumTSS + sumTSSi
      }
      .logThis(sumTSS, paste0(uniqGroups[z], " : SumTSS : ", chr[k]), logFile = logFile)

    }

    normBy <- mean(sumTSS[c(1:norm,(flank*2-norm+1):(flank*2+1))])

    df <- DataFrame(
      group = uniqGroups[z],
      x = seq_along(sumTSS) - flank - 1, 
      value = sumTSS, 
      normValue = sumTSS / normBy,
      smoothValue = .centerRollMean(sumTSS/normBy, 11)
    )

    .logThis(df, paste0(uniqGroups[z], " : TSSDf"), logFile = logFile)

    .logDiffTime(paste0(uniqGroups[z], " Finished Computing TSS (",z," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    df

  }, threads = threads) %>% Reduce("rbind", .)

  .logThis(dfTSS, paste0("All : TSSDf"), logFile = logFile)

  .endLogging(logFile = logFile)

  if(threads > 1){
    h5enableFileLocking()
  }
  
  if(returnDF){
    
    return(dfTSS)

  }else{

    plotDF <- data.frame(x=dfTSS$x,v=dfTSS$smoothValue,group=dfTSS$group)
    plotDF <- plotDF[sort(unique(c(1,seq(1,nrow(plotDF),11),nrow(plotDF)))), , drop = FALSE]
    
    if(is.null(pal)){
      pal <- paletteDiscrete(values=unique(plotDF$group))
    }

    p <- ggplot(plotDF, aes(x,v,color=group)) +
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("Distance From Center (bp)") +
      ylab("Normalized Insertion Profile") +
      scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0))

    p

  }

}


#' Plot the fragment size distribution for each sample
#' 
#' This function will plot a fragment size distribution for each sample. Only cells in the `ArchRProject` are used when making this plot.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param chromSizes A GRanges object of the chromosome lengths. See `getChromSizes` for more info.
#' @param maxSize The maximum fragment size (in basepairs) to be included when plotting the fragment size distribution.
#' @param pal A color palette representing the groups from groupBy in fragment size plot.
#' @param returnDF A boolean value that indicates whether to return a `data.frame` containing the plot information
#' instead of plotting the fragment size distribution.
#' @param threads An integer specifying the number of threads to use for calculation. By default this uses the number of threads set by `addArchRThreads()`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotFragmentSizes <- function(
  ArchRProj = NULL,
  groupBy = "Sample",
  chromSizes = getChromSizes(ArchRProj),
  maxSize = 750,
  pal = NULL,
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotFragmentSizes")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = maxSize, name = "maxSize", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotFragmentSizes Input-Parameters", logFile = logFile)

  chr <- paste0(seqnames(chromSizes))
  groups <- getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = FALSE)
  uniqGroups <- gtools::mixedsort(unique(groups[,1]))

  if(threads > 1){
     h5disableFileLocking()
  }

  dfFS <- .safelapply(seq_along(uniqGroups), function(x){

    .logDiffTime(paste0(uniqGroups[x], " Computing FragmentSizes (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    cellx <- rownames(groups)[which(paste0(groups[,1]) == uniqGroups[x])]

    for(i in seq_along(chr)){
      if(i == 1){
        fsi <- unlist(suppressMessages(getFragmentsFromProject(
          ArchRProj = ArchRProj,
          subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[i]],
          cellNames = cellx,
          logFile = logFile
        )), use.names=FALSE) %>% width %>% tabulate(nbins = maxSize)
      }else{
        fsi <- fsi + unlist(suppressMessages(getFragmentsFromProject(
          ArchRProj = ArchRProj,
          subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% chr[i]],
          cellNames = cellx,
          logFile = logFile
        )), use.names=FALSE) %>% width %>% tabulate(nbins = maxSize)
      }
      .logThis(fsi, paste0(uniqGroups[x], " : FragSizes : ", chr[i]), logFile = logFile)
    }

    df <- DataFrame(
      group = uniqGroups[x],
      fragmentSize = seq_along(fsi), 
      fragmentPercent = round(100*fsi/sum(fsi),4)
    )

    .logThis(df, paste0(uniqGroups[x], " : Frag DF"), logFile = logFile)

    .logDiffTime(paste0(uniqGroups[x], " Finished Computing FragmentSizes (",x," of ",length(uniqGroups),")!"), t1 = tstart, logFile = logFile)

    df

  }, threads = threads) %>% Reduce("rbind", .)

  .logThis(dfFS, paste0("All : FragSizes DF"), logFile = logFile)

  .endLogging(logFile = logFile)

  if(threads > 1){
    h5enableFileLocking()
  }

  if(returnDF){
    
    return(dfFS)

  }else{

    plotDF <- data.frame(dfFS)
    
    if(is.null(pal)){
      pal <- paletteDiscrete(values=unique(plotDF$group))
    }

    p <- ggplot(plotDF, aes(fragmentSize, fragmentPercent,color=group)) + 
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("ATAC-seq Fragment Size (bp)") +
      ylab("Percentage of Fragments") +
      scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$fragmentPercent)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$fragmentSize), max(plotDF$fragmentSize)), expand = c(0,0))

    p

  }

}




#' Plot fragment size distribution for each sample JJJ
#' 
#' This function will plot a TSS enrichment plot for each sample. Cells in `ArchRProject` are the only ones
#' used when making this plot.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param TSS A GRanges object containing stranded TSS ranges. By default will try to get this info from geneAnnotation stored in `ArchRProj`.
#' @param flank A numeric specifying how far in bp (+/-) to extend TSS for plotting.
#' @param norm A numeric specifying the number of base pairs from the ends of the flanks to be used for normalization. 
#' For example if `flank=2000` and `norm=100`, the TSS insertions will be normalized by +/- 1900-2000 bp from the TSS.
#' @param smooth A numeric describing the smoothing window to be applied to TSS plot.
#' @param returnDF A boolean for return plot as a `data.frame` instead.
#' @param threads An integer specifying the number of threads to use for calculation. By default this uses threads set by `addArchRThreads`.
#' @export
plotTSSEnrichment <- function(
  ArchRProj = NULL,
  TSS = getTSS(ArchRProj),
  flank = 2000, 
  norm = 100, 
  smooth = 11,
  returnDF = FALSE,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = TSS, name = "TSS", valid = c("granges"))
  .validInput(input = flank, name = "flank", valid = c("integer"))
  .validInput(input = norm, name = "norm", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  tstart <- Sys.time()

  ArrowFiles <- getArrowFiles(ArchRProj)
  chr <- .availableChr(ArrowFiles)
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  cellNames <- ArchRProj$cellNames
  TSS <- sort(sortSeqlevels(TSS))
  splitTSS <- split(resize(TSS,1,"start"), seqnames(TSS))[chr]
  window <- 2 * flank + 1

  dfTSS <- .safelapply(seq_along(ArrowFiles), function(x){

    .messageDiffTime(paste0(names(ArrowFiles)[x], " Computing TSS (",x," of ",length(ArrowFiles),")!"), tstart)

    for(i in seq_along(chr)){

      TSSi <- splitTSS[[chr[i]]]

        covi <- .getFragsFromArrow(
          ArrowFile = ArrowFiles[x], 
          chr = chr[i], 
          out = "IRanges", 
          cellNames = cellNames
        ) %>% {coverage(IRanges(c(start(.), end(.)), width = 1))}

        if(i == 1){
          sumTSS <- ArchR:::rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
        }else{
          sumTSS <- sumTSS + ArchR:::rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
        }
        
    }

    normBy <- mean(sumTSS[c(1:norm,(flank*2-norm+1):(flank*2+1))])

    df <- DataFrame(
      sampleName = names(ArrowFiles)[x],
      x = seq_along(sumTSS) - flank - 1, 
      value = sumTSS, 
      normValue = sumTSS/normBy,
      smoothValue = ArchR:::.centerRollMean(sumTSS/normBy, 11)
    )

    .messageDiffTime(paste0(names(ArrowFiles)[x], " Finished Computing TSS (",x," of ",length(ArrowFiles),")!"), tstart)

    df

  }, threads = threads) %>% Reduce("rbind", .)



  if(returnDF){
    
    return(dfTSS)

  }else{

    plotDF <- data.frame(x=dfTSS$x,v=dfTSS$smoothValue,sampleName=dfTSS$sampleName)
    
    p <- ggplot(plotDF, aes(x,v,color=sampleName)) + 
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("Distance From Center (bp)") +
      ylab("Normalized Insertion Profile") +
      scale_color_manual(values=paletteDiscrete(values=unique(plotDF$sampleName))) +
      scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0))

    p

  }

}


#' Plot fragment size distribution for each sample JJJ
#' 
#' This function will plot a fragment size distribution for each sample. Cells in `ArchRProject` are the only ones
#' used when making this plot.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param nbp Number of basepairs to plot fragment size distribution.
#' @param returnDF A boolean for return plot as a `data.frame` instead.
#' @param threads An integer specifying the number of threads to use for calculation. By default this uses threads set by `addArchRThreads`.
#' @export
plotFragmentSizes <- function(
  ArchRProj = NULL,
  nbp = 750,
  returnDF = FALSE,
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = nbp, name = "nbp", valid = c("integer"))
  .validInput(input = returnDF, name = "returnDF", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))

  tstart <- Sys.time()

  ArrowFiles <- getArrowFiles(ArchRProj)
  chr <- gtools::mixedsort(.availableChr(ArrowFiles))
  cellNames <- ArchRProj$cellNames

  dfFS <- .safelapply(seq_along(ArrowFiles), function(x){

    .messageDiffTime(paste0(names(ArrowFiles)[x], " Computing FragmentSizes (",x," of ",length(ArrowFiles),")!"), tstart)

    for(i in seq_along(chr)){

      if(i == 1){
          fsi <- .getFragsFromArrow(
            ArrowFile = ArrowFiles[x], 
            chr = chr[i], 
            out = "IRanges", 
            cellNames = cellNames
          ) %>% width %>% tabulate(nbins = nbp)
      }else{
          fsi <- fsi + .getFragsFromArrow(
            ArrowFile = ArrowFiles[x], 
            chr = chr[i], 
            out = "IRanges", 
            cellNames = cellNames
          ) %>% width %>% tabulate(nbins = nbp)
      }

    }

    df <- DataFrame(
      sampleName = names(ArrowFiles)[x],
      fragmentSize = seq_along(fsi), 
      fragmentPercent = round(100*fsi/sum(fsi),4)
    )

    .messageDiffTime(paste0(names(ArrowFiles)[x], " Finished Computing FragmentSizes (",x," of ",length(ArrowFiles),")!"), tstart)

    df

  }, threads = threads) %>% Reduce("rbind", .)



  if(returnDF){
    
    return(dfFS)

  }else{

    plotDF <- data.frame(dfFS)
    
    p <- ggplot(plotDF, aes(fragmentSize, fragmentPercent,color=sampleName)) + 
      geom_line(size = 1) +
      theme_ArchR() +
      xlab("ATAC-seq Fragment Size (bp)") +
      ylab("Percentage of Fragments") +
      scale_color_manual(values=paletteDiscrete(values=unique(plotDF$sampleName))) +
      scale_y_continuous(limits = c(0, max(plotDF$fragmentPercent)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$fragmentSize), max(plotDF$fragmentSize)), expand = c(0,0))

    p

  }

}


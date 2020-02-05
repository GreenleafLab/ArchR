#To be added and documented JJJ


.plotTSSEnrichment <- function(
  ArchRProj = NULL,
  TSS = getTSS(ArchRProj),
  chromSizes = getChromSizes(ArchRProj),
  flank = 2000, 
  norm = 100, 
  smooth = 11,
  returnDF = FALSE,
  threads = getArchRThreads()
  ){

  tstart <- Sys.time()

  ArrowFiles <- getArrowFiles(ArchRProj)
  chr <- .availableChr(ArrowFiles)
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  cellNames <- proj$cellNames
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


.plotFragmentSizes <- function(
  ArchRProj = NULL,
  nbins = 750,
  returnDF = FALSE,
  threads = getArchRThreads()
  ){

  tstart <- Sys.time()

  ArrowFiles <- getArrowFiles(ArchRProj)
  chr <- .availableChr(ArrowFiles)
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  cellNames <- proj$cellNames

  dfFS <- .safelapply(seq_along(ArrowFiles), function(x){

    .messageDiffTime(paste0(names(ArrowFiles)[x], " Computing FragmentSizes (",x," of ",length(ArrowFiles),")!"), tstart)

    for(i in seq_along(chr)){

      if(i == 1){
          fsi <- .getFragsFromArrow(
            ArrowFile = ArrowFiles[x], 
            chr = chr[i], 
            out = "IRanges", 
            cellNames = cellNames
          ) %>% width %>% tabulate(nbins = nbins)
      }else{
          fsi <- fsi + .getFragsFromArrow(
            ArrowFile = ArrowFiles[x], 
            chr = chr[i], 
            out = "IRanges", 
            cellNames = cellNames
          ) %>% width %>% tabulate(nbins = nbins)
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


##########################################################################################
# Transcription Factor Footprinting Methods
##########################################################################################

#' Calculate footprints from an ArchRProject
#' 
#' This function will get footprints for all samples in a given ArchRProject or a properly-formatted Summarized Experiment
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param positions A `list` or `GenomicRangesList` of `GRanges` containing the positions to incorporate into the footprint. Each position should be stranded.
#' @param plotName The prefix to add to the file name for the output PDF file containing the footprint plots.
#' @param groupBy The name of the column in `cellColData` used in the `addGroupCoverages()` function for grouping multiple cells together.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in `cellColData`.
#' This limits the groups used to perform footprinting.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for plotting the lines corresponding to the footprints.
#' @param flank The number of basepairs from the position center (+/-) to consider as the flank.
#' @param flankNorm The number of basepairs to consider at the edge of the flank region (+/-) to be used for footprint normalization.
#' @param smoothWindow The size in basepairs of the sliding window to be used for smoothing of the footprint signal.
#' @param minCells The minimum number of cells required in a given cell group to permit footprint generation.
#' @param nTop The number of genomic regions to consider. Only the top `nTop` genomic regions based on the "score" column in the `GRanges`
#' object will be considered for the footprint.
#' @param normMethod The name of the normalization method to use to normalize the footprint relative to the Tn5 insertion bias. Options
#' include "none", "subtract", "divide". "Subtract" means subtracting the normalized Tn5 Bias. "Divide" means dividing the normalized Tn5 Bias.
#' @param inputSE Input a previous footprint Summarized Experiment (returned after running `plotFootprints`) to be plotted instead of
#' regenerating the footprinting information.
#' @param height The height in inches to be used for the output PDF.
#' @param width The width in inches to be used for the output PDF file.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for
#' preventing overwritting of old plots.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @export
getFootprints <- function(
  ArchRProj = NULL,
  positions = NULL,
  plotName = "Plot-Footprints",
  groupBy = "Clusters",
  useGroups = NULL,
  flank = 250,
  minCells = 25,
  nTop = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getFootprints")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = positions, name = "positions", valid = c("grangeslist"))
  .validInput(input = plotName, name = "plotName", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)

  #####################################################
  # Compute Kmer Frequency Table 
  #####################################################
  coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy, minCells = minCells)
  coverageParams <- .getCoverageParams(ArchRProj = ArchRProj, groupBy = groupBy)
  kmerLength <- coverageParams$kmerLength

  .logThis(coverageMetadata, "coverageMetadata", logFile = logFile)
  .logThis(coverageParams, "coverageParams", logFile = logFile)

  if(!is.null(useGroups)){
    if(sum(coverageMetadata[,1] %in% useGroups) == 0){
      stop("No Groups found matching useGroups!")
    }
    coverageMetadata <- coverageMetadata[coverageMetadata[,1] %in% useGroups,]
  }

  genome <- getGenome(ArchRProj)
  .requirePackage(genome)
  .requirePackage("Biostrings", source = "bioc")
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)

  .logDiffTime("Computing Kmer Bias Table", tstart, verbose = verbose, logFile = logFile)
  kmerTableList <- .kmerPositionFrequency(
    featureList = positions, 
    genome = BSgenome, 
    flank = flank,
    k = kmerLength, 
    threads = 1,
    verbose = FALSE,
    logFile = logFile
  )

  #####################################################
  # Compute Footprints
  #####################################################
  .logDiffTime("Computing Footprints", tstart, verbose = verbose, logFile = logFile)
  footprintList <- .computeFootprints(
    featureList = positions, 
    coverageFiles = coverageMetadata$File, 
    flank = flank, 
    threads = threads, 
    verbose = FALSE,
    logFile = logFile
  )
  
  #####################################################
  # Compute Bias For Footprints
  #####################################################
  .logDiffTime("Computing Footprints Bias", tstart, verbose = verbose, logFile = logFile)
  footprintBiasList <- .computeFootprintsBias(
    kmerTableList = kmerTableList, 
    coverageFiles = coverageMetadata$File, 
    threads = threads, 
    verbose = FALSE
  )

  #####################################################
  # Summarize into SE
  #####################################################
  .logDiffTime("Summarizing Footprints", tstart, verbose = verbose, logFile = logFile)
  footAssay <- lapply(seq_along(positions), function(x){
    footMat <- lapply(seq_along(footprintList), function(y){
      footprintList[[y]][,x]
    }) %>% Reduce("cbind", .)
    colnames(footMat) <- coverageMetadata$Name
    biasMat <- lapply(seq_along(footprintBiasList), function(y){
      footprintBiasList[[y]][,x]
    }) %>% Reduce("cbind", .)
    colnames(biasMat) <- coverageMetadata$Name
    rbind(footMat, biasMat)
  }) %>% SimpleList
  names(footAssay) <- names(positions)

  #Clean GC
  rm(footprintList, footprintBiasList)
  gc()

  rowData <- DataFrame(
    x = c(seq(-flank, flank), seq(-flank, flank)), 
    type = c(rep("footprint", flank*2+1),rep("bias", flank*2+1))
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = footAssay, 
    colData = coverageMetadata,
    rowData = rowData
  )

  metadata(se)$Params <- SimpleList(kmerLength=kmerLength,flank=flank,date=Sys.Date())

  return(se)


}


#####################################################################################################
# Helpers for get Footprints
#####################################################################################################

.computeFootprintsBias <- function(
  kmerTableList = NULL, 
  coverageFiles = NULL, 
  threads = 1, 
  verbose = TRUE,
  logFile = NULL
  ){
  tstart <- Sys.time()
  out <- .safelapply(seq_along(coverageFiles), function(i){
      .logDiffTime(sprintf("Computing Footprints Bias %s of %s:", i, length(coverageFiles)),tstart,verbose=verbose, logFile = logFile)
      .computeFootprintsBiasSingle(kmerTableList, coverageFiles[i])
  }, threads = threads) %>% SimpleList
  return(out)
}

.computeFootprintsBiasSingle <- function(
  kmerTableList = NULL, 
  coverageFile = NULL,
  logFile = NULL
  ){
  kmerTableList <- as(kmerTableList, "list")
  oe <- h5read(coverageFile, "KmerBias/ObservedKmers") / h5read(coverageFile, "KmerBias/ExpectedKmers")
  names(oe) <- h5read(coverageFile, "KmerBias/Kmer")
   biasDF <- lapply(seq_along(kmerTableList), function(x){
    bias <- colSums(as.matrix(kmerTableList[[x]]) * as.vector(oe[rownames(kmerTableList[[x]])]))
    bias <- bias / sum(bias)
    bias    
  }) %>% Reduce("cbind", .) %>% data.frame
  gc()
  biasDF    
}

.computeFootprints <- function(
  featureList = NULL, 
  coverageFiles = NULL, 
  flank = 250, 
  threads = 1, 
  verbose = TRUE,
  logFile = NULL
  ){
  tstart <- Sys.time()
  out <- .safelapply(seq_along(coverageFiles), function(i){
    .computeFootprintsSingle(featureList, coverageFiles[i], flank, gc = TRUE, 
      prefix = sprintf("Computing Footprints %s of %s:", i, length(coverageFiles)),
      tstart = tstart, 
      verbose = verbose
      )
  }, threads = threads) %>% SimpleList
  return(out)
}

.computeFootprintsSingle <- function(
  featureList = NULL, 
  coverageFile = NULL, 
  flank = 250, 
  gc = FALSE, 
  prefix = "", 
  tstart = NULL, 
  verbose = TRUE,
  logFile = NULL
  ){
  window <- 2 * flank + 1
  featureNames <- names(featureList)
  featureList <- as(featureList, "list")
  allChr <- lapply(featureList, function(x) unique(as.character(seqnames(x)))) %>% unlist %>% unique %>% sort
  cov <- .getCoverageRle(coverageFile, allChr)
  footprintDF <- lapply(seq_along(featureList), function(x){
    outx <- tryCatch({
      
      featurex <- split(resize(featureList[[x]],1,"center"), seqnames(featureList[[x]]))
      intSeq <- intersect(names(featurex), names(cov))
      if(length(intSeq)==0){
        .logMessage(paste0("No intersecting chromsomes for feature ", names(featureList)[x], "!"))
        stop("No intersecting chromsomes for feature ", names(featureList)[x], "!")
      }
      outx <- rleSumsStranded(cov[intSeq], featurex[intSeq], window, as.integer) #Rcpp
      if(x %% 25 == 0 & gc){
        gc()
      }
      if(length(featureList) > 10){
        if(x %% 5 == 0){
          .logDiffTime(sprintf("%s %s Percent Completed", prefix, round(100 * x / length(featureList)),1), tstart, verbose=verbose, logFile = logFile)
        }
      }else{
          if(x == 1 | x == length(featureList)){
            .logDiffTime(sprintf("%s %s Percent Completed", prefix, round(100 * x / length(featureList)),1), tstart, verbose=verbose, logFile = logFile)
          }
      }
      outx
    
    }, error = function(e){

      errorList <- list(
        x = x,
        window = window,
        namex = if(exists("featurex", inherits = FALSE)) names(featureList)[x] else "namex",
        featurex = if(exists("featurex", inherits = FALSE)) featurex else "featurex",
        intSeq = if(exists("intSeq", inherits = FALSE)) intSeq else "intSeq",
        cov = if(exists("cov", inherits = FALSE)) cov else "cov"
      )

      .logError(e, fn = ".computeFootprintsSingle", info = basename(coverageFile), errorList = errorList, logFile = logFile)

    })
    outx
  }) %>% Reduce("cbind",.) %>% data.frame
  gc()
  footprintDF   
}

.getCoverageRle <- function(
  coverageFile = NULL, 
  allChr = NULL
  ){
  cov <- lapply(seq_along(allChr), function(x){
    Rle(
      lengths = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Lengths")), 
      values = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Values"))
    )
  }) %>% {as(.,"RleList")}
  names(cov) <- allChr
  cov
}

.kmerPositionFrequency <- function(
  featureList = NULL, 
  genome = NULL, 
  flank = 250, 
  k = 6, 
  threads = 1, 
  verbose = TRUE,
  logFile = NULL
  ){
  
  tstart <- Sys.time()
  genome <- validBSgenome(genome)
  window <- 2*flank + 1

  kmerList <- .safelapply(seq_along(featureList), function(i){
    .logDiffTime(sprintf("Computing Kmer Tables for %s of %s features", i, length(featureList)), tstart, verbose=verbose, logFile = logFile)
    bsv <- BSgenomeViews(genome , resize(featureList[[i]], window + k, "center"))
    bsv <- bsv[width(bsv) == window + k] #none that are trimmed!
    #BSgenome is already stranded
    #kmerPositionFrequencyCpp is Rcpp export for getting kmer position frequencies from strings
    kmerTable <- kmerPositionFrequencyCpp(as.character(bsv), rep(1L,length(bsv)), window, k, .getKmers(k)) #Rcpp
    return(kmerTable)
  }, threads = threads) %>% SimpleList
  names(kmerList) <- names(featureList)
  
  .logDiffTime("Finished Computing Kmer Tables", tstart)

  return(kmerList)
}

.getKmers <-function(
  k = NULL, 
  letters = c('A','C','G','T')
  ){
  kmers = ''
  for (i in seq_len(k)) {
    kmers <- unlist(lapply(kmers, function(x) paste0(x, letters)))
  }  
  return(kmers)
}


#####################################################################################################
# Plot Footprints
#####################################################################################################

#' Plot Footprints
#' 
#' This function will get footprints for all samples in a given ArchRProject or a properly-formatted Summarized Experiment
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param positions A `list` or `GenomicRangesList` of `GRanges` containing the positions to incorporate into the footprint. Each position should be stranded.
#' @param plotName The prefix to add to the file name for the output PDF file containing the footprint plots.
#' @param groupBy The name of the column in `cellColData` used in the `addGroupCoverages()` function for grouping multiple cells together.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in `cellColData`.
#' This limits the groups used to perform footprinting.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for plotting the lines corresponding to the footprints.
#' @param flank The number of basepairs from the position center (+/-) to consider as the flank.
#' @param flankNorm The number of basepairs to consider at the edge of the flank region (+/-) to be used for footprint normalization.
#' @param smoothWindow The size in basepairs of the sliding window to be used for smoothing of the footprint signal.
#' @param minCells The minimum number of cells required in a given cell group to permit footprint generation.
#' @param nTop The number of genomic regions to consider. Only the top `nTop` genomic regions based on the "score" column in the `GRanges`
#' object will be considered for the footprint.
#' @param normMethod The name of the normalization method to use to normalize the footprint relative to the Tn5 insertion bias. Options
#' include "none", "subtract", "divide". "Subtract" means subtracting the normalized Tn5 Bias. "Divide" means dividing the normalized Tn5 Bias.
#' @param inputSE Input a previous footprint Summarized Experiment (returned after running `plotFootprints`) to be plotted instead of
#' regenerating the footprinting information.
#' @param height The height in inches to be used for the output PDF.
#' @param width The width in inches to be used for the output PDF file.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for
#' preventing overwritting of old plots.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @export
plotFootprints <- function(
  seFoot = NULL,
  names = NULL,
  pal = NULL,
  flank = 250,
  flankNorm = 50,
  normMethod = "Subtract",
  smoothWindow = NULL,
  baseSize = 6,
  plot = TRUE,
  ArchRProj = NULL,
  plotName = paste0("Plot-Footprints-", normMethod),
  height = 6,
  width = 4,
  addDOC = TRUE,
  logFile = createLogFile("plotFootprints")
  ){

  tstart <- Sys.time()
  .startLogging(logFile = logFile)

  if(is.null(names)){
    names <- names(assays(seFoot))
  }

  if(length(names) > 25){
    if(!plot){
      if(force){
        .logMessage("Plotting more than 25 footprints can create large storage of ggplots. Continuing since force = TRUE", logFile = logFile)
        message("Plotting more than 25 footprints can create large storage of ggplots. Continuing since force = TRUE")
      }else{
        .logMessage("Plotting more than 25 footprints can create large storage of ggplots. Stopping since force = FALSE", logFile = logFile)
        stop("Plotting more than 25 footprints can create large storage of ggplots. Stopping since force = FALSE")
      }
    }
  }

  if(plot){

    name <- gsub("\\.pdf", "", plotName)
    if(is.null(ArchRProj)){
      outDir <- "Plots"
    }else{
      ArchRProj <- .validArchRProject(ArchRProj)
      outDir <- file.path(getOutputDirectory(ArchRProj), "Plots")
    }

    dir.create(outDir, showWarnings = FALSE)
    if(addDOC){
      doc <- gsub(":","-",stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2])
      filename <- file.path(outDir, paste0(name, "_Date-", Sys.Date(), "_Time-", doc, ".pdf"))
    }else{
      filename <- file.path(outDir, paste0(name, ".pdf"))
    }

    pdf(filename, width = width, height = height, useDingbats = FALSE)

  }

  ggList <- lapply(seq_along(names), function(x){

    .logDiffTime(sprintf("Plotting Footprint : %s (%s of %s)", names[x], x, length(names)), tstart, logFile = logFile)

    gg <- .ggFootprint(
      seFoot = seFoot,
      name = names[x],
      pal = pal,
      smoothWindow = smoothWindow,
      flank = flank,
      flankNorm = flankNorm,
      baseSize = baseSize,
      normMethod = normMethod,
      logFile = logFile
    )

    if(plot){
      if(x != 1){
        grid::grid.newpage()
      }
      grid::grid.draw(gg)
      return(0)
    }else{
      return(gg)
    }

  })
  
  .endLogging(logFile = logFile)

  if(!plot){
    names(ggList) <- names
    ggList
  }else{
    dev.off()
    return(invisible(0))
  }

}

.ggFootprint <- function(
  seFoot = NULL,
  name = NULL,
  pal = NULL,
  smoothWindow = NULL,
  flank = NULL,
  flankNorm = NULL,
  baseSize = 6,
  normMethod = NULL,
  logFile = NULL
  ){

  errorList <- list()

  #Get Footprint Info
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="footprint"),], name)
  biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="bias"),], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
  biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

  errorList$footMat <- footMat
  errorList$biasMat <- biasMat
  errorList$footDF <- footDF
  errorList$biasDF <- biasDF

  #Smooth Foot and Bias Mat because of sparsity
  if(!is.null(smoothWindow)){
    .logMessage("Applying smoothing window to footprint", logFile = logFile)
    footMat <- apply(footMat, 2, function(x) .centerRollMean(x, smoothWindow))
    biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, smoothWindow))
  }

  #Normalize Foot and Bias Mat
  .logMessage("Normalizing by flanking regions", logFile = logFile)
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat) / colMeans(footMat[idx, ,drop=FALSE]))
  biasMat <- t(t(biasMat) / colMeans(biasMat[idx, ,drop=FALSE]))

  errorList$footMatNorm <- footMat
  errorList$biasMatNorm <- footMat

  #Norm Foot By Bias
  if(tolower(normMethod) == "none"){
    title <- ""
  }else if(tolower(normMethod) == "subtract"){
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }else if(tolower(normMethod) == "divide"){
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat / biasMat
  }else{
    stop("normMethod not recognized!")
  }
  .logMessage(paste0("NormMethod = ", normMethod), logFile = logFile)

  #Get Mean and SD for each Assay
  footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x, 11)))

  errorList$footMatMean <- footMatMean
  errorList$footMatSd <- footMatSd
  errorList$biasMatMean <- biasMatMean
  errorList$biasMatSd <- biasMatSd
  errorList$smoothFoot <- smoothFoot

  #Create Plot Data Frames
  plotIdx <- seq_len(nrow(footMatMean)) #sort(unique(c(1, seq(1, nrow(footMatMean), smoothWindow), nrow(footMatMean))))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x){
    data.frame(
      x = footDF$x, 
      mean = footMatMean[,x], 
      sd = footMatSd[,x], 
      group = colnames(footMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))

  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x){
    data.frame(
      x = biasDF$x, 
      mean = biasMatMean[,x], 
      sd = biasMatSd[,x], 
      group = colnames(biasMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))

  errorList$plotFootDF <- plotFootDF
  errorList$plotBiasDF <- plotBiasDF

  out <- tryCatch({

    #Plot GG
    if(is.null(pal)){
      pal <- paletteDiscrete(values=gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
    }

    plotMax <- plotFootDF[order(plotFootDF$mean,decreasing=TRUE),]
    plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 50, ] #<= flank - flankNorm,]
    plotMax <- plotMax[!duplicated(plotMax$group),]
    plotMax <- plotMax[seq_len(ceiling(nrow(plotMax) / 4)), ]
    plotMax$x <- 25

    ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) + 
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() + 
      scale_color_manual(values = pal) + 
      scale_fill_manual(values = pal) + 
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE, 
        ylim = c(quantile(plotFootDF$mean, 0.0001), 1.15*quantile(smoothFoot, 0.999)), 
        xlim = c(min(plotFootDF$x),max(plotFootDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ggtitle(name) +
      guides(fill = FALSE) + 
      guides(color = FALSE) + ylab(paste0(title,"Normalized Insertions")) +
      ggrepel::geom_label_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA))

    ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) + 
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() + 
      scale_color_manual(values = pal) + 
      scale_fill_manual(values = pal) + 
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE, 
        ylim = c(quantile(plotBiasDF$mean, 0.0001), 1.05*quantile(plotBiasDF$mean, 0.999)), 
        xlim = c(min(plotBiasDF$x),max(plotBiasDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ylab("Tn5-Bias Normalized Insertions") + 
      theme(legend.position = "bottom", legend.box.background = element_rect(color = NA)) 
         
    ggAlignPlots(ggFoot, .ggSmallLegend(ggBias), sizes=c(2,1), draw = FALSE)

  }, error = function(e){

    .logError(e, fn = ".ggFootprint", info = name, errorList = errorList, logFile = logFile)

  })

  out

}

.ggSmallLegend <- function(
  gg = NULL,
  pointSize = 2,
  baseSize = 5,
  spaceLegend = 0.1
  ) {
    #https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
    gg +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = baseSize), 
              legend.text  = element_text(size = baseSize),
              legend.key.size = unit(spaceLegend, "lines"))
}

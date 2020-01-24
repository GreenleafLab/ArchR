##########################################################################################
# Transcription Factor Footprinting Methods
##########################################################################################

#' Plot footprints for an ArchRProject
#' 
#' This function will plot footprints for all samples in a given ArchRProject or a properly-formatted Summarized Experiment
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param positions A `list` or `GenomicRangesList` of `GRanges` containing the positions to incorporate into the footprint. Each position should be stranded.
#' @param plotName The prefix to add to the file name for the output PDF file containing the footprint plots.
#' @param groupBy The name of the column in `cellColData` used in `addGroupCoverages` for grouping multiple cells together (see `addGroupCoverages`).
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in `cellColData`. This limits the groups used to perform footprinting.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for plotting the lines corresponding to the footprints.
#' @param flank The number of basepairs from the position center (+/-) to consider as the flank.
#' @param flankNorm The number of basepairs to consider at the edge of the flank region (+/-) to be used for footprint normalization.
#' @param smoothWindow The size in basepairs of the sliding window to be used for smoothing of the footprint signal.
#' @param minCells The minimum number of cells required in a given cell group to permit footprint generation.
#' @param nTop The number of genomic regions to consider. Only the top `nTop` genomic regions based on the "score" column in the `GRanges` object will be considered for the footprint.
#' @param normMethod The name of the normalization method to use to normalize the footprint relative to the Tn5 insertion bias. Options include "none", "subtract", "divide". "Subtract" means subtracting the normalized Tn5 Bias. "Divide" means dividing the normalized Tn5 Bias.
#' @param inputSE Input a previous footprint Summarized Experiment (returned after running `plotFootprints`) to be plotted instead of regenerating the footprinting information.
#' @param height The height in inches to be used for the output PDF.
#' @param width The width in inches to be used for the output PDF file.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for preventing overwritting of old plots.
#' @param useSink A boolean value that indicates whether the `sink` function from base R should be used to hide messages during plotting.
#' @param threads The number of threads to be used for parallel computing.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... additional args
#' @export
plotFootprints <- function(
  ArchRProj = NULL,
  positions = NULL,
  plotName = "Plot-Footprints",
  groupBy = "Clusters",
  useGroups = NULL,
  pal = NULL,
  flank = 250,
  flankNorm = 50,
  smoothWindow = 10,
  minCells = 25,
  nTop = NULL,
  normMethod = "none",
  inputSE = NULL,
  height = 6,
  width = 4,
  addDOC = TRUE,
  useSink = TRUE,
  threads = getArchRThreads(),
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = positions, name = "positions", valid = c("grangeslist"))
  .validInput(input = plotName, name = "plotName", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  .validInput(input = flank, name = "flank", valid = c("integer"))
  .validInput(input = flankNorm, name = "flankNorm", valid = c("integer"))
  .validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = nTop, name = "nTop", valid = c("integer", "null"))
  .validInput(input = normMethod, name = "normMethod", valid = c("character"))
  .validInput(input = inputSE, name = "inputSE", valid = c("summarizedexperiment", "null"))
  .validInput(input = height, name = "height", valid = "integer")
  .validInput(input = width, name = "width", valid = "integer")
  .validInput(input = addDOC, name = "addDOC", valid = "boolean")
  .validInput(input = useSink, name = "useSink", valid = "boolean")
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
  .validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))

  tstart <- Sys.time()

  if(is.null(inputSE)){
    
    #Validate Positions
    if(!inherits(positions, "GenomicRangesList") & !inherits(positions, "list") & !inherits(positions, "SimpleList")){
      stop("Positions is not a list!")
    }
    positions <- as(positions, "list")
    valid <- lapply(positions, function(x) inherits(x, "GRanges")) %>% unlist %>% all
    if(!valid){
      stop("Positions is not a list of GenomicRanges!")
    }
    
    #If wanted can subset top positions
    if(!is.null(nTop)){
      posNames <- names(positions)
      positions <- lapply(seq_along(positions), function(x){
        positions[[x]][head(order(mcols(positions[[x]])$score, decreasing=TRUE), nTop)]
      })
      names(positions) <- posNames
    }

    #Get Footprints
    .messageDiffTime("Summarizing Footprints", tstart, addHeader = verboseAll)
    seFoot <- .summarizeFootprints(
      ArchRProj = ArchRProj, 
      positions = positions,
      groupBy = groupBy,
      useGroups = useGroups,
      minCells = minCells,
      flank = flank,
      threads = threads,
      verboseHeader = verboseHeader,
      verboseAll = verboseAll
    )

  }else{
    
    if(inherits(inputSE, "SummarizedExperiment")){
      seFoot <- inputSE
      rm(inputSE)
      gc()
      if(!is.null(useGroups)){
        if(sum(SummarizedExperiment::colData(seFoot)[,1] %in% useGroups) == 0){
          stop("No Groups found matching useGroups!")
        }
        seFoot <- seFoot[,SummarizedExperiment::colData(seFoot)[,1] %in% useGroups]
      }
    }else{
      stop("inputSE must be a footprint summarized experiment!")
    }

  }

  ############################################################################################
  # Plot Helper
  ############################################################################################

  .messageDiffTime("Plotting Footprints", tstart, addHeader = verboseAll)

  o <- tryCatch({
    if(useSink){
      tmpFile <- .tempfile()
      sink(tmpFile)
    }

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

    for(i in seq_along(seFoot@assays)){
      print(
        grid::grid.draw(.ggFootprint(
          seFoot = seFoot, 
          name = names(seFoot@assays)[i], 
          pal = pal, 
          smoothWindow = smoothWindow, 
          flank = flank, 
          flankNorm = flankNorm, 
          normMethod = normMethod
        )
      ))
      if(i != length(seFoot@assays)){
        grid::grid.newpage()
      }
    }
    dev.off()

    if(useSink){
      sink()
      file.remove(tmpFile)
    }

  }, error = function(x){

    suppressWarnings(sink())
    message(x)

  })

  seFoot

}

.ggFootprint <- function(seFoot, name, pal, smoothWindow, flank, flankNorm, baseSize = 6, normMethod, ...){

  #Get Footprint Info
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="footprint"),], name)
  biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="bias"),], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
  biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

  #Smooth Foot and Bias Mat because of sparsity
  footMat <- apply(footMat, 2, function(x) .centerRollMean(x, smoothWindow))
  biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, smoothWindow))

  #Normalize Foot and Bias Mat
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat) / colMeans(footMat[idx, ,drop=FALSE]))
  biasMat <- t(t(biasMat) / colMeans(biasMat[idx, ,drop=FALSE]))

  #Norm Foot By Bias
  if(tolower(normMethod) == "none"){
  }else if(tolower(normMethod) == "subtract"){
    footMat <- footMat - biasMat
  }else if(tolower(normMethod) == "divide"){
    footMat <- footMat / biasMat
  }else{
    stop("normMethod not recognized!")
  }

  #Get Mean and SD for each Assay
  footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)

  #Create Plot Data Frames
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x){
    data.frame(
      x = footDF$x, 
      mean = footMatMean[,x], 
      sd = footMatSd[,x], 
      group = colnames(footMatMean)[x]
      )
  }) %>% Reduce("rbind",. )
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))

  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x){
    data.frame(
      x = biasDF$x, 
      mean = biasMatMean[,x], 
      sd = biasMatSd[,x], 
      group = colnames(biasMatMean)[x]
      )
  }) %>% Reduce("rbind",. )
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))

  #Plot GG
  if(is.null(pal)){
    pal <- paletteDiscrete(values=SummarizedExperiment::colData(seFoot)$Group)
  }

  plotMax <- plotFootDF[order(plotFootDF$mean,decreasing=TRUE),]
  plotMax <- plotMax[abs(plotMax$x) <= flank - flankNorm,]
  plotMax <- plotMax[!duplicated(plotMax$group),]
  plotMax <- plotMax[seq_len(ceiling(nrow(plotMax) / 4)), ]
  plotMax$x <- 25

  ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) + 
    geom_line() + 
    scale_color_manual(values = pal) + 
    scale_fill_manual(values = pal) + 
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
    xlab("Distance to motif center (BP)") +
    coord_cartesian(
      expand = FALSE, 
      ylim = c(quantile(plotFootDF$mean, 0.0001), 1.15*quantile(plotFootDF$mean, 0.999)), 
      xlim = c(min(plotFootDF$x),max(plotFootDF$x))
    ) + theme_ArchR(baseSize = baseSize) + ggtitle(name) +
    guides(fill = FALSE) + 
    guides(color = FALSE) + ylab("Footprint Normalized Mean") +
    ggrepel::geom_label_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA))

  ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) + 
    geom_line() + 
    scale_color_manual(values = pal) + 
    scale_fill_manual(values = pal) + 
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
    xlab("Distance to motif center (BP)") +
    coord_cartesian(
      expand = FALSE, 
      ylim = c(quantile(plotBiasDF$mean, 0.0001), 1.05*quantile(plotBiasDF$mean, 0.999)), 
      xlim = c(min(plotBiasDF$x),max(plotBiasDF$x))
    ) + theme_ArchR(baseSize = baseSize) + ylab("Bias Normalized Mean") + 
    theme(legend.position = "bottom", legend.box.background = element_rect(color = NA)) 
       
  ggAlignPlots(ggFoot, .ggSmallLegend(ggBias), sizes=c(2,1), draw = FALSE)

}

.ggSmallLegend <- function(gg, pointSize = 2, baseSize = 5, spaceLegend = 0.1) {
    #https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
    gg +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = baseSize), 
              legend.text  = element_text(size = baseSize),
              legend.key.size = unit(spaceLegend, "lines"))
}

#####################################################################################################
# Summarize Footprints into a Summarized Experiment for Plotting 
#####################################################################################################

#' @export
.summarizeFootprints <- function(
  ArchRProj = NULL,
  positions,
  groupBy = "Clusters",
  useGroups = NULL,
  minCells = 25,
  flank = 250,
  threads = 16,
  force = FALSE,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  if(verboseAll){
    verboseHeader <- TRUE
  }

  tstart <- Sys.time()

  #####################################################
  # Compute Kmer Frequency Table 
  #####################################################
  coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy, minCells = minCells)
  coverageParams <- .getCoverageParams(ArchRProj = ArchRProj, groupBy = groupBy)
  kmerLength <- coverageParams$kmerLength

  if(!is.null(useGroups)){
    if(sum(coverageMetadata[,1] %in% useGroups) == 0){
      stop("No Groups found matching useGroups!")
    }
    coverageMetadata <- coverageMetadata[coverageMetadata[,1] %in% useGroups,]
  }

  genome <- getGenome(ArchRProj)
  .requirePackage(genome)
  .requirePackage("Biostrings")
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)

  .messageDiffTime("Computing Kmer Bias Table", tstart)
  kmerTableList <- .kmerPositionFrequency(positions, genome = BSgenome, flank = flank, k = kmerLength, threads = 1, verbose = verboseAll)

  #####################################################
  # Compute Footprints
  #####################################################
  .messageDiffTime("Computing Footprint", tstart)
  footprintList <- .computeFootprints(positions, coverageMetadata$File, flank = flank, threads = threads, verbose = verboseAll)

  #####################################################
  # Compute Bias For Footprints
  #####################################################
  .messageDiffTime("Computing Footprint Bias", tstart)
  footprintBiasList <- .computeFootprintsBias(kmerTableList, coverageMetadata$File, threads = threads, verbose = verboseAll)

  #####################################################
  # Summarize into SE
  #####################################################
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

.computeFootprintsBias <- function(kmerTableList, coverageFiles, threads = 8, verbose = TRUE){
  tstart <- Sys.time()
  out <- .safelapply(seq_along(coverageFiles), function(i){
      .messageDiffTime(sprintf("Computing Footprints Bias %s of %s:", i, length(coverageFiles)),tstart,verbose=verbose)
      .computeFootprintsBiasSingle(kmerTableList, coverageFiles[i])
  }, threads = threads) %>% SimpleList
  return(out)
}

.computeFootprintsBiasSingle <- function(kmerTableList, coverageFile){
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

.computeFootprints <- function(featureList, coverageFiles, flank = 250, threads = 8, verbose = TRUE){
  tstart <- Sys.time()
  out <- .safelapply(seq_along(coverageFiles), function(i){
    .computeFootprintsSingle(featureList, coverageFiles[i], flank, gc = TRUE, 
      pre = sprintf("Computing Footprints %s of %s:", i, length(coverageFiles)),
      tstart = tstart, 
      verbose = verbose
      )
  }, threads = threads) %>% SimpleList
  return(out)
}

.computeFootprintsSingle <- function(featureList, coverageFile, flank = 250, gc = FALSE, pre = "", tstart, verbose = TRUE){
  window <- 2 * flank + 1
  featureNames <- names(featureList)
  featureList <- as(featureList, "list")
  allChr <- lapply(featureList, function(x) unique(as.character(seqnames(x)))) %>% unlist %>% unique %>% sort
  cov <- .getCoverageRle(coverageFile, allChr)
  footprintDF <- lapply(seq_along(featureList), function(x){
    featurex <- split(resize(featureList[[x]],1,"center"), seqnames(featureList[[x]]))
    outx <- ArchR:::rleSumsStranded(cov, featurex, window, as.integer) #Rcpp
    if(x %% 25 == 0 & gc){
      gc()
    }
    if(length(featureList) > 10){
      if(x %% max(floor(length(featureList) * .25) ,1) == 0){
        .messageDiffTime(sprintf("%s %s Percent Completed", pre, floor(x / floor(length(featureList)/25) * 25)), tstart, verbose=verbose) 
      }
    }else{
        if(x == 1 | x == length(featureList)){
          .messageDiffTime(sprintf("%s %s Percent Completed", pre, round(100 * x / length(featureList)),1), tstart, verbose=verbose)  
        }
    }
    outx
  }) %>% Reduce("cbind",.) %>% data.frame
  gc()
  footprintDF   
}

.getCoverageRle <- function(coverageFile, allChr){
  cov <- lapply(seq_along(allChr), function(x){
    Rle(
      lengths = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Lengths")), 
      values = h5read(coverageFile, paste0("Coverage/",allChr[x],"/Values"))
    )
  }) %>% {as(.,"RleList")}
  names(cov) <- allChr
  cov
}

.kmerPositionFrequency <- function(featureList, genome, flank = 250, k = 6, threads = 8, verbose = TRUE){
  
  tstart <- Sys.time()
  genome <- validBSgenome(genome)
  window <- 2*flank + 1

  kmerList <- .safelapply(seq_along(featureList), function(i){
    .messageDiffTime(sprintf("Computing Kmer Tables for %s of %s features", i, length(featureList)), tstart, verbose=verbose)
    bsv <- BSgenomeViews(genome , resize(featureList[[i]], window + k, "center"))
    bsv <- bsv[width(bsv) == window + k] #none that are trimmed!
    #BSgenome is already stranded
    #kmerPositionFrequencyCpp is Rcpp export for getting kmer position frequencies from strings
    kmerTable <-  ArchR:::kmerPositionFrequencyCpp(as.character(bsv), rep(1L,length(bsv)), window, k, .getKmers(k)) #Rcpp
    return(kmerTable)
  }, threads = threads) %>% SimpleList
  names(kmerList) <- names(featureList)
  
  .messageDiffTime("Finished Computing Kmer Tables", tstart)

  return(kmerList)
}

.getKmers <-function(k, letters = c('A','C','G','T')){
  kmers = ''
  for (i in seq_len(k)) {
    kmers <- unlist(lapply(kmers, function(x) paste0(x, letters)))
  }  
  return(kmers)
}


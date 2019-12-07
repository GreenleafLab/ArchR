#' Create Arrow Files 
#' 
#' This function will create an Arrow Files from input files
#' for downstream analysis
#'
#' @param inputFiles input files (tabixFile, bamFile or textFile)
#' @param sampleNames sample names corresponding to input files
#' @param outputNames output names prefix (ie PBMC -> PBMC.arrow)
#' @param geneAnno geneAnnotation input for TSS Scores etc.
#' @param genomeAnno genomeAnnotation input for ChromSizes Nucleotide Information etc.
#' @param filterFrags min fragments per cell to be filtered for analyses such as tileMat etc.
#' @param filterTSS min TSS Score per cell to be filtered for analyses such as tileMat etc.
#' @param removeFilteredCells remove fragments corresponding to cells pass filterFrags and filterTSS
#' @param minFrags min fragments per cell to be immediately filtered
#' @param outDir out directory for QC information from sample to be plotted / saved
#' @param nucLength nucleosome length for id'ing fragments as sub-, mono-, or multi-nucleosome spanning
#' @param TSSParams TSS parameters for computing TSS scores
#' @param excludeChr exclude these chromosomes from analysis downstream (does not apply to fragments)
#' @param nChunk number of chunks per chromosome when reading in input files
#' @param bcTag barcode tag location in bam file (see ScanBam in Rsamtools)
#' @param bamFlag list of bam flags for reading in fragments from input files (see ScanBam in Rsamtools)
#' @param offsetPlus Tn5 offset of "+" stranded insertion (see Buenrostro 2013)
#' @param offsetMinus Tn5 offset of "-" stranded insertion (see Buenrostro 2013)
#' @param addTileMat addTileMatrix to ArrowFiles
#' @param TileMatParams additional parameters to pass to addTileMatrix (see addTileMatrix)
#' @param addGeneScoreMat addGeneScoreMatrix to ArrowFiles
#' @param GeneScoreMatParams additional parameters to pass to addGeneScoreMatrix (see addGeneScoreMatrix)
#' @param force force creation of arrow files if already exist
#' @param threads number threads for parallel execution
#' @param parallelParam parallel parameters for batch style execution
#' @param ... additional args
#' @export
createArrowFiles <- function(
  inputFiles = NULL, 
  sampleNames = NULL, 
  outputNames = paste0("./", sampleNames),
  validBaroces = NULL,
  geneAnno = NULL,
  genomeAnno = NULL,
  filterFrags = 1000,
  filterTSS = 4,
  removeFilteredCells = TRUE,
  minFrags = 500, 
  outDir = "QualityControl",
  nucLength = 147,
  TSSParams = list(),
  excludeChr = c("chrM", "chrY"),
  nChunk = 5,
  bcTag = "qname",
  gsubExpression = NULL,
  bamFlag = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  addTileMat = TRUE,
  TileMatParams = list(),
  addGeneScoreMat = TRUE,
  GeneScoreMatParams = list(),
  force = FALSE,
  threads = 1,
  parallelParam = NULL,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  ...
  ){

  dir.create(outDir, showWarnings = FALSE)

  #Add args to list
  args <- list()
  args <- append(args, mget(names(formals()),sys.frame(sys.nframe()))) #as.list(match.call())
  args$X <- seq_along(inputFiles)
  args$FUN <- .createArrow
  args$threads <- min(args$threads, length(inputFiles))
  args$registryDir <- file.path(outDir, "CreateArrowsRegistry")

  #Run With Parallel or lapply
  outArrows <- tryCatch({
    unlist(.batchlapply(args))
  },error = function(x){
    message("createArrowFiles has encountered an error, checking if any ArrowFiles completed...")
    for(i in seq_along(args$outputNames)){
      out <- paste0(args$outputNames[i],".arrow")
      if(file.exists(out)){
        o <- tryCatch({
          o <- h5read(out, "Metadata/Completed") #Check if completed
        },error = function(y){
          file.remove(out) #If not completed delete
        })
      }
    }
    print(paste0("Error Received : ", x))
    paste0(args$outputNames,".arrow")[file.exists(paste0(args$outputNames,".arrow"))]
  })

  return(outArrows)

}

#Main Function!
.createArrow <- function(
  i,
  inputFiles = NULL, 
  sampleNames = NULL, 
  outputNames = paste0("./", sampleName),
  nChunk = 3,
  offsetPlus = 4,
  offsetMinus = -5,
  geneAnno = NULL,
  genomeAnno = NULL,
  minFrags = 500, 
  removeFilteredCells = TRUE,
  filterFrags = 1000,
  filterTSS = 4,
  excludeChr = c("chrM", "chrY"),
  gsubExpression = NULL,
  bcTag = "qname",
  bamFlag = NULL,
  outDir = "QualityControl",
  nucLength = 147,
  TSSParams = list(),
  addTileMat = TRUE,
  TileMatParams = list(),
  addGeneScoreMat = TRUE,
  GeneScoreMatParams = list(),
  force = FALSE,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  tstart = NULL,
  ...
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  ArrowFiles <- paste0(outputNames, ".arrow")
  
  inputFile <- inputFiles[i]
  sampleName <- sampleNames[i] 
  outputName <- outputNames[i]
  ArrowFile <- ArrowFiles[i]
  outDir <- file.path(outDir, sampleName)
  dir.create(outDir, showWarnings = FALSE)

  prefix <- sprintf("(%s : %s of %s)", sampleName, i, length(inputFiles))

  .requirePackage("rhdf5")
  .requirePackage("Rsamtools")


  #Check if a completed file exists!
  if(file.exists(ArrowFile)){
    o <- tryCatch({
      o <- h5read(ArrowFile, "Metadata/Completed") #Check if completed
      if(force){
        .messageDiffTime(sprintf("%s Arrow Exists! Overriding since force = TRUE!", prefix), tstart, verbose = TRUE, addHeader = FALSE)
        file.remove(ArrowFile)
      }else{
        .messageDiffTime(sprintf("%s Arrow Exists! Marking as completed since force = FALSE!", prefix), tstart, verbose = TRUE, addHeader = FALSE)
        return(ArrowFile)
      }
    },error = function(y){
      .messageDiffTime(sprintf("%s Arrow Exists! Overriding since not completed!", prefix), tstart, verbose = TRUE, addHeader = FALSE)
      file.remove(ArrowFile) #If not completed delete
    })
  }

  #############################################################
  #Determine Arrow Construction Method
  #############################################################  
  fe <- .fileExtension(inputFile)
  if(fe == "gz" | fe == "bgz" | fe == "tsv" | fe == "txt"){
    if(.isTabix(inputFile)){
      readMethod <- "tabix"
    }else{
      stop("Error only supported inputs are tabix and bam!")
      readMethod <- "tsv"
    }
  }else if(fe == "bam"){
    readMethod <- "bam"
  }else{
    stop(sprintf("Read Method for %s Not Recognized!", fe))
  }

  .messageDiffTime(sprintf("%s Reading In Fragments from inputFiles (readMethod = %s)", prefix, readMethod), tstart, verbose = verboseHeader, addHeader = verboseAll)

  if(tolower(readMethod) == "tsv"){

    out <- .tsvToArrow(tsvFile = inputFile, outArrow = ArrowFile, 
            chromSizes = genomeAnno$chromSizes, genome = genomeAnno$genome, 
            minFrags = minFrags, sampleName = sampleName, prefix = prefix, 
            verboseHeader = verboseHeader, verboseAll = verboseAll, tstart = tstart, ...)

  }else if(tolower(readMethod)=="tabix"){
   
    tmp <- .tabixToTmp(tabixFile = inputFile, sampleName = sampleName, chromSizes = genomeAnno$chromSizes, nChunk = nChunk,
            gsubExpression = gsubExpression, prefix = prefix, verboseHeader = verboseHeader, 
            verboseAll = verboseAll, tstart = tstart,  ...)

    out <- .tmpToArrow(tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnno$genome, 
            minFrags = minFrags, sampleName = sampleName, prefix = prefix, 
            verboseHeader = verboseHeader, verboseAll = verboseAll, tstart = tstart, 
            chromSizes = genomeAnno$chromSizes, ...)
  
  }else if(tolower(readMethod)=="bam"){

    tmp <- .bamToTmp(bamFile = inputFile, sampleName = sampleName, 
            chromSizes = genomeAnno$chromSizes, bamFlag = bamFlag, 
            bcTag = bcTag, gsubExpression = gsubExpression, nChunk = nChunk, 
            offsetPlus = offsetPlus, offsetMinus = offsetMinus, prefix = prefix, 
            verboseHeader = verboseHeader, verboseAll = verboseAll, tstart = tstart, ...)

    out <- .tmpToArrow(tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnno$genome, 
            minFrags = minFrags, sampleName = sampleName, prefix = prefix, 
            verboseHeader = verboseHeader, verboseAll = verboseAll, tstart = tstart, 
            chromSizes = genomeAnno$chromSizes, ...)

  }else{
    
    stop(sprintf("Read Method : %s Not Recognized!", readMethod))

  }
  gc()
  
  #############################################################
  #Compute Fragment Information!
  #############################################################
  .messageDiffTime(sprintf("%s Adding Fragment Summary", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)

  fragSummary <- .fastFragmentInfo(ArrowFile = ArrowFile, cellNames = .availableCells(ArrowFile), nucLength = nucLength)
  Metadata <- fragSummary[[1]]
  plot <- tryCatch({
    
    tmpFile <- .tempfile()
    sink(tmpFile)

    dir.create(outDir, showWarnings = FALSE)
    pdf(file.path(outDir,paste0(sampleName,"-Fragment_Size_Distribution.pdf")),width=4,height=3,onefile=FALSE)
    plotDF <- data.frame(
      x = seq_along(fragSummary[[2]]), 
      percent = 100 * fragSummary[[2]]/sum(fragSummary[[2]])
    )
    gg <- ggplot(plotDF, aes(x = x, y = percent)) + theme_ArchR() + 
          geom_line(col = "darkblue", size = 0.25) + 
          coord_cartesian(xlim = c(1,750), ylim = c(0, max(plotDF$percent) * 1.1), expand = FALSE) + 
          xlab("Size of Fragments (bp) \n") + 
          ylab("Fragments (%)") + 
          ggtitle("Fragment Size Distribution")
    print(.fixPlotSize(gg, plotWidth = 4.5, plotHeight = 3.5, height = 3/4))
    dev.off()

    sink()
    file.remove(tmpFile)

  }, error = function(x){

      .messageDiffTime("Continuing through after error ggplot for Fragment Size Distribution", tstart)
      print(x)
      message("\n")

  })
  gc()

  #############################################################
  #Compute TSS Enrichment Scores Information!
  #############################################################
  .messageDiffTime(sprintf("%s Computing TSS Enrichment Scores", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)
  TSSParams$TSS <- geneAnno$TSS
  TSSParams$ArrowFile <- ArrowFile
  TSSParams$cellNames <- Metadata$cellNames
  TSSOut <- do.call(.fastTSSEnrichment, TSSParams)
  Metadata$TSSEnrichment <- TSSOut$tssScores
  Metadata$ReadsInTSS <- TSSOut$tssReads

  #Filter
  Metadata$Keep <- 1*(Metadata$nFrags >= filterFrags & Metadata$TSSEnrichment >= filterTSS)
  message(paste0(sampleName, " : Number of Cells Pass Filter = ", sum(Metadata$Keep)))
  message(paste0(sampleName, " : Median Frags = ", median(Metadata$nFrags[Metadata$Keep==1])))
  message(paste0(sampleName, " : Median TSS Enrichment = ", median(Metadata$TSSEnrichment[Metadata$Keep==1])))
  
  plot <- tryCatch({
    
    tmpFile <- .tempfile()
    sink(tmpFile)

    ggtitle <- sprintf("%s\n%s\n%s",
        paste0(sampleName, " : Number of Cells Pass Filter = ", sum(Metadata$Keep)),
        paste0("Median Frags = ", median(Metadata$nFrags[Metadata$Keep==1])),
        paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment[Metadata$Keep==1]))
      )

    pdf(file.path(outDir,paste0(sampleName,"-TSS_by_Unique_Frags.pdf")),width=4,height=4,onefile=FALSE)
    gg <- ggPoint(
      x = log10(Metadata$nFrags),
      y = Metadata$TSSEnrichment, 
      colorDensity = TRUE,
      continuousSet = "samba_night",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE) + 
      geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.5) +
      geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.5)
    print(.fixPlotSize(gg, plotWidth = 4, plotHeight = 4))
    dev.off()

    sink()
    file.remove(tmpFile)

  }, error = function(x) {

      .messageDiffTime("Continuing through after error ggplot for TSS by Frags", tstart)
      print(x)
      message("\n")

  })

  #Add To Metadata
  #Sanity Check Here to Make Sure!
  stopifnot(
    identical(
      paste0(h5read(ArrowFile, "Metadata/CellNames")), 
      paste0(stringr::str_split(Metadata$cellNames, pattern = "#", simplify = TRUE)[,2])
    )
  )
  o <- h5write(obj = Metadata$Keep, file = ArrowFile, name = "Metadata/PassQC")
  o <- h5write(obj = Metadata$nFrags, file = ArrowFile, name = "Metadata/nFrags")
  o <- h5write(obj = Metadata$nMonoFrags, file = ArrowFile, name = "Metadata/nMonoFrags")
  o <- h5write(obj = Metadata$nDiFrags, file = ArrowFile, name = "Metadata/nDiFrags")
  o <- h5write(obj = Metadata$nMultiFrags, file = ArrowFile, name = "Metadata/nMultiFrags")
  o <- h5write(obj = (Metadata$nDiFrags + Metadata$nMultiFrags) / Metadata$nMonoFrags, file = ArrowFile, name = "Metadata/NucleosomeRatio")
  o <- h5write(obj = Metadata$TSSEnrichment, file = ArrowFile, name = "Metadata/TSSEnrichment")
  o <- h5write(obj = Metadata$ReadsInTSS, file = ArrowFile, name = "Metadata/ReadsInTSS")
  gc()

  #############################################################
  # Remove Cells That Are Filtered?
  #############################################################
  
  if(removeFilteredCells){
    .messageDiffTime(sprintf("%s Removing Fragments from Filtered Cells", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)
    idx <- which(Metadata$Keep == 1)
    o <- .filterCellsFromArrow(inArrow = ArrowFile, cellNames = Metadata$cellNames[idx])
    o <- h5write(obj = Metadata$Keep[idx], file = ArrowFile, name = "Metadata/PassQC")
    o <- h5write(obj = Metadata$nFrags[idx], file = ArrowFile, name = "Metadata/nFrags")
    o <- h5write(obj = Metadata$nMonoFrags[idx], file = ArrowFile, name = "Metadata/nMonoFrags")
    o <- h5write(obj = Metadata$nDiFrags[idx], file = ArrowFile, name = "Metadata/nDiFrags")
    o <- h5write(obj = Metadata$nMultiFrags[idx], file = ArrowFile, name = "Metadata/nMultiFrags")
    o <- h5write(obj = ((Metadata$nDiFrags + Metadata$nMultiFrags) / Metadata$nMonoFrags)[idx], file = ArrowFile, name = "Metadata/NucleosomeRatio")
    o <- h5write(obj = Metadata$TSSEnrichment[idx], file = ArrowFile, name = "Metadata/TSSEnrichment")
    o <- h5write(obj = Metadata$ReadsInTSS[idx], file = ArrowFile, name = "Metadata/ReadsInTSS")
  }


  #############################################################
  # Create Tile Matrix
  #############################################################
  if(addTileMat){
    .messageDiffTime(sprintf("%s Adding TileMatrix", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)
    TileMatParams$i <- i
    TileMatParams$ArrowFiles <- ArrowFiles
    TileMatParams$cellNames <- Metadata$cellNames[idx]
    chromLengths <- end(genomeAnno$chromSizes)
    names(chromLengths) <- paste0(seqnames(genomeAnno$chromSizes))
    TileMatParams$chromLengths <- chromLengths
    TileMatParams$blacklist <- genomeAnno$blacklist
    TileMatParams$force <- TRUE
    TileMatParams$excludeChr <- excludeChr
    tileMat <- suppressMessages(do.call(.addTileMat, TileMatParams))
    gc()
  }

  #############################################################
  # Add Gene Score Matrix
  #############################################################
  if(addGeneScoreMat){
    .messageDiffTime(sprintf("%s Adding GeneScoreMatrix", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)
    GeneScoreMatParams$i <- i
    GeneScoreMatParams$ArrowFiles <- ArrowFiles
    GeneScoreMatParams$genes <- geneAnno$genes
    GeneScoreMatParams$cellNames <- Metadata$cellNames[which(Metadata$Keep==1)]
    GeneScoreMatParams$blacklist <- genomeAnno$blacklist
    GeneScoreMatParams$force <- TRUE
    GeneScoreMatParams$excludeChr <- excludeChr
    geneScoreMat <- suppressMessages(do.call(.addGeneScoreMat, GeneScoreMatParams))
    gc()
  }

  o <- h5closeAll()

  .messageDiffTime(sprintf("%s Finished Creating ArrowFile", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)
  o <- h5write(obj = "Finished", file = ArrowFile, name = "Metadata/Completed")

  ArrowFile <- paste0(outputName, ".arrow")  

  return(ArrowFile)

}

#########################################################################################################
# QC Methods 
#########################################################################################################

.fastFragmentInfo <- function(
  ArrowFile,
  cellNames = .availableCells(ArrowFile),
  nucLength = 147, 
  ...){

  #Info to get
  matNuc <- matrix(0, nrow = length(cellNames), ncol = 3)
  nFrags <- rep(0, length(cellNames))
  fragDist <- rep(0, 1000)

  chrArrow <- .availableChr(ArrowFile)
  for(x in seq_along(chrArrow)){

    #Read in frags
    fragx <- .getFragsFromArrow(ArrowFile = ArrowFile, chr = chrArrow[x], out = "IRanges", cellNames = cellNames)

    if(length(fragx) > 0){

      mcols(fragx)$RG@values <- S4Vectors::match(mcols(fragx)$RG@values, cellNames)
      nFrags <- nFrags + S4Vectors:::tabulate(mcols(fragx)$RG, nbins = length(cellNames))

      #Get Distributions
      fragDist <- fragDist + tabulate(width(fragx), nbins = 1000)
      w <- trunc(width(fragx)/nucLength) + 1
      w[w > 3] <- 3

      #Get Nuc Info
      matNuc <- matNuc + ArchR:::tabulate2dCpp(
        w, xmin = 1, xmax = 3, 
        as.integer(mcols(fragx)$RG), ymin = 1, ymax = length(cellNames)
      )   

    }

  }

  df <- DataFrame(matNuc)
  colnames(df) <- c("nMonoFrags", "nDiFrags", "nMultiFrags")
  df$cellNames <- cellNames
  df$nFrags <- nFrags
  df <- df[,c("cellNames","nFrags","nMonoFrags", "nDiFrags", "nMultiFrags")]

  out <- list(dfSummary = df, fragDistribution  = fragDist)
  return(out)

}

.fastTSSEnrichment <- function(
  TSS, 
  ArrowFile, 
  cellNames = NULL, 
  window = 101, 
  norm = 100, 
  flank = 2000, 
  minNorm = 1, 
  ...){

  tstart <- Sys.time()

  #Validate
  ArrowFile <- .validArrow(ArrowFile)
  TSS <- .validGRanges(TSS)

  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  #Create Window and Flank
  TSS <- resize(TSS, 1, fix = "start")
  strand(TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- resize(TSS, window, "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    #Positive Flank
    GRanges(seqnames(TSS), IRanges(end(TSS) + flank - norm + 1, end(TSS) + flank)),
    #Negative Flank
    GRanges(seqnames(TSS), IRanges(start(TSS) - flank, start(TSS) - flank + norm - 1))
  )
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)

  #Count
  #.messageDiffTime("Counting Around TSS!", tstart)
  
  countList <- .fastFeatureCounts(feature = tssFeatures, ArrowFile = ArrowFile, cellNames = cellNames)

  #Normalize per BP
  cWn <- countList$nWindow / window
  cFn <- countList$nFlank / norm

  #Compute scores
  tssScores <- 2 * cWn / (pmax(cFn[names(cWn)], minNorm)) #Multiply 2 because enrichment over average 2 flanks
  names(tssScores) <- cellNames
  tssScores <- round(tssScores, 3)

  #.messageDiffTime("Computed TSS Scores!", tstart)

  return(list(tssScores=tssScores, tssReads=cWn))

}

.fastFeatureCounts <- function(feature, ArrowFile, cellNames){
  
  tstart1 <- Sys.time()
  featureList <- split(feature, seqnames(feature))
  chrArrow <- .availableChr(ArrowFile)
  featureList <- featureList[chrArrow]
  if(length(featureList)==0){
    stop("Error No Overlap in Chromosomes and TSS Chromosomes!")
  }

  #Count Vector
  nWindow <- rep(0, length(cellNames))
  names(nWindow) <- cellNames

  nFlank <- rep(0, length(cellNames))
  names(nFlank) <- cellNames

  #Count
  #message(sprintf("Counting Insertions, %s minutes elapsed...", round(difftime(Sys.time(), tstart1, units = "mins"),3)))
  #pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  
  for(x in seq_along(featureList)){

    #setTxtProgressBar(pb,round(x*100/length(featureList),0))

    ###############################################################################
    # Get Fragments
    ###############################################################################
    featurex <- featureList[[x]]
    fragments <- .getFragsFromArrow(ArrowFile = ArrowFile, chr = names(featureList)[x], out = "IRanges", cellNames = cellNames)

    if(length(fragments) > 0){

      mcols(fragments)$RG@values <- match(mcols(fragments)$RG@values, cellNames)
      mcols(featurex)$typeIdx <- match(mcols(featurex)$type, c("window", "flank"))
      
      ###############################################################################
      # Count Each Insertion
      ###############################################################################
      for(y in seq_len(2)){
          
        if(y==1){
          temp <- IRanges(start(fragments), width=1)
        }else if(y==2){
          temp <- IRanges(end(fragments), width=1)
        }
        stopifnot(length(temp) == length(fragments))

        o <- findOverlaps(ranges(featurex), temp)
        remove(temp)
        gc()
        
        mat <- ArchR:::tabulate2dCpp(
              x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
              xmin = 1,
              xmax = length(cellNames),
              y = mcols(featurex)$typeIdx[queryHits(o)],
              ymin = 1,
              ymax = 2
          )

        #Add To
        nWindow <- nWindow + mat[1, ]
        nFlank <- nFlank + mat[2, ]
        rm(o, mat)

      }
      
      rm(fragments)
      gc()

    }

  }

  #message("\n")
  #.messageDiffTime("Finished Counting Insertions", tstart1)

  out <- list(nWindow = nWindow, nFlank = nFlank)

  return(out)

}

#########################################################################################################
# Methods to Turn Input File into a Temp File that can then be Efficiently converted to an Arrow!
#########################################################################################################
.isTabix <- function(file){
  tryCatch({
    TabixFile(file)
    TRUE
  }, error = function(x){
    tryCatch({
      message("Attempting to index ", file," as tabix...")
      indexTabix(file, format = "bed")
      TRUE
    }, error = function(y){
      FALSE
    })
  })
}

.tabixToTmp <- function(
  tabixFile, 
  sampleName,
  tmpFile = .tempfile(pattern = paste0("tmp-",sampleName,"-arrow"), fileext=".arrow"),
  chromSizes, 
  nChunk = 3,
  gsubExpression = NULL, 
  printEvery = 1,
  verboseHeader = TRUE,
  verboseAll = FALSE,
  prefix = "",
  tstart = NULL,
  ...
  ){

  .requirePackage("Rsamtools")

  #######################################################################################################
  # We will dump a chunked genome into an Hdf5 file in a memory efficient manner!
  #######################################################################################################

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  tstart2 <- Sys.time()
  
  if(verboseAll){
    printEvery <- 0.25
  }else{
    printEvery <- 1
  }

  nextPrint <- printEvery
  o <- h5closeAll()
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("Fragments"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = "tmp", file = tmpFile, name = "Metadata/Sample")

  tileChromSizes <- unlist(GenomicRanges::tile(chromSizes, nChunk))
  mcols(tileChromSizes)$chunkName <- paste0(seqnames(tileChromSizes),"#chunk",seq_along(tileChromSizes))
  for(x in seq_along(tileChromSizes)){

    if(as.numeric(difftime(Sys.time(),tstart2,units="mins")) > nextPrint){
      .messageDiffTime(sprintf("%s Reading TabixFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), tstart, 
        verbose = verboseHeader, addHeader = verboseAll)
      nextPrint <- nextPrint + printEvery
    }

    dt <- tryCatch({
      Rsamtools::scanTabix(tabixFile, param = tileChromSizes[x])[[1]] %>%
        textConnection %>% 
        {tryCatch(read.table(.), error = function(e) NULL)} %>% 
        {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}
    }, error = function(f){
      NULL
    })

    #Care for Break Points
    dt <- dt[dt$V2 >= start(tileChromSizes[x]),]

    if(all(!is.null(dt), nrow(dt) > 0)){
      if(!is.null(gsubExpression)){
        scanChunk$V4 <- gsub(gsubExpression, "", scanChunk$V4)
      }

      #Order by bc
      setkey(dt, V4)
      dt <- dt[order(V4)]
      RG <- Rle(paste0(dt$V4))

      chrTmp <- mcols(tileChromSizes)$chunkName[x]
      chrPos <- paste0("Fragments/",chrTmp,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chrTmp,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chrTmp,"/RGValues")
      lengthRG <- length(RG@lengths)
      o <- h5createGroup(tmpFile, paste0("Fragments/",chrTmp))
      o <- .suppressAll(h5createDataset(tmpFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
      o <- .suppressAll(h5createDataset(tmpFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
      o <- .suppressAll(h5createDataset(tmpFile, chrRGValues, storage.mode = "character", 
        dims = c(lengthRG, 1), level = 0, size = nchar(RG@values[1]) + 1))
      o <- h5write(obj = cbind(dt$V2,dt$V3-dt$V2), file = tmpFile, name = chrPos)
      o <- h5write(obj = RG@lengths, file = tmpFile, name = chrRGLengths)
      o <- h5write(obj = RG@values, file = tmpFile, name = chrRGValues)

      rm(dt, RG)
      gc()
    }
  }

  return(tmpFile)

}

.bamToTmp <- function(
  bamFile, 
  sampleName,
  tmpFile = .tempfile(pattern = paste0("tmp-",sampleName,"-arrow"), fileext=".arrow"), 
  chromSizes, 
  bamFlag = NULL,
  nChunk = 3,
  bcTag = "qname",
  gsubExpression = NULL, 
  offsetPlus = 4, 
  offsetMinus = -5, 
  verboseHeader = TRUE,
  verboseAll = FALSE,
  prefix = "",
  tstart = NULL,
  ...){

  .requirePackage("Rsamtools")

  #######################################################################################################
  # We will dump a chunked genome into an Hdf5 file in a memory efficient manner!
  #######################################################################################################
  if(is.null(bamFlag)){
    bamFlag <- scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE)
  }else if(inherits(bamFlag, "list") | inherits(bamFlag, "SimpleList")){
    bamFlag$isMinusStrand <- FALSE
    bamFlag$isProperPair <- TRUE
    bamFlag <- do.call(scanBamFlag, bamFlag)
  }else{
    stop("bamFlag must be a list or null!")
  }

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  tstart2 <- Sys.time()
 
  if(verboseAll){
    printEvery <- 0.25
  }else{
    printEvery <- 1
  }

  nextPrint <- printEvery

  o <- h5closeAll()
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("Fragments"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = "tmp", file = tmpFile, name = "Metadata/Sample")

  tileChromSizes <- unlist(tile(chromSizes, nChunk))
  mcols(tileChromSizes)$chunkName <- paste0(seqnames(tileChromSizes),"#chunk",seq_along(tileChromSizes))
  for(x in seq_along(tileChromSizes)){

    if(as.numeric(difftime(Sys.time(),tstart2,units="mins")) > nextPrint){
      .messageDiffTime(sprintf("%s Reading BamFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), tstart, 
        verbose = verboseHeader, addHeader = verboseAll)
      nextPrint <- nextPrint + printEvery
    }

    #If barcode is stored in read name use qname
    #Else look for barcode tag such as RG
    if(tolower(bcTag)=="qname"){
      
      dt <- tryCatch({

        scanChunk <- scanBam(bamFile,
                param = ScanBamParam(
                  flag = bamFlag,
                  what = c("qname", "pos", "isize"),
                  which = tileChromSizes[x]
               ))[[1]]

      if(!is.null(gsubExpression)){
        scanChunk$qname <- gsub(gsubExpression, "", scanChunk$qname)
      }

      #Create Data Table for faster indexing
      data.table(
          start = scanChunk$pos + offsetPlus,
          end = scanChunk$pos + abs(scanChunk$isize) - 1 + offsetMinus,
          RG = scanChunk$qname
        )

      }, error = function(e){

        NULL

      })


    }else{
      
      dt <- tryCatch({

        scanChunk <- scanBam(bamFile,
                param = ScanBamParam(
                  flag = bamFlag,
                  what = c("pos", "isize"),
                  tag = bcTag,
                  which = tileChromSizes[x]
               ))[[1]]

        if(!is.null(gsubExpression)){
          scanChunk$tag[[bcTag]] <- gsub(gsubExpression, "", scanChunk$tag[[bcTag]])
        }

        #Create Data Table for faster indexing
        dt <- data.table(
            start = scanChunk$pos + offsetPlus,
            end = scanChunk$pos + abs(scanChunk$isize) - 1 + offsetMinus,
            RG = scanChunk$tag[[bcTag]]
          )

      }, error = function(e){
        
        NULL

      })

    }
    
    #Clean Up Memory
    rm(scanChunk)

    #Care for Break Points
    dt <- dt[dt$start >= start(tileChromSizes[x]),]    

    if(all(!is.null(dt), nrow(dt) > 0)){

      #Order by bc
      setkey(dt, RG)
      dt <- dt[order(RG)]
      RG <- Rle(dt$RG)

      chrTmp <- mcols(tileChromSizes)$chunkName[x]
      chrPos <- paste0("Fragments/",chrTmp,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chrTmp,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chrTmp,"/RGValues")
      lengthRG <- length(RG@lengths)
      o <- h5createGroup(tmpFile, paste0("Fragments/",chrTmp))
      o <- .suppressAll(h5createDataset(tmpFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
      o <- .suppressAll(h5createDataset(tmpFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
      o <- .suppressAll(h5createDataset(tmpFile, chrRGValues, storage.mode = "character", 
                dims = c(lengthRG, 1), level = 0, size = nchar(RG@values[1]) + 1))
      o <- h5write(obj = cbind(dt$start,dt$end-dt$start), file = tmpFile, name = chrPos)
      o <- h5write(obj = RG@lengths, file = tmpFile, name = chrRGLengths)
      o <- h5write(obj = RG@values, file = tmpFile, name = chrRGValues)

      rm(dt, RG)
      gc()

    }

  }

  return(tmpFile)

}


#########################################################################################################
# Methods to temp file to arrow!
#########################################################################################################

.tmpToArrow <- function(
  tmpFile, 
  outArrow, 
  genome, 
  chromSizes,
  minFrags = 500, 
  sampleName, 
  verboseHeader = TRUE,
  verboseAll = FALSE,
  tstart = NULL,
  prefix = "",
  ...
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  .messageDiffTime(sprintf("%s Creating ArrowFile", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)

  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = sampleName, file = outArrow, name = "Metadata/Sample")
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5createGroup(outArrow, paste0("Fragments"))

  #Get Info
  chunkNames <- .availableChr(tmpFile)

  #######################################################################################################
  # First we will count the number of occurences per barcode!
  #######################################################################################################
  .messageDiffTime(sprintf("%s Counting Unique Barcodes", prefix), tstart, verbose = verboseAll)
  o <- h5closeAll()
  h5DF <- h5ls(tmpFile, recursive = TRUE)
  dtList <- lapply(seq_along(chunkNames), function(x){
    chrTmp <- chunkNames[x]
    nRG <- h5DF %>% 
      {.[.$group==paste0("/Fragments/",chrTmp) & .$name == "RGLengths",]$dim} %>% 
      {gsub(" x 1","",.)} %>% as.integer
    if(nRG > 0){
      dt <- data.table(
        values = h5read(tmpFile, paste0("Fragments/",chrTmp,"/RGValues")),
        lengths = h5read(tmpFile, paste0("Fragments/",chrTmp,"/RGLengths"))
        )
    }else{
      dt <- NULL
    }
    dt
  })
  names(dtList) <- chunkNames
  dt <- Reduce("rbind", dtList)
  dt <- dt[, sum(lengths.V1),by=list(values.V1)]
  
  #Order to reduce number of hyperslabs
  dt <- dt[order(V1,decreasing=TRUE)]
  bcPass <- BStringSet(dt$values.V1[dt$V1 >= minFrags])
  rm(dt)
  gc()

  #Add To Metadata
  o <- h5write(obj = as.character(bcPass), file = outArrow, name = "Metadata/CellNames")

  #######################################################################################################
  # Second we will dump the chunks into an Arrow File!
  #######################################################################################################
  chunkChr <- stringr::str_split(chunkNames, pattern = "#", simplify=TRUE)[,1]
  currentChunk <- 0
  uniqueChr <- unique(as.character(seqnames(chromSizes))) #sort(unique(chunkChr))

  for(x in seq_along(uniqueChr)){

    .messageDiffTime(sprintf("%s Adding Chromosome %s of %s", prefix, x, length(uniqueChr)), tstart, verbose = verboseAll)
    
    #Determine Ranges and RG Pre-Allocation
    chr <- uniqueChr[x]
    ix <- BiocGenerics::which(chunkChr == chr)

    if(length(ix) == 0){

      #HDF5 Write length 0
      chrPos <- paste0("Fragments/",chr,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chr,"/RGValues")
      o <- h5createGroup(outArrow, paste0("Fragments/",chr))
      o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(0, 2), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(0, 1), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(0, 1), level = 0, size = 4))

    }else{

      chunkNamex <- chunkNames[ix]
      dtListx <- dtList[ix] 

      #Read in Fragments!
      fragments <- lapply(seq_along(chunkNamex), function(x){
        .getFragsFromArrow(tmpFile, chr = chunkNamex[x], out = "IRanges")
      }) %>% Reduce("c", .)
      mcols(fragments)$RG@values <- stringr::str_split(mcols(fragments)$RG@values, pattern = "#", simplify=TRUE)[,2]

      #Order RG RLE based on bcPass
      fragments <- fragments[BiocGenerics::which(mcols(fragments)$RG %bcin% bcPass)]
      fragments <- fragments[order(S4Vectors::match(mcols(fragments)$RG, bcPass))]
      lengthRG <- length(mcols(fragments)$RG@lengths)

      #HDF5 Write
      chrPos <- paste0("Fragments/",chr,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chr,"/RGValues")
      o <- h5createGroup(outArrow, paste0("Fragments/",chr))
      o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(length(fragments), 2), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(lengthRG, 1), level = 0, 
              size = max(nchar(mcols(fragments)$RG@values)) + 1))
      o <- h5write(obj = cbind(start(fragments),width(fragments)), file = outArrow, name = chrPos)
      o <- h5write(obj = mcols(fragments)$RG@lengths, file = outArrow, name = chrRGLengths)
      o <- h5write(obj = mcols(fragments)$RG@values, file = outArrow, name = chrRGValues)

      #Free Some Memory!
      rm(fragments)
      gc()

    }

  }

  #Remove Tmp
  rmf <- file.remove(tmpFile)
  .messageDiffTime(sprintf("%s Finished Constructing ArrowFile", prefix), tstart, verbose = verboseHeader, addHeader = verboseAll)

  return(outArrow)

}

#########################################################################################################
# Methods to turn input file directly to arrow! These may not be memory friendly!
#########################################################################################################

.tsvToArrow <- function(
  tsvFile,
  outArrow,
  chromSizes,
  genome,
  minFrags = 500, 
  sampleName, 
  ...){

  tstart <- Sys.time()
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5write(obj = sampleName, file = outArrow, name = "Metadata/Sample")
  o <- h5createGroup(outArrow, paste0("Fragments"))

  #############################################################
  #Read in TSV File...
  #############################################################
  .messageDiffTime("Reading full inputTSV with data.table::fread", tstart, addHeader = TRUE)
  dt <- fread(tsvFile, sep = "\t", select = c(1,2,3,4))
  setkey(dt, V4) #Set Key
  dt <- dt[order(dt$V4),] #Sort Data.table
  dt <- DataFrame(chr = Rle(dt$V1), start = dt$V2, end = dt$V3, RG = Rle(dt$V4))
  
  #Order to reduce number of hyperslabs
  reOrderRG <- dt$RG@values[order(dt$RG@lengths, decreasing=TRUE)]
  dt <- dt[S4Vectors::match(dt$RG, reOrderRG),]
  gc()

  #############################################################
  #Filter Minimum because this would not be worth keeping at all!
  #############################################################
  idx <- BiocGenerics::which(dt$RG %bcin% dt$RG@values[dt$RG@lengths >= minFrags])
  .messageDiffTime(sprintf("Filtering Fragments less than %s Fragments (%s)", minFrags, 1 - round(length(idx) / nrow(dt),3)), tstart, addHeader = TRUE)
  dt <- dt[idx,]
  remove(idx)
  gc()

  #Add To Metadata
  o <- h5write(obj = as.character(dt$RG@values), file = outArrow, name = "Metadata/CellNames")

  #############################################################
  #Keep Those only in ChromSizes
  #############################################################
  uniqueChr <- unique(paste0(dt$chr@values))
  dt <- dt[BiocGenerics::which(dt$chr %bcin% paste0(seqnames(chromSizes))),]

  #############################################################
  #Check all chromSizes represented...
  #############################################################
  if(nrow(dt) == 0 | !all(paste0(seqnames(chromSizes)) %in% uniqueChr)){
    notIn <- paste0(seqnames(chromSizes)[BiocGenerics::which(seqnames(chromSizes) %bcni% uniqueChr)])
    stop(sprintf("Error no fragments in all seqnames of chromSizes (%s) are you sure this is the correct genome?",notIn))
  }

  #############################################################
  #Write Fragments
  #############################################################
  dt$start <- dt$start + 1
  expAll <- 0
  obsAll <- 0
  seqL <- 0
  for(i in seq_along(uniqueChr)){
    
    .messageDiffTime(sprintf("Writing Chromosome %s of %s to Arrow File!", i, length(uniqueChr)), tstart)
    chri <- uniqueChr[i]
    dti <- dt[BiocGenerics::which(dt$chr==chri),]
    chrPos <- paste0("Fragments/",chri,"/Ranges")
    chrRGLengths <- paste0("Fragments/",chri,"/RGLengths")
    chrRGValues <- paste0("Fragments/",chri,"/RGValues")
    lengthRG <- length(dti$RG@lengths)
    o <- h5createGroup(outArrow, paste0("Fragments/",chri))
    o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(nrow(dti), 2), level = 0))
    o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
    o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(lengthRG, 1), level = 0, size = nchar(dti$RG@values[1]) + 1))
    o <- h5write(obj = cbind(dti$start,dti$end-dti$start), file = outArrow, name = chrPos)
    o <- h5write(obj = dti$RG@lengths, file = outArrow, name = chrRGLengths)
    o <- h5write(obj = dti$RG@values, file = outArrow, name = chrRGValues)

    rm(dti)
    gc()

  }

  .messageDiffTime("Finished Constructing Arrow File!", tstart)

  #Clean Up
  rm(dt)
  gc()

  return(outArrow)

}


#########################################################################################################
# Filtering bad fragments!
#########################################################################################################

.filterCellsFromArrow <- function(inArrow, cellNames){

  tstart <- Sys.time()
  outArrow <- .tempfile(fileext = ".arrow")
  
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5write(obj = .sampleName(inArrow), file = outArrow, name = "Metadata/Sample")
  o <- h5write(obj = paste0(stringr::str_split(cellNames, pattern = "#", simplify = TRUE)[,2]), file = outArrow, name = "Metadata/CellNames")
  o <- h5createGroup(outArrow, paste0("Fragments"))

  allChr <- .availableChr(inArrow)
  
  for(i in seq_along(allChr)){
    
    chr <- allChr[i]
    fragments <- .getFragsFromArrow(inArrow, chr = chr)
    fragments <- fragments[BiocGenerics::which(mcols(fragments)$RG %bcin% cellNames)]

    if(length(fragments) == 0){

      #HDF5 Write length 0
      chrPos <- paste0("Fragments/",chr,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chr,"/RGValues")
      o <- h5createGroup(outArrow, paste0("Fragments/",chr))
      o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(0, 2), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(0, 1), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(0, 1), level = 0, size = 4))

    }else{

      mcols(fragments)$RG@values <- stringr::str_split(mcols(fragments)$RG@values, pattern = "#", simplify= TRUE)[,2]
      lengthRG <- length(mcols(fragments)$RG@lengths)
      
      #HDF5 Write
      chrPos <- paste0("Fragments/",chr,"/Ranges")
      chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
      chrRGValues <- paste0("Fragments/",chr,"/RGValues")
      o <- h5createGroup(outArrow, paste0("Fragments/",chr))
      o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(length(fragments), 2), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
      o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(lengthRG, 1), level = 0, 
              size = max(nchar(mcols(fragments)$RG@values)) + 1))
      o <- h5write(obj = cbind(start(fragments),width(fragments)), file = outArrow, name = chrPos)
      o <- h5write(obj = mcols(fragments)$RG@lengths, file = outArrow, name = chrRGLengths)
      o <- h5write(obj = mcols(fragments)$RG@values, file = outArrow, name = chrRGValues)

    }

    #Free Some Memory!
    rm(fragments)
    gc()

  }

  #Remove old Arrow
  rmf <- file.remove(inArrow)
  out <- .fileRename(from = outArrow, to = inArrow)

  .messageDiffTime("Finished Constructing Filtered Arrow File!", tstart)

  return(inArrow)

}

.fileRename <- function(from, to){

  if(!file.exists(from)){
    stop("Input file does not exist!")
  }
  
  tryCatch({
    
    .suppressAll(file.rename(from, to))
    
  }, error = function(x){

    tryCatch({

      system(paste0("mv ", from, " ", to))

      return(to)

    }, error = function(y){

      stop("File Moving/Renaming Failed!")

    })

  })

}







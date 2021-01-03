#' Create Arrow Files 
#' 
#' This function will create ArrowFiles from input files. These ArrowFiles are the main constituent for downstream analysis in ArchR.
#'
#' @param inputFiles A character vector containing the paths to the input files to use to generate the ArrowFiles.
#' These files can be in one of the following formats: (i) scATAC tabix files, (ii) fragment files, or (iii) bam files.
#' @param sampleNames A character vector containing the names to assign to the samples that correspond to the `inputFiles`.
#' Each input file should receive a unique sample name. This list should be in the same order as `inputFiles`.
#' @param outputNames The prefix to use for output files. Each input file should receive a unique output file name.
#' This list should be in the same order as "inputFiles". For example, if the predix is "PBMC" the output file will be named "PBMC.arrow"
#' @param validBarcodes A list of valid barcode strings to be used for filtering cells read from each input file
#' (see `getValidBarcodes()` for 10x fragment files).
#' @param geneAnnotation The geneAnnotation (see `createGeneAnnotation()`) to associate with the ArrowFiles. This is used downstream
#' to calculate TSS Enrichment Scores etc.
#' @param genomeAnnotation The genomeAnnotation (see `createGenomeAnnotation()`) to associate with the ArrowFiles. This is used
#' downstream to collect chromosome sizes and nucleotide information etc.
#' @param minTSS The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering for use
#' in downstream analyses. Cells with a TSS enrichment score greater than or equal to `minTSS` will be retained. TSS enrichment score
#' is a measurement of signal-to-background in ATAC-seq.
#' @param minFrags The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses.
#' Cells containing greater than or equal to `minFrags` total fragments wll be retained.
#' @param maxFrags The maximum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses.
#' Cells containing greater than or equal to `maxFrags` total fragments wll be retained.
#' @param QCDir The relative path to the output directory for QC-level information and plots for each sample/ArrowFile.
#' @param nucLength The length in basepairs that wraps around a nucleosome. This number is used for identifying fragments as
#' sub-nucleosome-spanning, mono-nucleosome-spanning, or multi-nucleosome-spanning.
#' @param promoterRegion A integer vector describing the number of basepairs upstream and downstream [c(upstream, downstream)] of the TSS to include 
#' as the promoter region for downstream calculation of things like the fraction of reads in promoters (FIP).
#' @param TSSParams A list of parameters for computing TSS Enrichment scores. This includes the `window` which is the size in basepairs
#' of the window centered at each TSS (default 101), the `flank` which is the size in basepairs of the flanking window (default 2000), 
#' and the `norm` which describes the size in basepairs of the flank window to be used for normalization of the TSS enrichment score (default 100). 
#' For example, given `window = 101, flank = 2000, norm = 100`, the accessibility within the 101-bp surrounding the TSS will be normalized
#' to the accessibility in the 100-bp bins from -2000 bp to -1901 bp and 1901:2000.
#' @param excludeChr A character vector containing the names of chromosomes to be excluded from downstream analyses. In most human/mouse
#' analyses, this includes the mitochondrial DNA (chrM) and the male sex chromosome (chrY). This does, however, not exclude the
#' corresponding fragments from being stored in the ArrowFile.
#' @param nChunk The number of chunks to divide each chromosome into to allow for low-memory parallelized reading of the `inputFiles`.
#' Higher numbers reduce memory usage but increase compute time.
#' @param bcTag The name of the field in the input bam file containing the barcode tag information. See `ScanBam` in Rsamtools.
#' @param gsubExpression A regular expression used to clean up the barcode tag string read in from a bam file. For example, if the
#' barcode is appended to the readname or qname field like for the mouse atlas data from Cusanovic* and Hill* et al. (2018), the
#' gsubExpression would be ":.*". This would retrieve the string after the colon as the barcode.
#' @param bamFlag A vector of bam flags to be used for reading in fragments from input bam files. Should be in the format of a
#' `scanBamFlag` passed to `ScanBam` in Rsamtools.
#' @param offsetPlus The numeric offset to apply to a "+" stranded Tn5 insertion to account for the precise Tn5 binding site.
#' See Buenrostro et al. Nature Methods 2013.
#' @param offsetMinus The numeric offset to apply to a "-" stranded Tn5 insertion to account for the precise Tn5 binding site.
#' See Buenrostro et al. Nature Methods 2013.
#' @param addTileMat A boolean value indicating whether to add a "Tile Matrix" to each ArrowFile. A Tile Matrix is a counts matrix that,
#' instead of using peaks, uses a fixed-width sliding window of bins across the whole genome. This matrix can be used in many downstream ArchR operations.
#' @param TileMatParams A list of parameters to pass to the `addTileMatrix()` function. See `addTileMatrix()` for options.
#' @param addGeneScoreMat A boolean value indicating whether to add a Gene-Score Matrix to each ArrowFile. A Gene-Score Matrix uses
#' ATAC-seq signal proximal to the TSS to estimate gene activity.
#' @param GeneScoreMatParams A list of parameters to pass to the `addGeneScoreMatrix()` function. See `addGeneScoreMatrix()` for options.
#' @param force A boolean value indicating whether to force ArrowFiles to be overwritten if they already exist.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param subThreading A boolean determining whether possible use threads within each multi-threaded subprocess if greater than the number of input samples.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param cleamTmp A boolean value that determines whether to clean temp folder of all intermediate ".arrow" files.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
#' 
createArrowFiles <- function(
  inputFiles = NULL, 
  sampleNames = names(inputFiles), 
  outputNames = sampleNames,
  validBarcodes = NULL,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  minTSS = 4,
  minFrags = 1000, 
  maxFrags = 100000,
  QCDir = "QualityControl",
  nucLength = 147,
  promoterRegion = c(2000, 100),
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
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  verbose = TRUE,
  cleanTmp = TRUE,
  logFile = createLogFile("createArrows"),
  filterFrags = NULL,
  filterTSS = NULL
  ){

  ################
  # NEW 
  ################

  #We have decided to force removal of filtered cells and thus we have now added messages describing this change
  #It is a simple change we just want to create a more consistent experience!
  removeFilteredCells <- TRUE
  if(!is.null(filterFrags)){
    message("filterFrags is no longer a valid input. Please use minFrags! Setting filterFrags value to minFrags!")
    minFrags <- filterFrags 
  }

  if(!is.null(filterTSS)){
    message("filterTSS is no longer a valid input. Please use minTSS! Setting filterTSS value to minTSS!")
    minTSS <- filterTSS
  }

  filterTSS <- minTSS
  filterFrags <- minFrags

  ################
  # NEW ^
  ################

  .validInput(input = inputFiles, name = "inputFiles", valid = c("character"))
  if(any(!file.exists(inputFiles))){
    stop("inputFiles do not exist :\n\t", paste0(inputFiles[!file.exists(inputFiles)], collapse = "\n\t"))
  }
  .validInput(input = sampleNames, name = "sampleNames", valid = c("character"))
  .validInput(input = outputNames, name = "outputNames", valid = c("character"))
  if(length(sampleNames) != length(inputFiles)){
    stop("sampleNames must be equal length to inputFiles")
  }
  if(length(outputNames) != length(inputFiles)){
    stop("outputNames must be equal length to inputFiles")
  }
  .validInput(input = validBarcodes, name = "validBarcodes", valid = c("list", "character", "null"))
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
  geneAnnotation <- .validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)
  .validInput(input = filterFrags, name = "filterFrags", valid = c("numeric"))
  .validInput(input = filterTSS, name = "filterTSS", valid = c("numeric"))
  .validInput(input = removeFilteredCells, name = "removeFilteredCells", valid = c("boolean"))
  .validInput(input = minFrags, name = "minFrags", valid = c("numeric"))
  .validInput(input = maxFrags, name = "maxFrags", valid = c("numeric"))
  .validInput(input = QCDir, name = "QCDir", valid = c("character"))
  .validInput(input = nucLength, name = "nucLength", valid = c("integer"))
  .validInput(input = promoterRegion, name = "promoterRegion", valid = c("integer"))
  .validInput(input = TSSParams, name = "TSSParams", valid = c("list"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = nChunk, name = "nChunk", valid = c("integer"))
  .validInput(input = bcTag, name = "bcTag", valid = c("character", "null"))
  .validInput(input = gsubExpression, name = "gsubExpression", valid = c("character", "null"))
  .validInput(input = bamFlag, name = "bamFlag", valid = c("integer", "list", "null"))
  .validInput(input = offsetPlus, name = "offsetPlus", valid = c("integer"))
  .validInput(input = offsetMinus, name = "offsetMinus", valid = c("integer"))
  .validInput(input = addTileMat, name = "addTileMat", valid = c("boolean"))
  .validInput(input = TileMatParams, name = "TileMatParams", valid = c("list"))
  .validInput(input = addGeneScoreMat, name = "addGeneScoreMat", valid = c("boolean"))
  .validInput(input = GeneScoreMatParams, name = "GeneScoreMatParams", valid = c("list"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam","null"))
  .validInput(input = subThreading, name = "subThreading", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = cleanTmp, name = "cleanTmp", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  dir.create(QCDir, showWarnings = FALSE)

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "createArrowFiles Input-Parameters", logFile = logFile)

  if(cleanTmp){
    .logMessage("Cleaning Temporary Files", logFile = logFile)
    o <- .suppressAll(file.remove(list.files("tmp", pattern = ".arrow", full.names = TRUE)))
  }

  #order inputFiles
  o <- tryCatch({
    order(file.info(inputFiles)$size, decreasing = TRUE)
    }, error=function(x){
    seq_along(inputFiles)
  })
  inputFiles <- inputFiles[o]
  sampleNames <- sampleNames[o]
  outputNames <- outputNames[o]

  #Add args to list
  args <- list()
  args <- append(args, mget(names(formals()),sys.frame(sys.nframe()))) #as.list(match.call())
  args$X <- seq_along(inputFiles)
  args$FUN <- .createArrow
  args$registryDir <- file.path(QCDir, "CreateArrowsRegistry")
  args$cleanTmp <- NULL

  if(subThreading){
    h5disableFileLocking()
  }else{
    args$threads <- length(inputFiles)
  }

  args$minTSS <- NULL

  #Run With Parallel or lapply
  outArrows <- tryCatch({
    unlist(.batchlapply(args))
  },error = function(x){
    .logMessage("createArrowFiles has encountered an error, checking if any ArrowFiles completed..", verbose = TRUE, logFile = logFile)
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
    paste0(args$outputNames,".arrow")[file.exists(paste0(args$outputNames,".arrow"))]
  })

  if(subThreading){
    h5enableFileLocking()
  }

  .endLogging(logFile = logFile)

  return(outArrows)

}

#Main Function!
.createArrow <- function(
  i = NULL,
  inputFiles = NULL, 
  sampleNames = NULL, 
  outputNames = paste0("./", sampleName),
  validBarcodes = NULL,
  nChunk = 3,
  offsetPlus = 4,
  offsetMinus = -5,
  geneAnnotation = NULL,
  genomeAnnotation = NULL,
  minFrags = 500, 
  maxFrags = 100000,
  removeFilteredCells = TRUE,
  filterFrags = 1000,
  filterTSS = 4,
  excludeChr = c("chrM", "chrY"),
  gsubExpression = NULL,
  bcTag = "qname",
  bamFlag = NULL,
  QCDir = "QualityControl",
  nucLength = 147,
  promoterRegion = c(2000, 100),
  TSSParams = list(),
  addTileMat = TRUE,
  TileMatParams = list(),
  addGeneScoreMat = TRUE,
  GeneScoreMatParams = list(),
  force = FALSE,
  verbose = TRUE,
  tstart = NULL,
  subThreads = 1,
  logFile = NULL
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  ArrowFiles <- paste0(outputNames, ".arrow")

  inputFile <- inputFiles[i]
  sampleName <- sampleNames[i] 
  outputName <- outputNames[i]
  ArrowFile <- ArrowFiles[i]
  prefix <- sprintf("(%s : %s of %s)", sampleName, i, length(inputFiles))

  .logHeader(sprintf("Creating Arrow File %s %s", ArrowFile, prefix), logFile)

  if(!is.null(validBarcodes)){
    if(length(inputFiles) == 1){
      if(inherits(validBarcodes, "list") | inherits(validBarcodes, "SimpleList")){
        if(length(validBarcodes) != length(inputFiles)){
          stop("validBarcodes must be same length as inputFiles!")
        }
        validBC <- validBarcodes[[sampleName]]
      }else{
        validBC <- validBarcodes
      }
    }else{
      if(length(validBarcodes) != length(inputFiles)){
        stop("validBarcodes must be same length as inputFiles!")
      }
      validBC <- validBarcodes[[sampleName]]
    }
  }else{
    validBC <- NULL
  }

  .logThis(validBC, name = "validBC", logFile)

  QCDir <- file.path(QCDir, sampleName)
  dir.create(QCDir, showWarnings = FALSE)

  .requirePackage("rhdf5", source = "bioc")
  .requirePackage("Rsamtools", source = "bioc")

  #Check if a completed file exists!
  if(file.exists(ArrowFile)){
    .logMessage(sprintf("%s Checking if completed file exists!", prefix), logFile = logFile)
    o <- tryCatch({
      o <- h5read(ArrowFile, "Metadata/Completed") #Check if completed
      if(force){
        .logDiffTime(sprintf("%s Arrow Exists! Overriding since force = TRUE!", prefix), t1 = tstart, verbose = TRUE, addHeader = FALSE, logFile = logFile)
        file.remove(ArrowFile)
      }else{
        .logDiffTime(sprintf("%s Arrow Exists! Marking as completed since force = FALSE!", prefix), t1 = tstart, verbose = TRUE, addHeader = FALSE, logFile = logFile)
        return(ArrowFile)
      }
    },error = function(y){
      .logDiffTime(sprintf("%s Arrow Exists! Overriding since not completed!", prefix), t1 = tstart, verbose = TRUE, addHeader = FALSE, logFile = logFile)
      file.remove(ArrowFile) #If not completed delete
    })
  }

  #############################################################
  #Determine Arrow Construction Method
  #############################################################  
  .logMessage(sprintf("%s Determining Arrow Method to use!", prefix), logFile = logFile)
  fe <- .fileExtension(inputFile)
  if(fe == "gz" | fe == "bgz" | fe == "tsv" | fe == "txt"){   
    if(.isTabix(inputFile)){
      readMethod <- "tabix"
    }else{
      stop("Error only supported inputs are tabix and bam!")
    }
  }else if(.isBam(inputFile)){
    readMethod <- "bam"
  }else{
    stop(sprintf("Read Method for %s Not Recognized!", fe))
  }

  .logDiffTime(sprintf("%s Reading In Fragments from inputFiles (readMethod = %s)", prefix, readMethod), t1 = tstart, verbose = verbose, logFile = logFile)

  if(tolower(readMethod)=="tabix"){

    tmp <- tryCatch({

      .tabixToTmp(tabixFile = inputFile, sampleName = sampleName, validBC = validBC,
                  chromSizes = genomeAnnotation$chromSizes, nChunk = nChunk, threads = subThreads,
                  gsubExpression = gsubExpression, prefix = prefix, verbose = verbose, 
                  tstart = tstart, logFile = logFile)
      
      }, error = function(e){

        errorList <- list(
          tabixFile = inputFile, sampleName = sampleName, validBC = validBC,
          chromSizes = genomeAnnotation$chromSizes, nChunk = nChunk, threads = subThreads,
          gsubExpression = gsubExpression, prefix = prefix, verbose = verbose, 
          tstart = tstart
        )

        .logError(e, fn = ".tabixToTmp", info = prefix, errorList = errorList, logFile = logFile)

    })

    out <- tryCatch({

      .tmpToArrow(tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnnotation$genome, 
                  minFrags = minFrags, maxFrags = maxFrags, sampleName = sampleName, prefix = prefix, threads = subThreads,
                  verbose = verbose, tstart = tstart, 
                  chromSizes = genomeAnnotation$chromSizes, removeFilteredCells = removeFilteredCells, logFile = logFile)
      
      }, error = function(e){
        
        errorList <- list(
          tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnnotation$genome, 
          minFrags = minFrags, maxFrags = maxFrags, sampleName = sampleName, prefix = prefix, threads = subThreads,
          verbose = verbose, tstart = tstart, 
          chromSizes = genomeAnnotation$chromSizes, removeFilteredCells = removeFilteredCells
        )

        .logError(e, fn = ".tmpToArrow", info = prefix, errorList = errorList, logFile = logFile)

    })
  
  }else if(tolower(readMethod)=="bam"){

    tmp <- tryCatch({

        .bamToTmp(bamFile = inputFile, sampleName = sampleName, validBC = validBC,
            chromSizes = genomeAnnotation$chromSizes, bamFlag = bamFlag, 
            bcTag = bcTag, gsubExpression = gsubExpression, nChunk = nChunk, threads = subThreads,
            offsetPlus = offsetPlus, offsetMinus = offsetMinus, prefix = prefix, 
            verbose = verbose, tstart = tstart, logFile = logFile)
      
      }, error = function(e){

        errorList <- list(
          bamFile = inputFile, sampleName = sampleName, validBC = validBC,
          chromSizes = genomeAnnotation$chromSizes, bamFlag = bamFlag, 
          bcTag = bcTag, gsubExpression = gsubExpression, nChunk = nChunk, threads = subThreads,
          offsetPlus = offsetPlus, offsetMinus = offsetMinus, prefix = prefix, 
          verbose = verbose, tstart = tstart
        )

        .logError(e, fn = ".bamToTmp", info = prefix, errorList = errorList, logFile = logFile)

    })

    out <- tryCatch({

      .tmpToArrow(tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnnotation$genome, 
                  minFrags = minFrags, maxFrags = maxFrags, sampleName = sampleName, prefix = prefix, threads = subThreads,
                  verbose = verbose, tstart = tstart, 
                  chromSizes = genomeAnnotation$chromSizes, removeFilteredCells = removeFilteredCells, logFile = logFile)

      }, error = function(e){
        
        errorList <- list(
          tmpFile = tmp, outArrow = ArrowFile, genome = genomeAnnotation$genome, 
          minFrags = minFrags, maxFrags = maxFrags, sampleName = sampleName, prefix = prefix, threads = subThreads,
          verbose = verbose, tstart = tstart, 
          chromSizes = genomeAnnotation$chromSizes, removeFilteredCells = removeFilteredCells
        )

        .logError(e, fn = ".tmpToArrow", info = prefix, errorList = errorList, logFile = logFile)

    })

  }else{
    
    .logStop(sprintf("Read Method : %s Not Recognized!", readMethod), logFile = logFile)

  }
  gc()
  
  #############################################################
  #Compute Fragment Information!
  #############################################################
  .logDiffTime(sprintf("%s Adding Fragment Summary", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)
  fragSummary <- .fastFragmentInfo(
    ArrowFile = ArrowFile, 
    cellNames = .availableCells(ArrowFile), 
    nucLength = nucLength,
    prefix = prefix,
    logFile = logFile
  )

  Metadata <- fragSummary[[1]]
  plot <- tryCatch({

    .logDiffTime(sprintf("%s Plotting Fragment Size Distribution", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)
    
    dir.create(QCDir, showWarnings = FALSE)
    pdf(file.path(QCDir,paste0(sampleName,"-Fragment_Size_Distribution.pdf")),width=4,height=3,onefile=FALSE)
    plotDF <- data.frame(
      x = seq_along(fragSummary[[2]]), 
      percent = 100 * fragSummary[[2]]/sum(fragSummary[[2]])
    )
    gg <- ggplot(plotDF, aes(x = x, y = percent)) + theme_ArchR(baseSize = 7) + 
          geom_line(col = "dodgerblue4", size = 0.5) + 
          coord_cartesian(xlim = c(0,750), ylim = c(0, max(plotDF$percent) * 1.1), expand = FALSE) + 
          xlab("Size of Fragments (bp) \n") + 
          ylab("Fragments (%)") + 
          ggtitle(paste0(sampleName,"\nnFrags = ", round(sum(Metadata[,2])/10^6, 2)," M\nFragment Size Distribution"))
    .fixPlotSize(gg, plotWidth = 4.5, plotHeight = 3.5, height = 3/4)
    dev.off()

  }, error = function(x){

      .logDiffTime("Continuing through after error ggplot for Fragment Size Distribution", t1 = tstart, logFile = logFile)
      #print(x)
      message("\n")

  })
  gc()

  #############################################################
  #Compute TSS Enrichment Scores Information!
  #############################################################
  .logDiffTime(sprintf("%s Computing TSS Enrichment Scores", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)
  TSSParams$TSS <- geneAnnotation$TSS
  TSSParams$ArrowFile <- ArrowFile
  TSSParams$cellNames <- Metadata$cellNames
  TSSParams$threads <- subThreads
  TSSParams$logFile <- logFile
  TSSParams$prefix <- prefix
  TSSOut <- do.call(.fastTSSEnrichment, TSSParams)
  Metadata$TSSEnrichment <- TSSOut$tssScores
  Metadata$ReadsInTSS <- TSSOut$tssReads

  #Filter
  Metadata$Keep <- 1*(Metadata$nFrags >= filterFrags & Metadata$TSSEnrichment >= filterTSS)

  if(sum(Metadata$Keep) < 3){
    mTSS <- round(median(Metadata$TSSEnrichment[Metadata$TSSEnrichment>0]), 2)
    mFrag <- round(median(Metadata$nFrags[Metadata$nFrags>0]), 2)
    .logStop(sprintf("Detected 2 or less cells pass filter (Non-Zero median TSS = %s, median Frags = %s) in file!\n       Check inputs such as 'filterFrags' or 'filterTSS' to keep cells! Exiting!", mTSS, mFrag), logFile = logFile)
  }

  .logDiffTime(sprintf("%s CellStats : Number of Cells Pass Filter = %s ", prefix, sum(Metadata$Keep)), t1 = tstart, verbose = verbose, logFile = logFile)
  .logDiffTime(sprintf("%s CellStats : Median Frags = %s ", prefix, median(Metadata$nFrags[Metadata$Keep==1])), t1 = tstart, verbose = verbose, logFile = logFile)
  .logDiffTime(sprintf("%s CellStats : Median TSS Enrichment = %s ", prefix, median(Metadata$TSSEnrichment[Metadata$Keep==1])), t1 = tstart, verbose = verbose, logFile = logFile)

  plot <- tryCatch({

   .logDiffTime(sprintf("%s Plotting TSS Enrichment Scores", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)
    
    ggtitle <- sprintf("%s\n%s\n%s",
        paste0(sampleName, "\nnCells Pass Filter = ", sum(Metadata$Keep)),
        paste0("Median Frags = ", median(Metadata$nFrags[Metadata$Keep==1])),
        paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment[Metadata$Keep==1]))
      )

    pdf(file.path(QCDir,paste0(sampleName,"-TSS_by_Unique_Frags.pdf")),width=4,height=4,onefile=FALSE)
    gg <- ggPoint(
      x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
      y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001), 
      colorDensity = TRUE,
      xlim = c(2.5, 5),
      ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
      baseSize = 6,
      continuousSet = "sambaNight",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE) + 
      geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
      geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)
    .fixPlotSize(gg, plotWidth = 4, plotHeight = 4)
    dev.off()

  }, error = function(x) {

      .logDiffTime("Continuing through after error ggplot for TSS by Frags", t1 = tstart, logFile = logFile)
      #message(x)
      message("\n")

  })

  #Add To Metadata
  .logDiffTime(sprintf("%s Adding Metadata to Fragments", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

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
  # Compute Other Feature Counts
  #############################################################

  .logDiffTime(sprintf("%s Adding Additional Feature Counts!", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  featureList <- list()
  featureList$Promoter <-  extendGR(
      gr = resize(geneAnnotation$genes, 1, "start"), 
      upstream = promoterRegion[1], 
      downstream = promoterRegion[2]
  )

  if(!is.null(genomeAnnotation$blacklist)){
    if(length(genomeAnnotation$blacklist) > 0){
      featureList$Blacklist <- genomeAnnotation$blacklist
    }
  }

  for(x in seq_along(featureList)){
    featurex <- featureList[[x]]
    mcols(featurex) <- NULL
    mcols(featurex)$FeatureName <- names(featureList)[x]
    if(x == 1){
      feature <- featurex
    }else{
      feature <- c(feature, featurex)
    }

  }

  countMat <- .fastFeatureCounts(
    feature = feature, 
    ArrowFile = ArrowFile, 
    cellNames = Metadata$cellNames, 
    threads = subThreads, 
    prefix = prefix, 
    logFile = logFile
  )

  for(j in seq_len(nrow(countMat))){

    featureName <- rownames(countMat)[j]
    Metadata[, paste0("ReadsIn", rownames(countMat)[j])] <- as.vector(countMat[j, ])
    Metadata[, paste0(rownames(countMat)[j], "Ratio")] <- as.vector(countMat[j, ]) / (Metadata$nFrags * 2)
    o <- h5write(obj = as.vector(countMat[j, ]), file = ArrowFile, name = paste0("Metadata/ReadsIn", rownames(countMat)[j]))
    o <- h5write(obj = as.vector(countMat[j, ]) / (Metadata$nFrags * 2), file = ArrowFile, name = paste0("Metadata/", rownames(countMat)[j]), "Ratio")

  }

  .logDiffTime(sprintf("%s Finished Adding Additional Feature Counts!", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

  #############################################################
  # Remove Cells That Are Filtered?
  #############################################################
  
  .logThis(Metadata, paste0(prefix, " Metadata"), logFile = logFile)
  saveRDS(Metadata, file.path(QCDir,paste0(sampleName,"-Pre-Filter-Metadata.rds")))

  if(removeFilteredCells){
    .logDiffTime(sprintf("%s Removing Fragments from Filtered Cells", prefix), t1 = tstart, verbose = verbose, logFile = logFile)
    idx <- which(Metadata$Keep == 1)
    o <- .filterCellsFromArrow(
      inArrow = ArrowFile, 
      cellNames = Metadata$cellNames[idx], 
      verbose = verbose, 
      prefix = prefix, 
      logFile = logFile,
      tstart = tstart
    )
    o <- h5write(obj = Metadata$Keep[idx], file = ArrowFile, name = "Metadata/PassQC")
    o <- h5write(obj = Metadata$nFrags[idx], file = ArrowFile, name = "Metadata/nFrags")
    o <- h5write(obj = Metadata$nMonoFrags[idx], file = ArrowFile, name = "Metadata/nMonoFrags")
    o <- h5write(obj = Metadata$nDiFrags[idx], file = ArrowFile, name = "Metadata/nDiFrags")
    o <- h5write(obj = Metadata$nMultiFrags[idx], file = ArrowFile, name = "Metadata/nMultiFrags")
    o <- h5write(obj = ((Metadata$nDiFrags + Metadata$nMultiFrags) / Metadata$nMonoFrags)[idx], file = ArrowFile, name = "Metadata/NucleosomeRatio")
    o <- h5write(obj = Metadata$TSSEnrichment[idx], file = ArrowFile, name = "Metadata/TSSEnrichment")
    o <- h5write(obj = Metadata$ReadsInTSS[idx], file = ArrowFile, name = "Metadata/ReadsInTSS")
    o <- h5write(obj = Metadata$ReadsInPromoter[idx], file = ArrowFile, name = "Metadata/ReadsInPromoter")
    o <- h5write(obj = Metadata$PromoterRatio[idx], file = ArrowFile, name = "Metadata/PromoterRatio")
    if(!is.null(genomeAnnotation$blacklist)){
      if(length(genomeAnnotation$blacklist) > 0){
        o <- h5write(obj = Metadata$ReadsInBlacklist[idx], file = ArrowFile, name = "Metadata/ReadsInBlacklist")
        o <- h5write(obj = Metadata$BlacklistRatio[idx], file = ArrowFile, name = "Metadata/BlacklistRatio")
      }
    }
    Metadata <- Metadata[idx, , drop = FALSE]
  }

  #############################################################
  # Create Tile Matrix
  #############################################################
  if(addTileMat){
    
    outT <- tryCatch({

      .logDiffTime(sprintf("%s Adding TileMatrix!", prefix), t1 = tstart, verbose = verbose, logFile = logFile)
      TileMatParams$i <- 1
      TileMatParams$ArrowFiles <- ArrowFile
      TileMatParams$cellNames <- Metadata$cellNames
      chromLengths <- end(genomeAnnotation$chromSizes)
      names(chromLengths) <- paste0(seqnames(genomeAnnotation$chromSizes))
      TileMatParams$chromLengths <- chromLengths
      TileMatParams$blacklist <- genomeAnnotation$blacklist
      TileMatParams$force <- TRUE
      TileMatParams$excludeChr <- excludeChr
      TileMatParams$logFile <- logFile
      tileMat <- suppressMessages(do.call(.addTileMat, TileMatParams))
      gc()
      .logDiffTime(sprintf("%s Finished Adding TileMatrix!", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

    }, error =function(e){

      .logError(e, fn = ".addTileMat", info = prefix, errorList = TileMatParams, logFile = logFile)

    })

  }

  #############################################################
  # Add Gene Score Matrix
  #############################################################
  if(addGeneScoreMat){

    outG <- tryCatch({

      .logDiffTime(sprintf("%s Adding GeneScoreMatrix!", prefix), t1 = tstart, verbose = verbose, logFile = logFile)
      GeneScoreMatParams$i <- 1
      GeneScoreMatParams$ArrowFiles <- ArrowFile
      GeneScoreMatParams$genes <- geneAnnotation$genes
      GeneScoreMatParams$cellNames <- Metadata$cellNames
      GeneScoreMatParams$blacklist <- genomeAnnotation$blacklist
      GeneScoreMatParams$force <- TRUE
      GeneScoreMatParams$excludeChr <- excludeChr
      GeneScoreMatParams$subThreads <- subThreads
      GeneScoreMatParams$logFile <- logFile
      geneScoreMat <- suppressWarnings(suppressMessages(do.call(.addGeneScoreMat, GeneScoreMatParams)))
      gc()
      .logDiffTime(sprintf("%s Finished Adding GeneScoreMatrix!", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)


    }, error = function(e){

      .logError(e, fn = ".addGeneScoreMat", info = prefix, errorList = GeneScoreMatParams, logFile = logFile)

    })
    
  }

  o <- h5closeAll()

  .logDiffTime(sprintf("%s Finished Creating Arrow File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)
  o <- h5write(obj = "Finished", file = ArrowFile, name = "Metadata/Completed")

  ArrowFile <- paste0(outputName, ".arrow")  

  return(ArrowFile)

}

#########################################################################################################
# QC Methods 
#########################################################################################################

.fastFragmentInfo <- function(
  ArrowFile = NULL,
  cellNames = .availableCells(ArrowFile),
  nucLength = 147,
  prefix = NULL,
  logFile = NULL
  ){

  .logDiffTime(paste0(prefix, " Computing fragment size info!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  #Info to get
  matNuc <- matrix(0, nrow = length(cellNames), ncol = 3)
  nFrags <- rep(0, length(cellNames))
  fragDist <- rep(0, 1000)

  chrArrow <- .availableChr(ArrowFile)
  for(x in seq_along(chrArrow)){

    countFrags <- tryCatch({

      #Read in frags
      fragx <- .getFragsFromArrow(ArrowFile = ArrowFile, chr = chrArrow[x], out = "IRanges", cellNames = cellNames)

      if(length(fragx) > 0){

        mcols(fragx)$RG@values <- S4Vectors::match(mcols(fragx)$RG@values, cellNames)
        nFragsx <- S4Vectors:::tabulate(mcols(fragx)$RG, nbins = length(cellNames))

        #Get Distributions
        fragDistx <- tabulate(width(fragx), nbins = 1000)
        w <- trunc(width(fragx)/nucLength) + 1
        w[w > 3] <- 3

        #Get Nuc Info
        matNucx <- tabulate2dCpp(
          w, xmin = 1, xmax = 3, 
          as.integer(mcols(fragx)$RG), ymin = 1, ymax = length(cellNames)
        )   

        list(nFrags = nFragsx, matNuc = matNucx, fragDist = fragDistx)

      }else{

        list(nFrags = NULL, matNuc = NULL, fragDist = NULL)

      }


    }, error = function(e){

      errorList <- list(
        x = x,
        chr = chrArrow[x],
        fragments = if(exists("fragx", inherits = FALSE)) fragx else "Error with Fragments!",
        nFrags = if(exists("nFragsx", inherits = FALSE)) nFragsx else "Error with Counting Fragment Barcodes!",
        fragDist = if(exists("fragDistx", inherits = FALSE)) fragDistx else "Error with Counting Fragment Width Numbers!",
        matNucx = if(exists("matNucx", inherits = FALSE)) matNucx else "Error with Counting Fragment Nucleosome Spanning Numbers!"
      )

      .logError(e, fn = ".fastFragmentInfo", info = prefix, errorList = errorList, logFile = logFile)

    })

    if(!is.null(countFrags[[1]])){
      nFrags <- nFrags + countFrags[[1]]
    }
    
    if(!is.null(countFrags[[2]])){
      matNuc <- matNuc + countFrags[[2]]
    }

    if(!is.null(countFrags[[3]])){
      fragDist <- fragDist + countFrags[[3]]
    }

  }

  df <- DataFrame(matNuc)
  colnames(df) <- c("nMonoFrags", "nDiFrags", "nMultiFrags")
  df$cellNames <- cellNames
  df$nFrags <- nFrags
  df <- df[,c("cellNames","nFrags","nMonoFrags", "nDiFrags", "nMultiFrags")]

  out <- list(dfSummary = df, fragDistribution  = fragDist)

  .logDiffTime(paste0(prefix, " Finished computing fragment size info!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  return(out)

}

.fastTSSEnrichment <- function(
  TSS = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL, 
  window = 101, 
  norm = 100, 
  flank = 2000, 
  minNorm = 0.2, #Handles low cell reads inflated TSS values
  maxFragSize = NULL,
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){

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
  #.logThis(tssFeatures, paste0(prefix, " tssFeatures"), logFile = logFile)

  #Counting
  countList <- .fastTSSCounts(
    feature = tssFeatures, 
    ArrowFile = ArrowFile, 
    cellNames = cellNames, 
    maxFragSize = maxFragSize,
    threads = threads, 
    prefix = prefix,
    logFile = logFile
  )

  #Normalize per BP
  cWn <- countList$nWindow / window
  cFn <- countList$nFlank / norm

  #Compute scores
  tssScores <- 2 * cWn / (pmax(cFn, minNorm))
  names(tssScores) <- cellNames
  tssScores <- round(tssScores, 3)

  .logDiffTime(paste0(prefix, " Computed TSS Scores!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  return(list(tssScores=tssScores, tssReads=cWn * window))

}

.fastTSSCounts <- function(
  feature = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL,
  maxFragSize = NULL, 
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){
  
  tstart1 <- Sys.time()
  featureList <- split(feature, seqnames(feature))
  chrArrow <- .availableChr(ArrowFile)
  featureList <- featureList[chrArrow]
  if(length(featureList)==0){
    stop("Error No Overlap in Chromosomes and TSS Chromosomes!")
  }

  #Count
  countDF <- .safelapply(seq_along(featureList), function(x){

    tryCatch({

      #Count Vector
      nWindow <- rep(0, length(cellNames))
      names(nWindow) <- cellNames

      nFlank <- rep(0, length(cellNames))
      names(nFlank) <- cellNames

      ###############################################################################
      # Get Fragments
      ###############################################################################
      featurex <- featureList[[x]]
      fragments <- .getFragsFromArrow(
        ArrowFile = ArrowFile, 
        chr = names(featureList)[x], 
        out = "IRanges", 
        cellNames = cellNames
      )

      if(length(fragments) > 0){
        if(!is.null(maxFragSize)){
          fragments <- fragments[width(fragments) <= maxFragSize]
        }
      }

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
          
          mat <- tabulate2dCpp(
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

      names(nWindow) <- NULL
      names(nFlank) <- NULL

      DataFrame(nWindow = nWindow, nFlank = nFlank)

    }, error = function(e){

      errorList <- list(
        x = x,
        chr = names(featureList)[x],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!",
        features = if(exists("featurex", inherits = FALSE)) featurex else "Error with Features!" 
      )

      .logError(e, fn = ".fastTSSCounts", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  for(i in seq_along(countDF)){
    if(i == 1){
      cDF <- countDF[[i]]
    }else{
      cDF$nWindow <- cDF$nWindow + countDF[[i]]$nWindow
      cDF$nFlank <- cDF$nFlank + countDF[[i]]$nFlank
    }
  }

  rm(countDF)
  gc()

  out <- list(nWindow = cDF[,1] , nFlank = cDF[, 2])

  return(out)

}

.fastFeatureCounts <- function(
  feature = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL, 
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){

  tstart1 <- Sys.time()
  featureNames <- unique(mcols(feature)$FeatureName)
  featureList <- split(feature, seqnames(feature))
  chrArrow <- .availableChr(ArrowFile)
  featureList <- featureList[chrArrow]
  if(length(featureList)==0){
    .logStop("Error No Overlap in Chromosomes and Feature Chromosomes!", logFile = logFile)
  }

  #Count
  countMat <- .safelapply(seq_along(featureList), function(x){

    tryCatch({

      m <- matrix(0, nrow = length(featureNames), ncol = length(cellNames))
      rownames(m) <- featureNames

      ###############################################################################
      # Get Fragments
      ###############################################################################
      featurex <- featureList[[x]]
      fragments <- .getFragsFromArrow(
        ArrowFile = ArrowFile, 
        chr = names(featureList)[x], 
        out = "IRanges", 
        cellNames = cellNames
      )

      if(length(fragments) > 0){

        mcols(fragments)$RG@values <- match(mcols(fragments)$RG@values, cellNames)
        mcols(featurex)$typeIdx <- match(mcols(featurex)$FeatureName, featureNames)
        
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
          
          m <- m + tabulate2dCpp(
                x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
                xmin = 1,
                xmax = length(cellNames),
                y = mcols(featurex)$typeIdx[queryHits(o)],
                ymin = 1,
                ymax = nrow(m)
            )

          rm(o)

        }
        
        rm(fragments)
        gc()

      }

      m

    }, error = function(e){

      errorList <- list(
        x = x,
        chr = names(featureList)[x],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!",
        features = if(exists("featurex", inherits = FALSE)) featurex else "Error with Features!" 
      )

      .logError(e, fn = ".fastFeatureCounts", info = prefix, errorList = errorList, logFile = logFile)

    })

  }, threads = threads)

  countMat <- Reduce("+", countMat)
  colnames(countMat) <- cellNames

  countMat

}

#########################################################################################################
# Validate Input Files!
#########################################################################################################
.isTabix <- function(file = NULL){
  tryCatch({
    TabixFile(file)
    TRUE
  }, error = function(x){
    tryCatch({
      message("Attempting to index ", file," as tabix..")
      indexTabix(file, format = "bed")
      TRUE
    }, error = function(y){
      FALSE
    })
  })
}

.isBam <- function(file = NULL){
  tryCatch({
    bf <- BamFile(file)
    if(file.exists(paste0(file, ".bai"))){
      return(TRUE) 
    }else{
      stop("Bam File Index Does Not Exist! Trying again...")
    }
  }, error = function(x){
    tryCatch({
      message("Attempting to index ", file," as bam...")
      indexBam(file)
      TRUE
    }, error = function(y){
      FALSE
    })
  })
}

#########################################################################################################
# Methods to Turn Input File into a Temp File that can then be Efficiently converted to an Arrow!
#########################################################################################################

.tabixToTmp <- function(
  tabixFile = NULL, 
  sampleName = NULL,
  validBC = NULL,
  tmpFile = .tempfile(pattern = paste0("tmp-",sampleName,"-arrow"), fileext=".arrow"),
  chromSizes = NULL, 
  nChunk = 3,
  gsubExpression = NULL, 
  verbose = TRUE,
  prefix = "",
  tstart = NULL,
  threads = 1,
  logFile = NULL
  ){

  .requirePackage("Rsamtools", source = "bioc")

  #######################################################################################################
  # We will dump a chunked genome into an Hdf5 file in a memory efficient manner!
  #######################################################################################################

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  tstart2 <- Sys.time()
  
  .logDiffTime(sprintf("%s Tabix Bed To Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  o <- h5closeAll()
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("Fragments"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = tmpFile, name = "ArchRVersion")
  o <- h5write(obj = "tmp", file = tmpFile, name = "Metadata/Sample")

  tileChromSizes <- unlist(GenomicRanges::tile(chromSizes, nChunk))
  mcols(tileChromSizes)$chunkName <- paste0(seqnames(tileChromSizes),"#chunk",seq_along(tileChromSizes))

  .logThis(tileChromSizes, name = "tileChromSizes", logFile = logFile)

  readTiledChrom <- .safelapply(seq_along(tileChromSizes), function(x){

    tryCatch({

      errorCheck <- 0

      if(threads == 1){
        if(x %% 10 == 0){
          .logDiffTime(sprintf("%s Reading TabixFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
            verbose = verbose,  logFile = logFile)
        }
      }else{
        if(x %% (2 * threads + 1) == 0){
          .logDiffTime(sprintf("%s Reading TabixFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
                    verbose = verbose, logFile = logFile)
        }
      }

      dt <- tryCatch({
        Rsamtools::scanTabix(tabixFile, param = tileChromSizes[x])[[1]] %>%
          textConnection %>% 
          {tryCatch(read.table(.), error = function(e) NULL)} %>% 
          {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}
      }, error = function(f){
          NULL
      })

      if(is.null(dt)){

        dt <- tryCatch({
          Rsamtools::scanTabix(tabixFile, param = .convertGRSeqnames(tileChromSizes[x], method = "remove"))[[1]] %>%
            textConnection %>% 
            {tryCatch(read.table(.), error = function(e) NULL)} %>% 
            {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}
        }, error = function(f){
            NULL
        })

        if(!is.null(dt)){
          .logMessage(msg = paste0(prefix, " found fragments when removed chromosome prefix : ", paste0(tileChromSizes[x])), logFile = logFile)          
        }

      }

      if(x == 1){
        .logThis(dt, name = paste0(prefix, " .tabixToTmp Fragments-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
        .logThis(unique(dt$V4), name = paste0(prefix, " .tabixToTmp Barcodes-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
      }

      if(is.null(dt)){
        return(list(tmpChrFile = NULL, errorCheck = errorCheck))
      }

      #Care for Break Points
      dt <- dt[dt$V2 >= start(tileChromSizes[x]),]

      if(!is.null(gsubExpression)){
        dt$V4 <- gsub(gsubExpression, "", dt$V4)
      }

      #Check for valid barcodes
      if(!is.null(validBC)){
        dt <- dt[dt$V4 %in% validBC, ]
      }

      if(all(!is.null(dt), nrow(dt) > 0)){

        errorCheck <- errorCheck + 1

        if(threads == 1){

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
            dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))
          o <- h5write(obj = cbind(dt$V2,dt$V3 - dt$V2 + 1), file = tmpFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = NULL, errorCheck = errorCheck))

        }else{

          chrTmp <- mcols(tileChromSizes)$chunkName[x]

          #Temporary File
          tmpChrFile <- paste0(gsub(".arrow", "", tmpFile), ".", chrTmp, ".arrow")
          if(file.exists(tmpChrFile)){
            file.remove(tmpChrFile)
          }

          o <- h5createFile(tmpChrFile)

          #Order by bc
          setkey(dt, V4)
          dt <- dt[order(V4)]
          RG <- Rle(paste0(dt$V4))

          chrPos <- paste0(chrTmp, "._.Ranges")
          chrRGLengths <- paste0(chrTmp, "._.RGLengths")
          chrRGValues <- paste0(chrTmp, "._.RGValues")
          lengthRG <- length(RG@lengths)

          o <- .suppressAll(h5createDataset(tmpChrFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGValues, storage.mode = "character", 
            dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))
          
          o <- h5write(obj = cbind(dt$V2,dt$V3 - dt$V2 + 1), file = tmpChrFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpChrFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpChrFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = tmpChrFile, errorCheck = errorCheck))

        }

      }

    }, error = function(e){

      errorList <- list(
        x = x,
        tileChromSizes = tileChromSizes[x],
        fragments = if(exists("dt", inherits = FALSE)) dt else "Error with Fragments!"
      )

      .logError(e, fn = ".tabixToTmp", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  if(threads > 1){

    #Parallel Linkage Hdf5
    .logDiffTime(sprintf("%s Parallel Hdf5 Linkage Temporary File", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

    file.remove(tmpFile)
    o <- h5closeAll()
    fid <- H5Fcreate(tmpFile)

    o <- h5createGroup(fid, paste0("Fragments"))
    o <- h5createGroup(fid, paste0("Metadata"))
    o <- h5write(obj = "Arrow", file = fid, name = "Class")
    o <- h5write(obj = paste0(packageVersion("ArchR")), file = fid, name = "ArchRVersion")
    o <- h5write(obj = "tmp", file = fid, name = "Metadata/Sample")

    tmpChrFiles <- lapply(readTiledChrom, function(x) x$tmpChrFile) %>% unlist

    for(i in seq_along(tmpChrFiles)){

      tmpChrFilei <- tmpChrFiles[i]
      splitNames <- stringr::str_split(h5ls(tmpChrFilei)$name, "._.", simplify=TRUE)
      chunkName <- splitNames[,1]
      group <- splitNames[,2]

      o <- h5createGroup(fid, paste0("Fragments/", chunkName[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[1],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[2],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[2]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[3],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[3]))

    }

    H5Fclose(fid)

  }

  errorCheck <- sum(lapply(readTiledChrom, function(x) x$errorCheck) %>% unlist)

  if(errorCheck == 0){
    if(!is.null(validBC)){
      .logStop("No fragments found! Possible error with validBarcodes!", logFile = logFile)
    }else{
      .logStop("No fragments found!", logFile = logFile)
    }
  }

  .logDiffTime(sprintf("%s Successful creation of Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  return(tmpFile)

}

.bamToTmp <- function(
  bamFile = NULL, 
  sampleName = NULL,
  tmpFile = .tempfile(pattern = paste0("tmp-",sampleName,"-arrow"), fileext=".arrow"), 
  validBC = NULL,
  chromSizes = NULL, 
  bamFlag = NULL,
  nChunk = 3,
  bcTag = "qname",
  gsubExpression = NULL, 
  offsetPlus = 4, 
  offsetMinus = -5, 
  verbose = TRUE,
  prefix = "",
  tstart = NULL,
  threads = 1,
  logFile = NULL
  ){

  .requirePackage("Rsamtools", source = "bioc")

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
  
  .logDiffTime(sprintf("%s Tabix Bam To Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  o <- h5closeAll()
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("Fragments"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = tmpFile, name = "ArchRVersion")
  o <- h5write(obj = "tmp", file = tmpFile, name = "Metadata/Sample")

  tileChromSizes <- unlist(tile(chromSizes, nChunk))
  mcols(tileChromSizes)$chunkName <- paste0(seqnames(tileChromSizes),"#chunk",seq_along(tileChromSizes))

  .logThis(tileChromSizes, name = "tileChromSizes", logFile = logFile)

  readTiledChrom <- .safelapply(seq_along(tileChromSizes), function(x){

    tryCatch({

      errorCheck <- 0

      if(threads == 1){
        if(x %% 10 == 0){
          .logDiffTime(sprintf("%s Reading BamFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
            verbose = verbose)
        }
      }else{
        if(x %% (2 * threads + 1) == 0){
          .logDiffTime(sprintf("%s Reading BamFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
                    verbose = verbose)
        }
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

        if(is.null(dt)){

          dt <- tryCatch({

            scanChunk <- scanBam(bamFile,
                    param = ScanBamParam(
                      flag = bamFlag,
                      what = c("qname", "pos", "isize"),
                      which = .convertGRSeqnames(tileChromSizes[x], method = "remove")
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

          if(!is.null(dt)){
            .logMessage(msg = paste0(prefix, " found fragments when removed chromosome prefix : ", paste0(tileChromSizes[x])), logFile = logFile)          
          }          

        }

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

        if(is.null(dt)){

          dt <- tryCatch({

            scanChunk <- scanBam(bamFile,
                    param = ScanBamParam(
                      flag = bamFlag,
                      what = c("pos", "isize"),
                      tag = bcTag,
                      which = .convertGRSeqnames(tileChromSizes[x], method = "remove")
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

          if(!is.null(dt)){
            .logMessage(msg = paste0(prefix, " found fragments when removed chromosome prefix : ", paste0(tileChromSizes[x])), logFile = logFile)          
          }

        }

      }
      
      #Clean Up Memory
      rm(scanChunk)

      if(is.null(dt)){
        return(list(tmpChrFile = NULL, errorCheck = errorCheck))
      }

      if(x == 1){
        .logThis(dt, name = paste0(prefix, " .bamToTmp Fragments-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
        .logThis(unique(dt$V4), name = paste0(prefix, " .bamToTmp Barcodes-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
      }

      #Care for Break Points
      dt <- dt[dt$start >= start(tileChromSizes[x]),] 
      dt <- dt[dt$end - dt$start >= 10, ] #Minimum Fragment Size

      #Check for valid barcodes
      if(!is.null(validBC)){
        dt <- dt[dt$RG %in% validBC, ]
      }

      if(all(!is.null(dt), nrow(dt) > 0)){

        errorCheck <- errorCheck + 1

        if(threads == 1){

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
                    dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))

          o <- h5write(obj = cbind(dt$start, dt$end - dt$start + 1), file = tmpFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = NULL, errorCheck = errorCheck))

        }else{

          chrTmp <- mcols(tileChromSizes)$chunkName[x]

          #Temporary File
          tmpChrFile <- paste0(gsub(".arrow", "", tmpFile), ".", chrTmp, ".arrow")
          if(file.exists(tmpChrFile)){
            file.remove(tmpChrFile)
          }

          o <- h5createFile(tmpChrFile)

          #Order by bc
          setkey(dt, RG)
          dt <- dt[order(RG)]
          RG <- Rle(dt$RG)

          chrPos <- paste0(chrTmp, "._.Ranges")
          chrRGLengths <- paste0(chrTmp, "._.RGLengths")
          chrRGValues <- paste0(chrTmp, "._.RGValues")
          lengthRG <- length(RG@lengths)

          o <- .suppressAll(h5createDataset(tmpChrFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGValues, storage.mode = "character", 
                    dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))

          o <- h5write(obj = cbind(dt$start, dt$end - dt$start + 1), file = tmpChrFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpChrFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpChrFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = tmpChrFile, errorCheck = errorCheck))

        }

      }

    }, error = function(e){

      errorList <- list(
        x = x,
        tileChromSizes = tileChromSizes[x],
        fragments = if(exists("dt", inherits = FALSE)) dt else "Error with Fragments!"
      )

      .logError(e, fn = ".bamToTmp", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  if(threads > 1){

    #Parallel Linkage Hdf5
    .logDiffTime(sprintf("%s Parallel Hdf5 Linkage Temporary File", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

    file.remove(tmpFile)
    o <- h5closeAll()
    fid <- H5Fcreate(tmpFile)

    o <- h5createGroup(fid, paste0("Fragments"))
    o <- h5createGroup(fid, paste0("Metadata"))
    o <- h5write(obj = "Arrow", file = fid, name = "Class")
    o <- h5write(obj = paste0(packageVersion("ArchR")), file = fid, name = "ArchRVersion")
    o <- h5write(obj = "tmp", file = fid, name = "Metadata/Sample")

    tmpChrFiles <- lapply(readTiledChrom, function(x) x$tmpChrFile) %>% unlist

    for(i in seq_along(tmpChrFiles)){

      tmpChrFilei <- tmpChrFiles[i]
      splitNames <- stringr::str_split(h5ls(tmpChrFilei)$name, "._.", simplify=TRUE)
      chunkName <- splitNames[,1]
      group <- splitNames[,2]

      o <- h5createGroup(fid, paste0("Fragments/", chunkName[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[1],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[2],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[2]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[3],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[3]))

    }

    H5Fclose(fid)

  }

  errorCheck <- sum(lapply(readTiledChrom, function(x) x$errorCheck) %>% unlist)


  if(errorCheck == 0){
    if(!is.null(validBC)){
      .logStop("No fragments found! Possible error with validBarcodes!", logFile = logFile)
    }else{
      .logStop("No fragments found!", logFile = logFile)
    }
  }

  .logDiffTime(sprintf("%s Successful creation of Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  return(tmpFile)

}

#########################################################################################################
# Methods to temp file to arrow!
#########################################################################################################

.tmpToArrow <- function(
  tmpFile = NULL, 
  outArrow = NULL, 
  genome = NULL, 
  chromSizes = NULL,
  minFrags = 500, 
  maxFrags = 100000, 
  sampleName = NULL, 
  verbose = TRUE,
  tstart = NULL,
  removeFilteredCells = FALSE,
  threads = 1,
  prefix = "",
  logFile = NULL
  ){

  if(!removeFilteredCells){
    threads <- 1 #there wont be a later filter step so we wont be linking here!
  }

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  .logDiffTime(sprintf("%s Creating ArrowFile From Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  o <- .suppressAll(file.remove(outArrow))
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = sampleName, file = outArrow, name = "Metadata/Sample")
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5createGroup(outArrow, paste0("Fragments"))

  #Get Info
  chunkNames <- .availableChr(tmpFile)
  .logThis(data.frame(chunkNames = as.character(chunkNames)), name = paste0(prefix, " chunkChrNames"), logFile = logFile)

  #######################################################################################################
  # First we will count the number of occurences per barcode!
  #######################################################################################################
  .logDiffTime(sprintf("%s Counting Unique Barcodes From Temporary File", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)
  o <- h5closeAll()
  dtList <- lapply(seq_along(chunkNames), function(x){
    chrTmp <- chunkNames[x]
    values  <- h5read(tmpFile, paste0("Fragments/", chrTmp, "/RGValues"))
    lengths <- h5read(tmpFile, paste0("Fragments/", chrTmp, "/RGLengths"))
    dt <- tryCatch({
      if(length(values) > 0 & length(lengths) > 0){
        d <- data.table(
          values = values,
          lengths = lengths
        ) 
      }else{
        d <- NULL
      }
      d
    }, error = function(x){
      NULL
    })
    dt
  })
  names(dtList) <- chunkNames
  dt <- Reduce("rbind", dtList)
  dt <- dt[, sum(lengths.V1),by=list(values.V1)]
  
  #Order to reduce number of hyperslabs
  dt <- dt[order(V1,decreasing=TRUE)]
  .logThis(dt, name = paste0(prefix, " BarcodeFrequencyTable"), logFile)
 
  bcPass <- BStringSet(dt$values.V1[dt$V1 >= minFrags & dt$V1 <= maxFrags])
  if(length(bcPass) < 3){
    .logStop(sprintf("Detected 2 or less cells (%s barcodes have greater than 50 fragments) in file!\n       Check inputs such as 'minFrags' or 'maxFrags' to keep cells! Exiting!", sum(dt$V1 > 50)), logFile = logFile)
  }
  .logThis(data.frame(bc = as.character(bcPass)), name = paste0(prefix, " BarcodesMinMaxFrags"), logFile = logFile)

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
  .logThis(data.frame(chr = as.character(uniqueChr)), name = paste0(prefix, " Unique Chromosomes"), logFile = logFile)

  if(threads > 1){
    dir.create("tmp", showWarnings = FALSE)
  }

  outList <- .safelapply(seq_along(uniqueChr), function(x){

    tryCatch({

      .logDiffTime(sprintf("%s Adding Chromosome %s of %s", prefix, x, length(uniqueChr)), t1 = tstart, verbose = FALSE, logFile = logFile)
      
      #Determine Ranges and RG Pre-Allocation
      chr <- uniqueChr[x]
      ix <- BiocGenerics::which(chunkChr == chr)

      if(threads == 1){

        if(length(ix) == 0){

          .logMessage(msg = paste0(prefix, " detected 0 Fragments for ", chr), logFile = logFile)

          #HDF5 Write length 0
          chrPos <- paste0("Fragments/",chr,"/Ranges")
          chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
          chrRGValues <- paste0("Fragments/",chr,"/RGValues")
          o <- h5createGroup(outArrow, paste0("Fragments/",chr))
          o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(0, 2), level = 0))
          o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(0, 1), level = 0))
          o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(0, 1), level = 0, size = 4))

          return(NULL)

        }else{

          chunkNamex <- chunkNames[ix]
          dtListx <- dtList[ix] 

          #Read in Fragments!
          fragments <- lapply(seq_along(chunkNamex), function(i){
            .getFragsFromArrow(tmpFile, chr = chunkNamex[i], out = "IRanges")
          }) %>% Reduce("c", .)
          mcols(fragments)$RG@values <- stringr::str_split(mcols(fragments)$RG@values, pattern = "#", simplify=TRUE)[,2]

          #Order RG RLE based on bcPass
          fragments <- fragments[BiocGenerics::which(mcols(fragments)$RG %bcin% bcPass)]
          fragments <- fragments[order(S4Vectors::match(mcols(fragments)$RG, bcPass))]
          lengthRG <- length(mcols(fragments)$RG@lengths)

          if(x == 1){
            .logThis(fragments, name = paste0(prefix, " .tmpToArrow Fragments-Chr-(",x," of ",length(uniqueChr),")-", uniqueChr[x]), logFile = logFile)
            .logThis(data.frame(bc = as.vector(mcols(fragments)$RG@values)), name = paste0(prefix, " .tmpToArrow Barcodes-Chr-(",x," of ",length(uniqueChr),")-", uniqueChr[x]), logFile = logFile)
          }

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

          return(NULL)

        }

      }else{

        #Temporary File
        tmpChrFile <- file.path("tmp", paste0(gsub(".arrow", "", outArrow), "#", chr, ".arrow"))

        if(file.exists(tmpChrFile)){
          file.remove(tmpChrFile)
        }
        o <- h5createFile(tmpChrFile)

        if(length(ix) == 0){

          #HDF5 Write length 0
          chrPos <- paste0(chr, "._.Ranges")
          chrRGLengths <- paste0(chr, "._.RGLengths")
          chrRGValues <- paste0(chr, "._.RGValues")

          o <- .suppressAll(h5createDataset(tmpChrFile, chrPos, storage.mode = "integer", dims = c(0, 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGLengths, storage.mode = "integer", dims = c(0, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGValues, storage.mode = "character", dims = c(0, 1), level = 0, size = 4))

          return(tmpChrFile)

        }else{

          chunkNamex <- chunkNames[ix]
          dtListx <- dtList[ix] 

          #Read in Fragments!
          fragments <- lapply(seq_along(chunkNamex), function(i){
            .getFragsFromArrow(tmpFile, chr = chunkNamex[i], out = "IRanges")
          }) %>% Reduce("c", .)
          mcols(fragments)$RG@values <- stringr::str_split(mcols(fragments)$RG@values, pattern = "#", simplify=TRUE)[,2]

          #Order RG RLE based on bcPass
          fragments <- fragments[BiocGenerics::which(mcols(fragments)$RG %bcin% bcPass)]
          fragments <- fragments[order(S4Vectors::match(mcols(fragments)$RG, bcPass))]
          lengthRG <- length(mcols(fragments)$RG@lengths)

          if(x == 1){
            .logThis(fragments, name = paste0(prefix, " .tmpToArrow Fragments-Chr-(",x," of ",length(uniqueChr),")-", uniqueChr[x]), logFile = logFile)
            .logThis(data.frame(bc = as.vector(mcols(fragments)$RG@values)), name = paste0(prefix, " .tmpToArrow Barcodes-Chr-(",x," of ",length(uniqueChr),")-", uniqueChr[x]), logFile = logFile)
          }

          chrPos <- paste0(chr, "._.Ranges")
          chrRGLengths <- paste0(chr, "._.RGLengths")
          chrRGValues <- paste0(chr, "._.RGValues")

          #HDF5 Write
          o <- .suppressAll(h5createDataset(tmpChrFile, chrPos, storage.mode = "integer", dims = c(length(fragments), 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGValues, storage.mode = "character", dims = c(lengthRG, 1), level = 0, 
                  size = max(nchar(mcols(fragments)$RG@values)) + 1))

          o <- h5write(obj = cbind(start(fragments),width(fragments)), file = tmpChrFile, name = chrPos)
          o <- h5write(obj = mcols(fragments)$RG@lengths, file = tmpChrFile, name = chrRGLengths)
          o <- h5write(obj = mcols(fragments)$RG@values, file = tmpChrFile, name = chrRGValues)

          #Free Some Memory!
          rm(fragments)
          gc()

          return(tmpChrFile)

        }

      }

    }, error = function(e){

      errorList <- list(
        x = x,
        chromosome = uniqueChr[x],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!"
      )

      .logError(e, fn = ".tmpToArrow", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  if(threads > 1){

    #Parallel Linkage Hdf5
    .logDiffTime(sprintf("%s Parallel Hdf5 Linkage Arrow File", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

    file.remove(outArrow)
    o <- h5closeAll()
    fid <- H5Fcreate(outArrow)

    o <- h5createGroup(fid, paste0("Fragments"))
    o <- h5createGroup(fid, paste0("Metadata"))
    o <- h5write(obj = "Arrow", file = fid, name = "Class")
    o <- h5write(obj = paste0(packageVersion("ArchR")), file = fid, name = "ArchRVersion")
    o <- h5write(obj = sampleName, file = fid, name = "Metadata/Sample")
    o <- h5write(obj = paste0(Sys.Date()), file = fid, name = "Metadata/Date")
    o <- h5write(obj = as.character(bcPass), file = fid, name = "Metadata/CellNames")

    tmpChrFiles <- unlist(outList)

    for(i in seq_along(tmpChrFiles)){

      tmpChrFilei <- tmpChrFiles[i]
      splitNames <- stringr::str_split(h5ls(tmpChrFilei)$name, "._.", simplify=TRUE)
      chunkName <- splitNames[,1]
      group <- splitNames[,2]

      o <- h5createGroup(fid, paste0("Fragments/", chunkName[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[1],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[2],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[2]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[3],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[3]))

    }

    H5Fclose(fid)

    #Remove all sub-linked tmp files
    rmf <- file.remove(list.files(dirname(tmpFile), pattern = paste0(gsub(".arrow", "", basename(tmpFile)), ".chr"), full.names=TRUE))

  }

  #Remove tmp file
  rmf <- file.remove(tmpFile)
  .logDiffTime(sprintf("%s Successful creation of Arrow File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  return(outArrow)

}

#########################################################################################################
# Filtering bad fragments!
#########################################################################################################

.filterCellsFromArrow <- function(
  inArrow = NULL, 
  cellNames  = NULL,
  logFile = NULL,
  prefix = NULL,
  verbose = TRUE,
  tstart = NULL
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  outArrow <- .tempfile(fileext = ".arrow")
  
  .logDiffTime(sprintf("%s Creating Filtered Arrow File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5write(obj = .sampleName(inArrow), file = outArrow, name = "Metadata/Sample")
  o <- h5write(obj = paste0(stringr::str_split(cellNames, pattern = "#", simplify = TRUE)[,2]), file = outArrow, name = "Metadata/CellNames")
  o <- h5createGroup(outArrow, paste0("Fragments"))

  allChr <- .availableChr(inArrow)
  
  for(i in seq_along(allChr)){

    tryCatch({

      chr <- allChr[i]
      fragments <- .getFragsFromArrow(inArrow, chr = chr)
      fragments <- fragments[BiocGenerics::which(mcols(fragments)$RG %bcin% cellNames)]

      if(length(fragments) == 0){

        .logMessage(msg = paste0(prefix, " detected 0 Fragments for ", chr), logFile = logFile)

        #HDF5 Write length 0
        chrPos <- paste0("Fragments/",chr,"/Ranges")
        chrRGLengths <- paste0("Fragments/",chr,"/RGLengths")
        chrRGValues <- paste0("Fragments/",chr,"/RGValues")
        o <- h5createGroup(outArrow, paste0("Fragments/",chr))
        o <- .suppressAll(h5createDataset(outArrow, chrPos, storage.mode = "integer", dims = c(0, 2), level = 0))
        o <- .suppressAll(h5createDataset(outArrow, chrRGLengths, storage.mode = "integer", dims = c(0, 1), level = 0))
        o <- .suppressAll(h5createDataset(outArrow, chrRGValues, storage.mode = "character", dims = c(0, 1), level = 0, size = 4))

      }else{

        if(i == 1){
          .logThis(fragments, name = paste0(prefix, " .filterCellsFromArrow Fragments-Chr-(",i," of ",length(allChr),")-", allChr[i]), logFile = logFile)
          .logThis(data.frame(bc = as.vector(mcols(fragments)$RG@values)), name = paste0(prefix, " .filterCellsFromArrow Barcodes-Chr-(",i," of ",length(allChr),")-", allChr[i]), logFile = logFile)
        }

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

    }, error = function(e){

      errorList <- list(
        x = i,
        chromosome = allChr[i],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!"
      )

      .logError(e, fn = ".filterCellsFromArrow", info = prefix, errorList = errorList, logFile = logFile)

    })

  }

  #Remove old Arrow
  rmf <- file.remove(inArrow)
  rmf <- .suppressAll(file.remove(list.files("tmp", pattern = paste0(gsub(".arrow", "", basename(inArrow)), "#chr"), full.names=TRUE)))
  out <- .fileRename(from = outArrow, to = inArrow)

  .logDiffTime(paste0(prefix, " Finished Constructing Filtered Arrow File!"), t1 = tstart, verbose = verbose, logFile = logFile)

  return(inArrow)

}

.fileRename <- function(from = NULL, to = NULL){

  if(!file.exists(from)){
    stop("Input file does not exist!")
  }
  
  tryCatch({
    
    o <- .suppressAll(file.rename(from, to))

    if(!o){
      stop("retry with mv")
    }
    
  }, error = function(x){

    tryCatch({

      system(paste0("mv ", from, " ", to))

      return(to)

    }, error = function(y){

      stop("File Moving/Renaming Failed!")

    })

  })

}

.convertGRSeqnames <- function(gr, method = "remove"){
  
  if(tolower(method) == "remove"){
  
    gr2 <- GRanges(seqnames = gsub("chr","",seqnames(gr)), ranges = ranges(gr), strand = strand(gr))
    mcols(gr2) <- mcols(gr)
  
  }else{
  
    idx <- grep("chr", seqnames(gr))
    if(length(idx)==0){
    
      seqnames2 <- paste0("chr", seqnames(gr))
    
    }else{
    
      seqnames2 <- seqnames(gr)
      seqnames2[idx] <- paste0("chr", seqnames(gr)[idx])
    
    }
    
    gr2 <- GRanges(seqnames = seqnames2, ranges = ranges(gr), strand = strand(gr))
    mcols(gr2) <- mcols(gr)
  
  }
  
  gr2

}














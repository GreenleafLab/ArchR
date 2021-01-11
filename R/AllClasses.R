#' @useDynLib ArchR
#' @importFrom Rcpp sourceCpp
NULL

setClassUnion("characterOrNull", c("character", "NULL"))
setClassUnion("GRangesOrNull", c("GRanges", "NULL"))

setClass("ArchRProject", 
  representation(
    projectMetadata = "SimpleList",
    projectSummary = "SimpleList",
    sampleColData = "DataFrame",
    sampleMetadata = "SimpleList",
    cellColData = "DataFrame", 
    cellMetadata = "SimpleList", 
    reducedDims = "SimpleList",
    embeddings = "SimpleList",
    peakSet = "GRangesOrNull",
    peakAnnotation = "SimpleList",
    geneAnnotation = "SimpleList",
    genomeAnnotation = "SimpleList",
    imputeWeights = "SimpleList"
  )
)

.validArrowFiles <- function(object){
  errors <- c()
  fe <- file.exists(object@sampleColData$ArrowFiles)
  if(any(!fe)){
    msg <- paste0("\nArrowFiles :\n  ", paste0(object@sampleColData$ArrowFiles[!fe], collapse=",\n  "), "\nDo not exist!")
    errors <- c(errors, msg)    
  }
  if (length(errors) == 0) TRUE else errors
}

setValidity("ArchRProject", .validArrowFiles)

setMethod("show", "ArchRProject",
  
  function(object) {
    scat <- function(fmt, vals=character(), exdent=2, n = 5, ...){
            vals <- ifelse(nzchar(vals), vals, "''")
            lbls <- paste(S4Vectors:::selectSome(vals, maxToShow = n), collapse=" ")
            txt <- sprintf(fmt, length(vals), lbls)
            cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    .ArchRLogo(ascii = "Package")
    cat("class:", class(object), "\n")
    cat("outputDirectory:", object@projectMetadata$outputDirectory, "\n")

    o <- tryCatch({
      object@cellColData$Sample
    }, error = function(x){
      stop(paste0("\nError accessing sample info from ArchRProject.",
        "\nThis is most likely the issue with saving the ArchRProject as an RDS",
        "\nand not with save/loadArchRProject. This bug has mostly been attributed",
        "\nto bioconductors DataFrame saving cross-compatability. We added a fix to this.",
        "\nPlease Try:",
        "\n\trecoverArchRProject(ArchRProj)",
        "\n\nIf that does not work please report to Github: https://github.com/GreenleafLab/ArchR/issues"
      ))
    })

    scat("samples(%d): %s\n", rownames(object@sampleColData))
    scat("sampleColData names(%d): %s\n", names(object@sampleColData))
    scat("cellColData names(%d): %s\n", names(object@cellColData))
    scat("numberOfCells(%d): %s\n", nrow(object@cellColData))
    scat("medianTSS(%d): %s\n", median(object@cellColData$TSSEnrichment))
    scat("medianFrags(%d): %s\n", median(object@cellColData$nFrags))

  }

)

#' Create ArchRProject from ArrowFiles
#' 
#' This function will create an ArchRProject from the provided ArrowFiles.
#'
#' @param ArrowFiles A character vector containing the relative paths to the ArrowFiles to be used.
#' @param outputDirectory A name for the relative path of the outputDirectory for ArchR results. Relative to the current working directory.
#' @param copyArrows A boolean value indicating whether ArrowFiles should be copied into `outputDirectory`.
#' @param geneAnnotation The `geneAnnotation` object (see `createGeneAnnotation()`) to be used for downstream analyses such as calculating
#' TSS Enrichment Scores, Gene Scores, etc.
#' @param genomeAnnotation The `genomeAnnotation` object (see `createGenomeAnnotation()`) to be used for downstream analyses requiring
#' genome information such as nucleotide information or chromosome sizes.
#' @param showLogo A boolean value indicating whether to show the ascii ArchR logo after successful creation of an `ArchRProject`.
#' @param threads The number of threads to use for parallel execution.
#' @export
ArchRProject <- function(
  ArrowFiles = NULL, 
  outputDirectory = "ArchROutput", 
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo = TRUE,
  threads = getArchRThreads()
  ){

  .validInput(input = ArrowFiles, name = "ArrowFiles", valid = "character")
  .validInput(input = outputDirectory, name = "outputDirectory", valid = "character")
  .validInput(input = copyArrows, name = "copyArrows", valid = "boolean")
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
  geneAnnotation <- .validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)
  .validInput(input = showLogo, name = "showLogo", valid = "boolean")
  .validInput(input = threads, name = "threads", valid = c("integer"))

  if(grepl(" ", outputDirectory)){
    stop("outputDirectory cannot have a space in the path! Path : ", outputDirectory)
  }
  dir.create(outputDirectory,showWarnings=FALSE)
  if(grepl(" ", normalizePath(outputDirectory))){
    stop("outputDirectory cannot have a space in the full path! Full path : ", normalizePath(outputDirectory))
  }
  sampleDirectory <- file.path(normalizePath(outputDirectory), "ArrowFiles")
  dir.create(sampleDirectory,showWarnings=FALSE)

  if(is.null(ArrowFiles)){
    stop("Need to Provide Arrow Files!")
  }

  threads <- min(threads, length(ArrowFiles))

  #Validate
  message("Validating Arrows...")
  if(any(!file.exists(ArrowFiles))){
    stop(paste0("Could not find ArrowFiles :\n", paste0(ArrowFiles[!file.exists(ArrowFiles)], collapse="\n")))
  }
  ArrowFiles <- unlist(lapply(ArrowFiles, .validArrow))

  message("Getting SampleNames...")
  sampleNames <- unlist(.safelapply(seq_along(ArrowFiles), function(x){
    if(getArchRVerbose()) message(x, " ", appendLF = FALSE)
    .sampleName(ArrowFiles[x])
  }, threads = threads))
  message("")

  if(any(duplicated(sampleNames))){
    stop("Error cannot have duplicate sampleNames, please add sampleNames that will overwrite the current sample name in Arrow file!")
  }

  if(length(sampleNames) != length(ArrowFiles)) stop("Samples is not equal to input ArrowFiles!")

  if(copyArrows){
    message("Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE")
    for(i in seq_along(ArrowFiles)){
      message(i, " ", appendLF = FALSE)
      cf <- file.copy(ArrowFiles[i], file.path(sampleDirectory, paste0(sampleNames[i], ".arrow")), overwrite = TRUE)
    }
    message("")
    ArrowFiles <- file.path(sampleDirectory, paste0(sampleNames, ".arrow"))
  }

  #Sample Information
  sampleColData <- DataFrame(row.names = sampleNames, ArrowFiles = ArrowFiles)
  sampleMetadata <- SimpleList(lapply(sampleNames, function(x) SimpleList()))
  names(sampleMetadata) <- sampleNames

  #Cell Information
  message("Getting Cell Metadata...")
  metadataList <- .safelapply(seq_along(ArrowFiles), function(x){
    if(getArchRVerbose()) message(x, " ", appendLF = FALSE)
    .getMetadata(ArrowFiles[x])
  }, threads = threads)
  message("")
  message("Merging Cell Metadata...")
  allCols <- unique(c("Sample",rev(sort(unique(unlist(lapply(metadataList,colnames)))))))
  cellColData <- lapply(seq_along(metadataList), function(x){
    mdx <- metadataList[[x]]
    idx <- which(allCols %ni% colnames(mdx))
    if(length(idx) > 0){
      for(i in seq_along(idx)){
        mdx[,allCols[idx]] <- NA 
      }
    }
    mdx[, allCols, drop = FALSE]
  }) %>% Reduce("rbind", .) %>% DataFrame

  message("Initializing ArchRProject...")
  AProj <- new("ArchRProject", 
    projectMetadata = SimpleList(outputDirectory = normalizePath(outputDirectory)),
    projectSummary = SimpleList(),
    sampleColData = sampleColData,
    sampleMetadata = sampleMetadata,
    cellColData = cellColData,
    cellMetadata = SimpleList(),
    reducedDims = SimpleList(),
    embeddings = SimpleList(),
    peakSet = NULL,
    peakAnnotation = SimpleList(),
    geneAnnotation = .validGeneAnnotation(geneAnnotation),
    genomeAnnotation = .validGenomeAnnotation(genomeAnnotation)
  )
  
  if(showLogo){
    .ArchRLogo(ascii = "Logo") 
  }

  AProj <- addProjectSummary(AProj, name = "DateOfCreation", summary = c("Date" = Sys.time()))

  AProj

}

#' Recover ArchRProject if broken sampleColData/cellColData
#' 
#' This function will recover an ArchRProject if it has broken sampleColData or cellColData due to different versions of bioconductor s4vectors.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
recoverArchRProject <- function(ArchRProj){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")

  if(!inherits(ArchRProj@cellColData, "DataFrame")){
    if(inherits(ArchRProj@cellColData, "DFrame")){
      ArchRProj@cellColData <- .recoverDataFrame(ArchRProj@cellColData)
    }else{
      stop("Unrecognized object for DataFrame in cellColData")
    }
  }

  if(!inherits(ArchRProj@sampleColData, "DataFrame")){
    if(inherits(ArchRProj@sampleColData, "DFrame")){
      ArchRProj@sampleColData <- .recoverDataFrame(ArchRProj@sampleColData)
    }else{
      stop("Unrecognized object for DataFrame in sampleColData")
    }
  }

  if(inherits(ArchRProj@peakSet, "GRanges")){

    peakSet <- tryCatch({
   
      ArchRProj@peakSet
   
    }, error = function(x){
      
      pSet <- ArchRProj@peakSet
      pSet@elementMetadata <- .recoverDataFrame(pSet@elementMetadata)
      mdata <- pSet@metadata
      mdata <- lapply(seq_along(mdata), function(x){
        if(inherits(mdata[[x]], "DFrame")){
          .recoverDataFrame(mdata[[x]])
        }else{
          mdata[[x]]
        }
      })
      names(mdata) <- names(pSet@metadata)
      pSet@metadata <- mdata
      pSet

    })

    ArchRProj@peakSet <- peakSet

  }

  ArchRProj

}

.recoverDataFrame <- function(DF){
  
  DFO <- DF

  rnNull <- (attr(DF, "rownames") == "\001NULL\001")[1]
  
  if(!rnNull){
    rn <- attr(DF, "rownames")
    DF <- DataFrame(row.names = attr(DF, "rownames"), attr(DF,"listData"))
  }else{
    DF <- DataFrame(attr(DF,"listData"))
  }
  
  if(length(attr(DFO, "metadata")) != 0){
    
    mdata <- attr(DFO, "metadata")

    mdata <- lapply(seq_along(mdata), function(x){
      
      mx <- mdata[[x]]
      
      if(inherits(mx, "DFrame")){
        rnNullx <- (attr(mx, "rownames") == "\001NULL\001")[1]
        if(!rnNull){
          rnx <- attr(mx, "rownames")
          mx <- DataFrame(row.names = attr(mx, "rownames"), attr(mx,"listData"))
        }else{
          mx <- DataFrame(attr(mx,"listData"))
        }
      }

      if(inherits(mx, "GRanges")){
        mx <- .recoverGRanges(mx)
      }

      mx

    })

    names(mdata) <- names(attr(DFO, "metadata"))
    metadata(DF) <- mdata

  }

  DF

}

.recoverGRanges <- function(GR){

  GRO <- tryCatch({
  
    GR[1]

    GR
  
  }, error = function(x){

    GR@elementMetadata <- .recoverDataFrame(GR@elementMetadata)
    mdata <- GR@metadata
    mdata <- lapply(seq_along(mdata), function(x){
      if(inherits(mdata[[x]], "DFrame")){
        .recoverDataFrame(mdata[[x]])
      }else{
        mdata[[x]]
      }
    })
    names(mdata) <- names(GR@metadata)
    GR@metadata <- mdata    

    GR

  })

  GR <- GRanges(seqnames = GRO@seqnames, ranges = GRO@ranges, strand = GRO@strand)
  metadata(GR) <- GRO@metadata
  if(nrow(GRO@elementMetadata) > 0){
    mcols(GR) <- GRO@elementMetadata
  }
  
  GR

}

#' Load Previous ArchRProject into R
#' 
#' This function will load a previously saved ArchRProject and re-normalize paths for usage.
#' 
#' @param path A character path to an `ArchRProject` directory that was previously saved using `saveArchRProject()`.
#' @param force A boolean value indicating whether missing optional `ArchRProject` components (i.e. peak annotations /
#' background peaks) should be ignored when re-normalizing file paths. If set to `FALSE` loading of the `ArchRProject`
#' will fail unless all components can be found.
#' @param showLogo A boolean value indicating whether to show the ascii ArchR logo after successful creation of an `ArchRProject`.
#' @export
loadArchRProject <- function(
  path = "./", 
  force = FALSE, 
  showLogo = TRUE
  ){

  .validInput(input = path, name = "path", valid = "character")
  .validInput(input = force, name = "force", valid = "boolean")
  .validInput(input = showLogo, name = "showLogo", valid = "boolean")

  path2Proj <- file.path(path, "Save-ArchR-Project.rds")
  
  if(!file.exists(path2Proj)){
    stop("Could not find previously saved ArchRProject in the path specified!")
  }

  ArchRProj <- recoverArchRProject(readRDS(path2Proj))
  outputDir <- getOutputDirectory(ArchRProj)
  outputDirNew <- normalizePath(path)

  #1. Arrows Paths
  ArrowFilesNew <- file.path(outputDirNew, "ArrowFiles", basename(ArchRProj@sampleColData$ArrowFiles))
  if(!all(file.exists(ArrowFilesNew))){
    stop("ArrowFiles do not exist in saved ArchRProject!")
  }
  ArchRProj@sampleColData$ArrowFiles <- ArrowFilesNew

  #2. Annotations Paths

  if(length(ArchRProj@peakAnnotation) > 0){
    
    keepAnno <- rep(TRUE, length(ArchRProj@peakAnnotation))

    for(i in seq_along(ArchRProj@peakAnnotation)){
      #Postions
      if(!is.null(ArchRProj@peakAnnotation[[i]]$Positions)){

        if(tolower(ArchRProj@peakAnnotation[[i]]$Positions) != "none"){

          PositionsNew <- gsub(outputDir, outputDirNew, ArchRProj@peakAnnotation[[i]]$Positions)
          if(!all(file.exists(PositionsNew))){
            if(force){
              keepAnno[i] <- FALSE
              message("Positions for peakAnnotation do not exist in saved ArchRProject!")
            }else{
              stop("Positions for peakAnnotation do not exist in saved ArchRProject!")
            }
          }
          ArchRProj@peakAnnotation[[i]]$Positions <- PositionsNew

        }

      }

      #Matches
      if(!is.null(ArchRProj@peakAnnotation[[i]]$Matches)){

        MatchesNew <- gsub(outputDir, outputDirNew, ArchRProj@peakAnnotation[[i]]$Matches)
        if(!all(file.exists(MatchesNew))){
          if(force){
            message("Matches for peakAnnotation do not exist in saved ArchRProject!")
            keepAnno[i] <- FALSE
          }else{
            stop("Matches for peakAnnotation do not exist in saved ArchRProject!")
          }
        }
        ArchRProj@peakAnnotation[[i]]$Matches <- MatchesNew

      }

    }

    ArchRProj@peakAnnotation <- ArchRProj@peakAnnotation[keepAnno]

  }


  #3. Background Peaks Paths
  if(!is.null(getPeakSet(ArchRProj))){

    if(!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks)){

      bgdPeaksNew <- gsub(outputDir, outputDirNew, metadata(getPeakSet(ArchRProj))$bgdPeaks)

      if(!all(file.exists(bgdPeaksNew))){
        
        if(force){
          message("BackgroundPeaks do not exist in saved ArchRProject!")
          metadata(ArchRProj@peakSet)$bgdPeaks <- NULL
        }else{
          stop("BackgroundPeaks do not exist in saved ArchRProject!")
        }

      }else{

        metadata(ArchRProj@peakSet)$bgdPeaks <- bgdPeaksNew

      }    

    }

  }

  #4. Set Output Directory 

  ArchRProj@projectMetadata$outputDirectory <- outputDirNew

  message("Successfully loaded ArchRProject!")
  if(showLogo){
      .ArchRLogo(ascii = "Logo")
  }  

  ArchRProj

}

#' Save ArchRProject for Later Usage
#' 
#' This function will organize arrows and project output into a directory and save the ArchRProject for later usage.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param outputDirectory A directory path to save all ArchR output and `ArchRProject` to. Default is outputDirectory of the `ArchRProject`.
#' @param overwrite When writing to outputDirectory, overwrite existing files with new files.
#' @param dropCells A boolean indicating whether to drop cells that are not in `ArchRProject` from corresponding Arrow Files.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param threads The number of threads to use for parallel execution.
#' @export
saveArchRProject <- function(
  ArchRProj = NULL,
  outputDirectory = getOutputDirectory(ArchRProj),
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = outputDirectory, name = "outputDirectory", valid = "character")
  .validInput(input = overwrite, name = "overwrite", valid = "boolean")
  .validInput(input = load, name = "load", valid = "boolean")

  if(grepl(" ", outputDirectory)){
    stop("outputDirectory cannot have a space in the path! Path : ", outputDirectory)
  }

  dir.create(outputDirectory, showWarnings=FALSE)
  outputDirectory <- normalizePath(outputDirectory)
  outDirOld <- normalizePath(getOutputDirectory(ArchRProj))
  
  newProj <- ArchRProj
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(ArrowFiles) %in% unique(newProj$Sample)]

  oldFiles <- list.files(outDirOld)
  oldFiles <- oldFiles[oldFiles %ni% c("ArrowFiles", "ImputeWeights", "Save-ArchR-Project.rds")]

  dir.create(file.path(outputDirectory, "ArrowFiles"), showWarnings=FALSE)
  ArrowFilesNew <- file.path(outputDirectory, "ArrowFiles", basename(ArrowFiles))
  names(ArrowFilesNew) <- names(ArrowFiles)

  if(outputDirectory != outDirOld){
    message("Copying ArchRProject to new outputDirectory : ", normalizePath(outputDirectory))
  }

  if(!identical(paste0(ArrowFiles), paste0(ArrowFilesNew))){

    #Copy Arrow Files
    message("Copying Arrow Files...")
    if(dropCells){
      cf <- .copyArrows(
        inArrows = ArrowFiles, 
        outArrows = ArrowFilesNew, 
        cellsKeep = ArchRProj$cellNames, 
        logFile = logFile, 
        threads = threads
      )
    }else{
      for(i in seq_along(ArrowFiles)){
        message(sprintf("Copying Arrow Files (%s of %s)", i, length(ArrowFiles)))
        cf <- file.copy(ArrowFiles[i], ArrowFilesNew[i], overwrite = overwrite)
      }
    }

  }else{

    if(dropCells){
      for(i in seq_along(ArrowFiles)){
        message(sprintf("Moving Arrow Files (%s of %s)", i, length(ArrowFiles)))
        cf <- .fileRename(ArrowFiles[i], paste0(ArrowFiles[i], "-old"))
      }
      cf <- .copyArrows(
        inArrows = paste0(ArrowFiles, "-old"), 
        outArrows = ArrowFilesNew, 
        cellsKeep = ArchRProj$cellNames, 
        logFile = logFile, 
        threads = threads
      )
      fe <- all(file.exists(ArrowFilesNew)) 
      fe2 <- all(file.exists(ArrowFiles)) 
      if(fe & fe2){
        rmf <- file.remove(paste0(ArrowFiles, "-old"))
      }
    }

  }

  if(outputDirectory != outDirOld){

    #Empty Impute Weights If Changing Directory Because This Could Be A Different Set of Cells
    if(!is.null(getImputeWeights(newProj))){
      message("Dropping ImputeWeights...")
      newProj@imputeWeights <- SimpleList()
    }

    #Copy Other Folders 2 layers nested
    message("Copying Other Files...")
    for(i in seq_along(oldFiles)){
      
      fin <- file.path(outDirOld, oldFiles[i])
      fout <- file.path(outputDirectory, oldFiles[i])
      message(sprintf("Copying Other Files (%s of %s): %s", i, length(oldFiles), basename(fin)))
      
      if(dir.exists(fin)){
      
        dir.create(file.path(outputDirectory, basename(fin)), showWarnings=FALSE)
        fin2 <- list.files(fin, full.names = TRUE)
      
        for(j in seq_along(fin2)){
      
          if(dir.exists(fin2[j])){
      
            dir.create(file.path(outputDirectory, basename(fin), basename(fin2)[j]), showWarnings=FALSE)
            fin3 <- list.files(fin2[j], full.names = TRUE)
      
            for(k in seq_along(fin3)){
      
              cf <- file.copy(fin3[k], file.path(fout, basename(fin3[k])), overwrite = overwrite)
      
            }
      
          }else{
      
            cf <- file.copy(fin2[j], file.path(fout, basename(fin2[j])), overwrite = overwrite)
      
          }
      
        }
      
      }else{
      
        cf <- file.copy(fin, fout, overwrite = overwrite)
      
      }

    }

    newProj@sampleColData <- newProj@sampleColData[names(ArrowFilesNew), , drop = FALSE]
    newProj@sampleColData$ArrowFiles <- ArrowFilesNew[rownames(newProj@sampleColData)]
  
  }

  message("Saving ArchRProject...")
  .safeSaveRDS(newProj, file.path(outputDirectory, "Save-ArchR-Project.rds"))
  
  if(load){
    message("Loading ArchRProject...")
    loadArchRProject(path = outputDirectory)
  }

}

#' Subset an ArchRProject for downstream analysis
#' 
#' This function will subset and ArchRProject by cells and save the output to a new directory and re-load the subsetted ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param cells A vector of cells to subset `ArchRProject` by. Alternatively can provide a subset `ArchRProject`.
#' @param outputDirectory A directory path to save all ArchR output and the subsetted `ArchRProject` to.
#' @param dropCells A boolean indicating whether to drop cells that are not in `ArchRProject` from corresponding Arrow Files.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param threads The number of threads to use for parallel execution. 
#' @param force If output directory exists overwrite.
#' @export
subsetArchRProject <- function(
  ArchRProj = NULL,
  cells = getCellNames(ArchRProj),
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = cells, name = "cells", valid = "character")
  .validInput(input = outputDirectory, name = "outputDirectory", valid = "character")

  outDirOld <- getOutputDirectory(ArchRProj)

  if(dir.exists(outputDirectory)){
    if(!force){
      stop("outputDirectory exists! Please set force = TRUE to overwrite existing directory!")
    }
  }

  if(outputDirectory == outDirOld){
    stop("outputDirectory must be different than ArchRProj outputDirectory to properly subset!")
  }

  saveArchRProject(
    ArchRProj = ArchRProj[cells, ], 
    outputDirectory = outputDirectory,
    load = TRUE,
    dropCells = dropCells,
    logFile = logFile,
    threads = threads
  )

}

#Accessor methods adapted from Seurat 
#https://github.com/satijalab/seurat/blob/87e2454817ed1d5d5aa2e9c949b9231f2231802f/R/objects.R

#'Accessing cellColData directly from dollar.sign accessor
#' 
#' This function will allow direct access to cellColData with a `$` accessor.
#'
#' @export
#'
".DollarNames.ArchRProject" <- function(x, pattern = ''){
  cn <- as.list(c("cellNames",colnames(x@cellColData)))
  names(cn) <- c("cellNames",colnames(x@cellColData))
  return(.DollarNames(x = cn, pattern = pattern))
}

#'Accessing cellColData directly from dollar.sign accessor
#' 
#' This function will allow direct access to cellColData with a `$` accessor.
#'
#' @export
#'
"$.ArchRProject" <- function(x, i){
  if(i=="cellNames"){
    return(rownames(x@cellColData))
  }else{
    val <- x@cellColData[[i, drop = TRUE]]
    return(as.vector(val))
  }
}

#' Add directly to cellColData directly from dollar.sign accessor
#' 
#' This function will allow adding directly to cellColData with a `$` accessor.
#'
#' @export
#'
"$<-.ArchRProject" <- function(x, i, value){
  if(i == "Sample"){
    stop("Sample is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(i == "cellNames"){
    stop("cellNames is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(i == "nFrags"){
    stop("nFrags is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(object.size(Rle(value)) < 2 * object.size(value)){ #Check if Rle is more efficient for storage purposes...
    value <- Rle(value)
  }
  if(!is.null(value)){
    if(length(value)==1){
      value <- Rle(value, lengths = nrow(x@cellColData))
    }
  }
  x@cellColData[[i]] <- value
  return(x)
}


#' Subset cells directly from ArchRProject
#' 
#' This function will allow adding directly to cellColData with a `$` accessor.
#'
#' @export
#'
"[.ArchRProject" <- function(x, i, j){
  cD <- x@cellColData
  
  if(missing(i)){
    return(x)
  }

  if(!missing(j)){
    message("Subsetting columns not supported this way to remove columns set them to NULL.\nEx. ArchRProj$Clusters <- NULL\nContinuing just with cell subsetting.")
  }
  
  if (is.logical(i)) {
    if (length(i) != nrow(cD)) {
      stop("Incorrect number of logical values provided to subset cells")
    }
    i <- rownames(cD)[i]
  }
  
  if (is.numeric(i)) {
    i <- rownames(cD)[i]
  }

  if(length(i) == 1){
    stop("Length of subsetting cells must be greater than 1!")
  }

  i <- unique(i)

  #First Subset CellColData
  x@cellColData <- cD[i, , drop=FALSE]
  cellsKeep <- rownames(x@cellColData)

  #Second Remove Impute Weights
  if(length(i) != nrow(cD)){
    if(length(x@imputeWeights) != 0){
    message("Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!")
    }
    x@imputeWeights <- SimpleList()
  }

  #Third Subset ReducedDims
  rD <- x@reducedDims
  rD2 <- lapply(seq_along(rD), function(x){
    rD[[x]][[1]] <- rD[[x]][[1]][cellsKeep, , drop = FALSE]
    rD[[x]]
  }) %>% SimpleList()
  names(rD2) <- names(rD)
  rD <- x@reducedDims
  rm(rD, rD2)

  #Fourth Subset Embeddings
  eD <- x@embeddings
  eD2 <- lapply(seq_along(eD), function(x){
    eD[[x]][[1]] <- eD[[x]][[1]][cellsKeep, , drop = FALSE]
    eD[[x]]
  }) %>% SimpleList()
  names(eD2) <- names(eD)
  x@embeddings <- eD2
  rm(eD, eD2)

  return(x)

}


setMethod(
  f = "colnames",
  signature = c("x" = "ArchRProject"),
  definition = function(x) {
    colnames(x@cellColData)
  }
)

setMethod(
  f = "rownames",
  signature = c("x" = "ArchRProject"),
  definition = function(x) {
    rownames(x@cellColData)
  }
)







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
#' @export
ArchRProject <- function(
  ArrowFiles = NULL, 
  outputDirectory = "ArchR_Output", 
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo = TRUE,
  threads = getArchRThreads()
  ){

  .validInput(input = ArrowFiles, name = "ArrowFiles", valid = "character")
  .validInput(input = outputDirectory, name = "outputDirectory", valid = "character")
  .validInput(input = copyArrows, name = "copyArrows", valid = "boolean")
  .validInput(input = geneAnnotation, name = "geneAnnotation", valid = c("list"))
  .validInput(input = genomeAnnotation, name = "genomeAnnotation", valid = c("list"))
  .validInput(input = showLogo, name = "showLogo", valid = "boolean")

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
    message(x, " ", appendLF = FALSE)
    .sampleName(ArrowFiles[x])
  }, threads = threads))
  message("")

  if(any(duplicated(sampleNames))){
    stop("Error cannot have duplicate sampleNames, please add sampleNames that will overwrite the current sample name in Arrow file!")
  }

  if(length(sampleNames) != length(ArrowFiles)) stop("Samples is not equal to input ArrowFiles!")

  dir.create(outputDirectory,showWarnings=FALSE)
  sampleDirectory <- file.path(normalizePath(outputDirectory), "ArrowFiles")
  dir.create(sampleDirectory,showWarnings=FALSE)

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
    message(x, " ", appendLF = FALSE)
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

#' Recover ArchRProject if Broken sampleColData/cellColData
#' 
#' This function will organize arrows and project output into a directory and save the ArchRProject for later usage.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param copyArrows A boolean indicating whether to copy (`TRUE`) or copy + remove (`FALSE`) original ArrowFiles prior to saving the `ArchRProject`.
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
   
      ArchRProj@peakSet[1]
   
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

#' Save ArchRProject for Later Usage
#' 
#' This function will organize arrows and project output into a directory and save the ArchRProject for later usage.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param copyArrows A boolean indicating whether to copy (`TRUE`) or copy + remove (`FALSE`) original ArrowFiles prior to saving the `ArchRProject`.
#' @export
saveArchRProject <- function(
  ArchRProj = NULL, 
  copyArrows = TRUE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = copyArrows, name = "copyArrows", valid = "boolean")

  outputDir <- getOutputDirectory(ArchRProj)
  
  #Set Up Arrow Files
  ArrowDir <- file.path(basename(outputDir), "ArrowFiles")
  dir.create(ArrowDir, showWarnings = FALSE)

  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFilesNew <- file.path(ArrowDir, basename(ArrowFiles))
  names(ArrowFilesNew) <- names(ArrowFiles)

  for(i in seq_along(ArrowFiles)){
    cf <- file.copy(ArrowFiles[i], ArrowFilesNew[i])
    if(!copyArrows){
      file.remove(ArrowFiles[i])
    }
  }

  ArchRProj@sampleColData$ArrowFiles <- ArrowFilesNew[rownames(ArchRProj@sampleColData)]

  saveRDS(ArchRProj, file.path(outputDir, "Save-ArchR-Project.rds"))

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
  ArrowFilesNew <- file.path(outputDirNew, gsub(paste0(basename(outputDir),"/"),"",ArchRProj@sampleColData$ArrowFiles))
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

  x@cellColData <- cD[i, , drop=FALSE]

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







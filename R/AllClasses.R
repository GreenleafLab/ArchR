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
#' @param geneAnnotation The `geneAnnotation` object (see `createGeneAnnotation()`) to be used for downstream analyses such as calculating TSS Enrichment Scores, Gene Scores, etc.
#' @param genomeAnnotation The `genomeAnnotation` object (see `createGenomeAnnotation()`) to be used for downstream analyses requiring genome information such as nucleotide information or chromosome sizes.
#' @param showLogo A boolean value indicating whether to show the ascii ArchR logo after successful creation of an `ArchRProject`.
#' @export
ArchRProject <- function(
  ArrowFiles = NULL, 
  outputDirectory = "ArchR_Output", 
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo = TRUE
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

  #Validate
  message("Validating Arrows...")
  ArrowFiles <- unlist(lapply(ArrowFiles, .validArrow))

  message("Getting SampleNames...")
  sampleNames <- unlist(lapply(seq_along(ArrowFiles), function(x) .sampleName(ArrowFiles[x])))

  if(any(duplicated(sampleNames))){
    stop("Error cannot have duplicate sampleNames, please add sampleNames that will overwrite the current sample name in Arrow file!")
  }

  if(length(sampleNames) != length(ArrowFiles)) stop("Samples is not equal to input ArrowFiles!")

  dir.create(outputDirectory,showWarnings=FALSE)
  sampleDirectory <- file.path(normalizePath(outputDirectory), "ArrowFiles")
  dir.create(sampleDirectory,showWarnings=FALSE)

  if(copyArrows){
    message("Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE")
    cf <- file.copy(ArrowFiles, file.path(sampleDirectory, paste0(sampleNames, ".arrow")), overwrite = TRUE)
    ArrowFiles <- file.path(sampleDirectory, paste0(sampleNames, ".arrow"))
  }

  #Sample Information
  sampleColData <- DataFrame(row.names = sampleNames, ArrowFiles = ArrowFiles)
  sampleMetadata <- SimpleList(lapply(sampleNames, function(x) SimpleList()))
  names(sampleMetadata) <- sampleNames

  #Cell Information
  metadataList <- lapply(ArrowFiles, .getMetadata)
  intCols <- Reduce("intersect",lapply(metadataList,colnames))
  cellColData <- lapply(metadataList, function(x) x[,intCols]) %>% Reduce("rbind",.)

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
#' @param force A boolean value indicating whether missing optional `ArchRProject` components (i.e. peak annotations / background peaks) should be ignored when re-normalizing file paths. If set to `FALSE` loading of the `ArchRProject` will fail unless all components can be found.
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

  ArchRProj <- readRDS(path2Proj)

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
    return(x@cellColData[[i, drop = TRUE]])
  }
}



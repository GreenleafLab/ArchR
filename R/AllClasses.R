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
#' This function will create an ArchRProject with given ArrowFiles.
#'
#' @param ArrowFiles A character vector containing the names of ArrowFiles to be used.
#' @param outputDirectory A name for the relative path of the outputDirectory for ArchR results 
#' @param copyArrows A boolean indicating whether ArrowFiles should be copied into outputDirectory
#' @param geneAnnotation 
#' @param genomeAnnotation 
#' @param showLogo A boolean indicating whether to show ArchR Logo after successful creation of an ArchRProject.
#' @export
ArchRProject <- function(
  ArrowFiles = NULL, 
  outputDirectory = "ArchR_Output", 
  copyArrows = FALSE,
  geneAnnotation = NULL,
  genomeAnnotation = NULL,
  showLogo = TRUE
  ){

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
    message("Copying ArrowFiles to Ouptut Directory!")
    cf <- file.copy(ArrowFiles, file.path(sampleDirectory, paste0(sampleNames, ".arrow")))
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

  proj <- new("ArchRProject", 
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
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )
  
  if(showLogo){
    .ArchRLogo(ascii = "Logo") 
  }

  proj <- addProjectSummary(proj, name = "DateOfCreation", summary = c("Date" = Sys.time()))

  proj

}

#Validity
.validArchRProject <- function(ArchRProj, ...){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

#' Save ArchRProject for Later Usage
#' 
#' This function will organize arrows and project output into a directory and save the ArchRProject for later usage.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param copyArrows A boolean indicating whether to copy or copy + remove original ArrowFiles prior to saving ArchRProject.
#' @export
saveArchRProject <- function(
  ArchRProj = NULL, 
  copyArrows = TRUE
  ){

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
#' @param path A character path to an ArchRProject directory that was previously saved.
#' @param force A boolean indicating when re-normalizing paths if an annotation/bdgPeaks is not found ignore and continue
#' @param showLogo show ArchRLogo upon completion.
#' @export
loadArchRProject <- function(
  path = "./", 
  force = FALSE, 
  showLogo = TRUE
  ){

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

  #4. Set Output Directory 

  ArchRProj@projectMetadata$outputDirectory <- outputDirNew

  message("Successfully loaded ArchRProject!")
  if(showLogo){
      .ArchRLogo(ascii = "Logo")
  }  

  ArchRProj

}

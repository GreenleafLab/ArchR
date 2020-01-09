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


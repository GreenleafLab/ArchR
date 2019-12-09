#' @useDynLib ArchR
#' @importFrom Rcpp sourceCpp
NULL

setClassUnion("characterOrNull", c("character", "NULL"))
setClassUnion("GRangesOrNull", c("GRanges", "NULL"))
setClassUnion("matrixOrNull",members = c("dgCMatrix","NULL"))

setClass("ArchRProject", 
  representation(
    projectMetadata = "SimpleList",
    projectSummary = "SimpleList",
    sampleColData = "DataFrame",
    sampleMetadata = "SimpleList",
    cellColData = "DataFrame", 
    cellMetadata = "SimpleList", #Where clustering output will go to
    reducedDims = "SimpleList", #Where clustering output will go to
    embeddings = "SimpleList", #Where clustering output will go to
    peakSet = "GRangesOrNull",
    annotations = "SimpleList", #MotifMatches ETC go here
    geneAnnotation = "SimpleList", #genes exons TSS
    genomeAnnotation = "SimpleList", #genome chromSizes BSgenome blacklist
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

#' @export
ArchRProject <- function(
  ArrowFiles=NULL, 
  sampleNames=NULL, 
  outputDirectory = "ArchR_Results", 
  copyArrows = FALSE,
  geneAnnotation = NULL,
  genomeAnnotation = NULL,
  showLogo = TRUE){

  if(is.null(ArrowFiles)){
    stop("Need to Provide Arrow Files!")
  }

  #Validate
  message("Validating Arrows...")
  ArrowFiles <- unlist(lapply(ArrowFiles, .validArrow))

  if(is.null(sampleNames)){
    message("Getting SampleNames...")
    sampleNames <- unlist(lapply(seq_along(ArrowFiles), function(x) .sampleName(ArrowFiles[x])))
  }

  if(any(duplicated(sampleNames))){
    stop("Error cannot have duplicate sampleNames, please add sampleNames that will overwrite the current sample name in Arrow file!")
  }

  if(length(sampleNames) != length(ArrowFiles)) stop("Samples is not equal to input ArrowFiles!")

  dir.create(outputDirectory,showWarnings=FALSE)
  sampleDirectory <- file.path(normalizePath(outputDirectory),"InputArrows")
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
    annotations = SimpleList(),
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation)
  if(showLogo){
    .ArchRLogo(ascii = "Logo") 
  }

  proj <- addProjectSummary(proj, name = "DateOfCreation", summary = c("Date" = Sys.time()))

  proj

}


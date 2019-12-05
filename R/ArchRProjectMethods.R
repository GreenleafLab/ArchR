##########################################################################################
# Validation Methods
##########################################################################################

.validArchRProject <- function(ArchRProj, ...){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

.validGeneAnnotation <- function(geneAnnotation, ...){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

.validGenomeAnnotation <- function(genomeAnnotation, ...){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

##########################################################################################
# Output Directory
##########################################################################################

#' Get outputDirectory from an ArchRProject
#' 
#' This function gets the outputDirectory from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getOutputDirectory <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  outDir <- ArchRProj@projectMetadata$outputDirectory
  return(outDir) 
}

##########################################################################################
# Sample Methods
##########################################################################################

#' Get ArrowFiles from an ArchRProject
#' 
#' This function gets the names of all ArrowFiles associated with a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getArrowFiles <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  af <- ArchRProj@sampleColData$ArrowFiles
  names(af) <- rownames(ArchRProj@sampleColData)
  return(af)
}

#' Get sampleNames from an ArchRProject
#' 
#' This function gets the sampleNames from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getSampleNames <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  snames <- rownames(ArchRProj@sampleColData)
  return(snames)
}

#' Get sampleColData from an ArchRProject
#' 
#' This function gets the sampleColData from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param select select a subset of column names from sampleColData
#' @param drop drop if selecting only one column name
#' @param ... additional args
#' @export
getSampleColData <- function(ArchRProj, select = NULL, drop = FALSE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  scd <- ArchRProj@sampleColData
  if(!is.null(select)){
    if(all(select %in% colnames(scd))){
      scd <- scd[,select,drop=drop]
    }else{
      stop("select Not Found in Colnames of sampleColData:\n", select[select %ni% colnames(scd)])
    }
  }
  return(scd)
}

#' Add information to sampleColData to an ArchRProject
#' 
#' This function adds new data to sampleColData in ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param data The data to add to sampleColData.
#' @param name The column header name to be used for this new data in sampleColData. If a column with this name already exists, you may set "force" equal to TRUE to overwrite the data in this column.
#' @param samples The names of the samples corresponding to data. Typically new data is added to all samples but you may use this argument to only add data to a subset of samples. Samples here data is not added are set to NA.
#' @param force A boolean (TRUE/FALSE) argument that indicates whether or not to overwrite data in a given column when the value passed to "name" already exists as a column name in sampleColData.
#' @param ... additional args
#' @export
addSampleColData <- function(ArchRProj, data = NULL, name = NULL, samples = rownames(sampleColData(ArchRProj)), force = FALSE){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(samples)){
    stop("Error samples must be provided")
  }
  if(is.null(data)){
    stop("Error data must be provided")
  }
  if(is.null(name)){
    stop("Error name is required for new column name!")
  }
  if(length(samples) != length(data)){
    stop("Error samples has to equal length of data!")
  }
  if(name %in% colnames(getSampleColData(ArchRProj))){
    if(force){
      message("Overriding previous entry for ", name)
    }else{
      message(paste0("Error previous entry for ", name, ", Set force = TRUE to override!"))
    }
  }
  ArchRProj@sampleColData[,name] <- NA
  ArchRProj@sampleColData[samples,name] <- data
  return(ArchRProj)
}

##########################################################################################
# Cell Methods
##########################################################################################

#' Get cellNames from an ArchRProject
#' 
#' This function gets the cellNames from a given ArchRProject object.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getCellNames <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  cnames <- rownames(ArchRProj@cellColData)
  return(cnames)
}

#' Get cellColData from an ArchRProject
#' 
#' This function gets the cellColData from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param select A character vector of column names to select from cellColData if you would like to subset the returned data.
#' @param drop A boolean argument to indicate whether additional data.frame information should be dropped if selecting only a single column name.
#' @param ... additional args
#' @export
getCellColData <- function(ArchRProj, select = NULL, drop = FALSE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  ccd <- data.frame(ArchRProj@cellColData)
  if(!is.null(select)){
    ccd2 <- lapply(seq_along(select), function(x){
      tryCatch({
        dplyr::mutate(ccd, tmpNewCol123=eval(parse(text=select[x])))[,"tmpNewCol123"]
      }, error = function(x){
        stop("select Not Found in Colnames of cellColData:\n",x)
      })
    }) %>% Reduce("cbind", .) %>% DataFrame
    colnames(ccd2) <- select
    rownames(ccd2) <- rownames(ccd)
    ccd <- ccd2
  }
  ccd <- DataFrame(ccd)
  if(drop){
    ccd <- ccd[,,drop=drop]
  }
  return(ccd)
}

#' Add information to cellColData in an ArchRProject
#' 
#' This function adds new data to cellColData in a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param data The data to add to cellColData.
#' @param name The column header name to be used for this new data in cellColData. If a column with this name already exists, you may set "force" equal to TRUE to overwrite the data in this column.
#' @param cells The names of the cells corresponding to "data". Typically new data is added to all cells but you may use this argument to only add data to a subset of cells. Cells where data is not added are set to NA.
#' @param force A boolean (TRUE/FALSE) argument that indicates whether or not to overwrite data in a given column when the value passed to "name" already exists as a column name in cellColData.
#' @param ... additional args
#' @export
addCellColData <- function(ArchRProj, data = NULL, name = NULL, cells = getCellNames(ArchRProj), force = FALSE, ...){

  ArchRProj <- .validArchRProject(ArchRProj)

  if(is.null(cells)){
    stop("Error cells must be provided")
  }

  if(is.null(data)){
    stop("Error data must be provided")
  }

  if(is.null(name)){
    stop("Error name is required for new column name!")
  }

  if(length(cells) != length(data)){
    stop("Error cells has to equal length of data!")
  }

  if(name %in% colnames(getCellColData(ArchRProj))){
    if(force){
      message("Overriding previous entry for ", name)
    }else{
      message(paste0("Error previous entry for ", name, ", Set force = TRUE to override!"))
    }
  }

  ArchRProj@cellColData[,name] <- NA
  ArchRProj@cellColData[cells,name] <- data

  return(ArchRProj)
}

##########################################################################################
# PeakSet Methods
##########################################################################################

#' Get the peak set from an ArchRProject
#' 
#' This function gets the peakSet as a GRanges object from an ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getPeakSet <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@peakSet)
}

#' Add a peak set to an ArchRProject
#' 
#' This function adds a peak set as a GRanges object to a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param peakSet A GRanges object containing the set of regions that define all peaks in the desired peak set.
#' @param force If a peakSet object has already been added to the given ArchRProject, the value of "force" determines whether or not to overwrite this peakSet.
#' @param ... additional args
#' @export
addPeakSet <- function(ArchRProj, peakSet, force = FALSE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(ArchRProj@peakSet) | force){
    #Index The Peak Set
    peakSet <- lapply(split(peakSet, seqnames(peakSet)), function(x){
      mcols(x)$idx <- seq_along(x)
      x
    }) %>% Reduce("c", .) %>% sortSeqlevels %>% sort
    ArchRProj@peakSet <- peakSet
  }else{
    stop("Error peakSet exists! Set force=TRUE to override!")
  }
  return(ArchRProj)
}

##########################################################################################
# Genome Annotation Methods
##########################################################################################

#' Get genomeAnnotation from an ArchRProject
#' 
#' This function gets the genomeAnnotation (in format QQQ) from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getGenomeAnnotation <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation)
}

#' Get blacklist from an ArchRProject
#' 
#' This function gets the blacklist (the regions to be excluded from analysis) as a GRanges from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getBlacklist <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$blacklist)
}

#' Get the genome used by an ArchRProject
#' 
#' This function gets the name of the genome from the genomeAnnotation used by a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getGenome <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$genome)
}

#' Get chromSizes from ArchRProject
#' 
#' This function gets the chromosome lengths as a GRanges onject from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getChromSizes <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$chromSizes)
}

#' Get chromLengths from ArchRProject
#' 
#' This function gets the chromosome lengths as a vector from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getChromLengths <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  cS <- ArchRProj@genomeAnnotation$chromSizes
  cL <- end(cS)
  names(cL) <- paste0(seqnames(cS))
  cL
  return(cL)
}

#' @export
.nullGenomeAnnotation <- function(){
  genome <- "none"
  chromSizes <- GRanges()
  blacklist <- GRanges()
  SimpleList(blacklist = blacklist, genome = genome, chromSizes = chromSizes)
}

##########################################################################################
# Gene Annotation Methods
##########################################################################################

#' Get geneAnnotation from an ArchRProject
#' 
#' This function gets the geneAnnotation (in format QQQ) from a given ArchRProject
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getGeneAnnotation <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@geneAnnotation)
}

#' Get the transcription start sites of all genes in an ArchRProject
#' 
#' This function gets the transcription start sites (TSSs) as a GRanges object of all genes from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getTSS <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@geneAnnotation$TSS)
}

#' Get the genes from an ArchRProject
#' 
#' This function gets the genes (in format QQQ) from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param symbols gene symbols to subset
#' @param ... additional args
#' @export
getGenes <- function(ArchRProj, symbols = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  genes <- ArchRProj@geneAnnotation$genes
  genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
  return(genes)
}

#' Get the exons from an ArchRProject
#' 
#' This function gets the exons (in format QQQ) from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param ... additional args
#' @export
getExons <- function(ArchRProj, symbols = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  exons <- ArchRProj@geneAnnotation$exons
  exons <- exons[which(tolower(exons$symbol) %in% tolower(symbols))]
  return(exons)
}

#' @export
.nullGeneAnnotation <- function(){
  genes <- GRanges("chr1", IRanges(1,1), symbol = "a")
  genes <- genes[-1]
  exons <- genes
  TSS <- genes
  SimpleList(genes = genes, exons = exons, TSS = TSS)
}

##########################################################################################
# Dimensionality Reduction / Embedding Methods
##########################################################################################

#' Get dimensionality reduction information stored in an ArchRProject
#' 
#' This function gets a dimensionality reduction object (i.e. UMAP, tSNE, etc) from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param reducedDims QQQ The name of the reducedDims object to retrieve from the designated ArchRProject. Options include QQQ.
#' @param return If set to "mat" or "matrix", the function will return the reducedDims object as a matrix. Otherwise, it will return the full reducedDims object.
#' @param ... additional args
#' @export
getReducedDims <- function(ArchRProj, reducedDims = "TileLSI", return = "matrix", ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(reducedDims %in% names(ArchRProj@reducedDims)){
    if(tolower(return)=="mat" | tolower(return)=="matrix"){
      out <- ArchRProj@reducedDims[[reducedDims]][[1]]
    }else{
      out <- ArchRProj@reducedDims[[reducedDims]]
    }
  }else{
    stop("reducedDims not in computed reduced dims!")
  }
  return(out)
}

#' Get embedding information stored in an ArchRProject
#' 
#' QQQ This function gets an embedding (i.e. QQQ) from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param embedding QQQ The name of the embedding object to retrieve from the designated ArchRProject. Options include QQQ.
#' @param return If set to "df", the function will return the embedding object as a data.frame. Otherwise, it will return the full embedding object.
#' @param ... additional args
#' @export
getEmbedding <- function(ArchRProj, embedding = "IterativeLSI", return = "df", ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(embedding %in% names(ArchRProj@embeddings)){
    if(tolower(return)=="df"){
      out <- ArchRProj@embeddings[[embedding]][[1]]
    }else{
      out <- ArchRProj@embeddings[[embedding]]
    }
  }else{
    stop("embedding not in computed embeddings!")
  }
  return(out)
}

##########################################################################################
# Annotation Methods
##########################################################################################

#' Get annotation from an ArchRProject
#' 
#' This function gets an annotation from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param name QQQ The name of the annotation object to retrieve from the designated ArchRProject. Options include QQQ.
#' @param ... additional args
#' @export
getAnnotation <- function(ArchRProj, name = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@annotations)){
      stop("Name is not in Annotations!")
    }
  }
  ArchRProj@annotations[[name]]
}

#' Get annotation positions from an ArchRProject
#' 
#' This function gets the annotation positions from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param name QQQ The name of the annotation object to retrieve from the designated ArchRProject. Options include QQQ.
#' @param annoName QQQ name to subset with annotations
#' @param ... additional args
#' @export
getPositions <- function(ArchRProj, name = NULL, annoName = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@annotations)){
      stop("Name is not in Annotations!")
    }
  }
  anno <- ArchRProj@annotations[[name]]
  idx <- grep("positions", names(anno), ignore.case=TRUE)
  if(length(idx)==0){
    stop("Annotation does not contain positions!")
  }
  positions <- readRDS(anno[[idx]])
  if(!is.null(annoName)){
    idx <- grep(annoName, names(positions), ignore.case=TRUE)
    if(length(idx)==0){
      stop("Positons do not contain annoName!")
    }
    positions <- positions[idx]
  }
  positions
}

#' Get annotation matches from an ArchRProject
#' 
#' This function gets annotation matches from a given ArchRProject.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param name name of annotations
#' @param annoName name to subset with annotations
#' @param ... additional args
#' @export
getMatches <- function(ArchRProj, name = NULL, annoName = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@annotations)){
      stop("Name is not in Annotations!")
    }
  }
  anno <- ArchRProj@annotations[[name]]
  idx <- grep("matches", names(anno), ignore.case=TRUE)
  if(length(idx)==0){
    stop("Annotation does not contain positions!")
  }
  matches <- readRDS(anno[[idx]])
  if(!is.null(annoName)){
    idx <- grep(annoName, colnames(matches), ignore.case=TRUE)
    if(length(idx)==0){
      stop("Matches do not contain annoName!")
    }
    matches <- matches[, idx, drop=FALSE]
  }
  matches
}

#' Add motif annotations to an ArchRProject
#' 
#' This function adds information about which peaks contain motifs to a given ArchRProject. For each peak, a binary value is stored indicating whether each motif is observed within the peak region.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param motifSet The motif set to be used for annotation. Options include: (i) "JASPAR2016", which gives the 2016 version of JASPAR motifs, (ii) "JASPAR2018", which gives the 2018 version of JASPAR motifs, or (iii) one of "human", "mouse", "encode", or "homer" which gives the corresponding motif sets from the chromVAR package. 
#' @param name QQQ of annotations to store as in ArchRProject
#' @param species QQQ The name of the species relevant to the supplied ArchRProject. This is used for QQQ. By default, this function will attempt to guess the species based on the value from getGenome.
#' @param collection QQQ If one of the JASPAR motif sets is used via "motifSet", this parameter allows you to indicate the JASPAR collection to be used. Possible options include "CORE", QQQ.
#' @param cutOff The p-value cutoff to be used for motif search (see the motimatchr package for more information).
#' @param w The width in basepairs to consider for motif matches (see the motimatchr package for more information).
#' @param ... additional args
#' @export
addMotifAnnotations <- function(
  ArchRProj = NULL,
  motifSet = "JASPAR2018",
  name = "Motif",
  species = NULL,
  collection = "CORE",
  cutOff = 5e-05, 
  w = 7,
  ...
  ){

  .requirePackage("motifmatchr", installInfo='BiocManager::install("motifmatchr")')
  ArchRProj <- .validArchRProject(ArchRProj)

  if(grepl("JASPAR",motifSet) & is.null(species)){
    if(grepl("hg19",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("hg38",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("mm9",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
    if(grepl("mm10",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
  }

  #############################################################
  # Get PWM List adapted from chromVAR!
  #############################################################
  tstart <- Sys.time()
  .messageDiffTime(paste0("Gettting Motif Set, Species : ", species), tstart)

  if(tolower(motifSet)=="jaspar2018"){
    .requirePackage("JASPAR2018",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else if(tolower(motifSet)=="jaspar2016"){
    .requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else if(tolower(motifSet)=="human"){
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("human_pwms_v2")
    motifs <- human_pwms_v2
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else if(tolower(motifSet)=="mouse"){
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("mouse_pwms_v2")
    motifs <- mouse_pwms_v2
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else if(tolower(motifSet)=="encode"){
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("encode_pwms")
    motifs <- encode_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else if(tolower(motifSet)=="homer"){
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("homer_pwms")
    motifs <- homer_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
  }else{
    stop("Error MotifSet Not Recognized!")
  }

  #############################################################
  # Get BSgenome Information!
  #############################################################
  genome <- ArchRProj@genomeAnnotation$genome
  .requirePackage(genome)
  BSgenome <- eval(parse(text = genome))
  BSgenome <- .validBSgenome(BSgenome)

  #############################################################
  # Calculate Motif Positions
  #############################################################
  .messageDiffTime("Finding Motif Positions with motifmatchr!", tstart)
  peakSet <- ArchRProj@peakSet
  motifPositions <- motifmatchr::matchMotifs(
      pwms = motifs,
      subject = peakSet,
      genome = BSgenome, 
      out = "positions", 
      p.cutoff = cutOff, 
      w = w
    )

  #############################################################
  # Motif Overlap Matrix
  #############################################################
  .messageDiffTime("Creating Motif Overlap Matrix", tstart)
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(peakSet, allPositions, ignore.strand=TRUE)
  motifMat <- Matrix::sparseMatrix(
    i = queryHits(overlapMotifs),
    j = match(names(allPositions),names(motifPositions))[subjectHits(overlapMotifs)],
    x = rep(TRUE, length(overlapMotifs)),
    dims = c(length(peakSet), length(motifPositions))
  )
  colnames(motifMat) <- names(motifPositions)
  motifMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = motifMat), rowRanges = peakSet)
  .messageDiffTime("Finished Getting Motif Info!", tstart)

  out <- SimpleList(
      motifSummary = motifSummary,
      motifMatches = motifMat,
      motifPositions = motifPositions,
      motifList = motifs,
      date = Sys.Date()
    )

  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings=FALSE)
  savePositions <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Positions-In-Peaks.rds"))
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Matches-In-Peaks.rds"))

  ArchRProj@annotations[[name]]$Name <- name
  ArchRProj@annotations[[name]]$motifs <- motifs
  ArchRProj@annotations[[name]]$motifSummary <- motifSummary
  ArchRProj@annotations[[name]]$Positions <- savePositions
  ArchRProj@annotations[[name]]$Matches <- saveMatches

  saveRDS(out, file.path(getOutputDirectory(ArchRProj),  "Annotations", paste0(name,"-In-Peaks-Summary.rds")), compress = FALSE)
  saveRDS(out$motifPositions, savePositions, compress = FALSE)
  saveRDS(out$motifMatches, saveMatches, compress = FALSE)

  return(ArchRProj)

}

.summarizeJASPARMotifs <- function(motifs){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)
  
}

.summarizeChromVARMotifs <- function(motifs){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      tags = motifs[[x]]@tags,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame

  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)

}

##########################################################################################
# Additional Methods
##########################################################################################

#' Return the available features that could be selected from a given data matrix within an ArchRProject
#' 
#' This function will identify available features from a given data matrix  (i.e. "GeneScoreMatrix", or "TileMatrix") and return them for downstream plotting utilities.
#' 
#' @param ArchRProj An ArchRProject object.
#' @param useMatrix QQQ The name of the data matrix as stored in the ArrowFiles of the ArchRProject. Options include "TileMatrix", "GeneScoreMatrix", QQQ.
#' @param select QQQ select a specific name with grep
#' @param ignore.case A boolean value indicating whether or not to ignore the case (upper-case / lower-case) when searching via grep for the string passed to "select".
#' @param ... additional args
#' @export
availableFeatures <- function(ArchRProj, useMatrix = "GeneScoreMatrix", select = NULL, ignore.case = TRUE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  fdf <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  if(is.null(select)){
    if(any(duplicated(paste0(fdf$name)))){
      paste0(fdf$seqnames,":",fdf$name)
    }else{
      fdf$name
    }
  }else{
    grepNames <- grep(select, fdf$name, value = TRUE, ignore.case = ignore.case)
    if(any(duplicated(grepNames))){
      grepIdx <- grep(select, fdf$name, ignore.case = ignore.case)
      grepNames <- paste0(fdf$seqnames[grepIdx],":",fdf$name[grepIdx])
    }
    if(all(c("deviations", "z") %in% unique(paste0(fdf$seqnames)))){
      grepNames <- rev(grepNames)
    }
    grepNames
  }

}

#' Plot PDF in outputDirectory of an ArchRProject
#' 
#' This function will save a plot or set of plots as a PDF file in the output directory of a given ArchRProject.
#' 
#' @param name The file name to be used for the output PDF file.
#' @param width The width in inches to be used for the output PDF file.
#' @param height The height in inches to be used for the output PDF.
#' @param ArchRProj An ArchRProject object.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for preventing overwritting of old plots.
#' @param useDingbats A boolean variable that determines wheter to use dingbats characters for plotting points.
#' @param ... additional args to pdf
#' @export
plotPDF <- function(name, width = 8, height = 8, ArchRProj = NULL, addDOC = TRUE, useDingbats = FALSE, ...){
  name <- gsub("\\.pdf", "", name)
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
  pdf(filename, width = width, height = height, useDingbats = useDingbats, ...)
}






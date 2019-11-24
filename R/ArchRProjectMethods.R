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

#' Get outputDirectory in ArchRProject
#' 
#' This function gets outputDirectory from ArchRProject
#' 
#' @param ArchRProj ArchRProject
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

#' Get ArrowFiles in ArchRProject
#' 
#' This function gets ArrowFiles in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getArrowFiles <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  af <- ArchRProj@sampleColData$ArrowFiles
  names(af) <- rownames(ArchRProj@sampleColData)
  return(af)
}

#' Get sampleNames in ArchRProject
#' 
#' This function gets sampleNames in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getSampleNames <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  snames <- rownames(ArchRProj@sampleColData)
  return(snames)
}

#' Get number of cells in ArchRProject
#' 
#' This function gets number of cells in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
nCells <- function(ArchRProj, ...){
  nrow(getCellColData(ArchRProj))
}

#' Get sampleColData in ArchRProject
#' 
#' This function gets sampleColData in ArchRProject
#' 
#' @param ArchRProj ArchRProject
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

#' Add information to sampleColData in ArchRProject
#' 
#' This function adds new data to sampleColData in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param data data to add to sampleColData
#' @param name new column name in sampleColData if already exists set force = TRUE to override
#' @param cells names of samples corresponding to data
#' @param force if name already exists in sampleColData set force = TRUE to override
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

#' Get cellNames in ArchRProject
#' 
#' This function gets cellNames in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getCellNames <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  cnames <- rownames(ArchRProj@cellColData)
  return(cnames)
}

#' Get cellColData in ArchRProject
#' 
#' This function gets sampleColData in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param select select a subset of column names from cellColData can put in a string function
#' @param drop drop if selecting only one column name
#' @param ... additional args
#' @export
getCellColData <- function(ArchRProj, select = NULL, drop = FALSE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  ccd <- data.frame(ArchRProj@cellColData, stringsAsFactors=FALSE)
  if(!is.null(select)){
    ccd2 <- lapply(seq_along(select), function(x){
      tryCatch({
        data.frame(dplyr::mutate(ccd, tmpNewCol123=eval(parse(text=select[x])))[,"tmpNewCol123"], stringsAsFactors=FALSE)
      }, error = function(x){
        stop("select Not Found in Colnames of cellColData:\n",x)
      })
    }) %>% Reduce("cbind", .) %>% {data.frame(.,stringsAsFactors=FALSE)}
    colnames(ccd2) <- select
    rownames(ccd2) <- rownames(ccd)
    ccd <- ccd2
  }
  ccd <- as(ccd, "DataFrame")
  if(drop){
    ccd <- ccd[,,drop=drop]
  }
  return(ccd)
}

#' Add information to cellColData in ArchRProject
#' 
#' This function adds new data to cellColData in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param data data to add to cellColData
#' @param name new column name in cellColData if already exists set force = TRUE to override
#' @param cells names of cells corresponding to data
#' @param force if name already exists in cellColData set force = TRUE to override
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

#' Get PeakSet from ArchRProject
#' 
#' This function gets peakSet from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getPeakSet <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@peakSet)
}

#' Add PeakSet to ArchRProject
#' 
#' This function adds a peakSet to an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param peakSet peakSet as a GRanges
#' @param force force overriding peakSet in ArchRProject
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

#' Get genomeAnnotation from ArchRProject
#' 
#' This function gets genomeAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getGenomeAnnotation <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation)
}

#' Get blacklist from ArchRProject
#' 
#' This function gets the blacklist as a GRanges from genomeAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getBlacklist <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$blacklist)
}

#' Get genome from ArchRProject
#' 
#' This function gets the genome from genomeAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getGenome <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$genome)
}

#' Get chromSizes from ArchRProject
#' 
#' This function gets chromosome lengths as GRanges from genomeAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getChromSizes <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@genomeAnnotation$chromSizes)
}

#' Get chromLengths from ArchRProject
#' 
#' This function gets chromosome lengths as a vector from genomeAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
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

#' Get geneAnnotation from ArchRProject
#' 
#' This function gets geneAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getGeneAnnotation <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@geneAnnotation)
}

#' Get TSS from ArchRProject
#' 
#' This function gets TSS from geneAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param ... additional args
#' @export
getTSS <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@geneAnnotation$TSS)
}

#' Get Genes from ArchRProject
#' 
#' This function gets genes from geneAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param symbols gene symbols to subset
#' @param ... additional args
#' @export
getGenes <- function(ArchRProj, symbols = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  genes <- ArchRProj@geneAnnotation$genes
  genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
  return(genes)
}

#' Get Exons from ArchRProject
#' 
#' This function gets exons from geneAnnotation in ArchRProject
#' 
#' @param ArchRProj ArchRProject
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

#' Get Reduced Dimensions from ArchRProject
#' 
#' This function gets an embedding from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param reducedDims reduced dimensions name in ArchRProject
#' @param return return reduced dimensions as matrix or all info
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

#' Get Embedding from ArchRProject
#' 
#' This function gets an embedding from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param embedding embedding name in ArchRProject
#' @param return return embedding as df or all info
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
# Project Summary
##########################################################################################

#' Get projectSummary from ArchRProject
#' 
#' This function prints the projectSummary from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param returnSummary return summary or just print
#' @param ... additional args
#' @export
getProjectSummary <- function(ArchRProj, returnSummary = FALSE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  pS <- ArchRProj@projectSummary
  o <- lapply(seq_along(pS), function(x){
    message(names(pS)[x], " :")
    p <- lapply(seq_along(pS[[x]]), function(y){
      message("\t", names(pS[[x]])[y], " : ", pS[[x]][y])
    })
    message("\n")
  })
  if(returnSummary){
    pS
  }else{
    0
  }
}

#' Add projectSummary tp ArchRProject
#' 
#' This function adds info to the projectSummary from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param name name of summary input
#' @param summary summary vector
#' @param ... additional args
#' @export
addProjectSummary <- function(ArchRProj, name, summary, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  pS <- ArchRProj@projectSummary
  name <- paste0(length(pS) + 1, "_", name)
  pS <- append(pS, SimpleList(summary))
  names(pS)[length(pS)] <- name
  ArchRProj@projectSummary <- pS
  ArchRProj
}

##########################################################################################
# Annotation Methods
##########################################################################################

#' Get Embedding from ArchRProject
#' 
#' This function gets an embedding from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param name name of annotations
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

#' Get Annotation Positions from ArchRProject
#' 
#' This function gets annotation positions from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param name name of annotations
#' @param annoName name to subset with annotations
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

#' Get Annotation Matches from ArchRProject
#' 
#' This function gets annotation matches from an ArchRProject
#' 
#' @param ArchRProj ArchRProject
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

#' Add Motif Annotations to ArchRProject
#' 
#' This function adds motif postions and matches to an ArchRProject
#' 
#' @param ArchRProj ArchRProject
#' @param motifSet motifSet JASPAR : JASPAR2016, JASPAR2018; chromVARmotifs : human, mouse, encode, homer 
#' @param name of annotations to store as in ArchRProject
#' @param species species relevant to dataset (default will guess based on getGenome)
#' @param collection JASPAR collection (default = CORE)
#' @param cutOff pvalue cutoff for motif search (see motimatchr)
#' @param w width to consider for motif (see motimatchr)
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

#' Return Available Features for a given Matrix in ArrowFiles within an ArchRProject
#' 
#' This function will identify available features for a matrix and return them for downstream
#' plotting utils.
#' 
#' @param ArchRProj ArchRProject
#' @param useMatrix Matrix Name as in Arrow Files (ie TileMatrix, GeneScoreMatrix, ...)
#' @param select select a specific name with grep
#' @param ignore.case ignore case when searching with select
#' @param ... additional args
#' @export
getFeatures <- function(ArchRProj, useMatrix = "GeneScoreMatrix", select = NULL, ignore.case = TRUE, ...){
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

#' Plot PDF in outputDirectory of ArchRProject
#' 
#' This function will plot PDF in output directory of an ArchRProject 
#' 
#' @param name name of PDF file
#' @param width width of PDF in inches
#' @param height height of PDF in inches
#' @param ArchRProj ArchRProject
#' @param addDOC add date of creation to end of plot file name
#' @param useDingbats use dingbats characters for plotting
#' @param ... additional args to pdf
#' @export
plotPDF <- function(..., name = "Plot", width = 6, height = 6, ArchRProj = NULL, addDOC = TRUE, 
  useDingbats = FALSE, plotList = NULL, useSink = TRUE){

  if(useSink){
    tmpFile <- .tempfile()
    sink(tmpFile)
  }

  if(is.null(plotList)){
    plotList <- list(...)
  }else{
    #plotList <- do.call(list, unlist(plotList, recursive=FALSE))
  }
  
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

  pdf(filename, width = width, height = height, useDingbats = useDingbats)
  for(i in seq_along(plotList)){
    
    if(inherits(plotList[[i]], "gg")){
      
      print("plotting ggplot!")

      print(.fixPlotSize(plotList[[i]], plotWidth = width, plotHeight = height, newPage = FALSE))
      if(i != length(plotList)){
        grid::grid.newpage()
      }
    
    }else if(inherits(plotList[[i]], "gtable")){

      print("plotting gtable!")
      
      print(grid::grid.draw(plotList[[i]]))
      if(i != length(plotList)){
        grid::grid.newpage()
      }

    }else if(attr(class(plotList[[i]]),"package") == "ComplexHeatmap"){
      
      print("plotting copmleheatmap!")

      padding <- 45
      draw(plotList[[i]], 
        padding = unit(c(padding, padding, padding, padding), "mm"), 
        heatmap_legend_side = "bot", 
        annotation_legend_side = "bot"
      )

    }else{

      print("plotting with print")
     
      print(plotList[[i]])

    }

  }
  dev.off()

  if(useSink){
    sink()
    file.remove(tmpFile)
  }

}






##########################################################################################
# Annotation Methods
##########################################################################################

#' Get peakAnnotation from an ArchRProject
#' 
#' This function gets a peakAnnotation from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the peakAnnotation object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @param ... additional args
#' @export
getPeakAnnotation <- function(ArchRProj, name = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@peakAnnotation)){
      stop("Name is not in peakAnnotation!")
    }
  }
  ArchRProj@peakAnnotation[[name]]
}

#' Get peakAnnotation positions from an ArchRProject
#' 
#' This function gets the peakAnnotation positions (i.e. Motifs) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the peakAnnotation object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @param annoName The name of a specific annotation to subset within the peakAnnotation.
#' @param ... additional args
#' @export
getPositions <- function(ArchRProj, name = NULL, annoName = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@peakAnnotation)){
      stop("Name is not in peakAnnotation!")
    }
  }
  anno <- ArchRProj@peakAnnotation[[name]]
  idx <- grep("positions", names(anno), ignore.case=TRUE)
  if(length(idx)==0){
    stop("peakAnnotation does not contain positions!")
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

#' Get peakAnnotation matches from an ArchRProject
#' 
#' This function gets peakAnnotation matches from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the annotation object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @param annoName The name of a specific annotation to subset within the peakAnnotation.
#' @param ... additional args
#' @export
getMatches <- function(ArchRProj, name = NULL, annoName = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@peakAnnotation)){
      stop("Name is not in peakAnnotation!")
    }
  }
  anno <- ArchRProj@peakAnnotation[[name]]
  idx <- grep("matches", names(anno), ignore.case=TRUE)
  if(length(idx)==0){
    stop("peakAnnotation does not contain positions!")
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
#' @param ArchRProj An `ArchRProject` object.
#' @param motifSet The motif set to be used for annotation. Options include: (i) "JASPAR2016", "JASPAR2018", "JASPAR2020" which gives the 2016, 2018 or 2020 version of JASPAR motifs or (ii) one of "cisbp", "encode", or "homer" which gives the corresponding motif sets from the chromVAR package. 
#' @param name The name of peakAnnotations to be stored as in `ArchRProject`
#' @param species The name of the species relevant to the supplied `ArchRProject`. This is used for identifying which motif to be used from CisBP/JASPAR. By default, this function will attempt to guess the species based on the value from `getGenome()`.
#' @param collection If one of the JASPAR motif sets is used via `motifSet`, this parameter allows you to indicate the JASPAR collection to be used. Possible options include "CORE", etc.
#' @param cutOff The  p-value cutoff to be used for motif search (see the `motimatchr` package for more information).
#' @param w The width in basepairs to consider for motif matches (see the `motimatchr` package for more information).
#' @param ... additional args
#' @export
addMotifAnnotations <- function(
  ArchRProj = NULL,
  motifSet = "cisbp",
  name = "Motif",
  species = NULL,
  collection = "CORE",
  cutOff = 5e-05, 
  w = 7,
  ...
  ){

  .requirePackage("motifmatchr", installInfo='BiocManager::install("motifmatchr")')
  ArchRProj <- .validArchRProject(ArchRProj)

  if(grepl("JASPAR|CISBP", motifSet, ignore.case = TRUE) & is.null(species)){
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

  if(tolower(motifSet)=="jaspar2020"){
    
    .requirePackage("JASPAR2020",installInfo='BiocManager::install("JASPAR2020")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
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

  }else if(tolower(motifSet)=="jaspar2016"){

    .requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="cisbp"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    if(tolower(species) == "mus musculus"){
      data("mouse_pwms_v2")
      motifs <- mouse_pwms_v2
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else if(tolower(species) == "homo sapiens"){
      data("human_pwms_v2")
      motifs <- human_pwms_v2
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }

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

  ArchRProj@peakAnnotation[[name]]$Name <- name
  ArchRProj@peakAnnotation[[name]]$motifs <- motifs
  ArchRProj@peakAnnotation[[name]]$motifSummary <- motifSummary
  ArchRProj@peakAnnotation[[name]]$Positions <- savePositions
  ArchRProj@peakAnnotation[[name]]$Matches <- saveMatches

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



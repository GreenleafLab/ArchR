
##########################################################################################
# Annotation Methods
##########################################################################################

#' Get peakAnnotation from an ArchRProject
#' 
#' This function gets a peakAnnotation from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the `peakAnnotation` object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @export
getPeakAnnotation <- function(ArchRProj = NULL, name = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  if(is.null(name)){
    name <- 1
  }else{
    if(name %ni% names(ArchRProj@peakAnnotation)){
      stop("Name is not in peakAnnotation!")
    }
  }
  ArchRProj@peakAnnotation[[name]]
}

#' Get peak annotation positions from an ArchRProject
#' 
#' This function gets the peak annotation positions (i.e. Motifs) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the `peakAnnotation` object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @param annoName The name of a specific annotation to subset within the `peakAnnotation`.
#' @export
getPositions <- function(ArchRProj = NULL, name = NULL, annoName = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  .validInput(input = annoName, name = "annoName", valid = c("character", "null"))
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

#' Get peak annotation matches from an ArchRProject
#' 
#' This function gets peak annotation matches from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the `peakAnnotation` object (i.e. Motifs) to retrieve from the designated `ArchRProject`.
#' @param annoName The name of a specific annotation to subset within the `peakAnnotation`.
#' @export
getMatches <- function(ArchRProj = NULL, name = NULL, annoName = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character", "null"))
  .validInput(input = annoName, name = "annoName", valid = c("character", "null"))
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

#' Add peak annotations to an ArchRProject
#' 
#' This function adds information about which peaks contain input regions to a given ArchRProject. For each peak, a binary value is stored indicating whether each region is observed within the peak region.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param regions A `list` of `GRanges` that are to be overlapped with the `peakSet` in the `ArchRProject`.
#' @param name The name of `peakAnnotation` object to be stored as in `ArchRProject`.
#' @param force A boolean value indicating whether to force the `peakAnnotation` object indicated by `name` to be overwritten if it already exist in the given `ArchRProject`.
#' @export
addPeakAnnotations <- function(
  ArchRProj = NULL,
  regions = NULL,
  name = "Region",
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = regions, name = "regions", valid = c("grangeslist", "list"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  tstart <- Sys.time()

  if(name %in% names(ArchRProj@peakAnnotation)){
    if(force){
      message("peakAnnotation name already exists! Overriding.")
    }else{
      stop("peakAnnotation name already exists! set force = TRUE to override!")
    }
  }

  if(inherits(regions, "GRanges")){

    regionPositions <- GRangesList(region = regions)

  }else{

    regionPositions <- lapply(seq_along(regions), function(x){
      
      if(inherits(regions[x], "GRanges")){

          gr <- .validGRanges(regions[x])

      }else if(is.character(regions[x])){
        
        gr <- tryCatch({
          .validGRanges(makeGRangesFromDataFrame(
            df = data.frame(data.table::fread(regions[x])), 
            keep.extra.columns = TRUE,
            seqnames.field = "V1",
            start.field = "V2",
            end.field = "V3"
          ))
        }, error = function(y){

          print(paste0("Could not successfully get region : ", regions[x]))

        })

      }else{
        
        stop("Unrecognized input in regions please input GRanges, GRangesList, or Paths to bed files!")
      
      }

      gr

    }) %>% GRangesList

    names(regionPositions) <- names(regions)

  }

  #############################################################
  # Peak Overlap Matrix
  #############################################################
  peakSet <- getPeakSet(ArchRProj)
  allPositions <- unlist(regionPositions)

  .messageDiffTime("Creating Peak Overlap Matrix", tstart)
  overlapRegions <- findOverlaps(peakSet, allPositions, ignore.strand=TRUE)
  regionMat <- Matrix::sparseMatrix(
    i = queryHits(overlapRegions),
    j = match(names(allPositions),names(regionPositions))[subjectHits(overlapRegions)],
    x = rep(TRUE, length(overlapRegions)),
    dims = c(length(peakSet), length(regionPositions))
  )
  colnames(regionMat) <- names(regionPositions)
  regionMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = regionMat), rowRanges = peakSet)

  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings=FALSE)
  savePositions <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Positions-In-Peaks.rds"))
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Matches-In-Peaks.rds"))

  out <- SimpleList(
      regionMatches = regionMat,
      regionPositions = regionPositions,
      date = Sys.Date()
    )

  ArchRProj@peakAnnotation[[name]]$Name <- name
  ArchRProj@peakAnnotation[[name]]$Positions <- savePositions
  ArchRProj@peakAnnotation[[name]]$Matches <- saveMatches

  saveRDS(out, file.path(getOutputDirectory(ArchRProj),  "Annotations", paste0(name,"-In-Peaks-Summary.rds")), compress = FALSE)
  saveRDS(out$regionPositions, savePositions, compress = FALSE)
  saveRDS(out$regionMatches, saveMatches, compress = FALSE)

  return(ArchRProj)

}

#' Add motif annotations to an ArchRProject
#' 
#' This function adds information about which peaks contain motifs to a given ArchRProject. For each peak, a binary value is stored indicating whether each motif is observed within the peak region.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param motifSet The motif set to be used for annotation. Options include: (i) "JASPAR2016", "JASPAR2018", "JASPAR2020" which gives the 2016, 2018 or 2020 version of JASPAR motifs or (ii) one of "cisbp", "encode", or "homer" which gives the corresponding motif sets from the `chromVAR` package. 
#' @param name The name of the `peakAnnotation` object to be stored in the provided `ArchRProject`
#' @param species The name of the species relevant to the supplied `ArchRProject`. This is used for identifying which motif to be used from CisBP/JASPAR. By default, this function will attempt to guess the species based on the value from `getGenome()`.
#' @param collection If one of the JASPAR motif sets is used via `motifSet`, this parameter allows you to indicate the JASPAR collection to be used. See `getMatrixSet()` from `TFBSTools` for all options to supply for collection.
#' @param cutOff The p-value cutoff to be used for motif search. The p-value is determined vs a background set of sequences (see `MOODS` for more details on this determination).
#' @param width The width in basepairs to consider for motif matches. See the `motimatchr` package for more information.
#' @param force A boolean value indicating whether to force the `peakAnnotation` object indicated by `name` to be overwritten if it already exist in the given `ArchRProject`.
#' @param ... Additional parameters to be passed to `TFBSTools::getMatrixSet` for getting a PWM object.
#' @export
addMotifAnnotations <- function(
  ArchRProj = NULL,
  motifSet = "cisbp",
  name = "Motif",
  species = NULL,
  collection = "CORE",
  cutOff = 5e-05, 
  width = 7,
  force = FALSE,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = motifSet, name = "motifSet", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = species, name = "species", valid = c("character", "null"))
  .validInput(input = collection, name = "collection", valid = c("character", "null"))
  .validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
  .validInput(input = w, name = "w", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  .requirePackage("motifmatchr", installInfo='BiocManager::install("motifmatchr")')

  if(name %in% names(ArchRProj@peakAnnotation)){
    if(force){
      message("peakAnnotation name already exists! Overriding.")
    }else{
      stop("peakAnnotation name already exists! set force = TRUE to override!")
    }
  }

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
  BSgenome <- validBSgenome(BSgenome)

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
      w = width
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

.summarizeJASPARMotifs <- function(motifs = NULL){

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

.summarizeChromVARMotifs <- function(motifs = NULL){

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

#' Add ArchR annotations to an ArchRProject
#' 
#' This function adds information about which peaks in the ArchR database contain input regions to a given ArchRProject. For each peak, a binary value is stored indicating whether each region is observed within the peak region.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param db A character indicating which database or a path to a database to use for peak annotation. Options include ArchR, LOLA, and a valid path to a file of class `ArchRAnno`.
#' @param collection A character indicating which collection within the database to collect for annotation. 
#' For ArchR, options JJJ include to be added.
#' For LOLA, options include "EncodeTFBS" "CistromeTFBS", "CistromeEpigenome", "Codex", or "SheffieldDnase".
#' If supplying a custom ArchRAnno please use a valid collection.
#' @param name The name of `peakAnnotation` object to be stored as in `ArchRProject`.
#' @param force A boolean value indicating whether to force the `peakAnnotation` object indicated by `name` to be overwritten if it already exist in the given `ArchRProject`.
#' @export
addArchRAnnotations <- function(
  ArchRProj = NULL,
  db = "LOLA",
  collection = "EncodeTFBS",
  name = collection,
  force = FALSE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = db, name = "db", valid = c("character"))
  .validInput(input = collection, name = "collection", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = force, name = "force", valid = c("boolean"))

  tstart <- Sys.time()

  if(name %in% names(ArchRProj@peakAnnotation)){
    if(force){
      message("peakAnnotation name already exists! Overriding.")
    }else{
      stop("peakAnnotation name already exists! set force = TRUE to override!")
    }
  }


  annoPath <- file.path(find.package("ArchR", NULL, quiet = TRUE), "data", "Annotations")
  dir.create(annoPath, showWarnings = FALSE)
  
  if(tolower(db) == "lola"){
  
    if(genome == "hg19"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/LOLA-Hg19-v1.Anno"
    }else if(genome == "hg38"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/LOLA-Hg38-v1.Anno"
    }else if(genome == "mm9"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/LOLA-Mm9-v1.Anno"
    }else if(genome == "mm10"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/LOLA-Mm10-v1.Anno"
    }
    
    #Download
    if(!file.exists(file.path(annoPath, basename(url)))){
      message("Annotation ", basename(url)," does not exist! Downloading...")
      download.file(
        url = url, 
        destfile = file.path(annoPath, basename(url)),
        quiet = FALSE
      )
    }
    AnnoFile <- file.path(annoPath, basename(url))

  }else if(tolower(db) == "archr"){

    if(genome == "hg19"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Hg19-v1.Anno"
    }else if(genome == "hg38"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Hg38-v1.Anno"
    }else if(genome == "mm9"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Mm9-v1.Anno"
    }else if(genome == "mm10"){
      url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Mm10-v1.Anno"
    }

    #Download
    if(!file.exists(file.path(annoPath, basename(url)))){
      message("Annotation ", basename(url)," does not exist! Downloading...")
      download.file(
        url = url, 
        destfile = file.path(annoPath, basename(url)),
        quiet = FALSE
      )
    }
    AnnoFile <- file.path(annoPath, basename(url))

  }else if(tolower(db) %ni% c("archr", "lola")){

    if(!file.exists(db)){
      stop("Database path does not exist! Please supply valid database for Annotations!")
    }
    AnnoFile <- db

  }else{

    stop("Please supply valid database for Annotations!")

  }

  #Check AnnoFile is ArchRAnno
  if (h5read(AnnoFile, "Class") != "ArchRAnno") {
    stop("Not Valid ArchRAnno!")
  }

  #Check if Collection is Valid
  h5d <- h5ls(AnnoFile)
  collections <- h5d[h5d$group == "/" & h5d$otype == "H5I_GROUP",]$name
  if(any(tolower(collections) == tolower(collection))){
    collection <- collections[tolower(collections) %in% tolower(collection)]
  }else{
    stop("Not Valid Collection in ArchRAnno Database!")
  }

  #############################################################
  # Peak Overlap Matrix
  #############################################################
  peakSet <- getPeakSet(ArchRProj)
  chr <- paste0(unique(seqnames(peakSet)))

  message("Annotating Chr: ", appendLF = FALSE)
  regionMat <- lapply(seq_along(chr), function(i){
    message(gsub("chr", " ", chr[i]), appendLF = FALSE)
    regions <- .getRegionsFromAnno(AnnoFile = AnnoFile, Group = collection, chr = chr[i])
    o <- DataFrame(findOverlaps(query = peakSet, subject = regions, ignore.strand = TRUE))
    o$ID <- as.integer(mcols(regions)[o$subjectHits, "ID"])
    o$subjectHits <- NULL
    o
  }) %>% Reduce("rbind", .)
  message("\n")

  regionMat <- Matrix::sparseMatrix(
    i = regionMat[, 1], 
    j = regionMat[, 2], 
    x = rep(TRUE, nrow(regionMat)), 
    dims = c(length(peakSet), length(unique(regionMat$ID)))
  )

  #Get Region Metadata
  regionMetadata <- DataFrame(h5read(AnnoFile, paste0(collection,"/Info")))
  rownames(regionMetadata) <- regionMetadata$name
  colnames(regionMat) <- rownames(regionMetadata)
  regionMetadata$name <- NULL

  #Matches Summarized Experiment
  regionMat <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(matches = regionMat), 
    rowRanges = peakSet, 
    colData = regionMetadata
  )

  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings=FALSE)
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Matches-In-Peaks.rds"))

  out <- SimpleList(
      regionMatches = regionMat,
      regionPositions = "None",
      date = Sys.Date()
    )

  ArchRProj@peakAnnotation[[name]]$Name <- name
  ArchRProj@peakAnnotation[[name]]$Positions <- "None"
  ArchRProj@peakAnnotation[[name]]$Matches <- saveMatches

  saveRDS(out, file.path(getOutputDirectory(ArchRProj),  "Annotations", paste0(name,"-In-Peaks-Summary.rds")), compress = FALSE)
  saveRDS(out$regionMatches, saveMatches, compress = FALSE)

  return(ArchRProj)

}

.getRegionsFromAnno <- function(
  AnnoFile = NULL, 
  Group = NULL,
  chr = NULL, 
  out = "GRanges", 
  method = "fast"
  ){

  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }

  o <- h5closeAll()
  if (h5read(AnnoFile, "Class") != "ArchRAnno") {
  stop("Not Valid ArchRAnno!")
  }
  
  if(chr %ni% .availableSeqnames(AnnoFile, Group)){
    stop("Error Chromosome not in AnnoFile!")
  }

  o <- h5closeAll()
  nRegions <- h5ls(AnnoFile, recursive = TRUE) %>% 
    {.[.$group==paste0("/",Group,"/",chr) & .$name == "Ranges",]$dim} %>% 
    {gsub(" x 2","",.)} %>% as.integer

  if(nRegions==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), ID = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$ID <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }


  if(tolower(method) == "fast"){
    
    output <- .h5read(AnnoFile, paste0(Group,"/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$ID <- Rle(
      values =  .h5read(AnnoFile, paste0(Group,"/",chr,"/IDValues"), method = method), 
      lengths = .h5read(AnnoFile, paste0(Group,"/",chr,"/IDLengths"), method = method)
    )

  }else{
    o <- h5closeAll()
  idRle <- Rle(h5read(AnnoFile, paste0(Group,"/",chr,"/IDValues")), h5read(AnnoFile, paste0(Group,"/",chr,"/IDLengths")))
  idx <- seq_along(idRle)
  if(length(idx) > 0){
    output <- h5read(AnnoFile, paste0(Group,"/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$ID <- idRle[idx]
  }else{
    output <- IRanges(start = 1, end = 1)
    mcols(output)$ID <- c("tmp")
    output <- output[-1,]
  }

  }
  
  o <- h5closeAll()

  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), ID = mcols(output)$ID)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$ID <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), ID = mcols(output)$ID)
      output <- output[-1,]
    }
  }

  return(output)

}

#' Peak Annotation Hypergeometric Enrichment in Marker Peaks.
#' 
#' This function will perform hypergeometric enrichment of a given peak annotation within the defined marker peaks.
#' 
#' @param seMarker  A `SummarizedExperiment` object returned by `markerFeatures()`.
#' @param ArchRProj An `ArchRProject` object.
#' @param peakAnnotation A `peakAnnotation` object in the provided `ArchRProject` to be used for hypergeometric testing.
#' @param matches A custom `peakAnnotation` matches object used as input for the hypergeometric test. See `motifmatchr::matchmotifs()` for additional information.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` to use. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param background A string that indicates whether to use a background set of matched peaks to compare against ("bgdPeaks") or all peaks ("all").
#' @export
peakAnnoEnrichment <- function(
  seMarker = NULL,
  ArchRProj = NULL,
  peakAnnotation = NULL,
  matches = NULL,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  background = "all"
  ){

  .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = peakAnnotation, name = "peakAnnotation", valid = c("character", "null"))
  .validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", "null"))
  .validInput(input = cutOff, name = "cutOff", valid = c("character"))
  .validInput(input = background, name = "background", valid = c("character"))

  tstart <- Sys.time()
  if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
    stop("Only markers identified from PeakMatrix can be used!")
  }

  if(is.null(matches)){
    matches <- getMatches(ArchRProj, peakAnnotation)
  }
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  mcols(r1) <- NULL

  r2 <- getPeakSet(ArchRProj)
  mcols(r2) <- NULL

  if(length(which(paste0(seqnames(r1),start(r1),end(r1), sep = "_") %ni% paste0(seqnames(r2),start(r2),end(r2),sep="_"))) != 0){
    stop("Peaks from matches do not match peakSet in ArchRProj!")
  }

  r3 <- GRanges(rowData(seMarker)$seqnames,IRanges(rowData(seMarker)$start, rowData(seMarker)$end))
  mcols(r3) <- NULL
  rownames(matches) <- paste0(seqnames(matches),start(matches),end(matches),sep="_")
  matches <- matches[paste0(seqnames(r3),start(r3),end(r3), sep = "_"), ]

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(tolower(background) %in% c("backgroundpeaks", "bgdpeaks", "background", "bgd")){
    method <- "bgd"
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  }else{
    method <- "all"
  }

  enrichList <- lapply(seq_len(ncol(seMarker)), function(x){
    .messageDiffTime(sprintf("Computing Enrichments %s of %s",x,ncol(seMarker)),tstart)
    idx <- which(passMat[, x])
    if(method == "bgd"){
      .computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx,])))
    }else{
      .computeEnrichment(matches, idx, seq_len(nrow(matches)))
    }
  }) %>% SimpleList
  names(enrichList) <- colnames(seMarker)

  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
      enrichList[[y]][colnames(matches),x,drop=FALSE]
    }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- colnames(enrichList[[1]])
  assays <- rev(assays)
  out <- SummarizedExperiment::SummarizedExperiment(assays=assays)

  out

}

.computeEnrichment <- function(matches = NULL, compare = NULL, background = NULL){

  matches <- .getAssay(matches,  grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
     pOut$BackgroundFrequency[x], # Number of all successes in background
     pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
     pOut$nCompare[x], # Number that were drawn
     lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)

  #Minus Log10 FDR
  pOut$mlog10FDR <- -log10(p.adjust(matrixStats::rowMaxs(cbind(10^-pOut$mlog10p, 4.940656e-324)), method = "fdr"))
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

  pOut

}

#' Heatmap of Peak Annotation Hypergeometric Enrichment in Marker Peaks.
#' 
#' This function will plot a heatmap of hypergeometric enrichment of peakAnnotation within the defined marker peaks.
#' 
#' @param seMarker  A `SummarizedExperiment` object returned by `peakAnnoEnrichment()`.
#' @param pal A custom continuous palette (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param n The number of top enriched peakAnnotations per column from the `seMarker` to display in the heatmap. This number can be lowered to improve visibility of the heatmap.
#' @param clusterCols A boolean indicating whether or not to cluster columns in the heatmap.
#' @param clusterRows A boolean indicating whether or not to cluster rows in the heatmap.
#' @param labelRows A boolean indicating whether or not to label all rows in the heatmap.
#' @export
enrichHeatmap <- function(
  seEnrich = NULL,
  pal = paletteContinuous(set = "white_blue_purple", n = 100),
  limits = c(0, 60),
  n = 10,
  clusterCols = TRUE,
  clusterRows = TRUE,
  labelRows = TRUE
  ){

  .validInput(input = seEnrich, name = "seEnrich", valid = c("SummarizedExperiment"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  .validInput(input = limits, name = "limits", valid = c("numeric"))
  .validInput(input = n, name = "n", valid = c("integer"))
  .validInput(input = clusterCols, name = "clusterCols", valid = c("boolean"))
  .validInput(input = clusterRows, name = "clusterRows", valid = c("boolean"))
  .validInput(input = labelRows, name = "labelRows", valid = c("boolean"))

  mat <- assays(seEnrich)[["mlog10FDR"]]
  mat[mat < min(limits)] <- min(limits)
  mat[mat > max(limits)] <- max(max(limits), 25)

  keep <- lapply(seq_len(ncol(mat)), function(x){
    rownames(mat)[head(order(mat[, x], decreasing = TRUE), n)]
  }) %>% unlist %>% unique

  ht <- .ArchRHeatmap(
    mat = as.matrix(mat[keep, ,drop = FALSE]),
    scale = FALSE,
    limits = c(min(mat), max(mat)),
    color = pal, 
    clusterCols = clusterCols, 
    clusterRows = clusterRows,
    labelRows = labelRows,
    labelCols = TRUE,
    showColDendrogram = TRUE,
    draw = FALSE,
    name = "Enrichment -log10(FDR)"
  )

  return(ht)

}


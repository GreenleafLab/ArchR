##########################################################################################
# Parallel Information
##########################################################################################

#' Add global number of threads for default parallel computing.
#' 
#' This function will set the global number of threads to be used for ArchR functions.
#' 
#' @param threads default number of threads to be used for parallel execution in ArchR functions by default.
#' @export
addArchRThreads <- function(threads = floor(parallel::detectCores()/ 2)){
  if(tolower(.Platform$OS.type) == "windows"){
    message("Detected windows OS, setting threads to 1.")
    threads <- 1
  }
  message("Setting default number of Parallel threads to ", threads, ".")
  assign("ArchRThreads", as.integer(threads), envir = .GlobalEnv)
}

#' Get global number of threads for default parallel computing.
#' 
#' This function will get the global number of threads to be used for ArchR functions.
#' 
#' @export
getArchRThreads <- function(){
  if(exists("ArchRThreads")){
    if(!is.integer(ArchRThreads)){
      1
    }else{
      ArchRThreads
    }
  }else{
    1
  }
}

##########################################################################################
# Create Gene/Genome Annotation
##########################################################################################

#' Create Genome Annotation for ArchR
#' 
#' This function will create a genome annotation that can be used for createArrowFiles, ArchRProject, etc.
#' 
#' @param genome A string that points to a BSgenome or a BSgenome object (ie hg38, BSgenome.Hsapiens.UCSC.hg38).
#' @param chromSizes A GRanges of chromosome start and end coordinates.
#' @param blacklist A GRanges of regions that should be excluded from analyses due to unwanted biases.
#' @param filter A boolean indicating whether non-normal chromosome scaffolds should be excluded.
#' @export
createGenomeAnnotation <- function(
  genome = NULL,
  chromSizes = NULL,
  blacklist = NULL,
  filter = TRUE
  ){

  if(is.null(genome) | is.null(blacklist) | is.null(chromSizes)){

    ##################
    message("Getting genome...")
    bsg <- .validBSgenome(genome)
    genome <- bsg@pkgname

    ##################
    message("Getting chromSizes...")
    chromSizes <- GRanges(names(seqlengths(bsg)), IRanges(1, seqlengths(bsg)))
    if(filter){
        chromSizes <- filterChrGR(chromSizes)
    }
    seqlengths(chromSizes) <- end(chromSizes)

    ##################
    message("Getting blacklist...")
    blacklist <- .getBlacklist(genome = bsg@provider_version)

  }else{

    bsg <- .validBSgenome(genome)
    genome <- bsg@pkgname
    
    chromSizes <- .validGRanges(chromSizes)
    
    blacklist <- .validGRanges(blacklist)

  }

  SimpleList(genome = genome, chromSizes = chromSizes, blacklist = blacklist)

}

#' Create Gene Annotation for ArchR
#' 
#' This function will create a gene annotation that can be used for createArrowFiles, ArchRProject, etc.
#' 
#' @param genome A string that specifies the genome (ie hg38, hg19, mm10, mm9).
#' @param TxDb A transcript database from Bioconductor which contains information for gene/transcript coordinates.
#' @param OrgDb An organism database from Bioconductor which contains information for gene/transcript symbols from ids.
#' @param genes A GRanges of gene coordinates (start to end). Needs to have a symbols column matching the exons symbols column.
#' @param exons A GRanges of gene exon coordinates. Needs to have a symbols column matching the genes symbols column
#' @param TSS A GRanges of transcription start sites (stranded) for computing TSS enrichment scores downstream.
#' @export
createGeneAnnnotation <- function(
  genome = NULL,
  TxDb = NULL,
  OrgDb = NULL,
  genes = NULL,
  exons = NULL,
  TSS = NULL
  ){

  if(is.null(genes) | is.null(exons) | is.null(TSS)){

    inGenes <- genes
    inExons <- exons
    inTSS <- TSS

    .requirePackage("GenomicFeatures")

    if(is.null(genome)) {
      if (is.null(TxDb) | is.null(OrgDb)) {
          stop("If no provided genome then you need TxDb and OrgDb!")
      }
    }

    if(!is.null(genome)){
      TxDb <- .getTxDb(genome)
      OrgDb <- .getOrgDb(genome)
    }

    ###########################
    message("Getting Genes...")
    genes <- GenomicFeatures::genes(TxDb)
    mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, 
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
    names(genes) <- NULL
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ###########################
    message("Getting Exons...")
    exons <- unlist(GenomicFeatures::exonsBy(TxDb, by = "tx"))
    exons$tx_id <- names(exons)
    mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(TxDb, keys = paste0(mcols(exons)$tx_id), 
        column = "GENEID", keytype = "TXID")[, "GENEID"])
    exons <- exons[!is.na(mcols(exons)$gene_id), ]
    mcols(exons)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(exons)$gene_id, 
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
    names(exons) <- NULL
    mcols(exons)$exon_id <- NULL
    mcols(exons)$exon_name <- NULL
    mcols(exons)$exon_rank <- NULL
    mcols(exons)$tx_id <- NULL
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

    ###########################
    message("Getting TSS...")
    TSS <- unique(resize(GenomicFeatures::transcripts(TxDb), width = 1, fix = "start"))

    if(!is.null(inGenes)){
      genes <- .validGRanges(inGenes)
    }

    if(!is.null(inExons)){
      exons <- .validGRanges(inExons)
    }

    if(!is.null(inTSS)){
      TSS <- .validGRanges(inTSS)
    }

  }else{

    genes <- .validGRanges(genes)
    exons <- .validGRanges(exons)
    TSS <- unique(.validGRanges(TSS))

  }

  SimpleList(genes = genes, exons = exons, TSS = TSS)

}

.getBlacklist <- function(genome){
  
  encodeBL <- c(
    "hg19" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz",
    "hg38" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
    "mm10" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
    "mm9" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
    "ce10" = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz",
    "dm3" = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz"
  )

  if(tolower(genome) %in% names(encodeBL)){
    bl <- tryCatch({
      blacklist <- import.bed(encodeBL[tolower(genome)])
    }, error = function(x){
      message("Blacklist not downloaded! Continuing without, be careful for downstream biases...")
      GRanges()
    })
  }else{
    message("Blacklist not downloaded! Continuing without, be careful for downstream biases...")
    bl <- GRanges()
  }

  bl

}

.getTxDb <- function(genome, filter = TRUE, install = TRUE){

  if(toupper(genome) == "HG19"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg19.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg19.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(toupper(genome) == "HG38"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg38.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg38.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(toupper(genome) == "MM9"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm9.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm9.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  }else if(toupper(genome) == "MM10"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm10.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm10.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(toupper(genome) == "SACCER3"){
    if(suppressWarnings(!require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", update=FALSE)
      }else{
        stop("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene is not installed!")
      }
    }
    library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
    txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
  }else if(toupper(genome) == "RHEMAC8"){
    if(suppressWarnings(!require(TxDb.Mmulatta.UCSC.rheMac8.refGene))){
      if(install){
        message("Package does not exist, now trying bioconductor...")
        BiocManager::install("TxDb.Mmulatta.UCSC.rheMac8.refGene", update=FALSE)
      }else{
        stop("TxDb.Mmulatta.UCSC.rheMac8.refGene is not installed!")
      }
    }
    library(TxDb.Mmulatta.UCSC.rheMac8.refGene)
    txdb <- TxDb.Mmulatta.UCSC.rheMac8.refGene
  }else{
    stop("Genome not recognized!")
  }

  if(filter){
    txdb <- filterChrGR(txdb)
  }

  return(txdb)

}

.getOrgDb <- function(genome){

  if(toupper(genome) == "HG19" | toupper(genome) == "HG38"){
    if(suppressWarnings(!require(org.Hs.eg.db))){
      message("Package does not exist, now trying bioconductor...")
      BiocManager::install("org.Hs.eg.db", update=FALSE)
    }
    library(org.Hs.eg.db)
    annodb <- org.Hs.eg.db
  }else if(toupper(genome) == "MM9" | toupper(genome) == "MM10"){
    if(suppressWarnings(!require(org.Mm.eg.db))){
      message("Package does not exist, now trying bioconductor...")
      BiocManager::install("org.Mm.eg.db", update=FALSE)
    }
    library(org.Mm.eg.db)
    annodb <- org.Mm.eg.db
  }else{
    stop("Genome not recognized!")
  }
  return(annodb)

}

.validGeneAnnotation <- function(geneAnnotation, ...){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){

    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])

  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }

  gA

}

.validGenomeAnnotation <- function(genomeAnnotation, ...){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){

    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    bsg <- .validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
    gA$genome <- bsg@pkgname
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])

  }else{

    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
  
  }

  gA

}

##########################################################################################
# Output Directory
##########################################################################################

#' Get outputDirectory from an ArchRProject
#' 
#' This function gets the outputDirectory from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param ... additional args
#' @export
getSampleNames <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  snames <- rownames(ArchRProj@sampleColData)
  return(snames)
}

#' Get number of cells from ArchRProject/ArrowFile
#' 
#' This function gets number of cells in ArchRProject/ArrowFile
#' 
#' @param input An `ArchRProject` object or ArrowFile.
#' @param ... additional args
#' @export
nCells <- function(input, ...){
  if(inherits(input, "ArchRProject")){
    nrow(getCellColData(input))
  }else if(inherits(input, "character")){
    if(file.exists(input)){
      length(.availableCells(input))
    }else{
      stop("File does not exist!")
    }
  }else{
    stop("Provide ArchRProj or ArrowFiles")
  }
}

#' Get sampleColData from an ArchRProject
#' 
#' This function gets the `sampleColData` from a given `ArchRProject`.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param select A character vector containing the column names to select from sampleColData.
#' @param drop A boolean value that indicates whether to drop the `dataframe` structure and convert to a vector if selecting only one column.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param data The data to add to `sampleColData`.
#' @param name The column header name to be used for this new data in `sampleColData`. If a column with this name already exists, you may set `force` equal to TRUE to overwrite the data in this column.
#' @param samples The names of the samples corresponding to `data`. Typically new data is added to all samples but you may use this argument to only add data to a subset of samples. Samples where `data` is not added are set to `NA`.
#' @param force A boolean value that indicates whether or not to overwrite data in a given column when the value passed to `name` already exists as a column name in `sampleColData`.
#' @param ... additional args
#' @export
addSampleColData <- function(ArchRProj, data = NULL, name = NULL, samples = NULL, force = FALSE){
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
      stop(paste0("Error previous entry for ", name, ", Set force = TRUE to override!"))
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
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param select A character vector of column names to select from `cellColData` if you would like to subset the returned data.
#' @param drop A boolean value that indicates whether to drop the `dataframe` structure and convert to a vector if selecting only one column.
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

#' Add information to cellColData in an ArchRProject
#' 
#' This function adds new data to cellColData in a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param data The data to add to `cellColData`.
#' @param name The column header name to be used for this new data in `cellColData`. If a column with this name already exists, you may set `force` equal to TRUE to overwrite the data in this column.
#' @param cells The names of the cells corresponding to `data`. Typically new data is added to all cells but you may use this argument to only add data to a subset of cells. Cells where `data` is not added are set to `NA`.
#' @param force A boolean value indicating whether or not to overwrite data in a given column when the value passed to `name` already exists as a column name in `cellColData`.
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
      stop(paste0("Error previous entry for ", name, ", Set force = TRUE to override!"))
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
#' This function gets the peak set as a GRanges object from an ArchRProject.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param peakSet A `GRanges` object containing the set of regions that define all peaks in the desired peak set.
#' @param force If a `peakSet` object has already been added to the given `ArchRProject`, the value of `force` determines whether or not to overwrite this `peakSet`.
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
#' This function gets the genomeAnnotation (see createGenomeAnnotation) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
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
#' This function gets the geneAnnotation (see createGeneAnnotation) from a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param ... additional args
#' @export
getTSS <- function(ArchRProj, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  return(ArchRProj@geneAnnotation$TSS)
}

#' Get the genes from an ArchRProject
#' 
#' This function gets the genes start to end coordinates as a GRanges from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param symbols A character vector containing the gene symbols to subset from the `geneAnnotation`.
#' @param ... additional args
#' @export
getGenes <- function(ArchRProj, symbols = NULL, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  genes <- ArchRProj@geneAnnotation$genes
  if(!is.null(symbols)){
    genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
  }
  return(genes)
}

#' Get the exons from an ArchRProject
#' 
#' This function gets the exons coordinates as a GRanges from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param symbols A character vector containing the gene symbols for the genes where exons should be extracted.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. IterativeLSI) to retrieve from the designated `ArchRProject`.
#' @param returnMatrix If set to "mat" or "matrix", the function will return the `reducedDims` object as a matrix with entries for each individual cell. Otherwise, it will return the full `reducedDims` object.
#' @param dimsToUse A vector containing the dimensions (i.e. 1:25) to return from the `reducedDims` object.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded.
#' @param ... additional args
#' @export
getReducedDims <- function(
  ArchRProj, 
  reducedDims = "IterativeLSI", 
  returnMatrix = TRUE, 
  dimsToUse = NULL,
  corCutOff = 0.75,
  ...
  ){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(reducedDims %in% names(ArchRProj@reducedDims)){
    corToDepth <- ArchRProj@reducedDims[[reducedDims]]$corToDepth
    if(!is.null(dimsToUse)){
      corToUse <- dimsToUse
    }else{
      corToUse <- seq_along(corToDepth)
    }
    idx <- which(abs(corToDepth[corToUse]) <= corCutOff)
    if(length(idx) != length(corToUse)){
      message("Filtering ", length(corToUse) - length(idx), " dims correlated > ", corCutOff, " to log10(depth + 1)")
    }
    if(returnMatrix){
      out <- ArchRProj@reducedDims[[reducedDims]][[1]][,corToUse[idx],drop=FALSE]
    }else{
      out <- ArchRProj@reducedDims[[reducedDims]]
      out[[1]] <- out[[1]][,corToUse[idx],drop=FALSE]
    }
  }else{
    stop("reducedDims not in computed reduced dims!")
  }
  return(out)
}

#' Get embedding information stored in an ArchRProject
#' 
#' This function gets an embedding (i.e. UMAP) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the `embeddings` object (i.e. UMAP, TSNE see embeddingOut of addEmbeddings) to retrieve from the designated `ArchRProject`.
#' @param returnDF A boolean value indicating whether to return the embedding object as a `data.frame`. Otherwise, it will return the full embedding object.
#' @param ... additional args
#' @export
getEmbedding <- function(ArchRProj, embedding = "UMAP", returnDF = TRUE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  if(embedding %in% names(ArchRProj@embeddings)){
    if(returnDF){
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
#' @param ArchRProj An `ArchRProject` object.
#' @param returnSummary A boolean value indicating whether to return a summary of the `ArchRProject` or to just print the summary.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the summary information to add to the `ArchRProject` object.
#' @param summary A vector to add as summary information to the `ArchRProject` object.
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
# Additional Methods
##########################################################################################

#' Get the features that could be selected from a given data matrix within an ArchRProject
#' 
#' This function will identify available features from a given data matrix  (i.e. "GeneScoreMatrix", or "TileMatrix") and return them for downstream plotting utilities.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix as stored in the ArrowFiles of the `ArchRProject`. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param select A string specifying a specific featureName (or rowname) found with grep
#' @param ignore.case A boolean value indicating whether or not to ignore the case (upper-case / lower-case) when searching via grep for the string passed to `select`.
#' @param ... additional args
#' @export
getFeatures <- function(ArchRProj, useMatrix = "GeneScoreMatrix", select = NULL, ignore.case = TRUE, ...){
  ArchRProj <- .validArchRProject(ArchRProj)
  fdf <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class"))
  if(is.null(select)){
    if(any(duplicated(paste0(fdf$name))) | matrixClass == "Sparse.Assays.Matrix"){
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
#' @param ... vector of plots to be plotted (if input is a list use plotList instead)
#' @param name The file name to be used for the output PDF file.
#' @param width The width in inches to be used for the output PDF file.
#' @param height The height in inches to be used for the output PDF.
#' @param ArchRProj An `ArchRProject` object to be used for getting plotDirectory in outputDirectory.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for preventing overwritting of old plots.
#' @param useDingbats A boolean variable that determines wheter to use dingbats characters for plotting points.
#' @param plotList A `list` of plots to be printed to the output PDF file. Each element of `plotList` should be a printable plot formatted object (ggplot2, plot, heatmap, etc).
#' @param useSink use sink to hide messages from plotting
#' @export
plotPDF <- function(..., name = "Plot", width = 6, 
  height = 6, ArchRProj = NULL, addDOC = TRUE, 
  useDingbats = FALSE, plotList = NULL, useSink = TRUE){

  if(is.null(plotList)){
    plotList <- list(...)
  }else{
    plotList2 <- list()
    for(i in seq_along(plotList)){
      if(inherits(plotList[[i]], "list")){
        for(j in seq_along(plotList[[i]])){
          plotList2[[length(plotList2) + 1]] <- plotList[[i]][[j]]
        }
      }else{
        plotList2[[length(plotList2) + 1]] <- plotList[[i]]
      }
    }
    plotList <- plotList2
    rm(plotList2)
    gc()
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

  if(useSink){
    tmpFile <- .tempfile()
    sink(tmpFile)
  }

  o <- tryCatch({

    pdf(filename, width = width, height = height, useDingbats = useDingbats)
    for(i in seq_along(plotList)){
      
      if(inherits(plotList[[i]], "gg")){
        
        print("plotting ggplot!")

        if(!is.null(attr(plotList[[i]], "ratioYX"))){
          print(.fixPlotSize(plotList[[i]], plotWidth = width, plotHeight = height, height = attr(plotList[[i]], "ratioYX"), newPage = FALSE))
        }else{
          print(.fixPlotSize(plotList[[i]], plotWidth = width, plotHeight = height, newPage = FALSE))
        }

        if(i != length(plotList)){
          grid::grid.newpage()
        }
      
      }else if(inherits(plotList[[i]], "gtable")){
        
        print(grid::grid.draw(plotList[[i]]))
        if(i != length(plotList)){
          grid::grid.newpage()
        }
      }else if(inherits(plotList[[i]], "HeatmapList") | inherits(plotList[[i]], "Heatmap") ){ 
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

  }, error = function(x){

    suppressWarnings(sink())
    message(x)

  })

  return(0)

}

#' Get Tutorial Data For ArchR
#' 
#' This function will download data for a given tutorial and return the input files required for ArchR
#' 
#' @param tutorial The name of the available tutorial for which to retreive the tutorial data. Options are "Hematopoiesis", "PBMC", "FreshFrozen". "Hematopoiesis" refers to hematopoieitic scATAC hierarchy. "PBMC" refers to a small standard PBMC scATAC dataset. "FreshFrozen" refers to a PBMC fresh and frozen scATAC dataset.
#' @param ... additional args
#' @export
getTutorialData <- function(tutorial = "hematopoiesis", ...){
  if(tolower(tutorial) %in% c("heme","hematopoiesis")){
    
    if(!dir.exists("Heme_Fragments")){
      download.file(
        url = "https://jeffgranja.s3.amazonaws.com/ArchR-Tutorial-Data/Heme/Heme_Fragments.zip", 
        destfile = "Heme_Fragments.zip"
      )
      unzip("Heme_Fragments.zip")
      if(dir.exists("Heme_Fragments")){
        file.remove("Heme_Fragments.zip")
      }else{
        stop("Download May Not Have Worked!")
      }
    }
    pathFragments <- "Heme_Fragments"

  }else if(tolower(tutorial) %in% c("pbmc")){

    if(!dir.exists("Pbmc_Fragments")){
      download.file(
        url = "https://jeffgranja.s3.amazonaws.com/ArchR-Tutorial-Data/Heme/Pbmc_Fragments.zip", 
        destfile = "Pbmc_Fragments.zip"
      )
      unzip("Pbmc_Fragments")
      if(dir.exists("Pbmc_Fragments")){
        file.remove("Pbmc_Fragments")
      }else{
        stop("Download May Not Have Worked!")
      }
    }
    pathFragments <- "Pbmc_Fragments"

  }else if(tolower(tutorial) %in% c("fresh_vs_frozen", "freshfrozen", "batch")){

    if(!dir.exists("Fresh_Frozen_Fragments")){
      download.file(
        url = "https://jeffgranja.s3.amazonaws.com/ArchR-Tutorial-Data/Heme/Fresh_Frozen_Fragments.zip", 
        destfile = "Fresh_Frozen_Fragments.zip"
      )
      unzip("Fresh_Frozen_Fragments")
      if(dir.exists("Fresh_Frozen_Fragments")){
        file.remove("Fresh_Frozen_Fragments")
      }else{
        stop("Download May Not Have Worked!")
      }
    }
    pathFragments <- "Fresh_Frozen_Fragments"

  }else{
  
    stop("There is no tutorial data for : ", tutorial)
  
  }

  inputFiles <- list.files(pathFragments, pattern = ".gz", full.names = TRUE)
  names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments, pattern = ".gz"))
  inputFiles <- inputFiles[!grepl(".tbi", inputFiles)]
  inputFiles

}

#' Get Input Files from paths to create arrows
#' 
#' This function will look for fragment files and bam files in the input paths and return the full path and sample names
#' 
#' @param paths A character vector of paths to search for usable input files.
#' @param ... additional args
#' @export
getInputFiles <- function(paths, ...){ 
  v <- lapply(paths, function(x){
    
    #Fragments
    inputFrags <- list.files(x, pattern = ".fragments.tsv.gz", full.names = TRUE)
    names(inputFrags) <- gsub(".fragments.tsv.gz", "", list.files(x, pattern = ".fragments.tsv.gz"))
    inputFrags <- inputFrags[!grepl(".tbi", inputFrags)]
    
    #Bams
    inputBams <- list.files(x, pattern = ".bam", full.names = TRUE)
    names(inputBams) <- gsub(".bam", "", list.files(x, pattern = ".bam"))
    inputBams <- inputBams[!grepl(".bai", inputBams)]
    
    c(inputFrags, inputBams)

  }) %>% unlist

  if(any(duplicated(names(v)))){
    names(v) <- paste0(names(v), "_", seq_along(v))
  }

  v

}

#' Get Valid Barcodes from 10x Cell Ranger output to pre-filter barcodes
#' 
#' This function will read in processed 10x cell ranger files and identify barcodes that are associated with a cell that passed QC.
#' 
#' @param csvFiles A character vector of file names to be read in for identification of valid cell barcodes.
#' @param sampleNames A character vector containing the sample names to be associated with each individual entry in `csvFiles`.
#' @param ... additional args
#' @export
getValidBarcodes <- function(csvFiles, sampleNames, ...){

  if(!all(file.exists(csvFiles))){
    stop("Not All csvFiles exists!")
  }

  barcodeList <- lapply(seq_along(csvFiles), function(x){
    df <- .suppressAll(data.frame(readr::read_csv(csvFiles[x])))
    if("cell_id" %ni% colnames(df)){
      stop("cell_id not in colnames of 10x singlecell.csv file! Are you sure inut is correct?")
    }
    as.character(df[which(paste0(df$cell_id) != "None"),]$barcode)
  }) %>% SimpleList
  names(barcodeList) <- sampleNames

  barcodeList

}


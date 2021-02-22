#' Create a genome annotation object for ArchR
#' 
#' This function will create a genome annotation object that can be used for creating ArrowFiles or an ArchRProject, etc.
#' 
#' @param genome Either (i) a string that is a valid `BSgenome` or (ii) a `BSgenome` object (ie "hg38" or "BSgenome.Hsapiens.UCSC.hg38").
#' @param chromSizes A `GRanges` object containing chromosome start and end coordinates.
#' @param blacklist A `GRanges` object containing regions that should be excluded from analyses due to unwanted biases.
#' @param filter A boolean value indicating whether non-standard chromosome scaffolds should be excluded.
#' These "non-standard" chromosomes are defined by `filterChrGR()`.
#' @param filterChr A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' If no manual removal is desired, `filterChr` should be set to `NULL`.
#' @export
createGenomeAnnotation <- function(
  genome = NULL,
  chromSizes = NULL,
  blacklist = NULL,
  filter = TRUE,
  filterChr = c("chrM")
  ){

  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = chromSizes, name = "chromSizes", valid = c("granges", "null"))
  .validInput(input = blacklist, name = "blacklist", valid = c("granges", "null"))
  .validInput(input = filter, name = "filter", valid = c("boolean"))
  .validInput(input = filterChr, name = "filterChr", valid = c("character", "null"))

  if(is.null(genome) | is.null(blacklist) | is.null(chromSizes)){

    ##################
    message("Getting genome..")
    bsg <- validBSgenome(genome)
    genome <- bsg@pkgname

    ##################
    message("Getting chromSizes..")
    chromSizes <- GRanges(names(seqlengths(bsg)), IRanges(1, seqlengths(bsg)))
    if(filter){
        chromSizes <- filterChrGR(chromSizes, remove = filterChr)
    }
    seqlengths(chromSizes) <- end(chromSizes)

    ##################
    message("Getting blacklist..")

    genomeName <- tryCatch({
      bsg@provider_version
    }, error = function(e){
      strsplit(bsg@pkgname,"\\.")[[1]][4]
    })

    blacklist <- .getBlacklist(genome = genomeName)

  }else{

    bsg <- validBSgenome(genome)
    genome <- bsg@pkgname
    
    chromSizes <- .validGRanges(chromSizes)
    
    blacklist <- .validGRanges(blacklist)

  }

  SimpleList(genome = genome, chromSizes = chromSizes, blacklist = blacklist)

}

#' Create a gene annotation object for ArchR
#' 
#' This function will create a gene annotation object that can be used for creating ArrowFiles or an ArchRProject, etc.
#' 
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9"). If `genome` is not supplied,
#' `TxDb` and `OrgDb` are required. If genome is supplied, `TxDb` and `OrgDb` will be ignored.
#' @param TxDb A `TxDb` object (transcript database) from Bioconductor which contains information for gene/transcript coordinates.
#' For example, from `txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene`.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids.
#' For example, from `orgdb <- org.Hs.eg.db`.
#' @param genes A `GRanges` object containing gene coordinates (start to end). Must have a symbols column matching the symbols column of `exons`.
#' @param exons A `GRanges` object containing gene exon coordinates. Must have a symbols column matching the symbols column of `genes`.
#' @param TSS A `GRanges` object containing standed transcription start site coordinates for computing TSS enrichment scores downstream.
#' @param annoStyle annotation style to map between gene names and various gene identifiers e.g. "ENTREZID", "ENSEMBL".
#' @export
createGeneAnnotation <- function(
  genome = NULL,
  TxDb = NULL,
  OrgDb = NULL,
  genes = NULL,
  exons = NULL,
  TSS = NULL,
  annoStyle = NULL
  ){

  .validInput(input = genome, name = "genome", valid = c("character", "null"))
  .validInput(input = TxDb, name = "TxDb", valid = c("txdb", "character", "null"))
  .validInput(input = OrgDb, name = "OrgDb", valid = c("orgdb", "character", "null"))
  .validInput(input = genes, name = "genes", valid = c("GRanges", "null"))
  .validInput(input = exons, name = "exons", valid = c("GRanges", "null"))
  .validInput(input = TSS, name = "TSS", valid = c("GRanges", "null"))
  .validInput(input = annoStyle, name = "annoStyle", valid = c("character", "null"))

  if(is.null(genes) | is.null(exons) | is.null(TSS)){

    inGenes <- genes
    inExons <- exons
    inTSS <- TSS

    .requirePackage("GenomicFeatures", source = "bioc")

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
    message("Getting Genes..")
    genes <- GenomicFeatures::genes(TxDb)

    if(is.null(annoStyle)){
      isEntrez <- mcols(genes)$symbol <- tryCatch({
        suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
        TRUE
      }, error = function(x){
        FALSE
      })

      isEnsembl <- mcols(genes)$symbol <- tryCatch({
        suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
        TRUE
      }, error = function(x){
        FALSE
      })

      if(isEntrez){
        annoStyle <- "ENTREZID"
      }else if(isEnsembl){
        annoStyle <- "ENSEMBL"
      }else{
        stop("Could not identify keytype for annotation format!")
      }

    }
    annoStyle <- toupper(annoStyle)

    message("Determined Annotation Style = ", annoStyle)

    ###########################
    mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = annoStyle, multiVals = "first"))
    mcols(genes)$symbol[is.na(mcols(genes)$symbol)] <- paste0("NA_", mcols(genes)$gene_id)[is.na(mcols(genes)$symbol)]
    names(genes) <- NULL
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ###########################
    message("Getting Exons..")
    exons <- unlist(GenomicFeatures::exonsBy(TxDb, by = "tx"))
    exons$tx_id <- names(exons)
    mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(TxDb, keys = paste0(mcols(exons)$tx_id), column = "GENEID", keytype = "TXID")[, "GENEID"])
    exons <- exons[!is.na(mcols(exons)$gene_id), ]
    mcols(exons)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(exons)$gene_id, column = "SYMBOL", keytype = annoStyle, multiVals = "first"))
    mcols(exons)$symbol[is.na(mcols(exons)$symbol)] <- paste0("NA_", mcols(exons)$gene_id)[is.na(mcols(exons)$symbol)]
    names(exons) <- NULL
    mcols(exons)$exon_id <- NULL
    mcols(exons)$exon_name <- NULL
    mcols(exons)$exon_rank <- NULL
    mcols(exons)$tx_id <- NULL
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

    ###########################
    message("Getting TSS..")
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

.getBlacklist <- function(genome = NULL){
  
  .validInput(input = genome, name = "genome", valid = "character")

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
      message("Blacklist not downloaded! Continuing without, be careful for downstream biases..")
      GRanges()
    })
  }else{
    message("Blacklist not downloaded! Continuing without, be careful for downstream biases..")
    bl <- GRanges()
  }

  bl

}

.getTxDb <- function(genome = NULL, filter = TRUE, install = TRUE){

  .validInput(input = genome, name = "genome", valid = "character")
  .validInput(input = filter, name = "filter", valid = "boolean")
  .validInput(input = install, name = "install", valid = "boolean")

  if(toupper(genome) == "HG19"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg19.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
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
        message("Package does not exist, now trying bioconductor..")
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
        message("Package does not exist, now trying bioconductor..")
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
        message("Package does not exist, now trying bioconductor..")
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
        message("Package does not exist, now trying bioconductor..")
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
        message("Package does not exist, now trying bioconductor..")
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

.getOrgDb <- function(genome = NULL){

  .validInput(input = genome, name = "genome", valid = "character")

  if(toupper(genome) == "HG19" | toupper(genome) == "HG38"){
    if(suppressWarnings(!require(org.Hs.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Hs.eg.db", update=FALSE)
    }
    library(org.Hs.eg.db)
    annodb <- org.Hs.eg.db
  }else if(toupper(genome) == "MM9" | toupper(genome) == "MM10"){
    if(suppressWarnings(!require(org.Mm.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Mm.eg.db", update=FALSE)
    }
    library(org.Mm.eg.db)
    annodb <- org.Mm.eg.db
  }else{
    stop("Genome not recognized!")
  }
  return(annodb)

}






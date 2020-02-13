##########################################################################################
# Parallel Information
##########################################################################################

#' Add a globally-applied number of threads to use for parallel computing. JJJ
#' 
#' This function will set the number of threads to be used for parallel computing across all ArchR functions.
#' 
#' @param threads The default number of threads to be used for parallel execution across all ArchR functions. This value is stored as a global environment variable, not part of the `ArchRProject`. This can be overwritten on a per-function basis using the given function's `threads` parameter.
#' @param force If you requested for more than nCPU - 2 threads. If `force = FALSE` ArchR will set the threads maximum at nCPU - 2 threads. To bypass this set `force = TRUE`.
#' @export
addArchRThreads <- function(threads = floor(parallel::detectCores()/ 2), force = FALSE){
  
  .validInput(input = threads, name = "threads", valid = "integer")
  .validInput(input = force, name = "force", valid = "boolean")

  if(tolower(.Platform$OS.type) == "windows"){
    message("Detected windows OS, setting threads to 1.")
    threads <- 1
  }

  if(threads >= parallel::detectCores() - 1){
    if(force){
      message("Input threads is equal to or greater than ncores minus 1 (",parallel::detectCores()-1,")\nOverriding since force = TRUE.")
    }else{
      message("Input threads is equal to or greater than ncores minus 1 (",parallel::detectCores()-1,")\nSetting cores to ncores minus 2. Set force = TRUE to set above this number!")
      threads <- parallel::detectCores()-2
    }
  }
  
  message("Setting default number of Parallel threads to ", threads, ".")
  assign(".ArchRThreads", as.integer(threads), envir = .GlobalEnv)

}

#' Get globally-applied number of threads to use for parallel computing.
#' 
#' This function will get the number of threads to be used for parallel execution across all ArchR functions.
#' 
#' @export
getArchRThreads <- function(){
  if(exists(".ArchRThreads")){
    if(!is.integer(.ArchRThreads)){
      message(".ArchRThreads : ", .ArchRThreads, " is not an integer. \nDid you mistakenly set this to a value without addArchRThreads? Deleting .ArchRThreads from global environment.")
      rm(list=".ArchRThreads", envir = .GlobalEnv ) # Remove this.
      1
    }else{
      .ArchRThreads
    }
  }else{
    1
  }
}

##########################################################################################
# Create Gene/Genome Annotation
##########################################################################################

#' Add a globally defined genome to all ArchR functions. JJJ
#' 
#' This function will set the genome across all ArchR functions.
#' 
#' @param genome The default genome to be used for all ArchR functions. This value is stored as a global environment variable, not part of the `ArchRProject`. This can be overwritten on a per-function basis using the given function's `geneAnnotation` and  `genomeAnnotation` parameter.
#' @param install Install the BSgenome associated with the ArchRGenome is not currently installed. This is useful for helping reduce user download requirements.
#' @export
addArchRGenome <- function(genome = NULL, install = TRUE){
  
  .validInput(input = genome, name = "genome", valid = "character")

  supportedGenomes <- c("hg19","hg38","mm9","mm10")
  
  if(tolower(genome) %ni% supportedGenomes){
    
    message("Genome : ", genome, " is not currently supported by ArchR.")
    message("Currently supported genomes : ", paste0(supportedGenomes, collapse = ","))
    message("To continue try building a custom geneAnnotation with createGeneAnnnotation,\nand genomeAnnotation with createGenomeAnnotation!")
 
  }else{

    #Check if BSgenome exists!
    if(tolower(genome)=="hg19"){
      if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
        if(install){
          message("BSgenome for hg19 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        }else{
          stop("BSgenome for hg19 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
        }
      }
    }else if(tolower(genome)=="hg38"){
      if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
        if(install){
          message("BSgenome for hg38 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
        }else{
          stop("BSgenome for hg38 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
        }
      }
    }else if(tolower(genome)=="mm9"){
      if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE)){
        if(install){
          message("BSgenome for mm9 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
          BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
        }else{
          stop("BSgenome for mm9 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
        }
      }
    }else if(tolower(genome)=="mm10"){
      if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
        if(install){
          message("BSgenome for mm10 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
          BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
        }else{
          stop("BSgenome for mm10 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
        }
      }
    }

    genome <- paste(toupper(substr(genome, 1, 1)), substr(genome, 2, nchar(genome)), sep="")
    
    message("Setting default genome to ", genome, ".")
    assign(".ArchRGenome", genome, envir = .GlobalEnv)
    
  }

  invisible(0)

}

#' Get a globally defined genome to all ArchR functions. JJJ
#' 
#' This function will get the genome across all ArchR functions. Additionally can return geneAnnotation and genomeAnnotation associated with ArchRGenome if desired.
#' 
#' @param geneAnnotation Export geneAnnotation set with addArchRGenome. The geneAnnotation (see `createGeneAnnotation()`) to associate with the ArrowFiles. This is used downstream to calculate TSS Enrichment Scores etc.
#' @param genomeAnnotation Export genomeAnnotation set with addArchRGenome. The genomeAnnotation (see `createGenomeAnnotation()`) to associate with the ArrowFiles. This is used downstream to collect chromosome sizes and nucleotide information etc.
#' @export
getArchRGenome <- function(
  geneAnnotation=FALSE, 
  genomeAnnotation=FALSE
  ){

  supportedGenomes <- c("hg19","hg38","mm9","mm10")

  if(exists(".ArchRGenome")){

    ag <- .ArchRGenome
    
    if(!is.character(ag)){
    
      return(NULL)
    
    }else{
      
      if(tolower(ag) %in% supportedGenomes){
        
        genome <- paste(toupper(substr(ag, 1, 1)), substr(ag, 2, nchar(ag)), sep="")
        
        if(geneAnnotation){

          message("Using GeneAnnotation set by addArchRGenome!")

          geneAnno <- paste0("geneAnno", genome)
          eval(parse(text=paste0("data(geneAnno",genome,")")))
          return(eval(parse(text=geneAnno)))
        
        }else if(genomeAnnotation){

          message("Using GenomeAnnotation set by addArchRGenome!")

          genomeAnno <- paste0("genomeAnno", genome)
          eval(parse(text=paste0("data(genomeAnno",genome,")")))
          return(eval(parse(text=genomeAnno)))

        }else{

          return(genome)

        }

      }else{
        
        rm(list=".ArchRGenome", envir = .GlobalEnv ) # Remove this.
        stop(".ArchRGenome : ", ag, " is not currently supported by ArchR. \nDid you mistakenly set this to a value without addArchRGenome?")
      
      }
    }

  }else{
    
    return(NULL)

  }

}


#' Create a genome annotation object for ArchR
#' 
#' This function will create a genome annotation object that can be used for creating ArrowFiles or an ArchRProject, etc.
#' 
#' @param genome Either (i) a string that is a valid `BSgenome` or (ii) a `BSgenome` object (ie "hg38" or "BSgenome.Hsapiens.UCSC.hg38").
#' @param chromSizes A `GRanges` object containing chromosome start and end coordinates.
#' @param blacklist A `GRanges` object containing regions that should be excluded from analyses due to unwanted biases.
#' @param filter A boolean value indicating whether non-standard chromosome scaffolds should be excluded. These "non-standard" chromosomes are defined by `filterChrGR()`.
#' @param filterChr A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels. If no manual removal is desired, `filterChr` should be set to `NULL`.
#' @export
createGenomeAnnotation <- function(
  genome = NULL,
  chromSizes = NULL,
  blacklist = NULL,
  filter = TRUE,
  filterChr = c("chrM")
  ){

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
    blacklist <- .getBlacklist(genome = bsg@provider_version)

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
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9"). If `genome` is not supplied, `TxDb` and `OrgDb` are required. If genome is supplied, `TxDb` and `OrgDb` will be ignored.
#' @param TxDb A `TxDb` object (transcript database) from Bioconductor which contains information for gene/transcript coordinates. For example, from `txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene`.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids. For example, from `orgdb <- org.Hs.eg.db`.
#' @param genes A `GRanges` object containing gene coordinates (start to end). Must have a symbols column matching the symbols column of `exons`.
#' @param exons A `GRanges` object containing gene exon coordinates. Must have a symbols column matching the symbols column of `genes`.
#' @param TSS A `GRanges` object containing standed transcription start site coordinates for computing TSS enrichment scores downstream.
#' @export
createGeneAnnnotation <- function(
  genome = NULL,
  TxDb = NULL,
  OrgDb = NULL,
  genes = NULL,
  exons = NULL,
  TSS = NULL
  ){

  .validInput(input = genome, name = "genome", valid = c("character", "null"))
  .validInput(input = TxDb, name = "TxDb", valid = c("txdb", "character", "null"))
  .validInput(input = OrgDb, name = "OrgDb", valid = c("orgdb", "character", "null"))
  .validInput(input = genes, name = "genes", valid = c("GRanges", "null"))
  .validInput(input = exons, name = "exons", valid = c("GRanges", "null"))
  .validInput(input = TSS, name = "TSS", valid = c("GRanges", "null"))

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
    mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, 
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
    names(genes) <- NULL
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ###########################
    message("Getting Exons..")
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


##########################################################################################
# Output Directory
##########################################################################################

#' Get outputDirectory from an ArchRProject
#' 
#' This function gets the outputDirectory from a given ArchRProject. If null this returns "QualityControl" directory.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getOutputDirectory <- function(ArchRProj = NULL){

  if(is.null(ArchRProj) | is.character(ArchRProj)){
    return("QualityControl")
  }

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")

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
#' @export
getArrowFiles <- function(ArchRProj = NULL){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  
  af <- ArchRProj@sampleColData$ArrowFiles
  
  names(af) <- rownames(ArchRProj@sampleColData)
  
  return(af)

}

#' Get the sample names from an ArchRProject
#' 
#' This function gets the names of all samples from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getSampleNames <- function(ArchRProj = NULL){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  
  snames <- rownames(ArchRProj@sampleColData)
  
  return(snames)

}

#' Get the number of cells from an ArchRProject/ArrowFile
#' 
#' This function gets number of cells from an ArchRProject or ArrowFile
#' 
#' @param input An `ArchRProject` object or the path to an ArrowFile.
#' @export
nCells <- function(input = NULL){
  
  .validInput(input = input, name = "input", valid = c("ArchRProject", "character"))

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
#' This function gets the sampleColData from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param select A character vector containing the column names to select from `sampleColData`.
#' @param drop A boolean value that indicates whether to drop the `dataframe` structure and convert to a vector if selecting only one column.
#' @export
getSampleColData <- function(ArchRProj = NULL, select = NULL, drop = FALSE){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = select, name = "select", valid = c("character", "null"))
  .validInput(input = drop, name = "drop", valid = "boolean")

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

#' Add information to sampleColData in an ArchRProject
#' 
#' This function adds new data to sampleColData in an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param data A vector containing the data to be added to `sampleColData`.
#' @param name The column header name to be used for this new data in `sampleColData`. If a column with this name already exists, you may set `force` equal to `TRUE` to overwrite the data in this column.
#' @param samples The names of the samples corresponding to `data`. Typically new data is added to all samples but you may use this argument to only add data to a subset of samples. Samples where `data` is not added are set to `NA`.
#' @param force A boolean value that indicates whether or not to overwrite data in a given column when the value passed to `name` already exists as a column name in `sampleColData`.
#' @export
addSampleColData <- function(ArchRProj = NULL, data = NULL, name = NULL, samples = NULL, force = FALSE){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = data, name = "data", valid = c("character", "integer", "numeric", "boolean"))
  .validInput(input = samples, name = "samlpes", valid = c("character"))
  .validInput(input = force, name = "force", valid = "boolean")
  .validInput(input = name, name = "name", valid = c("character"))

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
#' @export
getCellNames <- function(ArchRProj = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
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
#' @export
getCellColData <- function(ArchRProj = NULL, select = NULL, drop = FALSE){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = select, name = "select", valid = c("character", "null"))
  .validInput(input = drop, name = "drop", valid = "boolean")

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
#' @param name The column header name to be used for this new data in `cellColData`. If a column with this name already exists, you may set `force` equal to `TRUE` to overwrite the data in this column.
#' @param cells The names of the cells corresponding to `data`. Typically new data is added to all cells but you may use this argument to only add data to a subset of cells. Cells where `data` is not added are set to `NA`.
#' @param force A boolean value indicating whether or not to overwrite data in a given column when the value passed to `name` already exists as a column name in `cellColData`.
#' @export
addCellColData <- function(ArchRProj = NULL, data = NULL, name = NULL, cells =  NULL, force = FALSE){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = data, name = "data", valid = c("character", "integer", "numeric", "boolean"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = cells, name = "cells", valid = c("character"))
  .validInput(input = force, name = "force", valid = "boolean")

  if(name == "cellNames"){
    stop("cellNames is a protected column name in an ArchRProject!")
  }

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
#' @param ArchRProj An `ArchRProject` object.
#' @export
getPeakSet <- function(ArchRProj = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  return(ArchRProj@peakSet)
}

#' Add a peak set to an ArchRProject
#' 
#' This function adds a peak set as a GRanges object to a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param peakSet A `GRanges` object containing the set of regions that define all peaks in the desired peak set.
#' @param force If a `peakSet` object has already been added to the given `ArchRProject`, the value of `force` determines whether or not to overwrite this `peakSet`.
#' @export
addPeakSet <- function(ArchRProj = NULL, peakSet = NULL, force = FALSE){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = peakSet, name = "peakSet", valid = c("GRanges"))
  .validInput(input = force, name = "force", valid = c("boolean"))
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

#' Get the genomeAnnotation from an ArchRProject
#' 
#' This function gets the genomeAnnotation from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getGenomeAnnotation <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  return(ArchRProj@genomeAnnotation)
}

#' Get the blacklist from an ArchRProject
#' 
#' This function gets the blacklist (the regions to be excluded from analysis) as a GRanges object from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getBlacklist <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$blacklist)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  return(ArchRProj@genomeAnnotation$blacklist)
}

#' Get the genome used by an ArchRProject
#' 
#' This function gets the name of the genome from the genomeAnnotation used by a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getGenome <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$genome)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  return(ArchRProj@genomeAnnotation$genome)
}

#' Get chromSizes from ArchRProject
#' 
#' This function gets the chromosome lengths as a GRanges object from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getChromSizes <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$chromSizes)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  return(ArchRProj@genomeAnnotation$chromSizes)
}

#' Get chromLengths from ArchRProject
#' 
#' This function gets the chromosome lengths as a vector from the genomeAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getChromLengths <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      cS <- genomeAnnotation$chromSizes
      cL <- end(cS)
      names(cL) <- paste0(seqnames(cS))
      return(cL)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  cS <- ArchRProj@genomeAnnotation$chromSizes
  cL <- end(cS)
  names(cL) <- paste0(seqnames(cS))
  return(cL)
}

.nullGenomeAnnotation <- function(){
  genome <- "nullGenome"
  chromSizes <- GRanges()
  blacklist <- GRanges()
  SimpleList(blacklist = blacklist, genome = genome, chromSizes = chromSizes)
}

##########################################################################################
# Gene Annotation Methods
##########################################################################################

#' Get geneAnnotation from an ArchRProject
#' 
#' This function gets the geneAnnotation from a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getGeneAnnotation <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      return(geneAnnotation)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  ArchRProj@geneAnnotation
}

#' Get the transcription start sites of all genes in an ArchRProject
#' 
#' This function gets the transcription start sites (TSSs) as a GRanges object of all genes from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getTSS <- function(ArchRProj = NULL){
  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      return(geneAnnotation$TSS)
    }
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  ArchRProj@geneAnnotation$TSS
}

#' Get the genes from an ArchRProject
#' 
#' This function gets the genes start and end coordinates as a GRanges object from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param symbols A character vector containing the gene symbols to subset from the `geneAnnotation`.
#' @export
getGenes <- function(ArchRProj = NULL, symbols = NULL){
  
  .validInput(input = symbols, name = "symbols", valid = c("character", "null"))

  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      genes <- geneAnnotation$genes
      if(!is.null(symbols)){
        genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
      }
      return(genes)
    }
  }

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  
  genes <- ArchRProj@geneAnnotation$genes
  
  if(!is.null(symbols)){
    genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
  }

  genes

}

#' Get the exons from an ArchRProject
#' 
#' This function gets the exons coordinates as a GRanges object from the geneAnnotation of a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param symbols A character vector containing the gene symbols for the genes where exons should be extracted.
#' @export
getExons <- function(ArchRProj = NULL, symbols = NULL){

  .validInput(input = symbols, name = "symbols", valid = c("character", "null"))

  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      exons <- geneAnnotation$exons
      if(!is.null(symbols)){
        exons <- exons[which(tolower(exons$symbol) %in% tolower(symbols))]
      }
      return(exons)
    }
  }

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")

  exons <- ArchRProj@geneAnnotation$exons
  exons <- exons[which(tolower(exons$symbol) %in% tolower(symbols))]

  exons

}

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
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param returnMatrix If set to "mat" or "matrix", the function will return the `reducedDims` object as a matrix with entries for each individual cell. Otherwise, it will return the full `reducedDims` object.
#' @param dimsToUse A vector containing the dimensions (i.e. 1:30) to return from the `reducedDims` object.
#' @param scaleDims A boolean describing whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution of strong biases 
#' (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. 
#' If `NULL` this will scale the dimensions depending on if this were set true when the `reducedDims` were created by the dimensionality reduction method.
#' This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded.
#' @export
getReducedDims <- function(
  ArchRProj = NULL, 
  reducedDims = "IterativeLSI", 
  returnMatrix = TRUE, 
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75
  ){

  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = reducedDims, name = "reducedDims", valid = "character")
  .validInput(input = returnMatrix, name = "returnMatrix", valid = "boolean")
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  #########

  if(reducedDims %in% names(ArchRProj@reducedDims)){
    
    if(is.na(ArchRProj@reducedDims[[reducedDims]]$scaleDims[1])){
      scaleDims <- FALSE # if na this means dont scaleDims ever.
    }

    if(is.null(scaleDims)){
      scaleDims <- ArchRProj@reducedDims[[reducedDims]]$scaleDims
    }

    #Get Dimensions
    if(scaleDims){
      corToDepth <- ArchRProj@reducedDims[[reducedDims]]$corToDepth$scaled
      matDR <- .scaleDims(ArchRProj@reducedDims[[reducedDims]][[1]])
    }else{
      if(is.na(ArchRProj@reducedDims[[reducedDims]]$corToDepth[1])){
        corToDepth <- rep(0, ncol(ArchRProj@reducedDims[[reducedDims]][[1]]))
        matDR <- ArchRProj@reducedDims[[reducedDims]][[1]]
      }else{
        corToDepth <- ArchRProj@reducedDims[[reducedDims]]$corToDepth$none
        matDR <- ArchRProj@reducedDims[[reducedDims]][[1]]
      }
    }
    
    #Determine PCs to Keep
    if(!is.null(dimsToUse)){
      corToUse <- dimsToUse
    }else{
      corToUse <- seq_along(corToDepth)
    }
    idx <- which(abs(corToDepth[corToUse]) <= corCutOff)
    if(length(idx) != length(corToUse)){
      message("Filtering ", length(corToUse) - length(idx), " dims correlated > ", corCutOff, " to log10(depth + 1)")
    }

    #Return Values
    if(returnMatrix){

      cells <- rownames(matDR) %in% ArchRProj$cellNames
      
      return(matDR[cells,corToUse[idx],drop=FALSE])

    }else{
      
      cells <- rownames(matDR) %in% ArchRProj$cellNames

      out <- ArchRProj@reducedDims[[reducedDims]]
      out$dimsKept <- corToUse[idx]
      out[[1]] <- matDR[cells,corToUse[idx],drop=FALSE]
      return(out)

    }

  }else{

    stop(paste0("Reduced dimensions not in computed reduceDims, Current ones are : ", paste0(names(ArchRProj@reducedDims), collapse=",")))
  
  }

}

.scaleDims <- function(x, scaleMax = NULL){
  if(!is.null(scaleMax)){
    .rowZscores(m=x, min=-scaleMax, max = scaleMax, limit = TRUE)
  }else{
    .rowZscores(m=x)
  }
}

#' Get embedding information stored in an ArchRProject
#' 
#' This function gets an embedding (i.e. UMAP) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the `embeddings` object (i.e. UMAP, TSNE see `embeddingOut` of the `addEmbeddings()` function) to retrieve from the designated `ArchRProject`.
#' @param returnDF A boolean value indicating whether to return the embedding object as a `data.frame`. Otherwise, it will return the full embedding object.
#' @export
getEmbedding <- function(ArchRProj = NULL, embedding = "UMAP", returnDF = TRUE){

  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = embedding, name = "embedding", valid = "character")
  .validInput(input = returnDF, name = "returnDF", valid = "boolean")
  #########

  if(embedding %in% names(ArchRProj@embeddings)){
    if(returnDF){
      out <- ArchRProj@embeddings[[embedding]][[1]]
    }else{
      out <- ArchRProj@embeddings[[embedding]]
    }
  }else{
    stop(paste0("Embedding not in computed embeddings, Current ones are : ", paste0(names(ArchRProj@embeddings), collapse=",")))
  }
  return(out)
}

##########################################################################################
# Project Summary
##########################################################################################

#' Get projectSummary from an ArchRProject
#' 
#' This function prints the projectSummary from an ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param returnSummary A boolean value indicating whether to return a summary of the `ArchRProject` or to just print the summary.
#' @export
getProjectSummary <- function(ArchRProj = NULL, returnSummary = FALSE){

  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = returnSummary, name = "returnSummary", valid = "boolean")
  #########

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

#' Add projectSummary to an ArchRProject
#' 
#' This function adds info to the projectSummary of an ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name The name of the summary information to add to the `ArchRProject` object.
#' @param summary A vector to add as summary information to the `ArchRProject` object.
#' @export
addProjectSummary <- function(ArchRProj = NULL, name = NULL, summary = NULL){

  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = name, name = "name", valid = "character")
  #########

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
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
#' @param select A string specifying a specific feature name (or rowname) to be found with `grep`.
#' @param ignoreCase A boolean value indicating whether to ignore the case (upper-case / lower-case) when searching via grep for the string passed to `select`.
#' @export
getFeatures <- function(ArchRProj = NULL, useMatrix = "GeneScoreMatrix", select = NULL, ignoreCase = TRUE){

  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  .validInput(input = select, name = "select", valid = c("character", "null"))
  .validInput(input = ignoreCase, name = "ignoreCase", valid = "boolean")
  #########

  fdf <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class"))
  if(is.null(select)){
    if(any(duplicated(paste0(fdf$name))) | matrixClass == "Sparse.Assays.Matrix"){
      paste0(fdf$seqnames,":",fdf$name)
    }else{
      fdf$name
    }
  }else{
    grepNames <- grep(select, fdf$name, value = TRUE, ignore.case = ignoreCase)
    if(any(duplicated(grepNames))){
      grepIdx <- grep(select, fdf$name, ignore.case = ignoreCase)
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
#' This function will save a plot or set of plots as a PDF file in the outputDirectory of a given ArchRProject.
#' 
#' @param ... vector of plots to be plotted (if input is a list use plotList instead)
#' @param name The file name to be used for the output PDF file.
#' @param width The width in inches to be used for the output PDF file.
#' @param height The height in inches to be used for the output PDF.
#' @param ArchRProj An `ArchRProject` object to be used for retrieving the desired `outputDirectory` which will be used to store the output plots in a subfolder called "plots".
#' @param addDOC A boolean variable that determines whether to add the date of creation to the end of the PDF file name. This is useful for preventing overwritting of old plots.
#' @param useDingbats A boolean variable that determines wheter to use dingbats characters for plotting points.
#' @param plotList A `list` of plots to be printed to the output PDF file. Each element of `plotList` should be a printable plot formatted object (ggplot2, plot, heatmap, etc).
#' @param useSink A boolean value that indicates whether the `sink` function from base R should be used to hide messages during plotting.
#' @export
plotPDF <- function(..., name = "Plot", width = 6, 
  height = 6, ArchRProj = NULL, addDOC = TRUE, 
  useDingbats = FALSE, plotList = NULL, useSink = TRUE){

  #Validate
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = width, name = "width", valid = "numeric")
  .validInput(input = height, name = "height", valid = "numeric")
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject", "null"))
  .validInput(input = addDOC, name = "addDOC", valid = "boolean")
  .validInput(input = useDingbats, name = "useDingbats", valid = "boolean")
  .validInput(input = plotList, name = "plotList", valid = c("list","null"))
  .validInput(input = useSink, name = "useSink", valid = "boolean")
  #########

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
    .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
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

#' Get Relevant Data For ArchR Tutorials JJJ
#' 
#' This function will download data for a given tutorial and return the input files required for ArchR
#' 
#' @param tutorial The name of the available tutorial for which to retreive the tutorial data. Options are "Hematopoiesis" for now. "Hematopoiesis" refers to hematopoieitic scATAC hierarchy. "PBMC" refers to a small standard PBMC scATAC dataset. "FreshFrozen" refers to a PBMC fresh and frozen scATAC dataset.
#' @export
getTutorialData <- function(tutorial = "hematopoiesis"){

  #Validate
  .validInput(input = tutorial, name = "tutorial", valid = "character")
  #########

  if(tolower(tutorial) %in% c("heme","hematopoiesis")){
    
    if(!dir.exists("Heme_Fragments")){
      download.file(
        url = "https://jeffgranja.s3.amazonaws.com/ArchR-Tutorial-Data/Heme/HemeFragments.zip", 
        destfile = "HemeFragments.zip"
      )
      unzip("HemeFragments.zip")
      if(dir.exists("HemeFragments")){
        file.remove("HemeFragments.zip")
      }else{
        stop("Download May Not Have Worked!")
      }
    }
    pathFragments <- "HemeFragments"

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
#' @export
getInputFiles <- function(paths = NULL){ 

  #Validate
  .validInput(input = paths, name = "paths", valid = "character")
  #########

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
#' @param csvFiles A character vector of names from 10x CSV files to be read in for identification of valid cell barcodes.
#' @param sampleNames A character vector containing the sample names to be associated with each individual entry in `csvFiles`.
#' @export
getValidBarcodes <- function(csvFiles = NULL, sampleNames = NULL){

  #Validate
  .validInput(input = csvFiles, name = "csvFiles", valid = "character")
  .validInput(input = sampleNames, name = "sampleNames", valid = "character")
  #########

  if(length(sampleNames) != length(csvFiles)){
    stop("csvFiles and sampleNames must exist!")
  }

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


#' Get a list available matrices in the ArrowFiles storted in an ArchRProject JJJ
#' 
#' This function gets the available matrices from the ArrowFiles in a given ArchRProject object.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getAvailableMatrices <- function(ArchRProj = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .availableArrays(getArrowFiles(ArchRProj=ArchRProj))
}









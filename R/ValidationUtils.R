##########################################################################################
# Validation Methods
##########################################################################################

.validInput <- function(input = NULL, name = NULL, valid = NULL){

  valid <- unique(valid)

  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }

  if(!is.character(name)){
    stop("name must be a character!")
  }

  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }

  av <- FALSE

  for(i in seq_along(valid)){

    vi <- valid[i]

    if(vi == "integer" | vi == "wholenumber"){

      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }

    }else if(vi == "null"){

      cv <- is.null(input)

    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){

      cv <- is.logical(input)

    }else if(vi == "numeric"){

      cv <- is.numeric(input)

    }else if(vi == "vector"){

      cv <- is.vector(input)

    }else if(vi == "matrix"){

      cv <- is.matrix(input)

    }else if(vi == "sparsematrix"){

      cv <- is(input, "dgCMatrix")

    }else if(vi == "character"){

      cv <- is.character(input)

    }else if(vi == "factor"){

      cv <- is.factor(input)

    }else if(vi == "rlecharacter"){

      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }

    }else if(vi == "palette"){

      cv <- all(.isColor(input))

    }else if(vi == "timestamp"){

      cv <- is(input, "POSIXct")

    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){

      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)

    }else if(vi == "fileexists"){

      cv <- all(file.exists(input))

    }else if(vi == "direxists"){

      cv <- all(dir.exists(input))

    }else if(vi == "granges" | vi == "gr"){

      cv <- is(input, "GRanges")

    }else if(vi == "grangeslist" | vi == "grlist"){

      cv <- .isGRList(input)

    }else if(vi == "list" | vi == "simplelist"){

      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)

    }else if(vi == "bsgenome"){

      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)

    }else if(vi == "se" | vi == "summarizedexperiment"){

      cv <- is(input, "SummarizedExperiment")

    }else if(vi == "seurat" | vi == "seuratobject"){

      cv <- is(input, "Seurat")

    }else if(vi == "txdb"){

      cv <- is(input, "TxDb")

    }else if(vi == "orgdb"){

      cv <- is(input, "OrgDb")

    }else if(vi == "bsgenome"){

      cv <- is(input, "BSgenome")

    }else if(vi == "parallelparam"){

      cv <- is(input, "BatchtoolsParam")

    }else if(vi == "archrproj" | vi == "archrproject"){

      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ

    }else{

      stop("Validator is not currently supported by ArchR!")

    }

    if(cv){
      av <- TRUE
      break
    }   
     
  }

  if(av){

    return(invisible(TRUE))
  
  }else{

    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")

  }

}

#https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
.isWholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

#https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
.isColor <- function(x = NULL){
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
}

.isDiscrete <- function(x = NULL){
  is.factor(x) || is.character(x) || is.logical(x)
}

.isGRList <- function(x){
  isList <- grepl("list", class(x), ignore.case=TRUE)
  if(!isList){
    FALSE
  }else{
    allGR <- all(unlist(lapply(x, function(x) is(x, "GRanges") )))
    if(allGR){
      TRUE
    }else{
      FALSE
    }
  }
}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid ArchR-supported genome ("hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = masked, name = "masked", valid = c("boolean"))

  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}

.validTxDb <- function(TxDb = NULL){
  stopifnot(!is.null(TxDb))
  if(inherits(TxDb, "TxDb")){
    return(TxDb)
  }else if(is.character(TxDb)){
    return(getTxDb(TxDb)) #change
  }else{
    stop("Cannot validate TxDb options are a valid TxDb or character for getTxDb")
  }
}

.validOrgDb <- function(OrgDb = NULL){
  stopifnot(!is.null(OrgDb))
  if(inherits(OrgDb, "OrgDb")){
    return(OrgDb)
  }else if(is.character(OrgDb)){
    return(getOrgDb(OrgDb)) #change
  }else{
    stop("Cannot validate OrgDb options are a valid OrgDb or character for getOrgDb")
  }
}

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.validGeneAnnotation <- function(geneAnnotation = NULL){
  
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

.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
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
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])

  }else{

    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
  
  }

  gA

}

.validGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation){

  allSeqs <- unique(paste0(seqnames(genomeAnnotation$chromSizes)))

  geneSeqs <- unique(paste0(seqnames(geneAnnotation$genes)))
  if(!all(geneSeqs %in% allSeqs)){
    geneNotIn <- geneSeqs[which(geneSeqs %ni% allSeqs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(geneNotIn, collapse=","))
    geneAnnotation$genes <- .subsetSeqnamesGR(geneAnnotation$genes, names = allSeqs)
  }

  exonSeqs <- unique(paste0(seqnames(geneAnnotation$exons)))
  if(!all(exonSeqs %in% allSeqs)){
    exonNotIn <- exonSeqs[which(exonSeqs %ni% allSeqs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(exonNotIn, collapse=","))
    geneAnnotation$exons <- .subsetSeqnamesGR(geneAnnotation$exons, names = allSeqs)
  }

  TSSSeqs <- unique(paste0(seqnames(geneAnnotation$TSS)))
  if(!all(TSSSeqs %in% allSeqs)){
    TSSNotIn <- TSSSeqs[which(TSSSeqs %ni% allSeqs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(TSSNotIn, collapse=","))
    geneAnnotation$TSS <- .subsetSeqnamesGR(geneAnnotation$TSS, names = allSeqs)
  }

  geneAnnotation

}


.validArchRProject <- function(ArchRProj = NULL){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}




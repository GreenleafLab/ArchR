##########################################################################################
# Output Directory
##########################################################################################

#' Get outputDirectory from an ArchRProject
#' 
#' This function gets the outputDirectory from a given ArchRProject. If null this returns "QualityControl" directory.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getOutputDirectory <- function(
  ArchRProj = NULL
  ){

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
getArrowFiles <- function(
  ArchRProj = NULL
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  
  af <- ArchRProj@sampleColData$ArrowFiles
  
  names(af) <- rownames(ArchRProj@sampleColData)

  af <- af[unique(ArchRProj$Sample)]

  idx <- tryCatch({
      order(file.info(af)$size, decreasing = TRUE)
  }, error = function(x) {
      seq_along(af)
  })

  af <- af[idx]
  
  return(af)

}

#' Get the sample names from an ArchRProject
#' 
#' This function gets the names of all samples from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getSampleNames <- function(
  ArchRProj = NULL
  ){

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
nCells <- function(
  input = NULL
  ){
  
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

#' Get summary for Groups in ArchRProject
#' 
#' This function summarizes a numeric cellColData entry across groupings in a ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping multiple cells together for summarizing information.
#' @param select A character vector containing the column names to select from `cellColData`.
#' @param summary A character vector describing which method for summarizing across group. Options include "median", "mean", or "sum".
#' @param removeNA Remove NA's from summary method.
#' @export
getGroupSummary <- function(
  ArchRProj = NULL,
  groupBy = "Sample",
  select = "TSSEnrichment",
  summary = "median",
  removeNA = TRUE
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = select, name = "select", valid = c("character"))
  .validInput(input = summary, name = "summary", valid = c("character"))
  .validInput(input = removeNA, name = "removeNA", valid = c("boolean"))

  message("Getting ", summary, " of ", select, " for ", groupBy)

  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  select <- getCellColData(ArchRProj, select = select, drop = TRUE)

  if(!is.numeric(select)){
    stop("Select must be a numeric column!")
  }

  splitSelect <- split(select, groups)

  summarySelect <- lapply(seq_along(splitSelect), function(x){
    if(tolower(summary) == "median"){
      median(splitSelect[[x]], na.rm = removeNA)
    }else if(tolower(summary) == "mean"){
      mean(splitSelect[[x]], na.rm = removeNA)
    }else if(tolower(summary) == "sum"){
      sum(splitSelect[[x]], na.rm = removeNA)
    }else{
      stop("Summary Method Not Supported!")
    }
  }) %>% unlist

  names(summarySelect) <- names(splitSelect)
  summarySelect <- summarySelect[gtools::mixedsort(names(summarySelect))]

  summarySelect
  
}

#' Get sampleColData from an ArchRProject
#' 
#' This function gets the sampleColData from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param select A character vector containing the column names to select from `sampleColData`.
#' @param drop A boolean value that indicates whether to drop the `dataframe` structure and convert to a vector if selecting only one column.
#' @export
getSampleColData <- function(
  ArchRProj = NULL, 
  select = NULL, 
  drop = FALSE
  ){

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
#' @param name The column header name to be used for this new data in `sampleColData`.
#' If a column with this name already exists, you may set `force` equal to `TRUE` to overwrite the data in this column.
#' @param samples The names of the samples corresponding to `data`. Typically new data is added to all samples but you may
#' use this argument to only add data to a subset of samples. Samples where `data` is not added are set to `NA`.
#' @param force A boolean value that indicates whether or not to overwrite data in a given column when the value passed to `name`
#' already exists as a column name in `sampleColData`.
#' @export
addSampleColData <- function(ArchRProj = NULL, data = NULL, name = NULL, samples = NULL, force = FALSE){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = data, name = "data", valid = c("character", "integer", "numeric", "boolean"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = samples, name = "samlpes", valid = c("character"))
  .validInput(input = force, name = "force", valid = "boolean")


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
#' @param name The column header name to be used for this new data in `cellColData`. If a column with this name already exists,
#' you may set `force` equal to `TRUE` to overwrite the data in this column.
#' @param cells The names of the cells corresponding to `data`. Typically new data is added to all cells but you may use this
#' argument to only add data to a subset of cells. Cells where `data` is not added are set to `NA`.
#' @param force A boolean value indicating whether or not to overwrite data in a given column when the value passed to `name`
#' already exists as a column name in `cellColData`.
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
#' @param genomeAnnotation The genomeAnnotation (see `createGenomeAnnotation()`) to be used for generating peak metadata such as nucleotide
#' information (GC content) or chromosome sizes.
#' @param force If a `peakSet` object has already been added to the given `ArchRProject`, the value of `force` determines
#' whether or not to overwrite this `peakSet`.
#' @export
addPeakSet <- function(
  ArchRProj = NULL, 
  peakSet = NULL, 
  genomeAnnotation = getGenomeAnnotation(ArchRProj),
  force = FALSE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = peakSet, name = "peakSet", valid = c("GRanges"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
  
  if(is.null(ArchRProj@peakSet) | force){
   
    #Index The Peak Set
    peakSet <- lapply(split(peakSet, seqnames(peakSet)), function(x){
      mcols(x)$idx <- seq_along(x)
      x
    }) %>% Reduce("c", .) %>% sortSeqlevels %>% sort

    #Get NucleoTide Content
    peakSet <- tryCatch({
      .requirePackage("Biostrings",source="bioc")
      BSgenome <- eval(parse(text = genomeAnnotation$genome))
      BSgenome <- validBSgenome(BSgenome)
      nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peakSet))
      mcols(peakSet)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
      mcols(peakSet)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
      peakSet
    }, error = function(e){
      peakSet
    })

    #Add PeakSet
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation)
    }
    stop("getGenomeAnnotation : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  if(is.character(ArchRProj)){
    ArchRProj <- NULL
  }
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$blacklist)
    }
    stop("getBlacklist : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$genome)
    }
    stop("getGenome : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation$chromSizes)
    }
    stop("getChromSizes : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    genomeAnnotation <- getArchRGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      cS <- genomeAnnotation$chromSizes
      cL <- end(cS)
      names(cL) <- paste0(seqnames(cS))
      return(cL)
    }
    stop("getChromLengths : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      return(geneAnnotation)
    }
    stop("getGeneAnnotation : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
  if(is.null(ArchRProj)){
    geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
    if(!is.null(geneAnnotation)){
      return(geneAnnotation$TSS)
    }
    stop("getTSS : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
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
    stop("getGenes : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
  }

  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  
  genes <- ArchRProj@geneAnnotation$genes
  
  if(!is.null(symbols)){
    genes <- genes[which(tolower(genes$symbol) %in% tolower(symbols))]
  }

  if(inherits(mcols(genes)$symbol, "list") | inherits(mcols(genes)$symbol, "SimpleList")){
    stop("Found a list in genes symbol! This is an incorrect format. Please correct your genes!")
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
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject","null"))
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
    stop("getExons : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
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
#' @param returnMatrix If set to "mat" or "matrix", the function will return the `reducedDims` object as a matrix with entries for
#' each individual cell. Otherwise, it will return the full `reducedDims` object.
#' @param dimsToUse A vector containing the dimensions (i.e. 1:30) to return from the `reducedDims` object.
#' @param scaleDims A boolean describing whether to z-score the reduced dimensions for each cell. This is useful for minimizing the
#' contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If `NULL` this will scale the dimensions depending on if this were set true when the
#' `reducedDims` were created by the dimensionality reduction method. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation
#' to sequencing depth that is greater than the `corCutOff`, it will be excluded.
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
#' @param embedding The name of the `embeddings` object (i.e. UMAP, TSNE see `embeddingOut` of the `addEmbeddings()` function) to
#' retrieve from the designated `ArchRProject`.
#' @param returnDF A boolean value indicating whether to return the embedding object as a `data.frame`. Otherwise, it will return
#' the full embedding object.
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
  #.validInput(input = summary, name = "summary", valid = c("list", "vector", "character", "numeric"))
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
#' This function will identify available features from a given data matrix  (i.e. "GeneScoreMatrix", or "TileMatrix") and return
#' them for downstream plotting utilities.
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

#' Get the seqnames that could be selected from a given data matrix within an ArchRProject
#' 
#' This function will identify available seqnames from a given data matrix  (i.e. "GeneScoreMatrix", or "TileMatrix") and return
#' them for downstream plotting utilities.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix as stored in the ArrowFiles of the `ArchRProject`. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @export
getSeqnames <- function(ArchRProj = NULL, useMatrix = "GeneScoreMatrix"){
  #Validate
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  #########
  fdf <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  unique(paste0(fdf$seqnames))
}

#' Get a list available matrices in the ArrowFiles storted in an ArchRProject
#' 
#' This function gets the available matrices from the ArrowFiles in a given ArchRProject object.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @export
getAvailableMatrices <- function(ArchRProj = NULL){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .availableArrays(getArrowFiles(ArchRProj=ArchRProj))
}

#' This function will add total counts of scATAC cells in provided features into ArchRProject.
#' 
#' This function will add total counts of scATAC cells in provided features into ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param features A `GRanges` object of features to count scATAC-seq data in.
#' @param name A character defining the name of the features. "`name`Counts" and "`name`Ratio" will be added to the `ArchRProject`.
#' @param addRatio A boolean indicating whether to add the "`name`Ratio" to the `ArchRProject`.
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addFeatureCounts <- function(
  ArchRProj = NULL,
  features = NULL,
  name = NULL,
  addRatio = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("addFeatureCounts")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = features, name = "features", valid = c("granges"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = addRatio, name = "addRatio", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  ArrowFiles <- getArrowFiles(ArchRProj)
  cellNames <- ArchRProj$cellNames
  featuresList <- split(features, seqnames(features))

  h5disableFileLocking()

  countsDF <- .safelapply(seq_along(featuresList), function(i){

    chri <- names(featuresList)[i]
    cellTotal <- rep(0, length(cellNames))
    names(cellTotal) <- cellNames
    featuresi <- ranges(featuresList[[i]])
    .logDiffTime(paste0("Counting in ",chri," (", i, " of ", length(featuresList), ")"), t1 = tstart, logFile = logFile)

    for(j in seq_along(ArrowFiles)){

      fragmentsij <- .getFragsFromArrow(
        ArrowFile = ArrowFiles[j], 
        chr = chri, 
        out = "IRanges", 
        cellNames = cellNames
      )
      if(length(fragmentsij) > 0){
        
        #Set To Integers
        mcols(fragmentsij)$RG@values <- match(mcols(fragmentsij)$RG@values, cellNames)

        for(y in seq_len(2)){   
          if(y==1){
            temp <- IRanges(start(fragmentsij), width = 1)
          }else if(y==2){
            temp <- IRanges(end(fragmentsij), width = 1)
          }
          stopifnot(length(temp) == length(fragmentsij))
          tabSum <- S4Vectors:::tabulate(mcols(fragmentsij)$RG[queryHits(findOverlaps(temp, featuresi))])
          cellTotal[seq_along(tabSum)] <- cellTotal[seq_along(tabSum)] + tabSum
        }

      }

      if(j %% 3 == 0){
        gc()
      }

    }

    cellTotal

  }, threads = threads) %>% Reduce("rbind", .)

  totalCounts <- colSums(countsDF)
  countRatio <- totalCounts / (2 * ArchRProj$nFrags)

  .logDiffTime(sprintf("Adding %s to cellColData", paste0(name,"Counts")), t1 = tstart, verbose = TRUE, logFile = logFile)
  ArchRProj <- addCellColData(ArchRProj, data = totalCounts, cells = names(totalCounts),  name = paste0(name,"Counts"), force = TRUE)
 
  .logDiffTime(sprintf("Adding %s to cellColData", paste0(name,"Ratio")), t1 = tstart, verbose = TRUE, logFile = logFile)
  ArchRProj <- addCellColData(ArchRProj, data = countRatio, cells = names(totalCounts),  name = paste0(name,"Ratio"), force = TRUE)

  ArchRProj

}


# addColorPalette <- function(
#   ArchRProj = NULL,
#   pal = NULL
#   ){

# }

# getColorPalette <- function(
#   ArchRProj = NULL,
#   name = NULL
#   ){

# }


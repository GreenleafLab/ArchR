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
  # Need to put this code
}

.validGenomeAnnotation <- function(genomeAnnotation, ...){
  # Need to put this code
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
#' @param select QQQ A character vector containing the column names to select from sampleColData.
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
#' This function gets the genomeAnnotation (in format QQQ) from a given ArchRProject.
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
#' This function gets the geneAnnotation (in format QQQ) from a given ArchRProject
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
#' This function gets the genes (in format QQQ) from the geneAnnotation of a given ArchRProject.
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
#' This function gets the exons (in format QQQ) from the geneAnnotation of a given ArchRProject.
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
#' @param reducedDims QQQ The name of the `reducedDims` object to retrieve from the designated `ArchRProject`. Options include QQQ.
#' @param returnMatrix If set to "mat" or "matrix", the function will return the `reducedDims` object as a matrix with entries for each individual cell. Otherwise, it will return the full `reducedDims` object.
#' @param dimsToUse QQQ A vector containing the dimensions to return from the `reducedDims` object.
#' @param corCutOff QQQ A numeric cutoff for the correlation of each dimension to the sequencing depth.
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
#' QQQ This function gets an embedding (i.e. UMAP) from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding QQQ The name of the `embeddings` object to retrieve from the designated `ArchRProject`. Options include QQQ.
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
#' @param name QQQ The name of the summary information to add to the `ArchRProject` object.
#' @param summary QQQ A vector to add as summary information to the `ArchRProject` object.
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

#' Get annotation from an ArchRProject
#' 
#' This function gets an annotation from a given ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param name QQQ The name of the annotation object to retrieve from the designated `ArchRProject`. Options include QQQ.
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
#' @param ArchRProj An `ArchRProject` object.
#' @param name QQQ The name of the annotation object to retrieve from the designated `ArchRProject`. Options include QQQ.
#' @param annoName QQQ ??? name to subset with annotations
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
#' @param ArchRProj An `ArchRProject` object.
#' @param name QQQ The name of the annotation object to retrieve from the designated `ArchRProject`. Options include QQQ.
#' @param annoName QQQ ??? name to subset with annotations
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
#' @param ArchRProj An `ArchRProject` object.
#' @param motifSet The motif set to be used for annotation. Options include: (i) "JASPAR2016", which gives the 2016 version of JASPAR motifs, (ii) "JASPAR2018", which gives the 2018 version of JASPAR motifs, or (iii) one of "human", "mouse", "encode", or "homer" which gives the corresponding motif sets from the chromVAR package. 
#' @param name QQQ of annotations to store as in `ArchRProject`
#' @param species QQQ The name of the species relevant to the supplied `ArchRProject`. This is used for QQQ. By default, this function will attempt to guess the species based on the value from `getGenome()`.
#' @param collection QQQ If one of the JASPAR motif sets is used via `motifSet`, this parameter allows you to indicate the JASPAR collection to be used. Possible options include "CORE", QQQ.
#' @param cutOff QQQ The QQQhypergeometric p-value cutoff to be used for motif search (see the `motimatchr` package for more information).
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
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix QQQ The name of the data matrix as stored in the ArrowFiles of the `ArchRProject`. Options include "TileMatrix", "GeneScoreMatrix", QQQ.
#' @param select QQQ select a specific name with grep
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
#' @param name The file name to be used for the output PDF file.
#' @param width The width in inches to be used for the output PDF file.
#' @param height The height in inches to be used for the output PDF.
#' @param ArchRProj An `ArchRProject` object.
#' @param addDOC A boolean variable that determines whether to add the date of creation to end of the PDF file name. This is useful for preventing overwritting of old plots.
#' @param useDingbats A boolean variable that determines wheter to use dingbats characters for plotting points.
#' @param plotList QQQ A `list` of plots to be printed to the output PDF file. Each element of `plotList` should be a QQQ format object.
#' @param useSink QQQ ???
#' @param ... additional args to pdf
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

}

#' Get Tutorial Data For ArchR
#' 
#' This function will download data for a given tutorial and return the input files required for ArchR
#' 
#' @param tutorial The name of the available tutorial for which to retreive the tutorial data. Options are "Hematopoiesis", "PBMC", "FreshFrozen". "Hematopoiesis" refers to QQQ. "PBMC" refers to QQQ. "FreshFrozen" refers to QQQ data from Granja et al. Nature Biotechnology 2019.
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


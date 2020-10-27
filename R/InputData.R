#' Get Relevant Data For ArchR Tutorials
#' 
#' This function will download data for a given tutorial and return the input files required for ArchR.
#' 
#' @param tutorial The name of the available tutorial for which to retreive the tutorial data. Currently, the only available option is "Hematopoiesis".
#' "Hematopoiesis" is a small scATAC-seq dataset that spans the hematopoieitic hierarchy from stem cells to differentiated cells.
#' This dataset is made up of cells from peripheral blood, bone marrow, and CD34+ sorted bone marrow.
#' @param threads The number of threads to be used for parallel computing.
#' @export
getTutorialData <- function(
  tutorial = "hematopoiesis", 
  threads = getArchRThreads()
  ){

  #Validate
  .validInput(input = tutorial, name = "tutorial", valid = "character")
  .validInput(input = threads, name = "threads", valid = c("integer"))
  #########

  if(tolower(tutorial) %in% c("heme","hematopoiesis")){
    
    if(!dir.exists("HemeFragments")){

      filesUrl <- c(
        "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz",
        "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz",
        "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz"
      )

      dir.create("HemeFragments", showWarnings = FALSE)

      downloadFiles <- .safelapply(seq_along(filesUrl), function(x){
        download.file(
          url = filesUrl[x], 
          destfile = file.path("HemeFragments", basename(filesUrl[x]))
        )        
      }, threads = min(threads, length(filesUrl)))

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

#' Get PBMC Small Test Fragments
#' 
#' This function will download fragments for a small PBMC test dataset (2k Cells) spanning chr1 and 2 (~20MB).
#' 
#' @export
getTestFragments <- function(x){
  if(!file.exists("PBMCSmall.tsv.gz")){
    download.file(
      url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/PBMCSmall.tsv.gz",
      destfile = "PBMCSmall.tsv.gz"
    )
  }
  addArchRGenome("hg19test")
  c("PBMC" = "PBMCSmall.tsv.gz")
}

#' Get PBMC Small Test Project
#' 
#' This function will download an ArchRProject for a small PBMC test dataset (2k Cells) spanning chr1 and 2 (~2-300MB).
#' 
#' @export
getTestProject <- function(){
  if(!dir.exists("PBMCSmall")){
    download.file(
      url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/PBMCSmall.zip",
      destfile = "PBMCSmall.zip"
    )
    unzip("PBMCSmall.zip", exdir = getwd())
    file.remove("PBMCSmall.zip")
  }
  addArchRGenome("hg19test")
  loadArchRProject("PBMCSmall")
}

#' Get Input Files from paths to create arrows
#' 
#' This function will look for fragment files and bam files in the input paths and return the full path and sample names.
#' 
#' @param paths A character vector of paths to search for usable input files.
#' @export
getInputFiles <- function(
  paths = NULL
  ){ 

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
getValidBarcodes <- function(
  csvFiles = NULL, 
  sampleNames = NULL
  ){

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

##########################################################################################
# Cell Filtering Methods
##########################################################################################

#' Subset cells in an ArchRProject.
#' 
#' This function returns an ArchRProject object that contains a specified subset of cells.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param cellNames A character vector of `cellNames` that will be subsetted of the current `ArchRProject`.
#' @export
subsetCells <- function(
  ArchRProj = NULL, 
  cellNames = NULL
  ){  

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character"))

  ccd <- getCellColData(ArchRProj)
  if(!all(cellNames %in% rownames(ccd))){
    stop("Not all cellNames in ArchRProject! Please re-check your input!")
  }

  ArchRProj@cellColData <- ccd <- ccd[cellNames, , drop = FALSE]
  
  ArchRProj
  
}

#' Filter Doublets From an ArchRProject
#'
#' This function will filter doublets from an ArchRProject after addDoubletScores() has been run.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param cutEnrich The minimum numeric cutoff for `DoubletEnrichment`. This number is equivalent to the number of simulated
#' doublets identified as a nearest neighbor to the cell divided by the expected number given a random uniform distribution.
#' @param cutScore The minimum numeric cutoff for `DoubletScore` which represents the `-log10(binomial adjusted p-value)` for the `DoubletEnrichment`.
#' @param filterRatio The maximum ratio of predicted doublets to filter based on the number of pass-filter cells.
#' For example, if there are 5000 cells, the maximum would be `filterRatio * 5000^2 / (100000)` (which simplifies to `filterRatio * 5000 * 0.05`).
#' This `filterRatio` allows you to apply a consistent filter across multiple different samples that may have different
#' percentages of doublets because they were run with different cell loading concentrations.
#' The higher the `filterRatio`, the greater the number of cells potentially removed as doublets.
#' @export
filterDoublets <- function(ArchRProj = NULL, cutEnrich = 1, cutScore = -Inf, filterRatio = 1){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  .validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  .validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))

  if(any(grepl("filterDoublets", names(ArchRProj@projectSummary)))){
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }

  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))

  cellsFilter <- lapply(splitDF, function(y){

    x <- df[y, ,drop = FALSE]

    n <- nrow(x)

    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]
    
    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    } 
    
    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    } 

    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    }else{
      NULL
    }

  }) %>% unlist(use.names=FALSE)

  message("Filtering ", length(cellsFilter), " cells from ArchRProject!")
  tabRemove <- table(df[cellsFilter,]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for(i in seq_along(samples)){
    if(!is.na(tabRemove[samples[i]])){
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1),"%)")
    }else{
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }

  if(length(cellsFilter) > 0){
    
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]

  }
  
  ArchRProj <- addProjectSummary(ArchRProj = ArchRProj, name = "filterDoublets", 
    summary = c(cutEnrich = cutEnrich, cutScore = cutScore, filterRatio = filterRatio))

  ArchRProj

}


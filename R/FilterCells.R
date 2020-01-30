##########################################################################################
# Cell Filtering Methods
##########################################################################################

#' Filter cells in an ArchRProject
#' 
#' This function returns an ArchRProject object that has been filtered to remove cells that do not pass the filter critera suppplied in filterList.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param filterList A list of filters based on `cellColData` to apply when filtering cells. 
#' Format should be a named list where the name corresponds to a column name in `cellColData` and the value corresponds to the filter criteria. 
#' If numeric, a lower threshold is expected, below which cells are filtered (a higher threshold is optional and can be included using a numeric vector). 
#' If a chatacter vector, only the rows in `cellColData` with `rownames` corresponding to the supplied character values are kept. 
#' If a list or simpleList, the user can additionally supply filters that are applied to all samples or a subset of specified samples. 
#' For example, to apply a filter to all samples: list("TSSEnrichment" = c(4,25)).
#' Similarly, to aply a filter only to specific samples: list("TSSEnrichment" = list("Sample1" = c(4, 25), "Sample2" = c(5, 25))).
#' @param ... additional params
#' @export
filterCells <- function(
  ArchRProj = NULL, 
  filterList = NULL,
  ...
  ){  

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = filterList, name = "filterList", valid = c("list"))

  ccd <- getCellColData(ArchRProj)

  cellsPF <- lapply(seq_along(filterList), function(x){

    vx <- getCellColData(ArchRProj, names(filterList[x]), drop = TRUE)

    if(inherits(filterList[[x]], "numeric")){

      cutLow <- filterList[[x]][1]
      cutHigh <- if(is.na(filterList[[x]][2])) Inf else filterList[[x]][2]
      idx <- rownames(ccd)[vx >= cutLow & vx <= cutHigh]

    }else if(inherits(filterList[[x]], "character")){

      idx <- rownames(ccd)[vx %bcin% filterList[[x]]]
      
    }else if(inherits(filterList[[x]], "list") | inherits(filterList[[x]], "SimpleList")){

      idx <- lapply(seq_along(filterList[[x]]), function(y){

        if(names(filterList[[x]][y]) %ni% ccd$Sample){
          stop(names(filterList[[x]][y]), " is not in sampleNames of ArchR Project")
        }

        ccdy <- ccd[BiocGenerics::which(ccd$Sample == names(filterList[[x]][y])),]
        vy <- getCellColData(ArchRProj, names(filterList[x]), drop = TRUE)

        if(inherits(filterList[[x]][[y]], "numeric")){

          cutLow <- filterList[[x]][[y]][1]
          cutHigh <- if(is.na(filterList[[x]][[y]][2])) Inf else filterList[[x]][[y]][2]
          rownames(ccdy)[vy >= cutLow & vy <= cutHigh]

        }else if(inherits(filterList[[x]][[y]], "character")){

          rownames(ccdy)[vy %bcin% filterList[[x]][[y]]]

        }else{

          stop(names(filterList[x]), " is not a numeric, character or list!")

        }


      }) %>% Reduce("c", .)
      
    }else{

      stop(names(filterList[x]), " is not a numeric, character or list!")

    }

    idx

  }) %>% Reduce("intersect", .)


  if(length(cellsPF) == 0){
    stop("0 Cells passing filter please consider less stringent thresholds!")
  }

  ArchRProj@cellColData <- ccd[cellsPF,]
  ArchRProj
  
}

#' Filter plot for cells in an ArchRProject
#' 
#' This function plots a list of attributes with filter criteria to visualize which cells would pass filter.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param filterList A list of filters based on `cellColData` to apply when filtering cells. 
#' Format should be a named numeric list where the name corresponds to a column name in `cellColData` 
#' and the value corresponds to the filter criteria. A lower threshold is expected, below which 
#' cells are filtered (a higher threshold is optional and can be passed as a numeric vector along with the lower threshold).
#' Only the first 2 filters will be plotted.
#' @param sampleNames The sample names corresponding to the subset of samples to plot. If `NULL`, all samples are included.
#' @param ... additional params
#' @export
filterPlot <- function(
  ArchRProj = NULL, 
  filterList = NULL, 
  sampleNames = NULL,
  ...
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = filterList, name = "filterList", valid = c("list"))
  .validInput(input = sampleNames, name = "sampleNames", valid = c("character", "null"))

  ccd <- getCellColData(ArchRProj, unique(c("Sample", names(filterList))))
  if(!is.null(sampleNames)){
    ccd <- ccd[ccd$Sample %in% sampleNames,,drop=TRUE]
  }
  stopifnot(nrow(ccd)!=0)

  cutoffs <- lapply(seq_along(filterList), function(x){
      if(inherits(filterList[[x]], "numeric")){
        cutLow <- filterList[[x]][1]
        cutHigh <- if(is.na(filterList[[x]][2])) Inf else filterList[[x]][2]
      }else{
        stop(names(filterList[x]), " is not a numeric or list!")
      }
      data.frame(row.names=names(filterList)[x], cutLow = cutLow, cutHigh = cutHigh)
  }) %>% Reduce("rbind",.)

  if(nrow(cutoffs) == 1){

    ggViolin(
      x = ccd$Sample, 
      y = as.numeric(ccd[,rownames(cutoffs)[1]]), 
      points = TRUE, 
      ylabel = rownames(cutoffs)[1], 
      ...) + geom_hline(yintercept = c(cutoffs[1,1],cutoffs[1,2]), lty = "dashed")

  }else if(nrow(cutoffs) > 1){

    if(nrow(cutoffs) > 2){
      message("Only plotting first 2 cutoffs from filterList!")
    }

    ggPoint(
      x = as.numeric(ccd[,rownames(cutoffs)[1]]),
      y = as.numeric(ccd[,rownames(cutoffs)[2]]),
      xlabel = rownames(cutoffs)[1],
      ylabel = rownames(cutoffs)[2],
      colorDensity = TRUE,
      ...) + geom_vline(xintercept = c(cutoffs[1,1],cutoffs[1,2]), lty = "dashed") + 
      geom_hline(yintercept = c(cutoffs[2,1],cutoffs[2,2]), lty = "dashed")

  }else{

    stop("No Cutoffs Found!")

  }

}

#' Filter Doublets From an ArchRProject
#'
#' This function wil filter doublets from an ArchRProject after addDoubletScores has been run
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param cutEnrich The minimum numeric cutoff for `DoubletEnrichment`. This number is equivalent to the number of simulated doublets identified as a nearest neighbor to the cell divided by the expected number given a random uniform distribution.
#' @param cutScore The minimum numeric cutoff for `DoubletScore` which represents the -log10(binomial adjusted p-value) for the `DoubletEnrichment`.
#' @param filterRatio The maximum ratio of predicted doublets to filter based on the number of pass-filter cells.
#' For example, if there are 5,000 cells the maximum would be filterRatio * 5,000^2 / (100,000) (which simplifies to filterRatio * 5000 * 0.05).
#' This `filterRatio` allows you to apply a consistent filter across multiple different samples that
#' may have different percentages of doublets because they were run with different cell loading concentrations.
#' The higher the `filterRatio`, the greater the number of cells potentially removed as doublets.
#' @param ... additional params
#' @export
filterDoublets <- function(ArchRProj, cutEnrich = 1, cutScore = -Inf, filterRatio = 1, ...){

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
      x <- x[x$DoubletEnrichment >= cutEnrich, ]
    } 
    
    if(!is.null(cutScore)){
      x <- x[x$DoubletScore >= cutScore, ]
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
    summary = c("Date" = Sys.time(), cutEnrich = cutEnrich, cutScore = cutScore, filterRatio = filterRatio))

  ArchRProj

}


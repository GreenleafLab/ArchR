##########################################################################################
# Cell Filtering Methods
##########################################################################################

#' Filter cells in an ArchRProject
#' 
#' This function returns an ArchRProject object that has been filtered to remove cells that do not pass the filter critera suppplied in filterList.
#' 
#' @param ArchRProj An `ArchRProject` object
#' @param filterList QQQ A list of filters based on cellColData to apply when filtering cells. Format should be a named list where the name corresponds to a column name in cellColData and the value corresponds to the filter criteria. If numeric, a lower threshold is expected, below which cells are filtered (a higher threshold is optional). If a chatacter vector, rows in cellColData corresponding to the supplied character values are QQQkept. If a list or simpleList, QQQ.
#' @param ... additional params
#' @export
filterCells <- function(ArchRProj, filterList, ...){  
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
#' @param ArchRProj An `ArchRProject` object
#' @param filterList QQQ A list of filters based on cellColData to apply when filtering cells. Format should be a named numeric list where the name corresponds to a column name in cellColData and the value corresponds to the filter criteria. A lower threshold is expected, below which cells are filtered (a higher threshold is optional). Only the first 2 filters will be plotted.
#' @param sampleNames The sample names corresponding to the subset of samples to plot. If NULL, all samples are included.
#' @param ... additional params
#' @export
filterPlot <- function(ArchRProj, filterList, sampleNames = NULL, ...){

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
#' @param ArchRProj An `ArchRProject` object
#' @param cutEnrich QQQ The minimum numeric cutoff for `DoubletEnrichment`. This number is equivalent to the number of simulated doublets identified as a nearest neighbor to the cell divided by the expected number given a random uniform distribution.
#' @param cutScore QQQ The minimum numeric cutoff for `DoubletScore` which represents the -log10(binomial adjusted p-value) for QQQ(what is the pvalue from?).
#' @param filterRatio QQQ THIS IS VERY CONFUSING TO ME. WHY ISNT THIS JUST THE FRACTION OF YOUR CELLS THAT YOU ARE WILLING TO LOSE TO DOUBLETS? The maximum ratio of predicted doublets to filter based on the number of pass-filter cells. If there are 10,000 cells the maximum would be filterRatio * 10,000^2 / (100,000).
#' @param ... additional params
#' @export
filterDoublets <- function(ArchRProj, cutEnrich = 1, cutScore = -Inf, filterRatio = 1, ...){
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

  if(length(cellsFilter) > 0){
    
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]

  }
  
  ArchRProj

}


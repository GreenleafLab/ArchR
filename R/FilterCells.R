#' Filter cells in an ArchRProject
#' 
#' This function plots a list of filters to see which cells would pass filter.
#' 
#' @param ArchRProj an ArchR Project object
#' @param filterList list of filters for filtering cells from cellColData
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

#' Filter Plot for cells in an ArchRProject
#' 
#' This function plots a list of filters to see which cells would pass filter.
#' 
#' @param ArchRProj an ArchR Project object
#' @param filterList list of filters for filtering cells from cellColData (up to 2 will be plotted)
#' @param sampleNames specific samples to include
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


#' Filter Doublets From an ArchR Project
#'
#' This function wil filter doublets from an ArchRProject after addDoubletScores has been ran
#'
#' @param ArchRProj an ArchR Project object
#' @param cutEnrich minimum cutoff for doubletEnrichment which represents number of simulated doublets nearest a cell over the expected if uniform.
#' @param cutScore minimum cutoff for doubletScore which represents -log10 binomial adjusted p-value.
#' @param filterRatio filter ratio for max number of inferred doublets to remove based on the number of cells PF. If there are 10,000 cells the maximum would be filterRatio * 10,000^2 / (1000 * 100).
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








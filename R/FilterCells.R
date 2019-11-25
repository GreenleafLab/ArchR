#' Extend Filter then Normalize Scores for Summits
#' @param df dataframe where first column is sample names 2nd column is group information and 3rd column is MACS2 summit files
#' @param genome mm9, hg19 character or BSgenome object
#' @param blacklist regions to blacklist
#' @param extend how to extend summits (summit +- extend)
#' @param scorePerMillion normalized Score-per-million minimum to keep
#' @param selectionRules string with a formula containing n (majority = (n+1)/2, multiple samples = 2)
#' @export
filterCells <- function(ArchRProj, filterList){
  
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

#' Extend Filter then Normalize Scores for Summits
#' @param df dataframe where first column is sample names 2nd column is group information and 3rd column is MACS2 summit files
#' @param genome mm9, hg19 character or BSgenome object
#' @param blacklist regions to blacklist
#' @param extend how to extend summits (summit +- extend)
#' @param scorePerMillion normalized Score-per-million minimum to keep
#' @param selectionRules string with a formula containing n (majority = (n+1)/2, multiple samples = 2)
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


#' Extend Filter then Normalize Scores for Summits
#' @param df dataframe where first column is sample names 2nd column is group information and 3rd column is MACS2 summit files
#' @param genome mm9, hg19 character or BSgenome object
#' @param blacklist regions to blacklist
#' @param extend how to extend summits (summit +- extend)
#' @param scorePerMillion normalized Score-per-million minimum to keep
#' @param selectionRules string with a formula containing n (majority = (n+1)/2, multiple samples = 2)
#' @export
filterDoublets <- function(ArchRProj, cutEnrich = 1, cutScore = -Inf, fs = 1){


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

      head(rownames(x), fs * n * (n / 100000))

    }else{
      NULL
    }

  }) %>% unlist(use.names=FALSE)

  if(length(cellsFilter) > 0){
    
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]

  }
  
  ArchRProj

}








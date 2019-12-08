#' Extend Filter then Normalize Scores for Summits
#'
#'
#'
#' @param df dataframe where first column is sample names 2nd column is group information and 3rd column is MACS2 summit files
#' @param genome mm9, hg19 character or BSgenome object
#' @param blacklist regions to blacklist
#' @param extend how to extend summits (summit +- extend)
#' @param scorePerMillion normalized Score-per-million minimum to keep
#' @param selectionRules string with a formula containing n (majority = (n+1)/2, multiple samples = 2)
#' @export
FilterCells <- function(ArchRProj, filterList){
  
  ccd <- getCellColData(ArchRProj)

  cellsPF <- lapply(seq_along(filterList), function(x){

    if(names(filterList[x]) %ni% colnames(ccd)){
      stop(names(filterList[x]), " is not in colnames of cellColData")
    }

    vx <- ccd[,names(filterList[x])]

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
        vy <- ccdy[,names(filterList[x])]

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




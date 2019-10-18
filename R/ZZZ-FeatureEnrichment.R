#------------------------------------------------------
# Hypergeometric Testing
#------------------------------------------------------

#' Feature Matches Over Representation
#' 
#' This function takes in genes matches combined with a vector of peaks in target and peaks in bdg and returns hypergeometric minus log 10 pvalues
#' @param genes motif similarity matrix used for labeling family info default is null
#' @param compare vector of compare peaks idx which will be used for hypergeometric
#' @param background vector of background peaks idx which will be used for hypergeometric
#' @export
#'
featureEnrichment <- function(featureMatches, compare, background){
  
  suppressPackageStartupMessages(require(Matrix))
  
  #Prep
  stopifnot(length(grep("Matches", assayNames(featureMatches))) == 1)
  matches <- getAssay(featureMatches, grep("Matches", assayNames(featureMatches), value = TRUE))
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)
  
  pOut <- data.frame(feature = colnames(matches),
                      CompareFrequency = matchCompareTotal,
                      nCompare = nrow(matchCompare),
                      CompareProportion = matchCompareTotal/nrow(matchCompare),
                      BackgroundFrequency = matchBackgroundTotal,
                      nBackground = nrow(matchBackground),
                      BackgroundProporition = matchBackgroundTotal/nrow(matchBackground))
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #mlog10phyper
  pOut$mlog10phyper <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
                 pOut$BackgroundFrequency[x], # Number of all successes in background
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
                 pOut$nCompare[x], # Number that were drawn
                 lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(3)
  pOut$FDR <- p.adjust(10^-matrixStats::rowMins(as.matrix(data.frame(pOut$mlog10phyper, 250))), method = "BH")
  pOut$mlog10FDR <- -log10(pOut$FDR)
  pOut <- pOut[order(pOut$mlog10phyper, decreasing = TRUE), c(1, 10, 2:8, 9, 11)]

  return(pOut)
}





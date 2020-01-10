##########################################################################################
# Marker Feature Methods
##########################################################################################

#' Identify Marker Features for each cell grouping
#' 
#' This function will identify features that are definitional of each provided cell grouping where possible
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for marker feature identification.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in `cellColData`. This limits the groups used to perform marker feature identification.
#' @param bgdGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in `cellColData` to be used for background calculations in marker feature identification.
#' @param useMatrix The name of the matrix to be used for performing differential analyses. Options include "GeneScoreMatrix", "PeakMatrix", etc.
#' @param bias A character vector indicating the potential bias variables as a function (i.e. c("TSSEnrichment", "log10(nFrags)")) to account for in selecting a matched null group for marker feature identification. These should be column names from `cellColData`.
#' @param normBy The name of a numeric column in `cellColData` that should be normalized across cells (i.e. "ReadsInTSS") prior to performing marker feature identification.
#' @param testMethod The name of the pairwise test method to use in comparing cell groupings to the null cell grouping during marker feature identification. Valid options include "wilcoxon", "ttest", and "binomial".
#' @param maxCells The maximum number of cells to consider from a single cell group when performing marker feature identification.
#' @param scaleTo Normalization depth to center normalization to in normBy (default is 10,000).
#' @param threads The number of threads to be used for parallel computing.
#' @param k The number of nearby cells to use for selecting biased-matched background while accounting for bgdGroups proportions.
#' @param bufferRatio The buffering ratio of cells to enable best biased-matched background while accounting for bgdGroups proportions.
#' @param binarize A boolean value indicating whether to binarize the matrix prior to differential testing.
#' @param useSeqnames A character vector that indicates which seqnames should be used in marker feature identification. Features from seqnames that are not listed will be ignored. 
#' @param method The name of the method to be used for marker feature identification. Valid options are "ArchR" which will use the default ArchR method or "Venice" which will use the `Signac::VeniceMarker()` fucntion.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... additional args
#' @export
markerFeatures <- function(
    ArchRProj = NULL,
    groupBy = "Clusters",
    useGroups = NULL,
    bgdGroups = NULL,
    useMatrix = "GeneScoreMatrix",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    normBy = NULL,
    testMethod = "wilcoxon",
    maxCells = 500,
    scaleTo = 10^4,
    threads = 1,
    k = 100,
    bufferRatio = 0.8,
    binarize = FALSE,
    useSeqnames = NULL,
    method = "ArchR",
    verboseHeader = TRUE,
    verboseAll = FALSE,
    ...
    ){
  
  args <- append(args, mget(names(formals()),sys.frame(sys.nframe())))

  if(tolower(method) == "archr"){
 
    out <- do.call(.MarkersSC, args)
    
  }else if(tolower(method) == "venice"){

    .requirePackage("signac", installInfo = 'devtools::install_github("bioturing/signac")')

  }else{
  
    stop("Input Method Not Available!")
  
  }

  args$ArchRProj <- NULL
  metadata(out)$Params <- args

  return(out)

}

##################################################################################################
# Single Cell Implementation!
##################################################################################################
.MarkersSC <- function(
    ArchRProj = NULL,
    groupBy = "Clusters",
    useGroups = NULL,
    bgdGroups = NULL,
    normBy = NULL,
    maxCells = 500,
    scaleTo = 10^4,
    bufferRatio = 0.8,
    useSeqnames = NULL,
    bias = NULL,
    k = 100,
    threads = 8,
    binarize = FALSE,
    testMethod = "wilcoxon",
    useMatrix = "GeneScoreMatrix",
    markerParams = list(),
    verboseHeader = TRUE,
    verboseAll = FALSE,
    ...
  ){

    tstart <- Sys.time()
    
    #####################################################
    # Feature Info
    #####################################################
    ArrowFiles <- getArrowFiles(ArchRProj)
    featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
    matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))

    if(!is.null(useSeqnames)){
      if(length(useSeqnames) == 1){
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }else{
        if(matrixClass == "Sparse.Assays.Matrix"){
          stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames!\nPlease specify 1 seqname in useSeqnames to continue!\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available seqnames for input!")
        }
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }
    }else{
      if(matrixClass == "Sparse.Assays.Matrix"){
        stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames!\nPlease specify 1 seqname in useSeqnames to continue!\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available seqnames for input!")
      }
    }
    if(!(nrow(featureDF) > 1)){
      stop("Less than 1 feature is remaining in featureDF please check input!")
    }

    #####################################################
    # Match Bias Groups
    #####################################################
    .messageDiffTime("Matching Known Biases", tstart, addHeader = verboseAll)
    groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    colDat <- getCellColData(ArchRProj)

    matchObj <- .matchBiasCellGroups(
      input = colDat, 
      groups = groups,
      useGroups = useGroups,
      bgdGroups = bgdGroups,
      bias = bias,
      k = k,
      n = maxCells
    )

    #####################################################
    # Pairwise Test Per Seqnames
    #####################################################
    mColSums <- suppressMessages(.getColSums(ArrowFiles, seqnames = featureDF$seqnames@values, useMatrix = useMatrix, threads = threads))
    
    if(is.null(normBy)){
      if(tolower(testMethod) == "binomial"){
        normFactors <- NULL
      }else{
        if(tolower(useMatrix) %in% c("tilematrix", "peakmatrix")){
          normBy <- "ReadsInTSS"
          normFactors <- getCellColData(ArchRProj, normBy, drop=FALSE)
          normFactors[,1] <- median(normFactors[,1]) / normFactors[,1]
        }else{
          normFactors <- scaleTo / mColSums
          normFactors <- DataFrame(normFactors)
        }
      }
    }else{
      if(tolower(normBy) == "none"){
        normFactors <- NULL
      }else{
        normFactors <- scaleTo / mColSums
        normFactors <- DataFrame(normFactors)
      }
    }

    if(!is.null(normFactors)){
      normFactors[,1] <- normFactors[,1] * (scaleTo / median(normFactors[names(mColSums), 1] * mColSums))
    }
    
    diffList <- .safelapply(seq_along(matchObj[[1]]), function(x){
      .messageDiffTime(sprintf("Computing Pairwise Tests (%s of %s)", x, length(matchObj[[1]])), tstart, addHeader = verboseAll)
      .testMarkerSC(
          ArrowFiles = ArrowFiles,
          matchObj = matchObj, 
          group = names(matchObj[[1]])[x], 
          testMethod = testMethod, 
          threads = 1, 
          useMatrix = useMatrix,
          featureDF = featureDF,
          normFactors = normFactors,
          binarize = binarize
      )
    }, threads = threads)

    .messageDiffTime("Completed Pairwise Tests", tstart, addHeader = TRUE)

    #####################################################
    # Summarize Output
    #####################################################
    if(tolower(testMethod) == "wilcoxon"){
      pse <- SummarizedExperiment::SummarizedExperiment(
          assays = 
            SimpleList(
              Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
              Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
              FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
              AUC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$auc)) %>% Reduce("cbind",.),
              MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else if(tolower(testMethod) == "ttest"){
      pse <- SummarizedExperiment::SummarizedExperiment(
          assays = 
            SimpleList(
              Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
              Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
              Variance = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var1)) %>% Reduce("cbind",.),
              FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
              MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.),
              VarianceBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var2)) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else if(tolower(testMethod) == "binomial"){
      pse <- SummarizedExperiment::SummarizedExperiment(
          assays = 
            SimpleList(
              Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
              Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
              FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
              MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else{
      stop("Error Unrecognized Method!")
    }
    colnames(pse) <- names(matchObj[[1]])

    metadata(pse)$MatchInfo <- matchObj

    return(pse)

}

.testMarkerSC <- function(ArrowFiles, matchObj, group = NULL, testMethod = "ttest", useMatrix,
  threads = 1, featureDF, binarize = FALSE, normFactors = NULL){

  matchx <- matchObj[[1]][[group]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]
  
  if(!is.null(normFactors)){
    cellNF <- normFactors[cellsx,1]
    bgdNF <- normFactors[bgdx,1]
  }

  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  seqnames <- unique(featureDF$seqnames)

  pairwiseDF <- lapply(seq_along(seqnames), function(y){

    featureDFy <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames[y]), ]

    scMaty <- suppressMessages(.getPartialMatrix(
      ArrowFiles, 
      featureDF = featureDFy, 
      threads = threads, 
      useMatrix = useMatrix,
      cellNames = c(cellsx, bgdx),
      progress = FALSE
    ))
    rownames(scMaty) <- rownames(featureDFy)

    if(binarize){
      scMaty@x[scMaty@x > 0] <- 1
    }

    args <- list()
    if(!is.null(normFactors)){
      args$mat1 <- Matrix::t(Matrix::t(scMaty[, cellsx, drop = FALSE]) * cellNF)
      args$mat2 <- Matrix::t(Matrix::t(scMaty[, bgdx, drop = FALSE]) * bgdNF)
    }else{
      args$mat1 <- scMaty[, cellsx, drop = FALSE]
      args$mat2 <- scMaty[, bgdx, drop = FALSE]
    }

    if(tolower(testMethod) == "wilcoxon"){

      .suppressAll(do.call(.sparseMatWilcoxon, args))
    
    }else if(tolower(testMethod) == "ttest"){
    
      .suppressAll(do.call(.sparseMatTTest, args))
    
    }else if(tolower(testMethod) == "binomial"){

      if(!is.null(normFactors)){
        stop("Normfactors cannot be used with a binomial test!")
      }

      if(!binarize){
        stop("Binomial test requires binarization!")
      }
    
      .suppressAll(do.call(.sparseMatBinomTest, args))

    }else{
    
      stop("Error Unrecognized Method!")
    
    }

  }) %>% Reduce("rbind", .)

  idxFilter <- rowSums(pairwiseDF[,c("mean1","mean2")]) != 0
  pairwiseDF$fdr <- NA
  pairwiseDF$fdr[idxFilter] <- p.adjust(pairwiseDF$pval[idxFilter], method = "fdr")
  pairwiseDF <- pairwiseDF[rownames(featureDF), , drop = FALSE]
  pairwiseDF
  
}

#Wilcoxon Row-wise two matrices
.sparseMatWilcoxon <- function(mat1, mat2){
  
  .requirePackage("presto", installInfo = 'devtools::install_github("immunogenomics/presto")')
  df <- wilcoxauc(cbind(mat1,mat2), c(rep("Top", ncol(mat1)),rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]

  #Sparse Row Sums
  m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
  m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
  offset <- 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
  log2FC <- log2((m1 + offset)/(m2 + offset))
  log2Mean <- log2(((m1+offset) + (m2+offset)) / 2)

  out <- data.frame(
    log2Mean = log2Mean,
    log2FC = log2FC,
    fdr = df$padj, 
    pval = df$pval, 
    mean1 = Matrix::rowMeans(mat1, na.rm=TRUE), 
    mean2 = Matrix::rowMeans(mat2, na.rm=TRUE), 
    n = ncol(mat1),
    auc = df$auc
  )

  return(out)

}

#T-Test Row-wise two matrices
.sparseMatTTest <- function(mat1, mat2, m0 = 0){
    
    #Get Population Values
    n1 <- ncol(mat1)
    n2 <- ncol(mat2)
    n <- n1 + n2
    
    #Sparse Row Means
    m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
    m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
    
    #Sparse Row Variances
    v1 <- ArchR:::computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
    v2 <- ArchR:::computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
    
    #Calculate T Statistic
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
    tstat <- (m1-m2-m0)/se
    pvalue <- 2*pt(-abs(tstat), n - 2)
    fdr <- p.adjust(pvalue, method = "fdr")
    
    #Sparse Row Sums
    m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
    m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
    offset <- 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
    log2FC <- log2((m1 + offset)/(m2 + offset))
    log2Mean <- log2(((m1+offset) + (m2+offset)) / 2)

    out <- data.frame(
      log2Mean = log2Mean,
      log2FC = log2FC,
      fdr = fdr, 
      pval = pvalue, 
      mean1 = m1, 
      mean2 = m2, 
      var1 = v2,
      var2 = v2,
      n = n1
    )
    return(out)
}

#Binomial Test Row-wise two matrices
.sparseMatBinomTest <- function(mat1, mat2){
  
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  
  #Sparse Row Stats
  s1 <- Matrix::rowSums(mat1, na.rm=TRUE)
  m1 <- s1 / n1
  s2 <- Matrix::rowSums(mat2, na.rm=TRUE)
  m2 <- s2 / n2
  
  #Combute Binom.test 2-sided
  pval <- unlist(lapply(seq_along(s1), function(x){
    if(m1[x] >= m2[x]){
      p <- pbinom(q = s1[x]-1, size = n1, prob = (max(c(s2[x],1))) / n2, lower.tail=FALSE, log.p = TRUE)
    }else{
      p <- pbinom(q = s2[x]-1, size = n2, prob = (max(c(s1[x],1))) / n1, lower.tail=FALSE, log.p = TRUE)
    }
    p <- min(c(2 * exp(p), 1)) #handle 2-sided test
    p
  }))
  fdr <- p.adjust(pval, method = "fdr", length(pval))
  
  #Sparse Row Sums
  m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
  m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
  offset <- 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
  log2FC <- log2((m1 + offset)/(m2 + offset))
  log2Mean <- log2(((m1+offset) + (m2+offset)) / 2)

  out <- data.frame(
    log2Mean = log2Mean,
    log2FC = log2FC,
    fdr = fdr, 
    pval = pval, 
    mean1 = m1, 
    mean2 = m2, 
    n = n1
  )

  return(out)

}


.matchBiasCellGroups <- function(input, groups, useGroups, bgdGroups, bias, k = 100, n = 500, bufferRatio = 0.8){

  #Summary Function
  .summarizeColStats <- function(m, name = NULL){
    med <- apply(m, 2, median)
    mean <- colMeans(m)
    sd <- apply(m, 2, sd)
    loQ <- apply(m, 2, function(x) quantile(x, 0.25))
    hiQ <- apply(m, 2, function(x) quantile(x, 0.75))
    summaryDF <- t(data.frame(
      median = med,
      mean = mean,
      sd = sd,
      lowerQuartile = loQ,
      upperQuartile = hiQ
    )) %>% data.frame
    colnames(summaryDF) <- colnames(m)
    if(!is.null(name)){
      summaryDF$name <- name
    }
    summaryDF
  }

  #Set Seed
  set.seed(1)
  
  #Make sure input is dataframe
  input <- data.frame(input)
  
  #Norm using input string ie log10(nfrags)
  inputNorm <- lapply(seq_along(bias), function(x){
    plyr::mutate(input, o=eval(parse(text=bias[x])))$o
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNorm) <- rownames(input)

  #Quantile Normalization
  inputNormQ <- lapply(seq_len(ncol(inputNorm)), function(x){
    .getQuantiles(inputNorm[,x])
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNormQ) <- rownames(input)

  #Add Colnames
  colnames(inputNorm) <- bias
  colnames(inputNormQ) <- bias

  if(is.null(useGroups)){
    useGroups <- gtools::mixedsort(unique(paste0(groups)))
  }

  if(is.null(bgdGroups)){
    bgdGroups <- gtools::mixedsort(unique(paste0(groups)))
  }

  stopifnot(all(useGroups %in% unique(paste0(groups))))
  stopifnot(all(bgdGroups %in% unique(paste0(groups))))

  #Get proportion of each group
  prob <- table(groups) / length(groups)
  bgdProb <- prob[which(names(prob) %in% bgdGroups)] / sum(prob[which(names(prob) %in% bgdGroups)])

  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  matchList <- lapply(seq_along(useGroups), function(x){
    
    setTxtProgressBar(pb,round(x*100/length(useGroups),0))

    #############
    # Organize
    #############
    groupx <- useGroups[x]
    idx <- which(names(bgdProb) == groupx)
    if(length(idx) > 0 & length(idx) != length(bgdProb)){
      bgdProbx <- bgdProb[-idx]/sum(bgdProb[-idx])
    }else{
      bgdProbx <- bgdProb
    }

    idF <- which(groups == groupx)
    idB <- which(groups %in% names(bgdProbx))

    knnx <- .computeKNN(inputNormQ[idB, ], inputNormQ[idF, ], k = k)
    sx <- sample(seq_len(nrow(knnx)), nrow(knnx))

    minTotal <- min(n, length(sx) * bufferRatio)
    nx <- sort(floor(minTotal * bgdProbx))
    
    ###############
    # ID Matching
    ###############
    idX <- c()
    idY <- c()
    it <- 0
    
    if(any(nx <= 0)){
      nx[which(nx <= 0)] <- Inf
      nx <- sort(nx)
    }

    while(it < length(sx) & length(idX) < minTotal){
      
      it <- it + 1
      knnit <- knnx[sx[it],]
      groupit <- match(groups[idB][knnit],names(nx))
      selectUnique <- FALSE
      selectit <- 0
      oit <- order(groupit)
      
      while(!selectUnique){
        selectit <- selectit + 1
        itx <- which(oit==selectit)
        cellx <- knnit[itx]
        groupitx <- groupit[itx]
        if(is.infinite(nx[groupitx])){
          if(selectit == k){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }else{
          if(cellx %ni% idY){
            selectUnique <- TRUE
          }
          if(selectit == k){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }
      }
      
      if(!is.na(itx)){
        idX <- c(idX, sx[it])
        idY <- c(idY, cellx)
        nx[groupitx] <- nx[groupitx] - 1
        if(any(nx <= 0)){
          nx[which(nx <= 0)] <- Inf
          nx <- sort(nx)
        }
      }

      if(all(is.infinite(nx))){
        it <- length(sx)
      }

    }

    #####################
    # Convert Back to Normal Indexing
    #####################
    idX <- seq_len(nrow(inputNormQ))[idF][idX]
    idY <- seq_len(nrow(inputNormQ))[idB][idY]

    #####################
    # Matching Stats Groups
    #####################
    estbgd <- sort(floor(minTotal * bgdProbx))
    obsbgd <- rep(0, length(estbgd))
    names(obsbgd) <- names(estbgd)
    tabGroups <- table(groups[idY])
    obsbgd[names(tabGroups)] <- tabGroups
    estbgdP <- round(100 * estbgd / sum(estbgd),3)
    obsbgdP <- round(100 * obsbgd / sum(obsbgd),3)

    #####################
    # Matching Stats Bias Norm Values
    #####################
    forBias <- .summarizeColStats(inputNorm[idX,], name = "foreground")
    bgdBias <- .summarizeColStats(inputNorm[idY,], name = "background")

    out <- list(
        cells = idX, 
        bgd = idY, 
        summaryCells = forBias, 
        summaryBgd = bgdBias, 
        bgdGroups = rbind(estbgd, obsbgd),
        bgdGroupsProbs = rbind(estbgdP, obsbgdP),
        corbgdGroups = cor(estbgdP, obsbgdP),
        n = length(sx), 
        p = it / length(sx),
        group = groupx
      )

    return(out)

  }) %>% SimpleList
  names(matchList) <- useGroups
  
  message("\n")

  outList <- SimpleList(
    matchbgd = matchList,
    info = SimpleList(
        cells = rownames(input),
        groups = groups,
        biasNorm = inputNorm,
        biasNormQ = inputNormQ
      )
  )
  
  return(outList)

}


####################################################################################################
# Applications of Markers!
####################################################################################################

#' Plot a Heatmap of Identified Marker Features
#' 
#' This function will plot a heatmap of the results from markerFeatures
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted in the heatmap. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param log2Norm A boolean value indicating whether a log2 transformation whould be performed on the values in `seMarker` prior to plotting. Should be set to `TRUE` for counts-based assays (but not assays like `DeviationsMatrix`).
#' @param scaleTo scale to prior to log2 Normalization, if log2Norm is FALSE this does nothing
#' @param scaleRows A boolean value that indicates whether the heatmap should display row-wise z-scores instead of raw values.
#' @param limits A numeric vector of two numbers that represent the lower and upper color limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies rownames from `seMarker` to be excluded from the heatmap.
#' @param pal A custom continuous palette (see paletteContinuous) used to override the continuous palette for the heatmap.
#' @param binaryClusterRows A boolean value that indicates whether a binary sorting algorithm should be used for fast clustering of heatmap rows.
#' @param labelMarkers A character vector listing the `rownames` of `seMarker` that should be labeled on the side of the heatmap.
#' @param labelTop A boolean value that indicates whether the top features for each column in `seMarker` should be labeled on the side of the heatmap.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMat A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param invert A boolean value that indicates whether the heatmap will be inverted ie when looking for down-regulated markers (Log2FC < 0) instead of up-regulated markers (Log2FC > 0).
#' @param ... additional args
#' @export
markerHeatmap <- function(
  seMarker,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2,2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  labelMarkers = NULL,
  labelTop = NULL,
  labelRows = FALSE,
  returnMat = FALSE,
  invert = FALSE,
  ...
  ){

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  #Now Get Values
  if(plotLog2FC){
    mat <- SummarizedExperiment::assays(seMarker)[["Log2FC"]]
  }else{
    mat <- SummarizedExperiment::assays(seMarker)[["Mean"]]
    if(log2Norm){
      mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
    }
    if(scaleRows){
      mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }
  }
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  
  idx <- which(rowSums(passMat, na.rm = TRUE) > 0 & matrixStats::rowVars(mat) != 0)
  mat <- mat[idx,]
  passMat <- passMat[idx,]

  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }

  #add rownames
  rd <- SummarizedExperiment::rowData(seMarker)[idx,]
  if(is.null(rd$name)){
    rn <- paste0(rd$seqnames,":",rd$start,"-",rd$end)
  }else{
    if(sum(duplicated(rd$name)) > 0){
      rn <- paste0(rd$seqnames,":",rd$name)
    }else{
      rn <- rd$name
    }
  }
  rownames(mat) <- rn
  rownames(passMat) <- rn

  #identify to remove
  if(!is.null(grepExclude) & !is.null(rownames(mat))){
    idx2 <- which(!grepl(grepExclude, rownames(mat)))
    mat <- mat[idx2,]
  }

  if(nrow(mat)==0){
    stop("No Makers Found!")
  }

  if(!is.null(labelTop)){
    spmat <- passMat / rowSums(passMat)
    idx2 <- lapply(seq_len(ncol(spmat)), function(x){
      head(order(spmat[,x], decreasing = TRUE), labelTop)
    }) %>% unlist %>% unique %>% sort
    mat <- mat[idx2,]
    labelRows <- TRUE
  }

  if(binaryClusterRows){
    if(invert){
      bS <- .binarySort(-mat, lmat = passMat[rownames(mat), colnames(mat)])
      mat <- -bS[[1]][,colnames(mat)]
    }else{
      bS <- .binarySort(mat, lmat = passMat[rownames(mat), colnames(mat)])
      mat <- bS[[1]][,colnames(mat)]
    }
    clusterRows <- FALSE
    clusterCols <- bS[[2]]
  }else{
    clusterRows <- TRUE
    clusterCols <- TRUE
  }

  if(!is.null(labelMarkers)){
    mn <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
    mn <- mn[mn > 0]
  }else{
    mn <- NULL
  }

  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }

  message(sprintf("Identified %s markers!", nrow(mat)))

  if(is.null(pal)){
    if(is.null(metadata(seMarker)$Params$useMatrix)){
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }else if(tolower(metadata(seMarker)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blue_yellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solar_extra", n = 100)
    }
  }

  if(invert){
    pal <- rev(pal)
  }

  ht <- .ArchRHeatmap(
    mat = mat,
    scale = FALSE,
    limits = c(min(mat), max(mat)),
    color = pal, 
    clusterCols = clusterCols, 
    clusterRows = clusterRows,
    labelRows = labelRows,
    labelCols = TRUE,
    customRowLabel = mn,
    showColDendrogram = TRUE,
    draw = FALSE,
    ...
  )

  if(returnMat){
    return(mat)
  }else{
    return(ht)
  }

}

########################################################################################################
# Helpers for Nice Heatmap with Bioconductors ComplexHeamtap
########################################################################################################

.ArchRHeatmap <- function(
  mat, 
  scale = FALSE,
  limits = c(min(mat), max(mat)),
  colData = NULL, 
  color = paletteContinuous(set = "solar_extra", n = 100),
  clusterCols = TRUE,
  clusterRows = FALSE,
  labelCols = FALSE,
  labelRows = FALSE,
  colorMap = NULL,
  useRaster = TRUE,
  rasterQuality = 5,
  split = NULL,
  fontsize = 6,
  colAnnoPerRow = 4,
  showRowDendrogram = FALSE,
  showColDendrogram = FALSE,
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.75,
  rasterDevice = "png",
  padding = 45,
  borderColor = NA,
  draw = TRUE,
  name = ""
  ){
  
  #Packages
  .requirePackage("ComplexHeatmap")
  .requirePackage("circlize")
  
  #Z-score
  if (scale) {
    message("Scaling Matrix...")
    mat <- .rowZscores(mat, limit = FALSE)
    name <- paste0(name," Z-Scores")
  }
  
  #Get A Color map if null
  if (is.null(colorMap)) {
    colorMap <- .colorMapAnno(colData)
  }
  
  #Prepare ColorMap format for Complex Heatmap
  if (!is.null(colData)){
    colData = data.frame(colData)
    colorMap <- .colorMapForCH(colorMap, colData) #change
    showLegend <- .checkShowLegend(colorMap[match(names(colorMap), colnames(colData))]) #change
  }else {
    colorMap <- NULL
    showLegend <- NULL
  }
  
  #Prepare Limits if needed
  breaks <- NULL
  if (!is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    breaks <- seq(min(limits), max(limits), length.out = length(color))
    color <- circlize::colorRamp2(breaks, color)
  }
  
  if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
    anno_check_version_rows <- ComplexHeatmap::anno_mark
    anno_check_version_cols <- ComplexHeatmap::anno_mark
  }else{
    anno_check_version_rows <- ComplexHeatmap::row_anno_link
    anno_check_version_cols <- ComplexHeatmap::column_anno_link
  }

  #Annotation Heatmap
  if(!is.null(colData) & !is.null(customColLabel)){
    message("Adding Annotations...")
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customRowLabel]
    }
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        ),
      link = anno_check_version_cols(
        at = customColLabel, labels = customColLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs)

    )
  }else if(!is.null(colData)){
    message("Adding Annotations...")
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        )
    )
  }else if(is.null(colData) & !is.null(customColLabel)){
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customRowLabel]
    }
    message("Adding Annotations...")
    ht1Anno <- HeatmapAnnotation(
      link = anno_check_version_cols(
        at = customColLabel, labels = customColLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs)
    )
  }else{
    ht1Anno <- NULL
  }

  message("Preparing Main Heatmap...")
  ht1 <- Heatmap(
    
    #Main Stuff
    matrix = mat,
    name = name,
    col = color, 
    
    #Heatmap Legend
    heatmap_legend_param = 
      list(color_bar = "continuous", 
           legend_direction = "horizontal",
           legend_width = unit(5, "cm")
      ), 
    rect_gp = gpar(col = borderColor), 
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols, 
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontsize), 
    column_names_max_height = unit(100, "mm"),
    
    #Row Options
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontsize), 
    cluster_rows = clusterRows, 
    show_row_dend = showRowDendrogram, 
    clustering_method_rows = "ward.D2",
    split = split, 
    
    #Annotation
    top_annotation = ht1Anno, 

    #Raster Info
    use_raster = useRaster, 
    raster_device = rasterDevice, 
    raster_quality = rasterQuality
  )

  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    ht1 <- ht1 + rowAnnotation(link = 
        anno_check_version_rows(at = customRowLabel, labels = customRowLabelIDs),
        width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs))
  }

  if(draw){
    draw(ht1, 
      padding = unit(c(padding, padding, padding, padding), "mm"), 
      heatmap_legend_side = "bot", 
      annotation_legend_side = "bot")
  }else{
    ht1
  }

}

.colorMapForCH <- function(colorMap, colData){
  colorMap <- colorMap[which(names(colorMap) %in% colnames(colData))]
  colorMapCH <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete")){
      colorx <- colorMap[[x]]
    }else{
      vals <- colData[[names(colorMap)[x]]][!is.na(colData[[names(colorMap)[x]]])]
      s <-  seq(min(vals), max(vals), length.out = length(colorMap[[x]]))
      colorx <- circlize::colorRamp2(s, colorMap[[x]])
    }
    if(any(is.na(names(colorx)))){
      names(colorx)[is.na(names(colorx))] <- paste0("NA",seq_along(names(colorx)[is.na(names(colorx))]))
    }
    return(colorx)
  })
  names(colorMapCH) <- names(colorMap)
  return(colorMapCH)
}

.checkShowLegend <- function(colorMap, max_discrete = 30){
  show <- lapply(seq_along(colorMap), function(x){
      if(attr(colorMap[[x]],"discrete") && length(unique(colorMap[[x]])) > max_discrete){
        sl <- FALSE
      }else{
        sl <- TRUE
      }
      return(sl)
    }) %>% unlist
  names(show) <- names(colorMap)
  return(show)
}

.colorMapAnno <- function(colData, customAnno = NULL, discreteSet = "stallion", continuousSet = "solar_extra"){
  discreteCols <- sapply(colData,function(x) !is.numeric(x))
  if(!is.null(customAnno)){
    colorMap <- lapply(seq_along(discreteCols),function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      if(length(which(customAnno[,1] %in% names(discreteCols[x]))) > 0){
        if(length(which(customAnno[,2] %in% names(colors))) > 0){
          customAnnox <- customAnno[which(customAnno[,2] %in% names(colors)),]
          colors[which(names(colors) %in% customAnnox[,2])] <- paste0(customAnnox[match(names(colors),customAnnox[,2]),3])
        }
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }else{
    colorMap <- lapply(seq_along(discreteCols), function(x){
      if(discreteCols[x]){
       colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
       names(colors) <- unique(colData[[names(discreteCols[x])]])
       attr(colors, "discrete") <- TRUE
      }else{
       colors <- paletteContinuous(set = continuousSet)
       attr(colors, "discrete") <- FALSE
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }

}

.binarySort <- function(m, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){

  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }

  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  return(list(mat = m, hclust = hc))

}

#' Peak Annotation Hypergeometric Enrichment in Marker Peaks.
#' 
#' This function will perform hypergeometric enrichment of peakAnnotation within the defined Marker Peaks (see markerFeatures).
#' 
#' @param seMarker  A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param ArchRProj An `ArchRProject` object.
#' @param peakAnnotation A peakAnnotation in an `ArchRProject` to be used for hypergeometric test.
#' @param matches A custom peakAnnotations matches object used as input (see motifmatchr::matchmotifs).
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker`. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param background Whether to use background matched peaks "bgdPeaks" to compare against or all peaks "all" (see addBgdPeaks).
#' @param ... additional args
#' @export
markerAnnoEnrich <- function(
  seMarker = NULL,
  ArchRProj = NULL,
  peakAnnotation = NULL,
  matches = NULL,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  background = "bgdPeaks",
  ...
  ){

  tstart <- Sys.time()
  if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
    stop("Only markers identified from PeakMatrix can be used!")
  }

  if(is.null(matches)){
    matches <- getMatches(ArchRProj, peakAnnotation)
  }
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  mcols(r1) <- NULL

  r2 <- getPeakSet(ArchRProj)
  mcols(r2) <- NULL

  if(length(which(paste0(seqnames(r1),start(r1),end(r1), sep = "_") %ni% paste0(seqnames(r2),start(r2),end(r2),sep="_"))) != 0){
    stop("Peaks from matches do not match peakSet in ArchRProj!")
  }

  r3 <- GRanges(rowData(seMarker)$seqnames,IRanges(rowData(seMarker)$start, rowData(seMarker)$end))
  mcols(r3) <- NULL
  rownames(matches) <- paste0(seqnames(matches),start(matches),end(matches),sep="_")
  matches <- matches[paste0(seqnames(r3),start(r3),end(r3), sep = "_"), ]

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(tolower(background) %in% c("backgroundpeaks", "bgdpeaks", "background", "bgd")){
    method <- "bgd"
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  }else{
    method <- "all"
  }

  enrichList <- lapply(seq_len(ncol(seMarker)), function(x){
    .messageDiffTime(sprintf("Computing Enrichments %s of %s",x,ncol(seMarker)),tstart)
    idx <- which(passMat[, x])
    if(method == "bgd"){
      .computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx,])))
    }else{
      .computeEnrichment(matches, idx, seq_len(nrow(matches)))
    }
  }) %>% SimpleList
  names(enrichList) <- colnames(seMarker)

  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
      enrichList[[y]][colnames(matches),x,drop=FALSE]
    }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- colnames(enrichList[[1]])
  assays <- rev(assays)
  out <- SummarizedExperiment::SummarizedExperiment(assays=assays)

  out

}

.computeEnrichment <- function(matches, compare, background){

  matches <- .getAssay(matches,  grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
     pOut$BackgroundFrequency[x], # Number of all successes in background
     pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
     pOut$nCompare[x], # Number that were drawn
     lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)

  #Minus Log10 FDR
  pOut$mlog10FDR <- -log10(p.adjust(matrixStats::rowMaxs(cbind(10^-pOut$mlog10p, 4.940656e-324)), method = "fdr"))
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

  pOut

}

#' Identify Marker Feature Ranges
#' 
#' This function will identify Markers and return a GRangesList for each group of significant marker regions.
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker`. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param ... additional args
#' @export
markerGR <- function(
  seMarker,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  ...
  ){
  
  if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
    stop("Only markers identified from PeakMatrix can be used!")
  }

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  FDR <- assays(seMarker)[["FDR"]]

  rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, rowData(seMarker)$end))

  grL <- lapply(seq_len(ncol(passMat)), function(x){
    idx <- which(passMat[, x])
    rrx <- rr[idx]
    rrx$FDR <- FDR[idx, x]
    rrx
  }) %>% GenomicRangesList

  names(grL) <- colnames(seMarker)

  grL <- grL[gtools::mixedsort(names(grL))]

  grL

}

#' Plot Differential Markers
#' 
#' This function will plot one group/column of a differential markers se as a MA or Volcano plot.
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `ArchR::markerFeatures()`.
#' @param name column/group name of seMarker to be plotted.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param plotAs plot as "Volcano" or "MA" plot.
#' @param ... additional args
#' @export
markerPlot <- function(
  seMarker,
  name = NULL,
  cutOff = "FDR <= 0.01 & abs(Log2FC) >= 0.5",
  plotAs = "Volcano",
  scaleTo = 10^4,
  ...
  ){

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  passMat[is.na(passMat)] <- FALSE

  if(is.null(name)){
    name <- colnames(seMarker)[1]
  }
  
  LFC <- assays(seMarker[,name])$Log2FC
  LFC <- as.vector(as.matrix(LFC))

  FDR <- assays(seMarker[,name])$FDR
  FDR <- as.vector(as.matrix(FDR))
  FDR[is.na(FDR)] <- 1

  LM <- log2((assays(seMarker[,name])$Mean + assays(seMarker[,name])$MeanBGD)/2 + 1)
  LM <- as.vector(as.matrix(LM))

  color <- ifelse(passMat[, name], "Differential", "Not-Differential")
  color[color == "Differential"] <- ifelse(LFC[color == "Differential"] > 0, "Up-Regulated", "Down-Regulated")
  pal <- c("Up-Regulated" = "firebrick3", "Not-Differential" = "lightgrey", "Down-Regulated" = "dodgerblue3")

  idx <- c(which(!passMat[, name]), which(passMat[, name]))

  title <- sprintf("Number of features = %s\nNumber Up-Regulated = %s (%s Percent)\nNumber Down-Regulated = %s (%s Percent)",
    nrow(seMarker), sum(color=="Up-Regulated"), 
    round(100 * sum(color=="Up-Regulated") / nrow(seMarker), 2),
    sum(color=="Down-Regulated"), 
    round(100 * sum(color=="Down-Regulated") / nrow(seMarker), 2)
    )

  if(tolower(plotAs) == "ma"){
    ggPoint(
      x = LM[idx],
      y = LFC[idx], 
      color = color[idx], 
      rastr = TRUE, 
      pal = pal,
      xlabel = "Log2 Mean",
      ylabel = "Log2 Fold Change",
      title = title
    ) + geom_hline(yintercept = 0, lty = "dashed")
  }else if(tolower(plotAs) == "volcano"){
    ggPoint(
      x = LFC[idx],
      y = -log10(FDR[idx]), 
      color = color[idx], 
      rastr = TRUE, 
      pal = pal,
      xlabel = "Log2 Fold Change",
      ylabel = "-Log10 FDR",
      title = title
    ) + geom_vline(xintercept = 0, lty = "dashed")
  }else{
    stop("plotAs not recognized")
  }

}



#' Identify Marker Features for each Group
#' 
#' This function will identify a null set of cells that match biases per cell
#' while maintaining the input group proportions. Then it will compute a pairwise
#' test of the group vs the null set.
#' 
#' @param ArchRProj ArchR Project
#' @param groupBy group cells by this column in cellColData
#' @param useGroups use subset of groups in group column in cellColData for comparisons
#' @param bdgGroups use subset of groups in group column in cellColData for background
#' @param useMatrix matrix name in Arrow Files that will be used for identifying features
#' @param bias biases to account for in selecting null group using info from cellColData
#' @param normBy normalize by column in cellColData prior to test
#' @param testMethod pairwise test method group vs null
#' @param minCells minimum cells per group for testing
#' @param maxCells maximum cells per group for testing
#' @param k knn for matching cell biases
#' @param bufferRatio buffering ratio for matching cell biases
#' @param binarize binarize prior to testing
#' @param method marker identification method
#' @param useSeqnames specific seqnames to use only
#' @param verboseHeader verbose sections
#' @param verboseAll verbose sections and subsections
#' @param ... additional args
#' @export
markerFeatures <- function(
    ArchRProj = NULL,
    groupBy = "Clusters",
    useGroups = NULL,
    bdgGroups = NULL,
    useMatrix = "GeneScoreMatrix",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    normBy = NULL,
    testMethod = "wilcoxon",
    minCells = 50,
    maxCells = 500,
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
    bdgGroups = NULL,
    normBy = NULL,
    minCells = 50,
    maxCells = 500,
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
    if(!is.null(useSeqnames)){
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
    if(all(c("deviations","z") %in% unique(paste0(featureDF$seqnames)))){
      message("Detected using deviations matrix without using deviations or z!\nDefaulting to using z values!\nTo use deviations set useSeqnames='deviations'")
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% "z"),]
    }

    #####################################################
    # Match Bias Groups
    #####################################################
    .messageDiffTime("Matching Known Biases", tstart, addHeader = verboseAll)
    groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    
    if(!is.null(useGroups)){
      if(any(useGroups %ni% groups)){
        stop("Not all useGroups in Group names!")
      }
      groups <- groups[groups %in% useGroups]
    }

    if(!is.null(bdgGroups)){
      if(any(bdgGroups %ni% groups)){
        stop("Not all bdgGroups in Group names!")
      }
      groups <- groups[groups %in% bdgGroups]
    }

    matchObj <- .matchBiasCellGroups(
      input = getCellColData(ArchRProj), 
      groups = groups,
      bias = bias,
      k = k,
      n = maxCells
    )

    #####################################################
    # Pairwise Test Per Seqnames
    #####################################################
    .messageDiffTime("Computing Pairwise Tests", tstart, addHeader = verboseAll)
    if(is.null(normBy)){
      if(tolower(useMatrix) %in% c("tilematrix","peakmatrix")){
        normBy <- "ReadsInTSS"
        normFactors <- getCellColData(ArchRProj, normBy, drop=FALSE)
        normFactors[,1] <- median(normFactors[,1]) / normFactors[,1]
      }else{
        normFactors <- NULL
      }
    }else{
      normFactors <- NULL
    }

    diffList <- .safelapply(seq_along(matchObj[[1]]), function(x){
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

    #####################################################
    # Summarize Output
    #####################################################
    if(tolower(testMethod) == "wilcoxon"){
      pse <- SummarizedExperiment::SummarizedExperiment(
          assays = 
            SimpleList(
              Log2FC = lapply(seq_along(diffList), function(x) diffList[[x]]$log2FC) %>% Reduce("cbind",.),
              Mean = lapply(seq_along(diffList), function(x) diffList[[x]]$mean1) %>% Reduce("cbind",.),
              FDR = lapply(seq_along(diffList), function(x) diffList[[x]]$fdr) %>% Reduce("cbind",.),
              AUC = lapply(seq_along(diffList), function(x) diffList[[x]]$auc) %>% Reduce("cbind",.),
              MeanBDG = lapply(seq_along(diffList), function(x) diffList[[x]]$mean2) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else if(tolower(testMethod) == "ttest"){
      pse <- SummarizedExperiment::SummarizedExperiment(
          assays = 
            SimpleList(
              Log2FC = lapply(seq_along(diffList), function(x) diffList[[x]]$log2FC) %>% Reduce("cbind",.),
              Mean = lapply(seq_along(diffList), function(x) diffList[[x]]$mean1) %>% Reduce("cbind",.),
              Variance = lapply(seq_along(diffList), function(x) diffList[[x]]$var1) %>% Reduce("cbind",.),
              FDR = lapply(seq_along(diffList), function(x) diffList[[x]]$fdr) %>% Reduce("cbind",.),
              AUC = lapply(seq_along(diffList), function(x) diffList[[x]]$auc) %>% Reduce("cbind",.),
              MeanBDG = lapply(seq_along(diffList), function(x) diffList[[x]]$mean2) %>% Reduce("cbind",.),
              VarianceBDG = lapply(seq_along(diffList), function(x) diffList[[x]]$var2) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else{
      stop("Error Unrecognized Method!")
    }
    colnames(pse) <- names(matchObj[[1]])

    .messageDiffTime("Completed Pairwise Tests", tstart, addHeader = TRUE)

    return(pse)

}

.matchBiasCellGroups <- function(input, groups, bias, k = 100, n = 500, bufferRatio = 0.8){

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
  }) %>% Reduce("cbind", .)

  #Quantile Normalization
  inputNormQ <- lapply(seq_len(ncol(inputNorm)), function(x){
    .getQuantiles(inputNorm[,x])
  }) %>% Reduce("cbind", .)
  
  #Add Colnames
  colnames(inputNorm) <- bias
  colnames(inputNormQ) <- bias

  #Get proportion of each group
  prob <- table(groups) / length(groups)

  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  matchList <- lapply(seq_along(prob), function(x){
    
    setTxtProgressBar(pb,round(x*100/length(prob),0))

    #############
    # Organize
    #############
    probx <- prob[-x]/sum(prob[-x])
    id <- which(groups==names(prob)[x])
    knnx <- computeKNN(inputNormQ[-id,],inputNormQ[id,], k = k)
    sx <- sample(seq_len(nrow(knnx)), nrow(knnx))
    minTotal <- min(n, length(sx) * bufferRatio)
    nx <- sort(floor(minTotal * probx))
    
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
      groupit <- match(groups[-id][knnit],names(nx))
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
    idX <- seq_len(nrow(inputNormQ))[id][idX]
    idY <- seq_len(nrow(inputNormQ))[-id][idY]

    #####################
    # Matching Stats Groups
    #####################
    estBdg <- sort(floor(minTotal * probx))
    obsBdg <- rep(0, length(estBdg))
    names(obsBdg) <- names(estBdg)
    tabGroups <- table(groups[idY])
    obsBdg[names(tabGroups)] <- tabGroups
    estBdgP <- round(100 * estBdg / sum(estBdg),3)
    obsBdgP <- round(100 * obsBdg / sum(obsBdg),3)

    #####################
    # Matching Stats Bias Norm Values
    #####################
    forBias <- .summarizeColStats(inputNorm[idX,], name = "foreground")
    bdgBias <- .summarizeColStats(inputNorm[idY,], name = "background")

    out <- list(
        cells = idX, 
        bdg = idY, 
        summaryCells = forBias, 
        summaryBdg = bdgBias, 
        bdgGroups = rbind(estBdg, obsBdg),
        bdgGroupsProbs = rbind(estBdgP, obsBdgP),
        corBdgGroups = cor(estBdgP, obsBdgP),
        n = length(sx), 
        p = it / length(sx),
        group = names(prob)[x]
      )

    return(out)

  }) %>% SimpleList
  names(matchList) <- names(prob)
  
  message("\n")

  outList <- SimpleList(
    matchBdg = matchList,
    info = SimpleList(
        cells = rownames(input),
        groups = groups,
        biasNorm = inputNorm,
        biasNormQ = inputNormQ
      )
  )
  
  return(outList)

}

.testMarkerSC <- function(ArrowFiles, matchObj, group = NULL, testMethod = "ttest", useMatrix,
  threads = 1, featureDF, binarize = FALSE, normFactors = NULL){

  matchx <- matchObj[[1]][[group]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bdgx <- matchObj[[2]]$cells[matchx$bdg]
  
  if(!is.null(normFactors)){
    cellNF <- normFactors[cellsx,1]
    bdgNF <- normFactors[bdgx,1]
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
      cellNames = c(cellsx, bdgx),
      progress = FALSE
    ))
    rownames(scMaty) <- rownames(featureDFy)

    if(binarize){
      scMaty@x[scMaty@x > 0] <- 1
    }

    args <- list()
    args$mat1 <- scMaty[, cellsx, drop=FALSE]
    args$mat2 <- scMaty[, bdgx, drop=FALSE]

    if(!is.null(normFactors)){
      cellNF <- normFactors[cellsx,1]
      bdgNF <- normFactors[bdgx,1]
    }

    if(tolower(testMethod) == "wilcoxon"){

      .suppressAll(do.call(.sparseMatWilcoxon, args))
    
    }else if(tolower(testMethod) == "ttest"){
    
      .suppressAll(do.call(.sparseMatTTest, args))
    
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
  offset <- quantile(c(mat1@x,mat2@x), 0.99) * 10^-3
  .requirePackage("presto", installInfo = 'devtools::install_github("immunogenomics/presto")')
  df <- wilcoxauc(cbind(mat1,mat2), c(rep("Top", ncol(mat1)),rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]
  out <- data.frame(
    log2Mean = log2(df$avgExpr + offset),
    log2FC = df$logFC,
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
    offset <- quantile(c(mat1@x,mat2@x), 0.99) * 10^-3
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
    out <- data.frame(
      log2Mean = log2(((m1+offset) + (m2+offset)) / 2),
      log2FC = log2((m1+offset)/(m2+offset)),
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
  offset <- quantile(c(mat1@x,mat2@x), 0.99) * 10^-3
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  #Sparse Row Stats
  s1 <- Matrix::rowSums(mat1, na.rm=TRUE)
  m1 <- s1 / n1
  m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
  #Combute Binom.test
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  pval <- sapply(seq_along(s1), function(x){
    setTxtProgressBar(pb,round(x*100/length(s1),0))
    binom.test(s1[x], n1, m2[x], alternative="two.sided")$p.value
  })
  fdr <- p.adjust(pval, method = "fdr", length(pval))
  out <- data.frame(
    log2Mean = log2(((m1+offset) + (m2+offset)) / 2),
    log2FC = log2((m1+offset) / (m2+offset)),
    fdr = fdr, 
    pval = pval, 
    mean1 = m1, 
    mean2 = m2, 
    n = n1
  )
  return(out)
}






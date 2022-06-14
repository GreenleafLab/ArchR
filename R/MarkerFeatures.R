##########################################################################################
# Marker Feature Methods
##########################################################################################

#' @export
markerFeatures <- function(...){
    .Deprecated("getMarkerFeatures")
    getMarkerFeatures(...)
}

#' Identify Marker Features for each cell grouping
#' 
#' This function will identify features that are definitional of each provided cell grouping where possible
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for marker feature identification.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column
#' in `cellColData`. This limits the groups used to perform marker feature identification.
#' @param bgdGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column
#' in `cellColData` to be used for background calculations in marker feature identification.
#' @param useMatrix The name of the matrix to be used for performing differential analyses. Options include "GeneScoreMatrix", "PeakMatrix", etc.
#' @param bias A character vector indicating the potential bias variables (i.e. c("TSSEnrichment", "log10(nFrags)")) to account
#' for in selecting a matched null group for marker feature identification. These should be column names from `cellColData`.
#' @param normBy The name of a numeric column in `cellColData` that should be normalized across cells (i.e. "ReadsInTSS") prior
#' to performing marker feature identification.
#' @param testMethod The name of the pairwise test method to use in comparing cell groupings to the null cell grouping during
#' marker feature identification. Valid options include "wilcoxon", "ttest", and "binomial".
#' @param maxCells The maximum number of cells to consider from a single-cell group when performing marker feature identification.
#' @param scaleTo Each column in the matrix designated by `useMatrix` will be normalized to a column sum designated by `scaleTo`.
#' @param threads The number of threads to be used for parallel computing.
#' @param k The number of nearby cells to use for selecting a biased-matched background while accounting for `bgdGroups` proportions.
#' @param bufferRatio When generating optimal biased-matched background groups of cells to determine significance, it can be difficult
#' to find sufficient numbers of well-matched cells to create a background group made up of an equal number of cells. The `bufferRatio`
#' indicates the fraction of the total cells that must be obtained when creating the biased-matched group. For example to create a
#' biased-matched background for a group of 100 cells, when `bufferRatio` is set to 0.8 the biased-matched background group will be
#' composed of the 80 best-matched cells. This option provides flexibility in the generation of biased-matched background groups given
#' the stringency of also maintaining the group proportions from `bgdGroups`.
#' @param binarize A boolean value indicating whether to binarize the matrix prior to differential testing. This is useful when
#' `useMatrix` is an insertion counts-based matrix.
#' @param useSeqnames A character vector that indicates which `seqnames` should be plotted in the heatmap. Features from
#' `seqnames` that are not listed will be ignored. In the context of a `Sparse.Assays.Matrix`, such as a matrix containing chromVAR
#' deviations, the `seqnames` do not correspond to chromosomes, rather they correspond to the sub-portions of the matrix, for example
#' raw deviations ("deviations") or deviation z-scores ("z") for a chromVAR deviations matrix.
#' @param verbose A boolean value that determines whether standard output is printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMarkerFeatures <- function(
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
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = bgdGroups, name = "bgdGroups", valid = c("character", "null"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = bias, name = "bias", valid = c("character"))
  .validInput(input = normBy, name = "normBy", valid = c("character", "null"))
  .validInput(input = testMethod, name = "testMethod", valid = c("character", "null"))
  .validInput(input = maxCells, name = "maxCells", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = bufferRatio, name = "bufferRatio", valid = c("numeric"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character", "null"))

  args <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
  .startLogging(logFile = logFile)
  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Input-Parameters", logFile=logFile)
  out <- do.call(.MarkersSC, args)
  .endLogging(logFile = logFile)
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
  bias = NULL,
  k = 100,
  threads = 1,
  binarize = FALSE,
  useSeqnames = NULL,
  testMethod = "wilcoxon",
  useMatrix = "GeneScoreMatrix",
  markerParams = list(),
  verbose = TRUE,
  logFile = NULL
  ){

    tstart <- Sys.time()
    
    #####################################################
    # Feature Info
    #####################################################
    ArrowFiles <- getArrowFiles(ArchRProj)
    featureDF <- .getFeatureDF(head(ArrowFiles, 2), useMatrix)
    matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))

    .logThis(range(as.vector(table(paste0(featureDF$seqnames)))), "FeaturesPerSeqnames", logFile = logFile)

    isDeviations <- FALSE
    if(all(unique(paste0(featureDF$seqnames)) %in% c("z", "deviations"))){
      isDeviations <- TRUE
    }

    .logThis(featureDF, "FeatureDF", logFile=logFile)
    .logMessage(paste0("MatrixClass = ", matrixClass), logFile=logFile)

    seqnames <- unique(as.vector(featureDF$seqnames))
    useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
    if(length(useSeqnames)==0){
      useSeqnames <- NULL
    }

    if(!is.null(useSeqnames)){
      if(matrixClass == "Sparse.Assays.Matrix"){
        if(length(useSeqnames) == 1){
          featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
        }else{
          .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
            "Continuing with first seqname '", seqnames[1], "'!\n",
            "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
          useSeqnames <- seqnames[1]
          featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
        }
      }else{
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }
    }else{
      if(matrixClass == "Sparse.Assays.Matrix"){
        .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
          "Continuing with first seqname '", seqnames[1], "'!\n",
          "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
        useSeqnames <- seqnames[1]
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }
    }
    if(!(nrow(featureDF) > 1)){
      .logStop("Less than 1 feature is remaining in featureDF please check input!", logFile = logFile)
    }

    #####################################################
    # Match Bias Groups
    #####################################################
    .logDiffTime("Matching Known Biases", t1 = tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    colDat <- getCellColData(ArchRProj)

    matchObj <- .matchBiasCellGroups(
      input = colDat, 
      groups = groups,
      useGroups = useGroups,
      bgdGroups = bgdGroups,
      bias = bias,
      k = k,
      n = maxCells,
      bufferRatio = bufferRatio,
      logFile = logFile
    )

    #####################################################
    # Pairwise Test Per Seqnames
    #####################################################
    #ColSums
    mColSums <- tryCatch({
      suppressMessages(tmpColSum <- .getColSums(ArrowFiles, seqnames = featureDF$seqnames@values, useMatrix = useMatrix, threads = threads))
      tmpColSum[ArchRProj$cellNames]
    }, error = function(x){
      rep(1, nCells(ArchRProj))
    })

    if(all(mColSums==1) & is.null(normBy)){
      normBy <- "none"
    }
    
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
      }else if(normBy %in% colnames(ArchRProj@cellColData)) {
        normFactors <- getCellColData(ArchRProj, normBy, drop=FALSE)
        normFactors[,1] <- median(normFactors[,1]) / normFactors[,1]
      }else{
        .logMessage("Warning! Parameter 'normBy' was set to ", normBy," but no matching column was found in cellColData.\n",
                    "Continuing with normalization based on column sums of matrix!", verbose = verbose, logFile = logFile)
        normFactors <- scaleTo / mColSums
        normFactors <- DataFrame(normFactors)
      }
    }

    if(!is.null(normFactors)){
      normFactors[,1] <- normFactors[,1] * (scaleTo / median(normFactors[names(mColSums), 1] * mColSums))
    }
    
    diffList <- .safelapply(seq_along(matchObj[[1]]), function(x){
      .logDiffTime(sprintf("Computing Pairwise Tests (%s of %s)", x, length(matchObj[[1]])), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      .testMarkerSC(
          ArrowFiles = ArrowFiles,
          matchObj = matchObj, 
          group = names(matchObj[[1]])[x], 
          testMethod = testMethod, 
          threads = 1, 
          useMatrix = useMatrix,
          featureDF = featureDF,
          normFactors = normFactors,
          binarize = binarize,
          logFile = logFile
      )
    }, threads = threads)

    .logDiffTime("Completed Pairwise Tests", tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)

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
              Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
              MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
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
              Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
              MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
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
              Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
              MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
              MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
            ),
          rowData = featureDF
        )
    }else{
      stop("Error Unrecognized Method!")
    }
    colnames(pse) <- names(matchObj[[1]])

    metadata(pse)$MatchInfo <- matchObj

    if(isDeviations){
      assays(pse)[["Log2FC"]] <- NULL #This measure does not make sense with deviations matrices better to just remove
    }


    return(pse)

}

.testMarkerSC <- function(
  ArrowFiles = NULL,
  matchObj = NULL,
  group = NULL,
  testMethod = "ttest",
  useMatrix = NULL,
  threads = 1,
  featureDF,
  binarize = FALSE,
  normFactors = NULL,
  logFile = NULL
  ){

  matchx <- matchObj[[1]][[group]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]
  
  if(!is.null(normFactors)){
    cellNF <- as.numeric(normFactors[cellsx,1])
    bgdNF <- as.numeric(normFactors[bgdx,1])
  }

  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  seqnames <- unique(featureDF$seqnames)

  .logThis(cellsx, paste0(group, "_cellsx"), logFile = logFile)
  .logThis(bgdx, paste0(group, "_bgdx"), logFile = logFile)

  pairwiseDF <- lapply(seq_along(seqnames), function(y){

    .logMessage(sprintf("Pairwise Test %s : Seqnames %s", group, seqnames[y]), logFile = logFile)
    featureDFy <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames[y]), ]

    if(length(c(cellsx, bgdx)) == 0){
      stop(paste0("Cells in foreground and background are 0 for group = ", group))
    }

    scMaty <- suppressMessages(.getPartialMatrix(
      ArrowFiles, 
      featureDF = featureDFy, 
      threads = threads, 
      useMatrix = useMatrix,
      cellNames = c(cellsx, bgdx),
      progress = FALSE
    ))
    scMaty <- .checkSparseMatrix(scMaty, length(c(cellsx, bgdx)))
    .logThis(scMaty, paste0(group, "_", seqnames[y], "_scMaty"), logFile = logFile)
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

      o <- tryCatch({
        .suppressAll(do.call(.sparseMatWilcoxon, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatWilcoxon", info = seqnames[y], errorList = errorList, logFile = logFile)
      })
    
    }else if(tolower(testMethod) == "ttest"){
    
      o <- tryCatch({
        .suppressAll(do.call(.sparseMatTTest, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatTTest", info = seqnames[y], errorList = errorList, logFile = logFile)
      })
    
    }else if(tolower(testMethod) == "binomial"){

      if(!is.null(normFactors)){
        .logStop("Normfactors cannot be used with a binomial test!", logFile = logFile)
      }

      if(!binarize){
        .logStop("Binomial test requires binarization!", logFile = logFile)
      }
    
      o <- tryCatch({
        .suppressAll(do.call(.sparseMatBinomTest, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatBinomTest", info = seqnames[y], errorList = errorList, logFile = logFile)
      })

    }else{
    
      .logStop("Error Unrecognized Method!", logFile = logFile)
    
    }

    .logThis(o, paste0(group, "_", seqnames[y], "_diffResult"), logFile = logFile)

    o

  }) %>% Reduce("rbind", .)

  #Check for Mean being 0 for both Mean1 and Mean2
  idxFilter1 <- rowSums(pairwiseDF[,c("mean1","mean2")]) != 0

  #Check For NA in Either Mean1 Mean2
  idxFilter2 <- rowSums(is.na(pairwiseDF[,c("mean1","mean2")])) == 0
  
  #Combo Check
  idxFilter <- idxFilter1 & idxFilter2

  #FDR
  pairwiseDF$fdr <- NA
  pairwiseDF$fdr[idxFilter] <- p.adjust(pairwiseDF$pval[idxFilter], method = "fdr")
  pairwiseDF <- pairwiseDF[rownames(featureDF), , drop = FALSE]
  pairwiseDF
  
}

#Wilcoxon Row-wise two matrices
.sparseMatWilcoxon <- function(mat1 = NULL, mat2 = NULL){
  
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  stopifnot(n1==n2)

  .requirePackage("presto", installInfo = 'devtools::install_github("immunogenomics/presto")')
  df <- wilcoxauc(cbind(mat1,mat2), c(rep("Top", ncol(mat1)),rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]

  #Sparse Row Sums
  m1 <- Matrix::rowSums(mat1, na.rm=TRUE)
  m2 <- Matrix::rowSums(mat2, na.rm=TRUE)
  offset <- 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
  log2FC <- log2((m1 + offset) / (m2 + offset))
  log2Mean <- log2(((m1 + offset) + (m2 + offset)) / 2)

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
.sparseMatTTest <- function(mat1 = NULL, mat2 = NULL, m0 = 0){
    
    #Get Population Values
    n1 <- ncol(mat1)
    n2 <- ncol(mat2)
    stopifnot(n1==n2)
    n <- n1 + n2
    
    #Sparse Row Means
    m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
    m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
    
    #Sparse Row Variances
    v1 <- computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
    v2 <- computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
    
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
      mean1 = m1 / n1, 
      mean2 = m2 / n2, 
      var1 = v2,
      var2 = v2,
      n = n1
    )
    return(out)
}

#Binomial Test Row-wise two matrices
.sparseMatBinomTest <- function(mat1 = NULL, mat2 = NULL){
  
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  stopifnot(n1==n2)
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
  log2Mean <- log2(((m1 + offset) + (m2 + offset)) / 2)

  out <- data.frame(
    log2Mean = log2Mean,
    log2FC = log2FC,
    fdr = fdr, 
    pval = pval, 
    mean1 = m1 / n1, 
    mean2 = m2 / n2, 
    n = n1
  )

  return(out)

}

.matchBiasCellGroups <- function(
  input = NULL,
  groups = NULL,
  useGroups = NULL,
  bgdGroups = NULL,
  bias = NULL,
  k = 100,
  n = 500,
  seed = 1,
  bufferRatio = 0.8,
  logFile = NULL
  ){

  #Summary Function
  .summarizeColStats <- function(m = NULL, name = NULL){
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
  set.seed(seed)
  
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

  #pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  matchList <- lapply(seq_along(useGroups), function(x){
    
    #setTxtProgressBar(pb,round(x*100/length(useGroups),0))

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

    if(all(length(idF) * bgdProbx < 1)){
      if(length(idF) < length(bgdProbx)){
        bgdProbx <- bgdProbx[sample(names(bgdProbx), floor(length(idF) * bufferRatio))]
        bgdProbx[1:length(bgdProbx)] <- rep(1/length(bgdProbx), length(bgdProbx))
      }
    }

    idB <- which(groups %in% names(bgdProbx))

    if(k > length(idB)){
      .logMessage(paste0("Found less than 100 cells for background matching, Lowering k to ", length(idB)), verbose = TRUE, logFile = logFile)
      k2 <- length(idB)
    }else{
      k2 <- k
    }

    knnx <- .computeKNN(inputNormQ[idB, ,drop=FALSE], inputNormQ[idF, ,drop=FALSE], k = k2)
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
          if(selectit == k2){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }else{
          if(cellx %ni% idY){
            selectUnique <- TRUE
          }
          if(selectit == k2){
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
    forBias <- .summarizeColStats(inputNorm[idX,,drop=FALSE], name = "foreground")
    bgdBias <- .summarizeColStats(inputNorm[idY,,drop=FALSE], name = "background")

    out <- list(
        cells = idX, 
        bgd = idY, 
        summaryCells = forBias, 
        summaryBgd = bgdBias, 
        bgdGroups = rbind(estbgd, obsbgd),
        bgdGroupsProbs = rbind(estbgdP, obsbgdP),
        corbgdGroups = suppressWarnings(cor(estbgdP, obsbgdP)),
        n = length(sx), 
        p = it / length(sx),
        group = groupx,
        k = k2
      )

    .logThis(out, paste0("MatchSummary ", useGroups[x]), logFile = logFile)
    return(out)

  }) %>% SimpleList
  names(matchList) <- useGroups
  
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

#' @export
markerHeatmap <- function(...){
    .Deprecated("plotMarkerHeatmap")
    plotMarkerHeatmap(...)
}

#' Plot a Heatmap of Identified Marker Features
#' 
#' This function will plot a heatmap of the results from markerFeatures
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `getMarkerFeatures()`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted
#' in the heatmap. `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values in
#' `seMarker` prior to plotting. Should be set to `TRUE` for counts-based assays (but not assays like "DeviationsMatrix").
#' @param scaleTo Each column in the assay Mean from `seMarker` will be normalized to a column sum designated by `scaleTo`
#' prior to log2 normalization. If log2Norm is `FALSE` this option has no effect.
#' @param scaleRows A boolean value that indicates whether the heatmap should display row-wise z-scores instead of raw values.
#' @param limits A numeric vector of two numbers that represent the lower and upper limits of the heatmap color scheme.
#' @param grepExclude A character vector or string that indicates the `rownames` or a specific pattern that identifies
#' rownames from `seMarker` to be excluded from the heatmap.
#' @param pal A custom continuous palette from `ArchRPalettes` (see `paletteContinuous()`) used to override the default continuous palette for the heatmap.
#' @param binaryClusterRows A boolean value that indicates whether a binary sorting algorithm should be used for fast clustering of heatmap rows.
#' @param clusterCols A boolean value that indicates whether the columns of the marker heatmap should be clustered.
#' @param subsetMarkers A vector of rownames from seMarker to use for subsetting of seMarker to only plot specific features on the heatmap.
#' Note that these rownames are expected to be integers that come from `rownames(rowData(seMarker))`. If this parameter is used for
#' subsetting, then the values provided to `cutOff` are effectively ignored.
#' @param labelMarkers A character vector listing the `rownames` of `seMarker` that should be labeled on the side of the heatmap.
#' @param nLabel An integer value that indicates how many of the top `n` features for each column in `seMarker` should be labeled on the side of the heatmap.
#' To remove all feature labels, set `nLabel = 0`.
#' @param nPrint If provided `seMarker` is from "GeneScoreMatrix" print the top `n` genes for each group based on how uniquely up-regulated the gene is.
#' @param labelRows A boolean value that indicates whether all rows should be labeled on the side of the heatmap.
#' @param returnMatrix A boolean value that indicates whether the final heatmap matrix should be returned in lieu of plotting the actual heatmap.
#' @param transpose A boolean value that indicates whether the heatmap should be transposed prior to plotting or returning.
#' @param invert A boolean value that indicates whether the heatmap will display the features with the
#' lowest `log2(fold change)`. In this case, the heatmap will display features that are specifically lower in the given cell
#' group compared to all other cell groups. Additionally, the color palette is inverted for visualization. This is useful when
#' looking for down-regulated markers (`log2(fold change) < 0`) instead of up-regulated markers (`log2(fold change) > 0`). 
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotMarkerHeatmap <- function(
  seMarker = NULL,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  subsetMarkers = NULL,
  labelMarkers = NULL,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
  ){

  .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  .validInput(input = cutOff, name = "cutOff", valid = c("character"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
  .validInput(input = plotLog2FC, name = "plotLog2FC", valid = c("boolean"))
  .validInput(input = limits, name = "limits", valid = c("numeric"))
  .validInput(input = grepExclude, name = "grepExclude", valid = c("character", "null"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = binaryClusterRows, name = "binaryClusterRows", valid = c("boolean"))
  .validInput(input = clusterCols, name = "clusterCols", valid = c("boolean"))
  .validInput(input = subsetMarkers, name = "subsetMarkers", valid = c("integer", "null"))
  .validInput(input = labelMarkers, name = "labelMarkers", valid = c("character", "null"))
  .validInput(input = nLabel, name = "nLabel", valid = c("integer"))
  .validInput(input = nPrint, name = "nPrint", valid = c("integer"))
  .validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
  .validInput(input = returnMatrix, name = "returnMatrix", valid = c("boolean"))
  .validInput(input = transpose, name = "transpose", valid = c("boolean"))
  .validInput(input = invert, name = "invert", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "markerHeatmap Input-Parameters", logFile = logFile)

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  .logThis(passMat, "passMat", logFile = logFile)

  #Now Get Values
  if(ncol(seMarker) <= 2){
    if(!plotLog2FC){
      stop("Must use plotLog2FC = TRUE when ncol(seMarker) <= 2!")
    }
  }

  #Get Matrix
  if(plotLog2FC){
    mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Log2FC"]])
  }else{
    mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Mean"]])
    if(log2Norm){
      mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
    }
    if(scaleRows){
      mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }
  }
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  .logThis(mat, "mat", logFile = logFile) 

  if(ncol(mat) == 1){
    idx <- which(rowSums(passMat, na.rm = TRUE) > 0)
  }else{
    idx <- which(rowSums(passMat, na.rm = TRUE) > 0 & matrixStats::rowVars(mat) != 0 & !is.na(matrixStats::rowVars(mat)))
  }

  if(!is.null(subsetMarkers)) {
    if(length(which(subsetMarkers %ni% 1:nrow(mat))) == 0){
      idx <- subsetMarkers
    } else {
      stop("Rownames / indices provided to the subsetMarker parameter are outside of the boundaries of seMarker.")
    }
    
  }

  mat <- mat[idx,,drop=FALSE]
  passMat <- passMat[idx,,drop=FALSE]

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
    mat <- mat[idx2,,drop=FALSE]
  }

  if(nrow(mat)==0){
    stop("No Makers Found!")
  }

  spmat <- passMat / rowSums(passMat)
  #only print out identified marker genes if subsetMarkers is NULL
  if(is.null(subsetMarkers)) {
    if(metadata(seMarker)$Params$useMatrix == "GeneScoreMatrix"){
      message("Printing Top Marker Genes:")
      for(x in seq_len(ncol(spmat))){
        genes <- head(order(spmat[,x], decreasing = TRUE), nPrint)
        message(colnames(spmat)[x], ":")
        message("\t", paste(as.vector(rownames(mat)[genes]), collapse = ", "))
      }
    }
  }


  if(is.null(labelMarkers)){
    labelMarkers <- lapply(seq_len(ncol(spmat)), function(x){
      as.vector(rownames(mat)[head(order(spmat[,x], decreasing = TRUE), nLabel)])
    }) %>% unlist %>% unique
  }

  if(ncol(mat) == 1){
    binaryClusterRows <- FALSE
  }

  if(binaryClusterRows){
    if(invert){
      bS <- .binarySort(-mat, lmat = passMat[rownames(mat), colnames(mat),drop=FALSE], clusterCols = clusterCols)
      mat <- -bS[[1]][,colnames(mat),drop=FALSE]
    }else{
      bS <- .binarySort(mat, lmat = passMat[rownames(mat), colnames(mat),drop=FALSE], clusterCols = clusterCols)
      mat <- bS[[1]][,colnames(mat),drop=FALSE]
    }
    clusterRows <- FALSE
    if (clusterCols) {
      clusterCols <- bS[[2]]
    }
  }else{
    clusterRows <- TRUE
    clusterCols <- TRUE
  }

  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }

  message(sprintf("Identified %s markers!", nrow(mat)))

  if(is.null(pal)){
    if(is.null(metadata(seMarker)$Params$useMatrix)){
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }else if(tolower(metadata(seMarker)$Params$useMatrix)=="genescorematrix"){
      pal <- paletteContinuous(set = "blueYellow", n = 100)
    }else{
      pal <- paletteContinuous(set = "solarExtra", n = 100)
    }
  }

  if(invert){
    pal <- rev(pal)
  }

  print(labelMarkers)

  .logThis(mat, "mat-plot", logFile = logFile) 

  if(transpose){

    if(!is.null(clusterCols)){
      mat <- t(mat[seq_len(nrow(mat)), , drop = FALSE])
    }else{
      mat <- t(mat[seq_len(nrow(mat)), clusterCols$order, drop = FALSE])
    }

    if(!is.null(labelMarkers)){
      mn <- match(tolower(labelMarkers), tolower(colnames(mat)), nomatch = 0)
      mn <- mn[mn > 0]
    }else{
      mn <- NULL
    }

    if(returnMatrix){
      .endLogging(logFile = logFile)
      return(mat)
    }

    ht <- tryCatch({

      .ArchRHeatmap(
        mat = mat,
        scale = FALSE,
        limits = c(min(mat), max(mat)),
        color = pal, 
        clusterCols = clusterRows, 
        clusterRows = FALSE,
        labelRows = TRUE,
        labelCols = labelRows,
        customColLabel = mn,
        showRowDendrogram = TRUE,
        draw = FALSE,
        name = paste0("Column Z-Scores\n", ncol(mat), " features\n", metadata(seMarker)$Params$useMatrix)
      )

    }, error = function(e){

      errorList <- list(
        mat = mat,
        scale = FALSE,
        limits = c(min(mat), max(mat)),
        color = pal, 
        clusterCols = clusterRows, 
        clusterRows = FALSE,
        labelRows = TRUE,
        labelCols = labelRows,
        customColLabel = mn,
        showRowDendrogram = TRUE,
        draw = FALSE,
        name = paste0("Column Z-Scores\n", ncol(mat), " features\n", metadata(seMarker)$Params$useMatrix)
      )

    })

  }else{
    
    if(!is.null(labelMarkers)){
      mn <- match(tolower(labelMarkers), tolower(rownames(mat)), nomatch = 0)
      mn <- mn[mn > 0]
    }else{
      mn <- NULL
    }

    if(returnMatrix){
      .endLogging(logFile = logFile)
      return(mat)
    }

    ht <- tryCatch({

      .ArchRHeatmap(
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
        name = paste0("Row Z-Scores\n", nrow(mat), " features\n", metadata(seMarker)$Params$useMatrix)
      )

    }, error = function(e){

      errorList <- list(
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
        name = paste0("Row Z-Scores\n", nrow(mat), " features\n", metadata(seMarker)$Params$useMatrix)
      )

      .logError(e, fn = ".ArchRHeatmap", info = "", errorList = errorList, logFile = logFile)

    })


  }

  .endLogging(logFile = logFile)
  return(ht)

}

#' Get Marker Features from a marker summarized experiment
#' 
#' This function will identify Markers and return a List of Features or a GRangesList for each group of significant marker features.
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `getMarkerFeatures()`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker`. `cutoff` can contain any
#' of the `assayNames` from `seMarker`.
#' @param n An integer that indicates the maximum number of features to return per group.
#' @param returnGR A boolean indicating whether to return as a `GRanges` object. Only valid when `seMarker` is computed for a PeakMatrix.
#' @export
getMarkers <- function(
  seMarker = NULL,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  n = NULL,
  returnGR = FALSE
  ){

  .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  .validInput(input = cutOff, name = "cutOff", valid = c("character"))
  .validInput(input = n, name = "n", valid = c("integer", "null"))
  .validInput(input = returnGR, name = "returnGR", valid = c("boolean"))

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(returnGR){

    if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
      stop("Only markers can be returned as GRanges when PeakMatrix!")
    }

    rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, rowData(seMarker)$end))

    grL <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- rr[idx]
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$FDR),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList

    names(grL) <- colnames(seMarker)

    grL <- grL[gtools::mixedsort(names(grL))]

    return(grL)

  }else{

    markerList <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- SummarizedExperiment::rowData(seMarker[idx,])
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$FDR),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList

    names(markerList) <- colnames(seMarker)

    return(markerList)

  }

}

#' @export
markerPlot <- function(...){
    .Deprecated("plotMarkers")
    plotMarkers(...)
}

#' Plot Differential Markers
#' 
#' This function will plot one group/column of a differential markers as an MA or Volcano plot.
#' 
#' @param seMarker A `SummarizedExperiment` object returned by `getMarkerFeatures()`.
#' @param name The name of a column in `seMarker` (i.e. cell grouping in `groupBy` or `useGroups` for `getMarkerFeatures()`) to be plotted.
#' To see available options try `colnames(seMarker)`.
#' @param cutOff A valid-syntax logical statement that defines which marker features from `seMarker` will be plotted.
#' `cutoff` can contain any of the `assayNames` from `seMarker`.
#' @param plotAs A string indicating whether to plot a volcano plot ("Volcano") or an MA plot ("MA").
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize
#' lines and labels, just the internal portions of the plot.
#' @export
plotMarkers <- function(
  seMarker = NULL,
  name = NULL,
  cutOff = "FDR <= 0.01 & abs(Log2FC) >= 0.5",
  plotAs = "Volcano",
  scaleTo = 10^4,
  rastr = TRUE
  ){

  .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = cutOff, name = "cutOff", valid = c("character"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))

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
  
  FDR <- assays(seMarker[,name])$FDR
  FDR <- as.vector(as.matrix(FDR))
  FDR[is.na(FDR)] <- 1

  if(tolower(plotAs) == "volcanodiff"){
    Diff <- assays(seMarker[,name])$MeanDiff
    Diff <- as.vector(as.matrix(Diff))
    qDiff <- max(quantile(abs(Diff), probs = 0.999, na.rm=TRUE), 4) * 1.05
    color <- ifelse(passMat[, name], "Differential", "Not-Differential")
    color[color == "Differential"] <- ifelse(Diff[color == "Differential"] > 0, "Up-Regulated", "Down-Regulated")
  }else{
    LFC <- assays(seMarker[,name])$Log2FC
    LFC <- as.vector(as.matrix(LFC))
    qLFC <- max(quantile(abs(LFC), probs = 0.999, na.rm=TRUE), 4) * 1.05
    LM <- log2((assays(seMarker[,name])$Mean + assays(seMarker[,name])$MeanBGD)/2 + 1)
    LM <- as.vector(as.matrix(LM))
    color <- ifelse(passMat[, name], "Differential", "Not-Differential")
    color[color == "Differential"] <- ifelse(LFC[color == "Differential"] > 0, "Up-Regulated", "Down-Regulated")
  }

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
      ylim = c(-qLFC, qLFC),
      size = 1,
      extend = 0,
      rastr = rastr, 
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Log2 Mean",
      ylabel = "Log2 Fold Change",
      title = title
    ) + geom_hline(yintercept = 0, lty = "dashed") + 
    scale_y_continuous(breaks = seq(-100, 100, 2), limits = c(-qLFC, qLFC), expand = c(0,0))
  }else if(tolower(plotAs) == "volcano"){
    ggPoint(
      x = LFC[idx],
      y = -log10(FDR[idx]), 
      color = color[idx], 
      xlim = c(-qLFC, qLFC),
      extend = 0,
      size = 1,
      rastr = rastr, 
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Log2 Fold Change",
      ylabel = "-Log10 FDR",
      title = title
    ) + geom_vline(xintercept = 0, lty = "dashed") + 
    scale_x_continuous(breaks = seq(-100, 100, 2), limits = c(-qLFC, qLFC), expand = c(0,0))
  }else if(tolower(plotAs) == "volcanodiff"){
    ggPoint(
      x = Diff[idx],
      y = -log10(FDR[idx]), 
      color = color[idx], 
      xlim = c(-qDiff, qDiff),
      extend = 0,
      size = 1,
      rastr = rastr, 
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Mean Difference",
      ylabel = "-Log10 FDR",
      title = title
    ) + geom_vline(xintercept = 0, lty = "dashed")
  }else{
    stop("plotAs not recognized")
  }

}



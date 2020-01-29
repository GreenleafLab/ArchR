##########################################################################################
# Validation Methods
##########################################################################################
#' @export
.validInput <- function(input = NULL, name = NULL, valid = NULL){

  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }

  if(!is.character(name)){
    stop("name must be a character!")
  }

  av <- lapply(seq_along(valid), function(i){

    vi <- valid[i]

    if(vi == "integer" | vi == "wholenumber"){

      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1))) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }

    }else if(vi == "null"){

      cv <- is.null(input)

    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){

      cv <- is.logical(input)

    }else if(vi == "numeric"){

      cv <- is.numeric(input)

    }else if(vi == "matrix"){

      cv <- is.matrix(input)

    }else if(vi == "sparsematrix"){

      cv <- inherits(input, "dgCMatrix")

    }else if(vi == "character"){

      cv <- is.character(input)

    }else if(vi == "rlecharacter"){

      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }

    }else if(vi == "palette"){

      #https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
      isColor <- function(x){
       unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
      }
      cv <- all(isColor(input))

    }else if(vi == "timestamp"){

      cv <- inherits(input, "POSIXct")

    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){

      cv1 <- is.data.frame(input)
      cv2 <- inherits(input, "DataFrame")
      cv <- any(cv1, cv2)

    }else if(vi == "fileexists"){

      cv <- all(file.exists(input))

    }else if(vi == "direxists"){

      cv <- all(dir.exists(input))

    }else if(vi == "granges" | vi == "gr"){

      cv <- inherits(input, "GRanges")

    }else if(vi == "grangeslist" | vi == "grlist"){

      cv <- all(unlist(lapply(input, function(x) inherits(x, "GRanges"))))

    }else if(vi == "list" | vi == "simplelist"){

      cv1 <- is.list(input)
      cv2 <- inherits(input, "SimpleList")
      cv <- any(cv1, cv2)

    }else if(vi == "se" | vi == "summarizedexperiment"){

      cv <- inherits(input, "SummarizedExperiment")

    }else if(vi == "txdb"){

      cv <- inherits(input, "TxDb")

    }else if(vi == "orgdb"){

      cv <- inherits(input, "OrgDb")

    }else if(vi == "bsgenome"){

      cv <- inherits(input, "BSgenome")

    }else if(vi == "parallelparam"){

      cv <- inherits(input, "BatchtoolsParam")

    }else if(vi == "archrproj" | vi == "archrproject"){

      cv <- inherits(input, "ArchRProject")

    }else{

      stop("Validator is not currently supported by ArchR!")

    }

  }) %>% Reduce("c", .) %>% any

  if(av){

    return(invisible(TRUE))
  
  }else{

    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")

  }

}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome The name of a valid genome (for example "hg38", "hg19", or "mm10"), name of a BSgenome package (BSgenome.Hsapiens.UCSC.hg19) or a BSgenome object.
#' @param masked A boolean describing whether or not to access the masked genome versino. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}

#' @export
.validTxDb <- function(TxDb = NULL){
  stopifnot(!is.null(TxDb))
  if(inherits(TxDb, "TxDb")){
    return(TxDb)
  }else if(is.character(TxDb)){
    return(getTxDb(TxDb)) #change
  }else{
    stop("Cannot validate TxDb options are a valid TxDb or character for getTxDb")
  }
}

#' @export
.validOrgDb <- function(OrgDb = NULL){
  stopifnot(!is.null(OrgDb))
  if(inherits(OrgDb, "OrgDb")){
    return(OrgDb)
  }else if(is.character(OrgDb)){
    return(getOrgDb(OrgDb)) #change
  }else{
    stop("Cannot validate OrgDb options are a valid OrgDb or character for getOrgDb")
  }
}

#' @export
.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

#' @export
.validGeneAnnotation <- function(geneAnnotation, ...){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){

    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])

  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }

  gA

}

#' @export
.validGenomeAnnotation <- function(genomeAnnotation, ...){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){

    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])

  }else{

    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
  
  }

  gA

}

#' @export
.validArchRProject <- function(ArchRProj, ...){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

##########################################################################################
# S4Vectors/BiocGenerics Within Methods
##########################################################################################

#' Negated Value Matching
#'
#' This function is the reciprocal of %in%. See the match funciton in base R.
#'
#' @param x The value to search for in `table`.
#' @param table The set of values to serve as the base for the match function.
#' @export
"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)

#' Generic matching function for S4Vector objects
#'
#' This function provides a generic matching function for S4Vector objects primarily to avoid ambiguity.
#'
#' @param x An `S4Vector` object to search for in `table`.
#' @param table The set of `S4Vector` objects to serve as the base for the match function.
#' @export
'%bcin%' <- function(x, table) S4Vectors::match(x, table, nomatch = 0) > 0

#' Negated matching function for S4Vector objects
#'
#' This function provides the reciprocal of %bcin% for S4Vector objects primarily to avoid ambiguity.
#'
#' @param x An `S4Vector` object to search for in `table`.
#' @param table The set of `S4Vector` objects to serve as the base for the match function.
#' @export
'%bcni%' <- function(x, table) !(S4Vectors::match(x, table, nomatch = 0) > 0)

##########################################################################################
# Helper Intermediate Methods
##########################################################################################

#' @export
.mergeParams <- function(paramInput, paramDefault){
  for(i in seq_along(paramDefault)){
    if(!(names(paramDefault)[i] %in% names(paramInput))){
      paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
    }
  }
  return(paramInput)
}

#' @export
.requirePackage <- function(x, load = TRUE, installInfo = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

#' @export
.messageDiffTime <- function(main = "", t1 = NULL, verbose = TRUE, addHeader = FALSE, t2 = Sys.time(), units = "mins", header = "###########", tail = "elapsed...", precision = 3){
  if(verbose){
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),precision))
      if(addHeader){
        message(sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header))
      }else{
        message(sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail))
      }
    }, error = function(x){
      message("Time Error : ", x)
    })
  }
  return(0)
}

##########################################################################################
# Lapply Methods
##########################################################################################

#' @export
.safelapply <- function(..., threads = 1, preschedule = FALSE){

  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }

  if(threads > 1){

    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

    for(i in seq_along(o)){
      if(inherits(o[[i]], "try-error")){
        stop(o)
      }
    }

  }else{

    o <- lapply(...)

  }

  o

}

#' @export
.batchlapply <- function(args, sequential = FALSE){

  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }

  #Determine Parallel Backend
  if(inherits(args$parallelParam, "BatchtoolsParam")){

    .messageDiffTime("Batch Execution w/ BatchTools through BiocParallel!", args$tstart)

    require(BiocParallel)
    
    args$parallelParam <- btParam
    #Unlink registry Directory
    if(dir.exists(args$registryDir)){
      #Clean Up Registry
      unlink(args$registryDir, recursive = TRUE)# Delete registry directory
    }

    #Set Up Registry For Runnning
    args$parallelParam$registryargs <- batchtoolsRegistryargs(
      file.dir = args$registryDir,
      work.dir = getwd(),
      packages = character(0L),
      namespaces = character(0L),
      source = character(0L),
      load = character(0L)
    )

    #Register
    BPPARAM <- args$parallelParam
    register(BPPARAM)

    #Add To Args
    args$BPPARAM <- BPPARAM

    if("..." %in% names(args)){
      args["..."] <- NULL
    }

    #Run
    outlist <- do.call(bplapply, args)

  }else{

    .messageDiffTime("Batch Execution w/ safelapply!", args$tstart)
    if(sequential){
      args$subThreads <- args$threads
      args$threads <- 1
    }else{
      if(args$threads > length(args$X)){
        args$subThreads <- floor( args$threads / length(args$X) )
        args$threads <- length(args$X)
      }else{
        args$subThreads <- 1
      }
    }
    outlist <- do.call(.safelapply, args)

  }

  return(outlist)

}

##########################################################################################
# Stat/Summary Methods
##########################################################################################

#' @export
.rowZscores <- function(m,min=-2,max=2,limit=FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

#' @export
.computeROC <- function(labels, scores, name="ROC"){
  calcAUC <- function(TPR, FPR){
    # http://blog.revolutionanalytics.com/2016/11/calculating-auc.html
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    out <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2
    return(out)
  }
  labels <- labels[order(scores, decreasing=TRUE)]
  df <- data.frame(
    False_Positive_Rate = cumsum(!labels)/sum(!labels),
    True_Positive_Rate =  cumsum(labels)/sum(labels)
    )
  df$AUC <- round(calcAUC(df$True_Positive_Rate,df$False_Positive_Rate),3)
  df$name <- name
  return(df)
}

#' @export
.getQuantiles <- function(v, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

#' @export
.rowScale <- function(mat, min = NULL, max = NULL){
  if(!is.null(min)){
    rMin <- min
  }else{
    rMin <- matrixStats::rowMins(mat)
  }
  if(!is.null(max)){
    rMax <- max
  }else{
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff/rScale
  out <- list(mat=matScale, min=rMax, max=rMin)
  return(out)
}

#' @export
.quantileCut <- function(x, lo = 0.025, hi = 0.975, maxIf0 = TRUE){
  q <- quantile(x, probs = c(lo,hi))
  if(q[2] == 0){
    if(maxIf0){
      q[2] <- max(x)
    }
  }
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

#' @export
.normalizeCols <- function(mat, colSm = NULL, scaleTo = NULL){
    if(is.null(colSm)){
        colSm <- Matrix::colSums(mat)
    }
    if(!is.null(scaleTo)){
        mat@x <- scaleTo * mat@x / rep.int(colSm, Matrix::diff(mat@p))
    }else{
        mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
    }
    return(mat)
}

#' @export
.confusionMatrix <- function(i,j){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#' @export
.safeSubset <- function(mat, subsetRows = NULL, subsetCols = NULL){
  
  if(!is.null(subsetRows)){
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(length(idxNotIn), ncol = ncol(mat)))
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows,]
  }

  if(!is.null(subsetCols)){
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[,subsetCols]
  }

  mat

}

#' @export
.groupMeans <- function(mat, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

#' @export
.groupSums <- function(mat, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

#' @export
.groupSds <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gs <- lapply(unique(groups), function(x){
    if (sparse){
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }else{
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gs) <- unique(groups)
  return(gs)
}

#' @export
.centerRollMean <- function(v, k){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}

##########################################################################################
# Miscellaneous Methods
##########################################################################################

#' @export
.cleanParams <- function(params, maxSize = 0.05){
  lapply(params, function(x){
    os <- object.size(x)
    if(os > 10^6 * maxSize){
      NULL
    }else{
      x
    }
  })
}

#' @export
.splitEvery <- function(x, n){
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if(is.atomic(x)){
    split(x, ceiling(seq_along(x) / n))
  }else{
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}

#' @export
.suppressAll <- function(expr){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

#' @export
.getAssay <- function(se, assayName = NULL){
  assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(assayNames(se),collapse=", ")))
  }
  return(o)
}

#' @export
.fileExtension <- function (x){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

#' @export
.checkPath <- function(u = NULL, path = NULL, throwError = TRUE){
  if(is.null(u)){
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE){
    if (Sys.which(x) == "") {
      if(!is.null(path) && file.exists(file.path(path, x))){
        o <- TRUE
      }else{
        if(throwError){
          stop(x, " not found in path, please add ", x, " to path!")
        }else{
          o <- FALSE
        }
      }
    }else{
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

#' @export
.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){

  dir.create(tmpdir, showWarnings = FALSE)

  if(addDOC){
    doc <- paste0("Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }

  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))

}

#' @export
.ArchRLogo <- function(ascii = "Logo"){
  Ascii <- list(
    Package = c("
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    "),

    #modified from cyu@athena.mit.edu
    Logo = c("
                                                   / |
                                                 /    \\\
            .                                  /      |.
            \\\\\\                              /        |.
              \\\\\\                          /           `|.
                \\\\\\                      /              |.
                  \\\                    /                |\\\
                  \\\\#####\\\           /                  ||
                ==###########>      /                   ||
                 \\\\##==......\\\    /                     ||
            ______ =       =|__ /__                     ||      \\\\\\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \\\               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\\\####\\\\________________,--\\\\_##,/
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    ")
  )
  message(Ascii[[ascii]])
}


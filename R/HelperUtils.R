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
# Helper to try to reformat fragment files appropriately if a bug is found
##########################################################################################

#' Reformat Fragment Files to be Tabix and Chr Sorted JJJ
#'
#' This function provides help in reformatting Fragment Files for reading in createArrowFiles.
#' It will handle weird anomalies found that cause errors in reading tabix bgzip'd fragment files.
#'
#' @param fragmentFiles a character vector of paths to fragment files to be reformatted
#' @param seqnamesIsChr a boolean describing to check if seqnames containt "chr".
#' @export
reformatFragmentFiles <- function(
  fragmentFiles = NULL,
  seqnamesIsChr = TRUE
  ){
  options(scipen = 999)
  .requirePackage("data.table")
  .requirePackage("Rsamtools")
  for(i in seq_along(fragmentFiles)){
    message(i, " of ", length(fragmentFiles))
    dt <- data.table::fread(fragmentFiles[i])
    dt <- dt[order(dt$V1,dt$V2,dt$V3), ]
    if(seqnamesIsChr){
      idxRemove1 <- which(substr(dt$V1,1,3) != "chr")
    }else{
      idxRemove1 <- c()
    }
    idxRemove2 <- which(dt$V2 != as.integer(dt$V2))
    idxRemove3 <- which(dt$V3 != as.integer(dt$V3))
    #get all
    idxRemove <- unique(c(idxRemove1, idxRemove2, idxRemove3))
    if(length(idxRemove) > 0){
      dt <- dt[-idxRemove,]
    }
    fileNew <- gsub(".tsv.bgz|.tsv.gz", "-Reformat.tsv", fragmentFiles[i])
    data.table::fwrite(dt, fileNew, sep = "\t", col.names = FALSE)
    Rsamtools::bgzip(fileNew)
    ArchR:::.fileRename(paste0(fileNew, ".bgz"), paste0(fileNew, ".gz"))
  }
}

##########################################################################################
# Validation Methods
##########################################################################################
.validInput <- function(input = NULL, name = NULL, valid = NULL){

  valid <- unique(valid)

  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }

  if(!is.character(name)){
    stop("name must be a character!")
  }

  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }

  av <- FALSE

  for(i in seq_along(valid)){

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

      cv <- all(.isColor(input))

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

      #I think there could be a bug here JJJ
      cv <- all(unlist(lapply(input, function(x) inherits(x, "GRanges"))))

    }else if(vi == "list" | vi == "simplelist"){

      cv1 <- is.list(input)
      cv2 <- inherits(input, "SimpleList")
      cv <- any(cv1, cv2)

    }else if(vi == "se" | vi == "summarizedexperiment"){

      cv <- inherits(input, "SummarizedExperiment")

    }else if(vi == "seurat" | vi == "seuratobject"){

      cv <- inherits(input, "Seurat")

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
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! JJJ

    }else{

      stop("Validator is not currently supported by ArchR!")

    }

    if(cv){
      av <- TRUE
      break
    }   
     
  }

  if(av){

    return(invisible(TRUE))
  
  }else{

    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")

  }

}

#https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
.isColor <- function(x = NULL){
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid genome (for example "hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
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

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.validGeneAnnotation <- function(geneAnnotation = NULL){
  
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

.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
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

.validArchRProject <- function(ArchRProj = NULL){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}

##########################################################################################
# Helper Intermediate Methods
##########################################################################################

.mergeParams <- function(paramInput = NULL, paramDefault = NULL){
  for(i in seq_along(paramDefault)){
    if(!(names(paramDefault)[i] %in% names(paramInput))){
      paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
    }
  }
  return(paramInput)
}

.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

#' @export
#JJJ
mapLabels <- function(labels = NULL, newLabels = NULL, oldLabels = names(newLabels)){

  .validInput(input = labels, name = "labels", valid = c("character"))
  .validInput(input = newLabels, name = "newLabels", valid = c("character"))
  .validInput(input = oldLabels, name = "oldLabels", valid = c("character"))

  if(length(newLabels) != length(oldLabels)){
    stop("newLabels and oldLabels must be equal length!")
  }

  if(!requireNamespace("plyr", quietly = TRUE)){
    labels <- paste0(labels)
    oldLabels <- paste0(oldLabels)
    newLabels <- paste0(newLabels)
    labelsNew <- labels
    for(i in seq_along(oldLabels)){
        labelsNew[labels == oldLabels[i]] <- newLabels[i]
    }
    paste0(labelsNew)
  }else{
    paste0(plyr::mapvalues(x = labels, from = oldLabels, to = newLabels))
  }

}

##########################################################################################
# Stat/Summary Methods
##########################################################################################

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

.computeROC <- function(labels = NULL, scores = NULL, name="ROC"){
  .calcAUC <- function(TPR = NULL, FPR = NULL){
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
  df$AUC <- round(.calcAUC(df$True_Positive_Rate,df$False_Positive_Rate),3)
  df$name <- name
  return(df)
}

.getQuantiles <- function(v = NULL, len = length(v)){
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

.rowScale <- function(mat = NULL, min = NULL, max = NULL){
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
  out <- list(mat=matScale, min=rMin, max=rMax)
  return(out)
}

.quantileCut <- function(x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE){
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

.normalizeCols <- function(mat = NULL, colSm = NULL, scaleTo = NULL){
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

.confusionMatrix <- function(i = NULL, j = NULL){
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

.safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL){
  
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

.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
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

.groupSums <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
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

.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
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

.centerRollMean <- function(v = NULL, k = NULL){
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

.splitEvery <- function(x = NULL, n = NULL){
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if(is.atomic(x)){
    split(x, ceiling(seq_along(x) / n))
  }else{
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}

.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

.fileExtension <- function (x = NULL){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

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

.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){

  dir.create(tmpdir, showWarnings = FALSE)

  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }

  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))

}

.ArchRLogo <- function(ascii = "Logo", messageLogo = TRUE){
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

  if(messageLogo){
    message(Ascii[[ascii]])
  }else{
    Ascii[[ascii]]
  }

}

########
# Developer Utils
########

.devMode <- function(package = "ArchR"){
  # fn <- unclass(lsf.str(envir = asNamespace(package), all = TRUE))
  # for(i in seq_along(fn)){
  #   tryCatch({
  #     assign(fn[i], paste0(package,':::', fn[i]), envir=globalenv())
  #     #eval(parse(text=paste0(fn[i], paste0('<<-',package,':::'), fn[i])))
  #   }, error = function(x){
  #   })
  # }
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for(i in seq_along(fn)){
    tryCatch({
      eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
    }, error = function(x){
    })
  }
}

#https://stackoverflow.com/questions/42734547/generating-random-strings
.randomStr <- function(letters = 10, n = 1){
  a <- do.call(paste0, replicate(letters, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
.isWholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

.convertToPNG <- function(
  ArchRProj = NULL,
  paths = c("QualityControl"),
  recursive = TRUE,
  outDir = "Figures"
  ){

  #If error try
  #brew install fontconfig

  .requirePackage("pdftools", source = "cran")

  if(!is.null(ArchRProj)){
    paths <- c(paths, file.path(getOutputDirectory(ArchRProj), "Plots"))
  }
  
  pdfFiles <- lapply(seq_along(paths), function(i){
    if(recursive){
      dirs <- list.dirs(paths[i], recursive = FALSE, full.names = FALSE)
      if(length(dirs) > 0){
        pdfs <- lapply(seq_along(dirs), function(j){
          list.files(file.path(paths[i], dirs[j]), full.names = TRUE, pattern = "\\.pdf")
        }) %>% unlist
      }else{
        pdfs <- c()
      }
      pdfs <- c(list.files(paths[i], full.names = TRUE, pattern = "\\.pdf"), pdfs)
    }else{
      pdfs <- list.files(paths[i], full.names = TRUE, pattern = "\\.pdf")
    }
    pdfs
  }) %>% unlist

  dir.create(outDir, showWarnings = FALSE)

  for(i in seq_along(pdfFiles)){
    print(i)
    tryCatch({
      pdf_convert(
        pdfFiles[i], 
        format = "png", 
        pages = NULL, 
        filenames = file.path(outDir, gsub("\\.pdf", "_%d.png",basename(pdfFiles[i]))),
        dpi = 300, 
        opw = "", 
        upw = "", 
        verbose = TRUE
      )
    },error=function(x){
      0
    })
  }

}


##########################################################################################
# Batch Methods
##########################################################################################

.safelapply <- function(..., threads = 1, preschedule = FALSE){

  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }

  if(threads > 1){

    # o <- .modifiedMclapply(..., mc.cores = threads, mc.preschedule = preschedule, stop.on.error = TRUE)

    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

    for(i in seq_along(o)){ #JJJ this doesnt seem to work?
      if(inherits(o[[i]], "try-error")){
        stop(o)
      }
    }

  }else{

    o <- lapply(...)

  }

  o

}

.batchlapply <- function(args = NULL, sequential = FALSE){

  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }

  #Determine Parallel Backend
  if(inherits(args$parallelParam, "BatchtoolsParam")){

    stop("Batchtools not yet fully supported please use local parallel threading!")

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
    args <- args[names(args) %ni% c("threads", "parallelParam", "subThreading")]
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

    args <- args[names(args) %ni% c("registryDir", "parallelParam", "subThreading")]
    outlist <- do.call(.safelapply, args)

  }

  return(outlist)

}

#https://stackoverflow.com/questions/56353069/r-parallel-abort-all-mclapply-operations
#This is super useful for safe parallel apply erroring
.modifiedMclapply <- function(
    X, 
    FUN, 
    ..., 
    mc.preschedule = TRUE,
    mc.set.seed = TRUE, 
    mc.silent = FALSE,
    mc.cores = getOption("mc.cores", 2L), 
    mc.cleanup = TRUE, 
    mc.allow.recursive = TRUE,
    affinity.list = NULL, 
    stop.on.error = TRUE
    ){
    stop.on.error <- stop.on.error[1]        #!!
    stopifnot(is.logical(stop.on.error))     #!!
    cores <- as.integer(mc.cores)
    if ((is.na(cores) || cores < 1L) && is.null(affinity.list)) 
        stop("'mc.cores' must be >= 1")
    parallel:::.check_ncores(cores)
    if (parallel:::isChild() && !isTRUE(mc.allow.recursive)) 
        return(lapply(X = X, FUN = FUN, ...))
    if (!is.vector(X) || is.object(X)) 
        X <- as.list(X)
    if (!is.null(affinity.list) && length(affinity.list) < length(X)) 
        stop("affinity.list and X must have the same length")
    if (mc.set.seed) 
        mc.reset.stream()
    if (length(X) < 2) {
        old.aff <- mcaffinity()
        mcaffinity(affinity.list[[1]])
        res <- lapply(X = X, FUN = FUN, ...)
        mcaffinity(old.aff)
        return(res)
    }
    if (length(X) < cores) 
        cores <- length(X)
    if (cores < 2L && is.null(affinity.list)) 
        return(lapply(X = X, FUN = FUN, ...))
    jobs <- list()
    parallel:::prepareCleanup()
    on.exit(parallel:::cleanup(mc.cleanup))
    if (!mc.preschedule) {
        FUN <- match.fun(FUN)
        if (length(X) <= cores && is.null(affinity.list)) {
            jobs <- lapply(seq_along(X), function(i) mcparallel(FUN(X[[i]], 
                ...), name = names(X)[i], mc.set.seed = mc.set.seed, 
                silent = mc.silent))
            res <- mccollect(jobs)
            if (length(res) == length(X)) 
                names(res) <- names(X)
            has.errors <- sum(sapply(res, inherits, "try-error"))
        }
        else {
            sx <- seq_along(X)
            res <- vector("list", length(sx))
            names(res) <- names(X)
            fin <- rep(FALSE, length(X))
            if (!is.null(affinity.list)) {
                cores <- max(unlist(x = affinity.list, recursive = TRUE))
                d0 <- logical(cores)
                cpu.map <- lapply(sx, function(i) {
                  data <- d0
                  data[as.vector(affinity.list[[i]])] <- TRUE
                  data
                })
                ava <- do.call(rbind, cpu.map)
            }
            else {
                ava <- matrix(TRUE, nrow = length(X), ncol = cores)
            }
            jobid <- integer(cores)
            for (i in 1:cores) {
                jobid[i] <- match(TRUE, ava[, i])
                ava[jobid[i], ] <- FALSE
            }
            if (anyNA(jobid)) {
                unused <- which(is.na(jobid))
                jobid <- jobid[-unused]
                ava <- ava[, -unused, drop = FALSE]
            }
            jobs <- lapply(jobid, function(i) mcparallel(FUN(X[[i]], 
                ...), mc.set.seed = mc.set.seed, silent = mc.silent, 
                mc.affinity = affinity.list[[i]]))
            jobsp <- parallel:::processID(jobs)
            has.errors <- 0L
            delivered.result <- 0L
            while (!all(fin)) {
                s <- parallel:::selectChildren(jobs[!is.na(jobsp)], -1)
                if (is.null(s)) 
                  break
                if (is.integer(s)) 
                  for (ch in s) {
                    ji <- match(TRUE, jobsp == ch)
                    ci <- jobid[ji]
                    r <- parallel:::readChild(ch)
                    if (is.raw(r)) {
                      child.res <- unserialize(r)
                      if (inherits(child.res, "try-error")){
                          if(stop.on.error)                     #!!
                              stop("error in process X = ", ci, "\n", attr(child.res, "condition")$message) #!!
                          has.errors <- has.errors + 1L
                      }
                      if (!is.null(child.res)) 
                        res[[ci]] <- child.res
                      delivered.result <- delivered.result + 
                        1L
                    }
                    else {
                      fin[ci] <- TRUE
                      jobsp[ji] <- jobid[ji] <- NA
                      if (any(ava)) {
                        nexti <- which.max(ava[, ji])
                        if (!is.na(nexti)) {
                          jobid[ji] <- nexti
                          jobs[[ji]] <- mcparallel(FUN(X[[nexti]], 
                            ...), mc.set.seed = mc.set.seed, 
                            silent = mc.silent, mc.affinity = affinity.list[[nexti]])
                          jobsp[ji] <- parallel:::processID(jobs[[ji]])
                          ava[nexti, ] <- FALSE
                        }
                      }
                    }
                  }
            }
            nores <- length(X) - delivered.result
            if (nores > 0) 
                warning(sprintf(ngettext(nores, "%d parallel function call did not deliver a result", 
                  "%d parallel function calls did not deliver results"), 
                  nores), domain = NA)
        }
        if (has.errors) 
            warning(gettextf("%d function calls resulted in an error", 
                has.errors), domain = NA)
        return(res)
    }
    if (!is.null(affinity.list)) 
        warning("'mc.preschedule' must be false if 'affinity.list' is used")
    sindex <- lapply(seq_len(cores), function(i) seq(i, length(X), 
        by = cores))
    schedule <- lapply(seq_len(cores), function(i) X[seq(i, length(X), 
        by = cores)])
    ch <- list()
    res <- vector("list", length(X))
    names(res) <- names(X)
    cp <- rep(0L, cores)
    fin <- rep(FALSE, cores)
    dr <- rep(FALSE, cores)
    inner.do <- function(core) {
        S <- schedule[[core]]
        f <- parallel:::mcfork()
        if (isTRUE(mc.set.seed)) 
            parallel:::mc.advance.stream()
        if (inherits(f, "masterProcess")) {
            on.exit(mcexit(1L, structure("fatal error in wrapper code", 
                class = "try-error")))
            if (isTRUE(mc.set.seed)) 
                parallel:::mc.set.stream()
            if (isTRUE(mc.silent)) 
                closeStdout(TRUE)
            parallel:::sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
            parallel:::mcexit(0L)
        }
        jobs[[core]] <<- ch[[core]] <<- f
        cp[core] <<- parallel:::processID(f)
        NULL
    }
    job.res <- lapply(seq_len(cores), inner.do)
    ac <- cp[cp > 0]
    has.errors <- integer(0)
    while (!all(fin)) {
        s <- parallel:::selectChildren(ac[!fin], -1)
        if (is.null(s)) 
            break
        if (is.integer(s)) 
            for (ch in s) {
                a <- parallel:::readChild(ch)
                if (is.integer(a)) {
                  core <- which(cp == a)
                  fin[core] <- TRUE
                }
                else if (is.raw(a)) {
                  core <- which(cp == attr(a, "pid"))
                  job.res[[core]] <- ijr <- unserialize(a)
                  if (inherits(ijr, "try-error")){ 
                    has.errors <- c(has.errors, core)
                    if(stop.on.error)  #!!
                        stop("error in one of X = ", paste(schedule[[core]], collapse=", "), "\n", attr(ijr, "condition")$message) #!!
                  }
                  dr[core] <- TRUE
                }
                else if (is.null(a)) {
                  core <- which(cp == ch)
                  fin[core] <- TRUE
                }
            }
    }
    for (i in seq_len(cores)) {
        this <- job.res[[i]]
        if (inherits(this, "try-error")) {
            for (j in sindex[[i]]) res[[j]] <- this
        }
        else if (!is.null(this)) 
            res[sindex[[i]]] <- this
    }
    nores <- cores - sum(dr)
    if (nores > 0) 
        warning(sprintf(ngettext(nores, "scheduled core %s did not deliver a result, all values of the job will be affected", 
            "scheduled cores %s did not deliver results, all values of the jobs will be affected"), 
            paste(which(dr == FALSE), collapse = ", ")), domain = NA)
    if (length(has.errors)) {
        if (length(has.errors) == cores) 
            warning("all scheduled cores encountered errors in user code")
        else warning(sprintf(ngettext(has.errors, "scheduled core %s encountered error in user code, all values of the job will be affected", 
            "scheduled cores %s encountered errors in user code, all values of the jobs will be affected"), 
            paste(has.errors, collapse = ", ")), domain = NA)
    }
    res
}

.retryCatch <- function(expr, ..., maxAttempts = 3, warnAttempts = FALSE, nameFN = "FN", printInfo = NULL){
  currentAttempt <- 0
  completed <- FALSE
  while(!completed & currentAttempt <= maxAttempts){
    currentAttempt <- currentAttempt + 1
    if(currentAttempt > 1){
      message(nameFN, " : Error occured, attempting again (", currentAttempt - 1, " of ", maxAttempts, ")")
    }
    ###########################################################
    tryResult <- tryCatch({
      #########################################################
      #Try Catch Statement Here
      if(warnAttempts){
        out <- return(expr)
      }else{
        out <- suppressWarnings(return(expr))
      }
      #########################################################
      list(out = out, completed = TRUE)
    }, error = function(e){
      list(out = e, completed = FALSE)
    }, ...)
    ###########################################################
    completed <- tryResult$completed
  }
  if(!completed){
    message(nameFN, " : Error occured and could not be resolved after ", maxAttempts, " additional attempts!")
    if(!is.null(printInfo)){
      message("Error occured at ", printInfo)
    }
    print(tryResult[[1]])
    stop()
  }
  
  tryResult[[1]]

}



##########################################################################################
# SparseMatrix Utilities
##########################################################################################

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

.safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL){
  
  if(!is.null(subsetRows)){
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(length(idxNotIn), ncol = ncol(mat)))
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows,,drop=FALSE]
  }

  if(!is.null(subsetCols)){
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[,subsetCols,drop=FALSE]
  }

  mat

}

###################
# Binary
###################

.colBinarySums <- function(m){
  m@x <- 1 * (m@x > 0)
  Matrix::colSums(m)
}

.rowBinarySums <- function(m){
  m@x <- 1 * (m@x > 0)
  Matrix::rowSums(m)
}

###################
# rowVars
###################

.sparseRowVars <- function(m, rM = NULL){
  if(require("sparseMatrixStats", quietly = TRUE)){
    sparseMatrixStats::rowVars(m)
  }else{
    message("Using ArchR sparse rowVars recommended to install `sparseMatrixStats`!")
    if(is.null(rM)){
      rM <- Matrix::rowMeans(m)
    }
    computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  }
}

.sparesRowSds <- function(m){
  sqrt(.sparseRowVars(m=m))
}

###################
# Geo
###################

.expm1 <- function(x, cap = 100){
  x[x > cap] <- cap
  expm1(x)
}

.sparseRowGeoVars <- function(m, rM = NULL){
  m@x <- .expm1(m@x)
  .sparseRowVars(m, rM = rM)
}

.sparseRowGeoMeans <- function(m){
  m@x <- .expm1(m@x)
  Matrix::rowMeans(m)
}

###################
# Log2
###################

.log2p1 <- function(m){
  m@x <- log2(m@x + 1)
  m
}

.sparseRowLog2p1Means <- function(m){
  m@x <- log2(m@x + 1)
  Matrix::rowMeans(m)
}

.sparseColLog2p1Means <- function(m){
  m@x <- log2(m@x + 1)
  Matrix::rowMeans(m)  
}

.sparseRowLog2p1Vars <- function(m, rM = NULL){
  m@x <- log2(m@x + 1)
  .sparseRowVars(m, rM = rM)
}

###################
# Other
###################
.sparseRowVarsStd <- function(
  m = NULL, 
  rM = NULL, 
  expSd = NULL, 
  vmax = sqrt(ncol(m)),
  method = "Seurat"
  ){

  if(is.null(rM)){
    rM <- Matrix::rowMeans(m)
  }

  if(require("Seurat", quietly = TRUE) & tolower(method) == "seurat"){
    
    #Since these Seurat Functions Are Not Exported
    #We will have a backup which is slower but should give identical results
    o <- tryCatch({
    
      Seurat:::SparseRowVarStd(
        mat = m,
        mu = rM,
        sd = expSd,
        vmax = vmax,
        display_progress = FALSE
      )
    
    }, error = function(e){

      #Vals
      m@x <- m@x - rM[m@i + 1]
      m@x <- m@x / expSd[m@i + 1]
      m@x <- (pmin(m@x, vmax))^2

      #Compute Row Vars
      v <- Matrix::rowSums(m)

      #Determine 0s
      m@x <- rep(1, length(m@x))
      n <- (ncol(m) - Matrix::rowSums(m))
      v <- v + n * (rM / expSd)^2
      v[is.na(v)] <- 0
      o <- as.vector(v / (ncol(m)-1))
      o

    })

  }else{

    #Vals
    m@x <- m@x - rM[m@i + 1]
    m@x <- m@x / expSd[m@i + 1]
    m@x <- (pmin(m@x, vmax))^2

    #Compute Row Vars
    v <- Matrix::rowSums(m)

    #Determine 0s
    m@x <- rep(1, length(m@x))
    n <- (ncol(m) - Matrix::rowSums(m))
    v <- v + n * (rM / expSd)^2
    v[is.na(v)] <- 0
    o <- as.vector(v / (ncol(m)-1))

  }

  o

}

.sparseRowScale <- function(m, max = 10, method = "seurat", stats = FALSE){

  #Since these Seurat Functions Are Not Exported
  #We will have a backup which is slower but should give identical results
  if(stats){
    method <- "ArchR"
  }

  if(require("Seurat", quietly = TRUE) & tolower(method) == "seurat"){
    m <- tryCatch({
      rn <- rownames(m)
      cn <- colnames(m)
      m <- Seurat:::FastSparseRowScale(m, scale_max = max)
      colnames(m) <- cn
      rownames(m) <- rn
      m
    }, error = function(e){
      rM <- Matrix::rowMeans(m)
      rV <- .sparseRowVars(m, rM = rM)
      m <- m - rM
      m <- m / sqrt(rV)
      m@x[m@x > max] <- max
      m@x[m@x < -max] <- -max
      m
    })
  }else{
    rM <- Matrix::rowMeans(m)
    rV <- .sparseRowVars(m, rM = rM)
    m <- m - rM
    m <- m / sqrt(rV)
    m@x[m@x > max] <- max
    m@x[m@x < -max] <- -max
  }

  if(stats){
    list(m=m, rM = rM, rV=rV)
  }else{
    m
  }

}

.sparseRowGeoDisp <- function(
    m = NULL
  ){
  rM <- .sparseRowGeoMeans(m)
  rV <- .sparseRowGeoVars(m, rM = m)
  rD <- log(rV / rM)
  as.vector(rD)
}





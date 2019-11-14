####################################################################
# Reading fragments from Arrow Files
####################################################################

#' Read Fragments from Arrow
#' 
#' This function for each sample will independently compute counts for each feature
#' per cell in the Arrow File
#'
#' @param ArrowFile ArchRProject or ArrowFiles
#' @param chr GRanges to count for each cell
#' @param cellNames matrix output name in ArrowFiles cannot be a protected matrix name
#' @param method ceiling for the number of counts per feature
#' @param verbose binarize matrix
#' @param ... additional params
#' @export
getFragmentsFromArrow <- function(
  ArrowFile, 
  chr = NULL, 
  cellNames = NULL, 
  method = "fast",
  verbose = TRUE,
  ...){

  ArrowFile <- .validArrow(ArrowFile)

  if(is.null(chr)){
    chr <- .availableSeqnames(ArrowFile, subGroup = "Fragments")
  }

  if(any(chr %ni% .availableSeqnames(ArrowFile, subGroup = "Fragments"))){
    stop("Error Chromosome not in ArrowFile!")
  }
  
  tstart <- Sys.time()
  out <- lapply(seq_along(chr), function(x){
    .messageDiffTime(sprintf("Reading Chr %s of %s", x, length(chr)), tstart, verbose = verbose)
    .getFragsFromArrow(ArrowFile = ArrowFile, chr = chr[x], out = "GRanges", method = method)
  }) %>% GenomicRangesList

  .messageDiffTime("Merging", tstart, verbose = verbose)

  out <- .suppressAll(unlist(out))

  out

}

#' @export
.getFragsFromArrow <- function(
  ArrowFile, 
  chr = NULL, 
  out = "GRanges", 
  cellNames = NULL, 
  method = "fast",
  ...){

  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }

  o <- h5closeAll()
  ArrowFile <- .validArrow(ArrowFile)
  
  if(chr %ni% .availableSeqnames(ArrowFile)){
    stop("Error Chromosome not in ArrowFile!")
  }

  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)

  o <- h5closeAll()
  nFrags <- h5ls(ArrowFile, recursive = TRUE) %>% 
    {.[.$group==paste0("/Fragments/",chr) & .$name == "Ranges",]$dim} %>% 
    {gsub(" x 2","",.)} %>% as.integer

  if(nFrags==0){
    output <- IRanges(start = 1, end = 1)
    mcols(output)$RG <- c("tmp")
    output <- output[-1,]
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
    }
    return(output)
  }


  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }

  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){

      stop("None of input cellNames are in ArrowFile availableCells!")

    }else{

      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
          {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }

  }
  
  o <- h5closeAll()

  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }

  return(output)
}

####################################################################
# Reading Matrices/Arrays from Arrow Files
####################################################################

#' Read Fragments from Arrow
#' 
#' This function for each sample will independently compute counts for each feature
#' per cell in the Arrow File
#'
#' @param ArrowFile ArchRProject or ArrowFiles
#' @param useMatrix matrix name to get from Arrow
#' @param useSeqnames use a subset of seqnames for matrix
#' @param cellNames ceiling for the number of counts per feature
#' @param verbose binarize matrix
#' @param ... additional params
#' @export
getMatrixFromArrow <- function(
  ArrowFile, 
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  cellNames = NULL, 
  verbose = TRUE,
  ...){

  ArrowFile <- .validArrow(ArrowFile)

  seqnames <- .availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- .getFeatureDF(ArrowFile, subGroup = useMatrix)

  if(!is.null(useSeqnames)){
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }

  if(length(seqnames) == 0){
    stop("No seqnames available!")
  }

  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]

  mat <- .getMatFromArrow(
    ArrowFile = ArrowFile, 
    featureDF = featureDF, 
    cellNames = cellNames, 
    useMatrix = useMatrix,
    binarize = binarize,
    useIndex = FALSE
  )

  if(all(c("z", "deviations") %in% seqnames)){
    mat <- as(split(mat, featureDF$seqnames), "SimpleList")
    featureDF <- featureDF[featureDF$seqnames=="deviations", "name", drop = FALSE]
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix
  }

  colData <- .getMetadata(ArrowFile)

  if(useMatrix == "PeakMatrix"){
    se <- SummarizedExperiment(
      assays = mat,
      rowRanges = getPeakSet(ArchRProj),
      colData = colData[colnames(mat[[1]]),,drop=FALSE]
    )
  }else{
    se <- SummarizedExperiment(
      assays = mat,
      rowData = featureDF,
      colData = colData[colnames(mat[[1]]),,drop=FALSE]
    )
  }

  se

}

#' @export
.getMatFromArrow <- function(
  ArrowFile, 
  featureDF = NULL, 
  binarize = NULL, 
  cellNames = NULL,
  useMatrix = "TileMatrix", 
  useIndex = FALSE,
  threads = 1,
  ...
  ){

  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }

  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }

  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))

  o <- h5closeAll()

  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }

  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }

  seqnames <- unique(featureDF$seqnames)

  mat <- .safelapply(seq_along(seqnames), function(x){

    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx

    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
      )

    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]

    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }

    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)

    return(mat)

  }, threads = threads) %>% Reduce("rbind", .)

  o <- h5closeAll()

  colnames(mat) <- matColNames[idxCols]

  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL

  return(mat)

}


####################################################################
# Helper read functioning
####################################################################
#' @export
.getGroupMatrix <- function(
  ArrowFiles, 
  featureDF, 
  groupList,
  threads = 1, 
  useIndex = FALSE, 
  verbose = TRUE, 
  useMatrix = "TileMatrix",
  tstart = NULL,
  ...
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################
  seqnames <- unique(featureDF$seqnames)
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  cellNames <- unlist(groupList, use.names = FALSE)

  mat <- .safelapply(seq_along(seqnames), function(x){

    .messageDiffTime(sprintf("Constructing Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)

    #Construct Matrix
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex), ]

    matChr <- matrix(0, nrow = nrow(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDFx)

    for(y in seq_along(ArrowFiles)){
     
      maty <- .getMatFromArrow(
        ArrowFile = ArrowFiles[y], 
        useMatrix = useMatrix,
        featureDF = featureDFx, 
        cellNames = cellNames, 
        useIndex = useIndex
      )

      for(z in seq_along(groupList)){

        #Check Cells In Group
        cellsGroupz <- groupList[[z]]
        idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)

        #If In Group RowSums
        if(length(idx) > 0){
          matChr[,z] <- Matrix::rowSums(maty[,idx,drop=FALSE])
        }

      } 

    }

    matChr

  }, threads = threads) %>% Reduce("rbind", .)

  mat <- mat[rownames(featureDF), , drop = FALSE]
  
  .messageDiffTime("Successfully Created Group Matrix", tstart, verbose = verbose)

  return(mat)
  
}

#' @export
.getPartialMatrix <- function(
  ArrowFiles, 
  featureDF, 
  cellNames, 
  progress = TRUE, 
  threads = 1, 
  useMatrix = "TileMatrix",
  doSampleCells = FALSE, 
  sampledCellNames = NULL, 
  tmpPath = .tempfile(pattern = paste0("tmp-partial-mat")), 
  useIndex = FALSE,
  tstart = NULL,
  verbose = TRUE,
  ...
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################

  mat <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .messageDiffTime(sprintf("Getting Partial Matrix %s of %s", x, length(ArrowFiles)), tstart, verbose = verbose)

    o <- h5closeAll()
    matx <- .getMatFromArrow(
      ArrowFile = ArrowFiles[x], 
      featureDF = featureDF, 
      cellNames = cellNames,
      useMatrix = useMatrix, 
      useIndex = useIndex
    )
   
    if(doSampleCells){

      #Save Temporary Matrix
      outx <- paste0(tmpPath, "-", .sampleName(ArrowFiles[x]), ".rds")
      saveRDS(matx, outx, compress = FALSE)     

      #Sample Matrix 
      matx <- matx[, which(colnames(matx) %in% sampledCellNames),drop = FALSE]
      
      return(list(mat = matx, out = outx))

    }else{
      
      return(matx)

    }

  }, threads = threads)
  

  if(doSampleCells){

    matFiles <- lapply(mat, function(x) x[[2]]) %>% Reduce("c", .)
    mat <- lapply(mat, function(x) x[[1]]) %>% Reduce("cbind", .)
    mat <- mat[,sampledCellNames]

    .messageDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(list(mat = mat, matFiles = matFiles))

  }else{

    mat <- Reduce("cbind", mat)
    mat <- mat[,cellNames]
    
    .messageDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(mat)

  }


}

########################################################################
# Compute Summary Statistics!
########################################################################

#' @export 
.getRowSums <- function(ArrowFiles, seqnames, useMatrix, verbose = TRUE, tstart = NULL, filter0 = FALSE, threads = 1){
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  #Compute RowSums
  rowSumsDF <- .safelapply(seq_along(seqnames), function(x){
    o <- h5closeAll()
    chr <- seqnames[x]
    for(y in seq_along(ArrowFiles)){
      if(y == 1){
        sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", chr, "/rowSums"))
      }else{
        sumy1 <- h5read(ArrowFiles[y], paste0(useMatrix, "/", chr, "/rowSums"))
        #The way we designed sparse matrix holds true that the rows are in order every tile even in rS = 0!
        if(length(sumy1) > length(sumy)){
          sumy1[seq_along(sumy)] <- sumy1[seq_along(sumy)] + sumy
          sumy <- sumy1
        }else{
          sumy[seq_along(sumy1)] <- sumy[seq_along(sumy1)] + sumy1 
        }
      }
    }
    #Return Setup In Feature DF Format (seqnames, idx columns)
    DataFrame(seqnames = Rle(chr, lengths = length(sumy)), idx = seq_along(sumy), value = as.vector(sumy))
  }, threads = threads) %>% Reduce("rbind", .)
  if(filter0){
    rowSumsDF <- rowSumsDF[rowSumsDF$value > 0, ]
  }
  .messageDiffTime("Successfully Created RowSums DataFrame", tstart, verbose = verbose)
  return(rowSumsDF)
}


#' @export 
.getColSums <- function(ArrowFiles, seqnames, useMatrix, verbose = TRUE, tstart = NULL, threads = 1){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Compute ColSums
  cS <- .safelapply(seq_along(seqnames), function(x){
    
    lapply(seq_along(ArrowFiles), function(y){
      
      o <- h5closeAll()
      cSy <- h5read(ArrowFile, paste0(useMatrix, "/", seqnames[x], "/colSums"))
      rownames(cSy) <- .availableCells(ArrowFile, useMatrix)
      cSy
      
    }) %>% Reduce("rbind", .)
    
  }, threads = threads) %>% Reduce("cbind", .) %>% rowSums

  .messageDiffTime("Successfully Computed colSums", tstart, verbose = verbose)

  return(cS)

}

# h5read implementation for optimal reading
#' @export
.h5read <- function(file, name, method = "fast", index = NULL, start = NULL, block = NULL, count = NULL, ...){
  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count, ...)
  }
  o <- h5closeAll()
  return(res)
}




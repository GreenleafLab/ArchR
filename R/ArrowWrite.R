.initializeMat <- function(
  ArrowFile = NULL,
  Group = NULL,
  Class = "Double",
  Units = "none",
  cellNames = NULL,
  featureDF = NULL,
  params = NULL,
  date = Sys.Date(),
  force = FALSE
  ){

  #Add Group Entry of SparseMatrix Format
  #This Includes the following format
  #
  # Info 
  #   - Class - Sparse.Integer.Matrix = Sparse Matrix with Integer Entries
  #           - Sparse.Binary.Matrix = Sparse Matrix with Binary ie no x values
  #           - Sparse.Double.Matrix = Sparse Matrix with Double/Numeric Entries
  #           - Sparse.Assays.Matrix = Sparse Matrices with same feature names (same cell x feature)
  #                                    separated by different seqnames (ie assayNames)
  #   - CellNames ie Colnames
  #   - FeatureDF dataframe that describes the rows of each seqname
  #   - Params Params that are used for construction to be checked when comparing Arrows
  #   - Date Date of Creation
  # Chr1
  #   - i, j (as an Rle), x, and rowSums,colSums,rowVars,etc.
  # Chr2
  # Chr3
  # ..
  #

  if(!suppressMessages(h5createGroup(ArrowFile, paste0(Group)))){
    if(force){
      h5delete(ArrowFile, paste0(Group))
      h5createGroup(ArrowFile, paste0(Group))
    }else{
      stop("Matrix Group Already Exists! Set force = TRUE to overwrite!")
    }
  }
  o <- h5createGroup(ArrowFile, paste0(Group, "/Info"))

  if(tolower(Class)=="binary"){

    o <- h5write(obj = "Sparse.Binary.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))

  }else if(tolower(Class)=="integer"){

    o <- h5write(obj = "Sparse.Integer.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))

  }else if(tolower(Class)=="double"){

    o <- h5write(obj = "Sparse.Double.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))

  }else if(tolower(Class)=="assays"){

    o <- h5write(obj = "Sparse.Assays.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))

  }else{

    stop("Matrix Class Not Supported!")

  }

  ##########
  # Add Units To Class
  ##########
  if(!is.character(Units)){
    stop("Please provide Units as character when writing matrix to Arrow!")
  }
  o <- h5write(obj = Units, file = ArrowFile, name = paste0(Group, "/Info/Units"))  

  ##########
  # Cell Names in Arrow
  ##########
  splitNames <- stringr::str_split(cellNames, pattern = "#", simplify=TRUE)
  if(ncol(splitNames) > 2){
    stop("Found error with cell names containing multiple # characters!")
  }else{
    cellNames <- splitNames[,ncol(splitNames)]
  }
  o <- h5write(obj = cellNames, file = ArrowFile, name = paste0(Group,"/Info/CellNames"))

  ##########
  # FeatureDF in Arrow
  ##########
  df <- data.frame(featureDF, stringsAsFactors = FALSE)
  stopifnot(all(c("seqnames","idx") %in% colnames(featureDF)))
  o <- h5write(obj = df, file = ArrowFile, name = paste0(Group,"/Info/FeatureDF"))

  ##########
  # Parameters for Matrix for Validity in Arrow
  ##########  
  o <- h5write(obj = params, file = ArrowFile, name = paste0(Group,"/Info/Params"))

  ##########
  # Date of Creation
  ##########  
  o <- h5write(obj = paste0(date), file = ArrowFile, name = paste0(Group,"/Info/Date"))

  return(0)

}

.addMatToArrow <- function(
  mat = NULL, 
  ArrowFile = NULL,
  Group = NULL,
  binarize = FALSE, 
  addRowSums = FALSE,
  addColSums = FALSE,
  addRowMeans = FALSE,
  addRowVars = FALSE,
  addColMeansLog2 = FALSE,     #New
  addRowMeansLog2 = FALSE,     #New
  addRowVarsLog2 = FALSE,
  addRowMeansGeo = FALSE,      #New
  addRowVarsGeo = FALSE,       #New
  addRowSumsBinary = FALSE,    #New
  addColSumsBinary = FALSE,    #New
  addRowMeansLog2Norm = FALSE, #New
  addRowVarsLog2Norm = FALSE,  #New
  scaleTo = 10000,             #New
  colSm = NULL,                #New
  logFile = NULL,
  version = 2                  #New
  ){

  stopifnot(inherits(mat, "dgCMatrix"))

  checkCells <- .availableCells(ArrowFile, dirname(Group))
  if(!identical(paste0(colnames(mat)), paste0(checkCells))){
    .logThis(colnames(mat), "colnames", logFile=logFile)
    .logThis(checkCells, "cellNames", logFile=logFile)
    .logMessage(paste0("Mismatch = ", sum(colnames(mat) != checkCells)))
    .logMessage("CellNames in Matrix Group do not Match CellNames in Matrix Being Written!",logFile=logFile)
    stop("CellNames in Matrix Group do not Match CellNames in Matrix Being Written!")
  }

  #Create Group
  o <- h5closeAll()
  o <- h5createGroup(ArrowFile, Group)
  
  if(version == 1){

    #Convert Columns to Rle
    j <- Rle(findInterval(seq(mat@x)-1,mat@p[-1]) + 1)

    #Info
    lengthRle <- length(j@lengths)
    lengthI <- length(mat@i)

    #Create Data Set
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/i"), storage.mode = "integer", 
      dims = c(lengthI, 1), level = getArchRH5Level()))

    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jLengths"), storage.mode = "integer", 
      dims = c(lengthRle, 1), level = getArchRH5Level()))

    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jValues"), storage.mode = "integer", 
      dims = c(lengthRle, 1), level = getArchRH5Level()))
    
    #Write Data Set
    o <- .suppressAll(h5write(obj = mat@i + 1, file = ArrowFile, name = paste0(Group,"/i")))
    o <- .suppressAll(h5write(obj = j@lengths, file = ArrowFile, name = paste0(Group,"/jLengths")))
    o <- .suppressAll(h5write(obj = j@values, file = ArrowFile, name = paste0(Group,"/jValues")))
    
    #If binary dont store x
    if(!binarize){

      o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/x"), storage.mode = "double", 
        dims = c(lengthI, 1), level = getArchRH5Level()))

      o <- .suppressAll(h5write(obj = mat@x, file = ArrowFile, name = paste0(Group, "/x")))

    }else{

      mat@x[mat@x > 0] <- 1

    }

    rm(j)

  }else if(version == 2){

    #Info
    lengthP <- length(mat@p)
    lengthI <- length(mat@i)

    if(binarize){
      mat@x[mat@x > 0] <- 1
    }

    #Write Data Set
    o <- .suppressAll(h5write(obj = mat@i, file = ArrowFile, name = paste0(Group,"/indices")))
    o <- .suppressAll(h5write(obj = mat@p, file = ArrowFile, name = paste0(Group,"/indptr")))
    o <- .suppressAll(h5write(obj = mat@x, file = ArrowFile, name = paste0(Group, "/data")))
    o <- .suppressAll(h5write(obj = dim(mat), file = ArrowFile, name = paste0(Group, "/shape")))

  }else{

    stop("ArchR Write Method Does Not Exist! Only versions = 1 and 2 exist!")

  }

  ####################
  # Normal
  ####################

  cS <- 0
  rS <- 0
  rM <- 0
  rV <- 0

  if(addColSums){
    cS <- Matrix::colSums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/colSums"), storage.mode = "double", 
      dims = c(ncol(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = cS, file = ArrowFile, name = paste0(Group, "/colSums")))
  }

  if(addRowSums){
    rS <- Matrix::rowSums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowSums"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rS, file = ArrowFile, name = paste0(Group, "/rowSums")))
  }

  if(addRowMeans | addRowVars){
    rM <- Matrix::rowMeans(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeans"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rM, file = ArrowFile, name = paste0(Group, "/rowMeans")))
  }

  if(addRowVars){
    rV <- .sparseRowVars(m=mat, rM=rM)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVars"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rV, file = ArrowFile, name = paste0(Group, "/rowVars")))
  }

  rm(cS, rS, rM, rV)

  ####################
  # Binary
  ####################

  cSB <- 0
  rSB <- 0

  if(addColSumsBinary){
    cSB <- .colBinarySums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/colSumsBinary"), storage.mode = "double", 
      dims = c(ncol(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = cSB, file = ArrowFile, name = paste0(Group, "/colSumsBinary")))

  }

  if(addRowSumsBinary){
    rSB <- .rowBinarySums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowSumsBinary"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rSB, file = ArrowFile, name = paste0(Group, "/rowSumsBinary")))
  }

  rm(cSB, rSB)

  ####################
  # Log2
  ####################

  cMLog <- 0
  rMLog <- 0
  rVLog <- 0

  if(addColMeansLog2){
    cMLog <- .sparseColLog2p1Means(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/colMeansLog2"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = cMLog, file = ArrowFile, name = paste0(Group, "/colMeansLog2")))
  }

  if(addRowMeansLog2 | addRowVarsLog2){
    rMLog <- .sparseRowLog2p1Means(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeansLog2"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rMLog, file = ArrowFile, name = paste0(Group, "/rowMeansLog2")))
  }

  if(addRowVarsLog2){
    rVLog <- .sparseRowLog2p1Vars(mat, rM = rMLog)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVarsLog2"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rVLog, file = ArrowFile, name = paste0(Group, "/rowVarsLog2")))
  }

  rm(cMLog, rMLog, rVLog)

  ####################
  # Geo
  ####################

  rMGeo <- 0
  rVGeo <- 0

  if(addRowMeansGeo | addRowVarsGeo){
    rMGeo <- .sparseRowGeoMeans(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeansGeo"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rMGeo, file = ArrowFile, name = paste0(Group, "/rowMeansGeo")))
  }

  if(addRowVarsGeo){
    rVGeo <- .sparseRowGeoVars(mat, rM = rMGeo)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVarsGeo"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rVGeo, file = ArrowFile, name = paste0(Group, "/rowVarsGeo")))
  }

  rm(rMGeo, rVGeo)

  ####################
  # Log2 Depth Norm
  ####################

  rMLogNorm <- 0
  rVLogNorm <- 0

  if(addRowMeansLog2Norm | addRowVarsLog2Norm){
    mat <- .normalizeCols(mat, colSm = colSm, scaleTo = scaleTo)
    rMLogNorm <- .sparseRowLog2p1Means(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeansLog2Norm"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rMLogNorm, file = ArrowFile, name = paste0(Group, "/rowMeansLog2Norm")))
  }

  if(addRowVarsLog2Norm){
    rVLogNorm <- .sparseRowLog2p1Vars(mat, rM = rMLogNorm)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVarsLog2Norm"), storage.mode = "double", 
      dims = c(nrow(mat), 1), level = getArchRH5Level()))
    o <- .suppressAll(h5write(obj = rVLogNorm, file = ArrowFile, name = paste0(Group, "/rowVarsLog2Norm")))
  }

  rm(rMLogNorm, rVLogNorm)

  ####################
  #Clean Up Memorys
  ####################
  rm(mat)
  gc()

  o <- h5closeAll()

  return(0)

}


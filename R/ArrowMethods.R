####################################################################
# Hidden Helper Utils for Arrow Files
####################################################################

#' @export
.validArrow <- function(ArrowFile){
  o <- h5closeAll()
  if(h5read(ArrowFile,"Class")!="Arrow"){
    stop("Not Valid Arrow!")
  }
  o <- h5closeAll()
  return(ArrowFile)
}

#' @export
.isProtectedArray <- function(matrixName, exclude = NULL){
  protectedArrays <- tolower(c("peakmatrix", "tilematrix", "genescorematrix"))
  if(!is.null(exclude)){
    protectedArrays <- protectedArrays[protectedArrays %ni% tolower(exclude)]
  }
  if(tolower(matrixName) %in% protectedArrays){
    stop(sprintf("Error %s cannot be used as this conflicts with another predefined matrix function!", matrixName))
  }
  matrixName
}

#' @export
.availableSeqnames <- function(ArrowFiles, subGroup = "Fragments"){
  o <- h5closeAll()
  seqList <- lapply(seq_along(ArrowFiles), function(x){
    seqnames <- h5ls(ArrowFiles[x]) %>% {.[.$group==paste0("/",subGroup),]$name}
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  })
  if(!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],seqList[[1]]))))){
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(seqList[[1]])
}

#' @export
.availableChr <- function(ArrowFiles, subGroup = "Fragments"){
  seqnames <- .availableSeqnames(ArrowFiles, subGroup)
  seqnames <- seqnames[grep("chr", seqnames, ignore.case = TRUE)]
  if(length(seqnames) == 0){
    stop("No Chr Found in ArrowFiles!")
  }
  return(seqnames)
}

#' @export
.availableCells <- function(ArrowFile, subGroup = NULL, passQC = TRUE){
  if(is.null(subGroup)){
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, "Metadata/CellNames")
    if(passQC){
      passQC <- tryCatch({
        h5read(ArrowFile, "Metadata/PassQC")
      }, error = function(x){
        rep(1, length(cellNames))
      })
      cellNames <- cellNames[which(passQC==1)]
    }
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }else{
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, paste0(subGroup, "/Info/CellNames"))
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  return(paste0(sampleName,"#",cellNames))
}

#' @export
.sampleName <- function(ArrowFile){
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}

#' @export
.summarizeArrowContent <- function(ArrowFile){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  h5DF <- h5ls(ArrowFile)

  #Re-Organize Content Info
  h5DF <- h5DF[-which(h5DF$group == "/"),]
  groups <- stringr::str_split(h5DF$group, pattern = "/", simplify=TRUE)[,2]
  groupList <- split(h5DF, groups)

  #Split Nested Lists
  groupList2 <- lapply(seq_along(groupList), function(x){
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]),"",groupDFx$group)
    if(all(groupx=="")){
      groupDFx
    }else{
      subDF <- groupDFx[-which(groupx == ""),]
      split(subDF, stringr::str_split(subDF$group, pattern = "/", simplify=TRUE)[,3])
    }
  })
  names(groupList2) <- names(groupList)


  o <- h5closeAll()

  return(groupList2)

}

#' @export
.getMetadata <- function(ArrowFile){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- .summarizeArrowContent(ArrowFile)$Metadata

  #Which are same dimensions as cell names
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name=="CellNames"]),]

  #Load these into a S4 DataFrame
  md <- lapply(seq_len(nrow(arrowMD)), function(x){
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x],"/",arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)

  #Correct CellNames
  md$CellNames <- paste0(sampleName,"#",md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md)=="CellNames")]
  md <- md[,order(colnames(md))]

  o <- h5closeAll()

  return(md)
}

#' @export
.getFeatureDF <- function(ArrowFiles, subGroup = "TileMatrix"){
  
  .helpFeatureDF <- function(ArrowFile, subGroup){
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup,"/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }

  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)

  if(length(ArrowFiles) > 1){
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- lapply(seq_along(ArrowFiles), function(x){
        fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
        identical(fdfx, fdf)
    }) %>% unlist %>% all
    if(!checkIdentical){
      stop("Error not all FeatureDF for asssay is the same!")
    }
  }

  #Re-Order for Split Check!
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% {lapply(seq_along(.), function(x) .[[x]])} %>% Reduce("c", .)
  fdf[newOrder,]
  
}


# exportShiny function -----------------------------------------------------------
#' Export a Shiny App based on ArchRProj
#' 
#' Generate all files required for an autonomous Shiny app to display browser tracks and embeds.
#'
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outputDir The name of the directory for the Shiny App files. 
#' @param groupBy The name of the column in cellColData to use for grouping cells together for generating sequencing tracks. Only one cell grouping is allowed.
#'        defaults to "Clusters".
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportShinyArchR <- function(
  ArchRProj = NULL,
  outputDir = "Shiny",
  subOutputDir = "inputData",
  groupBy = "Clusters",
  embedding = "UMAP",
  tileSize = 100,
  threads = getArchRThreads(),
  logFile = createLogFile("exportShinyArchR")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outputDir, name = "outputDir", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')
  
  mainDir <- getOutputDirectory(ArchRProj)

  allMatrices <- getAvailableMatrices(ArchRProj)
  matrices <- list()
  imputeMatricesList <- list()
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProj)
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  
  if(!file.exists(file.path(outputDir, subOutputDir, "matrices.rds")) && !file.exists(file.path(outputDir, subOutputDir, "imputeMatricesList.rds"))){
    for(allmatrices in allMatrices){
      print(allmatrices)
      name <- paste0(allmatrices, "_names")
      result = assign(name, getFeatures(ArchRProj = ArchRProj, useMatrix = allmatrices))
      saveRDS(result, paste0("./", outputDir, "/", subOutputDir, "/", allmatrices,"_names.rds"))
    
      if(!is.null(result)){
        # nameColor <- paste0("colorMat", allmatrices)
        matrix = Matrix(.getMatrixValues(
          ArchRProj = ArchRProj, 
          name = result,
          matrixName = allmatrices,
          log2Norm = FALSE,
          threads = threads), sparse = TRUE)
        
        matrices[[allmatrices]] = matrix
        
        matList = matrix[,rownames(df), drop=FALSE]
        
        # assign(nameColor, matList)
        .logThis(matList, paste0(allmatrices,"-Before-Impute"), logFile = logFile)
        if(getArchRVerbose()) message("Imputing Matrix")
        
        colorImputeMat <- imputeMatrix(mat = as.matrix(matList), imputeWeights = imputeWeights, logFile = logFile)
        # assign(paste0(nameColor, "_Impute"), imputeMat)
        
        if(!inherits(colorImputeMat, "matrix")){
          colorImputeMat <- matrix(colorImputeMat, ncol = nrow(df))
          colnames(colorImputeMat) <- rownames(df)
        }
        imputeMatricesList[[allmatrices]] <- colorImputeMat
        
        
      }else{
        message(allmatrices, " is NULL.")
      }
    }
    saveRDS(matrices, paste0(file.path(mainDir, outputDir, subOutputDir, "matrices.rds")))
    saveRDS(imputeMatricesList, paste0(file.path(mainDir, outputDir, subOutputDir, "/imputeMatricesList.rds")))
  }else{
    
    message("matrices and imputeMatricesList already exist. reading from local files...")
    
    matrices <- readRDS(file.path(mainDir, outputDir, subOutputDir, "matrices.rds"))
    imputeMatricesList <- readRDS(file.path(mainDir, outputDir, subOutputDir, "/imputeMatricesList.rds"))
  }
  
  if(is.null(groupBy)){
    stop("groupBy must be provided")
  } else if(groupBy %ni% colnames(getCellColData(ArchRProj))){
    stop("groupBy must be a column in cellColData")
  }else{
    print(paste0("groupBy:", groupBy))
  }

  # Check that the embedding exists in ArchRProj@embeddings
  if(embedding %ni% names(ArchRProj@embeddings)){
    stop("embedding doesn't exist in ArchRProj@embeddings")
  }else{
    print(paste0("embedding:", embedding))
  }
  
  # Make directory for Shiny App 
  if(!dir.exists(outputDir)) {
    
    dir.create(file.path(mainDir, outputDir), showWarnings = TRUE)
    
    ## Check the links for the files
    filesUrl <- data.frame(
      fileUrl = c(
        "https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/app.R",
        "https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/global.R", 
        "https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/server.R", 
        "https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/ui.R"
      ),
      md5sum = c(
        "77502e1f195e21d2f7a4e8ac9c96e65e",
        "618613b486e4f8c0101f4c05c69723b0",
        "a8d5ae747841055ef230ba496bcfe937"
      ),
      stringsAsFactors = FALSE
    )
    
    .downloadFiles(filesUrl = filesUrl, pathDownload = file.path(mainDir, outputDir), threads = threads)
    
  }else{
    message("Using existing Shiny files...")
  }
  
  # Create a copy of the ArchRProj 
  ArchRProjShiny <- ArchRProj

  # Add metadata to ArchRProjShiny
  if (groupBy %ni% colnames(ArchRProjShiny@cellColData)) {
    stop("groupBy is not part of cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("Some entries in the column indicated by groupBy have NA values. Please subset your project using subsetArchRProject() to only contain cells with values for groupBy")
  } else {
    ArchRProjShiny@projectMetadata[["groupBy"]] <- groupBy
  }
  ArchRProjShiny@projectMetadata[["tileSize"]] <- tileSize
  units <- tryCatch({
    .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
  },error=function(e){
    "values"
  })
  ArchRProjShiny@projectMetadata[["units"]] <- units
  ArchrProjShiny <- saveArchRProject(ArchRProj = ArchRProjShiny, outputDirectory = 
  file.path(mainDir, outputDir, "Save-ArchRProjShiny"), dropCells = TRUE, overwrite = TRUE)
  
  # Create fragment files 
  if(length(list.files(file.path(outputDir, "fragments"))) == 0){
    .getGroupFragsFromProj(ArchRProj = ArchRProj, groupBy = groupBy, outDir = file.path(outputDir, "fragments"))
  }else{ 
    message("Fragment files already exist...")
  }
  
  # Create coverage objects
  if(length(list.files(file.path(outputDir, "coverage"))) == 0){
    .getClusterCoverage(ArchRProj = ArchRProj, tileSize = tileSize, groupBy = groupBy, outDir = file.path(outputDir, "coverage"))
  }else{
    message("Coverage files already exist...")
  }

  # Create directory to save everything that will be preprocessed 
  dir.create(file.path(mainDir, outputDir, subOutputDir), showWarnings = TRUE)
 
  if(!file.exists(file.path(mainDir, outputDir, subOutputDir, "features.rds"))){
    gene_names <- getFeatures(ArchRProj = ArchRProj)
    saveRDS(gene_names, file.path(mainDir, outputDir, subOutputDir, "features.rds"))
  }else{
    message("gene_names already exists...")
    gene_names <- readRDS(file.path(mainDir, outputDir, subOutputDir, "features.rds"))
  }
  
  if(!file.exists(file.path(mainDir, outputDir, subOutputDir, "/embeddingMaps.rds"))){  
    embeddingMaps <- list()
    embedNames <- colnames(ArchRProjShiny@cellColData)[][colnames(ArchRProjShiny@cellColData) %in% groupBy]
    
    embeddingMaps <- .safelapply(1:length(embedNames), function(x){
      print(embedNames[x])
      tryCatch({
        embed <- plotEmbedding(
          ArchRProj = ArchRProjShiny,
          baseSize=12,
          colorBy = "cellColData",
          name = embedNames[x],
          embedding = embedding,
          rastr = FALSE,
          size=0.5,
          matrices = matrices,
          imputeMatricesList = imputeMatricesList,
        )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12),
                                                          legend.title = element_text(size = 12),legend.text = element_text(size = 6))
        
        embeddingMaps[[embedNames[[x]]]] <- embed
      },
      error = function(e){
        print(e)
      })
    })
    
    saveRDS(embeddingMaps, paste0("./", outputDir, "/", subOutputDir,"/embeddingMaps.rds"))
    
    
  }else{
    message("embeddingMaps already exists...")
    embeddingMaps <- readRDS(paste0(outputDir, "/", subOutputDir,"/embeddingMaps.rds"))
  }
  

# Create an HDF5 containing the nativeRaster vectors for the main matrices
if (!file.exists(file.path(outputDir, subOutputDir, "mainEmbeds.h5"))) {
  
  if(groupBy %in% colnames(ArchRProjShiny@cellColData)){
  
  mainEmbed(ArchRProj = ArchRProj,
            outDirEmbed = file.path(outputDir, subOutputDir),
            names = groupBy,
            matrices =  matrices,
            imputeMatricesList = imputeMatricesList,
            Shiny = ShinyArchR
          )          
  }else{
    
    message(groupBy, "is not defined in ArchRProj...")
    
  }
  
} else{
  message("H5 for main embeds already exists...")
}


if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/plotBlank72.h5"))){
  
  matrixEmbeds(
    ArchRProj = ArchRProj,
    outputDirEmbeds = paste0(outputDir,"/", subOutputDir),
    embedding = embedding,
    matrices = matrices,
    imputeMatricesList = imputeMatricesList,
    threads = getArchRThreads(),
    verbose = TRUE,
    logFile = createLogFile("matrixEmbeds")
  )
  
}else{
  
  message("H5 file already exists...")
  
}  
## delete unnecessary files -----------------------------------------------------------------
unlink("./fragments", recursive = TRUE) 
unlink("./ArchRLogs", recursive = TRUE) 

## ready to launch ---------------------------------------------------------------
message("App created! To launch, 
          ArchRProj <- loadArchRProject('",getOutputDirectory(ArchRProj),"') and 
          run shiny::runApp('", outputDir, "') from parent directory")
#  runApp("myappdir")


}


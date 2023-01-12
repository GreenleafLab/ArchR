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
  if(length(list.files(file.path(mainDir, outputDir, "coverage"))) == 0){
    .getClusterCoverage(ArchRProj = ArchRProj, tileSize = tileSize, groupBy = groupBy, outDir = file.path(mainDir, outputDir, "coverage"))
  }else{
    message("Coverage files already exist...")
  }

  # Create directory to save input data to Shinyapps.io (everything that will be preprocessed)
  dir.create(file.path(mainDir, outputDir, subOutputDir), showWarnings = TRUE)
 
  # if(!file.exists(file.path(mainDir, outputDir, subOutputDir, "features.rds"))){
  #   gene_names <- getFeatures(ArchRProj = ArchRProj)
  #   saveRDS(gene_names, file.path(mainDir, outputDir, subOutputDir, "features.rds"))
  # }else{
  #   message("gene_names already exists...")
  #   gene_names <- readRDS(file.path(mainDir, outputDir, subOutputDir, "features.rds"))
  # }

  allMatrices <- getAvailableMatrices(ArchRProj)
  matrices <- list()
  imputeMatrices <- list()
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProj)
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  
  if(!file.exists(file.path(outputDir, subOutputDir, "matrices.rds")) && !file.exists(file.path(outputDir, subOutputDir, "imputeMatrices.rds"))){
    for(matName in allMatrices){
      matFeaturesNames <- paste0(matName, "_names")
      result = assign(matFeaturesNames, getFeatures(ArchRProj = ArchRProj, useMatrix = matName))
      saveRDS(result, file.path(outputDir, subOutputDir, matName, "_names.rds"))
    
      if(!is.null(result)){
        
         mat = Matrix(.getMatrixValues(
          ArchRProj = ArchRProj, 
          name = result,
          matrixName = mat,
          log2Norm = FALSE,
          threads = threads), sparse = TRUE)
        
        matrices[[matName]] = mat
        matList = mat[,rownames(df), drop=FALSE]
        .logThis(matList, paste0(matName,"-Before-Impute"), logFile = logFile)

        if(getArchRVerbose()) message("Imputing Matrix")
        imputeMat <- imputeMatrix(mat = as.matrix(matList), imputeWeights = imputeWeights, logFile = logFile)

        if(!inherits(imputeMat, "matrix")){
          imputeMat <- mat(imputeMat, ncol = nrow(df))
          colnames(imputeMat) <- rownames(df)
        }
        imputeMatrices[[matName]] <- imputeMat
        
        
      }else{
        message(matName, " is NULL.")
      }
    }
    matrices$allColorBy= .availableArrays(head(getArrowFiles(ArchRProj), 2))
    saveRDS(matrices, file.path(outputDir, subOutputDir, "matrices.rds"))
    saveRDS(imputeMatrices, file.path(outputDir, subOutputDir, "imputeMatrices.rds"))
  }else{
    
    message("matrices and imputeMatrices already exist. reading from local files...")
    
    matrices <- readRDS(file.path(mainDir, outputDir, subOutputDir, "matrices.rds"))
    imputeMatrices <- readRDS(file.path(mainDir, outputDir, subOutputDir, "imputeMatrices.rds"))
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
  
# mainEmbed will create an HDF5 containing the nativeRaster vectors for the main matrices
if (!file.exists(file.path(mainDir, outputDir, subOutputDir, "mainEmbeds.h5"))) {
  mainEmbed(ArchRProj = ArchRProj,
            outDirEmbed = file.path(mainDir, outputDir, subOutputDir),
            colorBy = "cellColData",
            names = groupBy,
            embeddingDF = df,
            matrices =  matrices,
            imputeMatrices = imputeMatrices,
            Shiny = TRUE
          )          
} else{
  message("H5 for main embeddings already exists...")
}

if(!file.exists(file.path(outputDir, subOutputDir, "plotBlank72.h5"))){
  
  matrixEmbeds(
    ArchRProj = ArchRProj,
    outputDirEmbed = file.path(mainDir, outputDir, subOutputDir),
    embedding = embedding,
    matrices = matrices,
    imputeMatrices = imputeMatrices,
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
          ArchRProj <- loadArchRProject('", mainDir,"') and 
          run shiny::runApp('", outputDir, "') from parent directory")
#  runApp("myappdir")

}


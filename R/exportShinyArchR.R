# exportShiny function -----------------------------------------------------------
#' Export a Shiny App based on ArchRProj
#' 
#' Generate all files required for an autonomous Shiny app to display browser tracks and UMAPs.
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
    
    dir.create(outputDir)

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
    
    .downloadFiles(filesUrl = filesUrl, pathDownload = outputDir, threads = threads)

  }else{
    message("Using existing Shiny files...")
  }
  
  # Create a copy of the ArchRProj object
  ArchRProjShiny <- ArchRProj
  # Add metadata to ArchRProjShiny
  if (groupBy %ni% colnames(ArchRProjShiny@cellColData)) {
    stop("groupBy is not part of cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("Some entries in the column indicated by groupBy have NA values. 
          This is not allowed. Please subset your project using subsetArchRProject() to only contain cells with values for groupBy")
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
  
  # The following gives error: Error in file.copy(oldPath, outputDirectory, recursive = TRUE, overwrite = overwrite) : 
  # attempt to copy a directory to itself (That is why I commented it out)
  # ArchrProjShiny <- saveArchRProject(ArchRProj = ArchRProjShiny, outputDirectory = "Save-ArchRProjShiny", dropCells = TRUE, overwrite = TRUE)
  
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

  dir.create(file.path(getOutputDirectory(ArchRProj), outputDir, subOutputDir),showWarnings = TRUE)

  ## colorMats without Impute Weights----------------------------------------------------------------
  
  # check which matrices are available
  allMatrices <- getAvailableMatrices(ArchRProj)

  # Get gene and motif names and save as RDS
  if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/gene_names_GSM.rds"))){
    if ("GeneScoreMatrix" %in% allMatrices){
      gene_names_GSM <- getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneScoreMatrix")
      saveRDS(gene_names_GSM, paste0("./", outputDir, "/", subOutputDir,"/gene_names_GSM.rds"))
    }else{
      message("GeneScoreMatrix does not exist...")
    }
  }else{
    message("gene_names_GSM already exists...")
    gene_names_GSM <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/gene_names_GSM.rds"))
  }
  
  
if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/gene_names_GIM.rds"))){
  if ("GeneIntegrationMatrix" %in% allMatrices){
    gene_names_GIM <- getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneIntegrationMatrix")
    saveRDS(gene_names_GSM, paste0("./", outputDir, "/", subOutputDir,"/gene_names_GIM.rds"))
  }else{
    
    message("GeneIntegrationMatrix does not exist...")
  }
}else{
    message("gene_names_GIM already exists...")
    gene_names_GIM <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/gene_names_GIM.rds"))
  }
  

  if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/motif_names.rds"))){
    if ("MotifMatrix" %in% allMatrices){
    motif_names <- getFeatures(ArchRProj = ArchRProj, useMatrix = "MotifMatrix") %>% 
      gsub(".*:", "", .) %>% unique(.)
    saveRDS(motif_names, paste0("./", outputDir, "/", subOutputDir,"/motif_names.rds"))
    }else{
      
      message("MotifMatrix does not exist...")
      
    }
  }else{
    
    message("motif_names already exists...")
    motif_names <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/motif_names.rds"))
  }
  
  if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/matrices.rds"))){
    matrices <- list()
    #GSM colorMat
    colorMatGSM <- Matrix(.getMatrixValues(
      ArchRProj = ArchRProj, 
      name = gene_names_GSM,
      matrixName = "GeneScoreMatrix",
      log2Norm = FALSE,
      threads = threads,
    ), sparse = TRUE)
    matrices$"GeneScoreMatrix" <- colorMatGSM 
    
    #GIM
    colorMatGIM <- Matrix(.getMatrixValues(
      ArchRProj = ArchRProj,
      name = gene_names_GIM,
      matrixName = "GeneIntegrationMatrix",
      log2Norm = FALSE,
      threads = threads
    ),sparse = TRUE)
    matrices$"GeneIntegrationMatrix" <- colorMatGIM
    
    colorMatMM <- Matrix(.getMatrixValues(
      ArchRProj = ArchRProj,
      #name = getFeatures(ArchRProj, "MotifMatrix")
      name =  paste0("deviations:", motif_names), #used deviations:
      matrixName = "MotifMatrix",
      log2Norm = FALSE,
      threads = threads
    ), sparse = TRUE)
    matrices$"MotifMatrix" <- colorMatMM
    
    matrices$allColorBy=c(.availableArrays(head(getArrowFiles(ArchRProj), 2)))
    saveRDS(matrices, paste0("./", outputDir, "/", subOutputDir,"/matrices.rds"))
    
  }else{
    message("matrices already exist...")
    matrices <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/matrices.rds"))
  }
  ## Impute Weights ------------------------------------------------------------
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProj)
  if(!is.null(imputeWeights)) {
    df <- getEmbedding(ArchRProj, embedding = "UMAP", returnDF = TRUE)
    
    if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/imputeMatricesList.rds"))){
      imputeMatricesList <- list()
      # colorMats for each colorBy
      
      # GSM
      # colorMatGSM <- matrices$"GeneScoreMatrix"
      # colorMatGSM <- colorMatGSM[,rownames(df), drop=FALSE]
      colorMatGSM <- matrices[["GeneScoreMatrix"]][,rownames(df), drop=FALSE]
      
      
      
      .logThis(colorMatGSM, "colorMatGSM-Before-Impute", logFile = logFile)
      
      if(getArchRVerbose()) message("Imputing Matrix")
      colorMatGSM_Impute <- imputeMatrix(mat = as.matrix(colorMatGSM), imputeWeights = imputeWeights, logFile = logFile)
      if(!inherits(colorMatGSM_Impute, "matrix")){
        colorMatGSM_Impute <- matrix(colorMatGSM_Impute, ncol = nrow(df))
        colnames(colorMatGSM_Impute) <- rownames(df)
      }
      
      # .logThis(colorMat_Impute, "colorMatGSM-After-Impute", logFile = logFile)
      
      imputeMatricesList$"GeneScoreMatrix" <- colorMatGSM_Impute
      
      # GIM
      # colorMatGIM <- matrices$"GeneIntegrationMatrix"
      # colorMatGIM <- colorMatGIM[,rownames(df), drop=FALSE]
      colorMatGIM <- matrices[["GeneIntegrationMatrix"]][,rownames(df), drop=FALSE]
      
      
      .logThis(colorMatGIM, "colorMatGIM-Before-Impute", logFile = logFile)
      
      if(getArchRVerbose()) message("Imputing Matrix")
      colorMatGIM_Impute <- imputeMatrix(mat = as.matrix(colorMatGIM), imputeWeights = imputeWeights, logFile = logFile)
      if(!inherits(colorMatGIM_Impute, "matrix")){
        colorMatGIM_Impute <- matrix(colorMatGIM_Impute, ncol = nrow(df))
        colnames(colorMatGIM_Impute) <- rownames(df)
      }
      
      .logThis(colorMatGIM_Impute, "colorMatGIM-After-Impute", logFile = logFile)
      
      imputeMatricesList$"GeneIntegrationMatrix" <- colorMatGIM_Impute
      
      # Motif Matrix
      # colorMatMM <- matrices$"MotifMatrix"
      # colorMatMM <- colorMatMM[,rownames(df), drop=FALSE]
      colorMatMM <- matrices[["MotifMatrix"]][,rownames(df), drop=FALSE]
      
      
      .logThis(colorMatMM, "colorMatMM-Before-Impute", logFile = logFile)
      
      if(getArchRVerbose()) message("Imputing Matrix")
      colorMatMM_Impute <- imputeMatrix(mat = as.matrix(colorMatMM), imputeWeights = imputeWeights, logFile = logFile)
      if(!inherits(colorMatMM_Impute, "matrix")){
        colorMatMM_Impute <- matrix(colorMatMM_Impute, ncol = nrow(df))
        colnames(colorMatMM_Impute) <- rownames(df)
      }
      
      .logThis(colorMatMM_Impute, "colorMatMM-After-Impute", logFile = logFile)
      
      imputeMatricesList$"MotifMatrix" <- colorMatMM_Impute
      
      saveRDS(imputeMatricesList,paste0("./", outputDir, "/", subOutputDir,"/imputeMatricesList.rds"))
    }else{
      message("imputeMatricesList already exists...")
      imputeMatricesList <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/imputeMatricesList.rds"))
    }
  }
    # Create an HDF5 containing the nativeRaster vectors for the main matrices
    if (!file.exists(file.path(outputDir, subOutputDir, "mainEmbeds.h5"))) {
      
      mainEmbed(ArchRProj = ArchRProj,
                outDirEmbed = file.path(outputDir, subOutputDir),
                names = colnames(ArchRProjShiny@cellColData),
                matrices =  matrices,
                imputeMatricesList = imputeMatricesList,
                Shiny = ShinyArchR
                )
    } else{
      message("H5 for main embeds already exists...")
    }
    
    
    if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/plotBlank72.h5"))){
      
      shinyRasterUMAPs(
        ArchRProj = ArchRProj,
        outputDirUmaps = paste0(outputDir,"/", subOutputDir),
        threads = getArchRThreads(),
        verbose = TRUE,
        logFile = createLogFile("ShinyRasterUMAPs")
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


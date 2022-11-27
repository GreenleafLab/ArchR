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
  
  # Check that all columns exist in cellColData 
  if(is.null(groupBy)){
    stop("groupBy must be provided")
  } else if(name %ni% colnames(getCellColData(ArchRProj))){
    stop("groupBy must be a column in cellColData")
  }
  # Check that the embedding exists in ArchRProj@embeddings
  if(name %in% names(ArchRProj@embeddings)){
    stop("embedding doesn't exist in ArchRProj@embeddings")
  }

  # Make directory for Shiny App 
  if(!dir.exists(outputDir)) {
    
    dir.create(outputDir)
    
    filesUrl <- data.frame(
      fileUrl = c(
      https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/app.R 
      https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/global.R 
      https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/server.R 
      https://jeffgranja.s3.amazonaws.com/ArchR/Shiny/ui.R
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
  ArchrProjShiny <- saveArchRProject(ArchRProj = ArchRProjShiny, outputDirectory = "Save-ArchRProjShiny", dropCells = TRUE, overwrite = FALSE)
  
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
  ## main umaps -----------------------------------------------------------------

  dir.create(file.path(getOutputDirectory(ArchRProj), outputDir, "inputData")"showWarnings = FALSE")

  # need arrowFiles to getFeatures so need to save genes as RDS
  # TODO change hardcoding of these paths: Should be file.path(getOutputDirectory(ArchRProj), outputDir, "inputData", "features.rds"
  if(!file.exists(file.path(getOutputDirectory(ArchRProj), outputDir, "inputData", "features.rds"))){
    gene_names <- getFeatures(ArchRProj = ArchRProj)
    saveRDS(gene_names, "./Shiny/data/inputData/gene_names.rds")
  }else{
    message("gene_names already exists...")
    gene_names <- readRDS("./Shiny/data/inputData/gene_names.rds")
    
  }
  
  if(!file.exists("./Shiny/data/inputData/umaps.rds")){  
    umaps <- list()
    cluster_umap <- plotEmbedding(
      ArchRProj = ArchRProjShiny,
      baseSize=12,
      colorBy = groupBy,
      name = "Clusters",
      embedding = embedding,
      rastr = FALSE,
      size=0.5,
    )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12),
                                                      legend.title = element_text(size = 12),legend.text = element_text(size = 6))
    umaps[["Clusters"]] <- cluster_umap
    
    sample_umap <- plotEmbedding(
      ArchRProj = ArchRProj,
      baseSize=12,
      colorBy = "cellColData",
      name = "Sample",
      embedding = "UMAP",
      rastr = FALSE,
      size=0.5
    )+ ggtitle("Colored by original identity")+theme(text=element_text(size=12),
                                                     legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    umaps[["Sample"]] <- sample_umap
    
    constrained_umap <- plotEmbedding(
      ArchRProj = ArchRProjShiny,
      colorBy = "cellColData",
      name = "predictedGroup_Co",
      rastr = FALSE,
      baseSize=12,
      size=0.5
    )+ggtitle("UMAP: constrained integration")+theme(text=element_text(size=12),
                                                     legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    umaps[["Constrained"]] <- constrained_umap
    
    unconstrained_umap <- plotEmbedding(
      ArchRProj = ArchRProjShiny,
      embedding = "UMAP",
      colorBy = "cellColData",
      name = "predictedGroup_Un",
      baseSize=12,
      rastr = FALSE,
      size=0.5
    )+ggtitle("UMAP: unconstrained integration")+theme(text=element_text(size=12),
                                                       legend.title = element_text(size = 12),legend.text = element_text(size = 6))
    # saveRDS(unconstrained_umap, "./UMAPs/unconstrained_umap.rds")
    umaps[["unconstrained"]] <- unconstrained_umap
    
    constrained_remapped_umap <- plotEmbedding(
      ArchRProj = ArchRProjShiny,
      colorBy = "cellColData",
      name = "Clusters2",
      rastr = FALSE,
    )+ggtitle("UMAP: Constrained remapped clusters")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    # saveRDS(constrained_remapped_umap, "./UMAPs/constrained_remapped_umap.rds")
    umaps[["Constrained remap"]] <- constrained_remapped_umap
    
    saveRDS(umaps, "./Shiny/data/inputData/umaps.rds")}else{
      message("umaps already exists...")
      umaps <- readRDS("./Shiny/data/inputData/umaps.rds")
    }
  
  ## colorMats without Impute Weights----------------------------------------------------------------
  
  #TODO check with matrices are available
  allMatrices <- getAvailableMatrices(ArchRProj)

  #Get gene and motif names and save as RDS
  if(!file.exists("./Shiny/data/inputData/gene_names_GSM.rds")){
    # TODO check if ArchRProj has GSM
    if ("GeneScoreMatrix" %in% allMatrices){
      gene_names_GSM <- getFeatures(ArchRProj = ArchRProj, matrix = "GeneScoreMatrix")
      saveRDS(gene_names_GSM, "./Shiny/data/inputData/gene_names_GSM.rds")
    {
      # else skip to checking if next matrix exists and output message saying no GSM
    }
  }else{
    message("gene_names_GSM already exists...")
    gene_names_GSM <- readRDS("./Shiny/data/inputData/gene_names_GSM.rds")
  }
  
  if(!file.exists("./Shiny/data/inputData/gene_names_GIM.rds")){
     # TODO check if ArchRProj has GIM
    gene_names_GIM <- getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneIntegrationMatrix")
    saveRDS(gene_names_GIM, "./Shiny/data/inputData/gene_names_GIM.rds")
  }else{
    message("gene_names_GIM already exists...")
    gene_names_GIM <- readRDS("./Shiny/data/inputData/gene_names_GIM.rds")
  }
  
  if(!file.exists("./Shiny/data/inputData/motif_names.rds")){
     # TODO check if ArchRProj has MM
    motif_names <- getFeatures(ArchRProj = ArchRProj, useMatrix = "MotifMatrix") %>% 
      gsub(".*:", "", .) %>% unique(.)
    saveRDS(motif_names, "./Shiny/data/inputData/motif_names.rds")
  }else{
    
    message("motif_names already exists...")
    motif_names <- readRDS("./Shiny/data/inputData/motif_names.rds")
  }
  
  if(!file.exists("./Shiny/data/inputData/matrices.rds")){
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
      name = getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneIntegrationMatrix"),
      matrixName = "GeneIntegrationMatrix",
      log2Norm = FALSE,
      threads = threads
    ),sparse = TRUE)
    matrices$"GeneIntegrationMatrix" <- colorMatGIM
    
    #colorMatMM has 1740 rows because in name = getFeatures() returns the 870 z: + the 870 deviations:
    colorMatMM <- Matrix(.getMatrixValues(
      ArchRProj = ArchRProj,
      #name = getFeatures(ArchRProj, "MotifMatrix")
      name = paste0("deviations:", motif_names), #used deviations:
      matrixName = "MotifMatrix",
      log2Norm = FALSE,
      threads = threads
    ), sparse = TRUE)
    matrices$"MotifMatrix" <- colorMatMM
    
    #TODO modify this so it only has the matrices we are actually supporting
    matrices$allColorBy=c("colData", "cellColData", .availableArrays(head(getArrowFiles(ArchRProj), 2)))
    # shouldn't save rds because it's too hefty for ShinyApps
    saveRDS(matrices,"./Shiny/data/inputData/matrices.rds")
    
  }else{
    message("matrices already exist...")
    matrices <- readRDS("./Shiny/data/inputData/matrices.rds")
  }
  ## Impute Weights ------------------------------------------------------------
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProj)
  if(!is.null(imputeWeights)) {
    df <- getEmbedding(ArchRProj, embedding = "UMAP", returnDF = TRUE)
    
    
    if(!file.exists("./Shiny/data/inputData/imputeMatricesList.rds")){
      imputeMatricesList <- list()
      # colorMats for each colorBy
      
      # GSM
      colorMatGSM <- matrices$"GeneScoreMatrix"
      colorMatGSM <- colorMatGSM[,rownames(df), drop=FALSE]
      
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
      colorMatGIM <- matrices$"GeneIntegrationMatrix"
      colorMatGIM <- colorMatGIM[,rownames(df), drop=FALSE]
      
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
      colorMatMM <- matrices$"MotifMatrix"
      colorMatMM <- colorMatMM[,rownames(df), drop=FALSE]
      
      .logThis(colorMatMM, "colorMatMM-Before-Impute", logFile = logFile)
      
      if(getArchRVerbose()) message("Imputing Matrix")
      colorMatMM_Impute <- imputeMatrix(mat = as.matrix(colorMatMM), imputeWeights = imputeWeights, logFile = logFile)
      if(!inherits(colorMatMM_Impute, "matrix")){
        colorMatMM_Impute <- matrix(colorMatMM_Impute, ncol = nrow(df))
        colnames(colorMatMM_Impute) <- rownames(df)
      }
      
      .logThis(colorMatMM_Impute, "colorMatMM-After-Impute", logFile = logFile)
      
      imputeMatricesList$"MotifMatrix" <- colorMatMM_Impute
      
      saveRDS(imputeMatricesList,"./Shiny/data/inputData/imputeMatricesList.rds")
    }else{
      message("imputeMatricesList already exists...")
      imputeMatricesList <- readRDS("./Shiny/data/inputData/imputeMatricesList.rds")
    }
    
    if(!file.exists("./Shiny/data/inputData/plotBlank72.h5")){
      
      rasterUmaps(
        ArchRProj = ArchRProj,
        FDR = 0.01,
        Log2FC = 1.25,
        outputDir = "Shiny",
        groupBy = "Clusters",
        tileSize = 100,
        markerList = NULL,
        threads = getArchRThreads(),
        verbose = TRUE,
        logFile = createLogFile("exportShinyArchR")
      )
      
    }else{
      
      message("H5 file already exists...")
      
    }  
    ## delete unnecessary files -----------------------------------------------------------------
    unlink("./fragments", recursive = TRUE) 
    unlink("./ArchRLogs", recursive = TRUE) 
    
    ## ready to launch ---------------------------------------------------------------
    message("App created! To launch, 
          ArchRProj <- loadArchRProject('path to ArchRProject/') and 
          run shiny::runApp('", outputDir, "') from parent directory")
    #  runApp("myappdir")
  }  
}


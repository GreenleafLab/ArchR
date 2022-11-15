#' Export a Shiny App based on ArchRProj
#'
#' Generate all files required for an autonomous Shiny app to display your browser tracks.
#'
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outputDir The name of the directory for the Shiny App files.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for summarizing.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportShinyArchR <- function(ArchRProj = NULL,
                             outputDir = "Shiny",
                             groupBy = "Clusters",
                             tileSize = 100,
                             threads = getArchRThreads(),
                             verbose = TRUE,
                             logFile = createLogFile("exportShinyArchR")) {
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outputDir, name = "outputDir", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')
  
  # Make directory for Shiny App 
  if(!dir.exists(outputDir)) {
    
    dir.create(outputDir)
    # if(length(dir(outDir,  all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
    #   stop("Please specify a new or empty directory")
    # }
    
    filesUrl <- c(
      "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/app.R",
      "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/global.R",
      "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/server.R",
      "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/ui.R"
    )
    
    downloadFiles <- lapply(seq_along(filesUrl), function(x){
      download.file(
        url = filesUrl[x], 
        destfile = file.path(outputDir, basename(filesUrl[x]))
      )
    })
    
  }else{
    message("Using existing Shiny files...")
  }
  
  # Create a copy of the ArchRProj object
  ArchRProjShiny <- ArchRProj
  # Add metadata to ArchRProjShiny
  if (is.na(paste0("ArchRProj$", groupBy))) {
    stop("groupBy is not part of cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("incomplete data. some NA observations for groupBy")
  } else {
    ArchRProjShiny@projectMetadata[["groupBy"]] <- groupBy
  }
  ArchRProjShiny@projectMetadata[["tileSize"]] <- tileSize
  saveArchRProject(ArchRProj = ArchRProj, outputDirectory = "Save-ArchRProjShiny")
  
  # Create fragment files 
  .getGroupFragsFromProj(ArchRProj = ArchRProjShiny, groupBy = groupBy, outDir = file.path(outputDir, "fragments"))
  
  # Create coverage objects
  .getClusterCoverage(ArchRProj = ArchRProjShiny, tileSize = tileSize, groupBy = groupBy, outDir = file.path(outputDir, "coverage"))
  
  ## main umaps -----------------------------------------------------------------
  dir.create("UMAPs") 
  
  units <- tryCatch({
    .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
  },error=function(e){
    "values"
  })
  
  ArchRProjShiny@projectMetadata[["units"]] <- units
  
  # need arrowFiles to getFeatures so need to save genes as RDS
  gene_names <- getFeatures(ArchRProj = ArchRProj)
  saveRDS(gene_names, "./inputData/gene_names.rds")
  
  umaps <- list()
  cluster_umap <- plotEmbedding(
    ArchRProj = ArchRProjShiny,
    baseSize=12,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP",
    rastr = FALSE,
    size=0.5,
  )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12), 
                                                    legend.title = element_text(size = 12),legend.text = element_text(size = 6))
  umaps[["Clusters"]] <- cluster_umap
  # saveRDS(cluster_umap, "./UMAPs/cluster_umap.rds")
  
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
  # saveRDS(sample_umap, "./UMAPs/sample_umap.rds")
  
  constrained_umap <- plotEmbedding(
    ArchRProj = ArchRProjShiny,
    colorBy = "cellColData",
    name = "predictedGroup_Co",
    rastr = FALSE,
    baseSize=12,
    size=0.5
  )+ggtitle("UMAP: constrained integration")+theme(text=element_text(size=12), 
                                                   legend.title = element_text( size = 12),legend.text = element_text(size = 6))
  # saveRDS(constrained_umap, "./UMAPs/constrained_umap.rds")
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
  
  saveRDS(umaps, "./inputData/umaps.rds")
  umaps <- readRDS("./inputData/umaps.rds")
  
  ## colorMats without Impute Weights----------------------------------------------------------------
  
  # Get gene and motif names and save as RDS
  gene_names_GSM <- getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneScoreMatrix")
  saveRDS(gene_names_GSM, file="./inputData/geneNamesGSM.rds")
  
  gene_names_GIM <- getFeatures(ArchRProj = ArchRProj, useMatrix = "GeneIntegrationMatrix")
  saveRDS(gene_names_GIM, file = "./inputData/geneNamesGIM.rds")
  
  motif_names <- getFeatures(ArchRProj = ArchRProj, useMatrix = "MotifMatrix") %>% 
    gsub(".*:", "", .) %>% unique(.)
  saveRDS(motif_names, "./inputData/markerListMM.rds")
  
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
    name = paste0("deviations:", markerListMM), #used deviations:
    matrixName = "MotifMatrix",
    log2Norm = FALSE,
    threads = threads
  ), sparse = TRUE)
  matrices$"MotifMatrix" <- colorMatMM
  
  # TODO modify this so it only has the matrices we are actually supporting
  matrices$allColorBy=c("colData", "cellColData", .availableArrays(head(getArrowFiles(ArchRProj), 2)))
  
  ## Impute Weights ------------------------------------------------------------
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProj)
  if(!is.null(imputeWeights)) {
    df <- getEmbedding(ArchRProj, embedding = "UMAP", returnDF = TRUE)
    
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
    
    .logThis(colorMat_Impute, "colorMatGSM-After-Impute", logFile = logFile)
    
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
  }
  
  saveRDS(imputeMatricesList,"~/tests/raster/inputData/imputeMatricesList.rds")
  imputeMatricesList <- readRDS("imputeMatricesList.rds")
  
  ## delete unnecessary files -----------------------------------------------------------------
  unlink("./fragments", recursive = TRUE) 
  unlink("./ArchRLogs", recursive = TRUE) 
  
  ## ready to launch ---------------------------------------------------------------
  message("App created! To launch, 
          ArchRProj <- loadArchRProject('path to ArchRProject/') and 
          run shiny::runApp('", outputDir, "') from parent directory")
  #  runApp("myappdir")
}



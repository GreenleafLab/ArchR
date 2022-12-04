# Setting up ----------------------------------------------------------------------

library(shinycssloaders)
library(hexbin)
library(magick)
library(gridExtra)
library(grid)
library(patchwork)
library(shinybusy)
library(cowplot)
library(ggpubr)
library(farver)
library(rhdf5)
library(plotfunctions)
library(raster)
library(jpeg)
library(sparseMatrixStats)
library(BiocManager)
# options(repos = BiocManager::repositories())
library(AnnotationDbi)
library(BSgenome)
library(Biobase)
library(BiocGenerics)
library(BiocParallel)
library(Biostrings)
library(CNEr)
library(ComplexHeatmap)
# options(download.file.method = "libcurl")
# devtools::install_github("selcukorkmaz/ArchR", ref = "dev")
library(ArchR)

# specify whether you use a local machine or the shiny app
ShinyArchR = TRUE

# specify desired number of threads 
addArchRThreads(threads = 1) 
# specify genome version. Default hg19 set
addArchRGenome("hg19")
set.seed(1)

ArchRProj=loadArchRProject(path = "Save-ProjHeme5/")
ArchRProj <- addImputeWeights(ArchRProj = ArchRProj)
setwd(getOutputDirectory(ArchRProj))

# myLoadArchRProject -----------------------------------
#' Load Previous ArchRProject into R
#'
#' This function will load a previously saved ArchRProject and re-normalize paths for usage.
#'
#' @param path A character path to an `ArchRProject` directory that was previously saved using `saveArchRProject()`.
#' @param force A boolean value indicating whether missing optional `ArchRProject` components (i.e. peak annotations /
#' background peaks) should be ignored when re-normalizing file paths. If set to `FALSE` loading of the `ArchRProject`
#' will fail unless all components can be found.
#' @param showLogo A boolean value indicating whether to show the ascii ArchR logo after successful creation of an `ArchRProject`.
#' @export
myLoadArchRProject <- function(path = "./",
                               force = FALSE,
                               showLogo = TRUE) {
  .validInput(input = path,
              name = "path",
              valid = "character")
  .validInput(input = force,
              name = "force",
              valid = "boolean")
  .validInput(input = showLogo,
              name = "showLogo",
              valid = "boolean")
  
  path2Proj <- file.path(path, "Save-ArchR-Project.rds")
  
  if (!file.exists(path2Proj)) {
    stop("Could not find previously saved ArchRProject in the path specified!")
  }
  
  ArchRProj <- recoverArchRProject(readRDS(path2Proj))
  outputDir <- getOutputDirectory(ArchRProj)
  outputDirNew <- normalizePath(path)
  
  
  ArchRProj@projectMetadata$outputDirectory <- outputDirNew
  
  message("Successfully loaded ArchRProject!")
  if (showLogo) {
    .ArchRLogo(ascii = "Logo")
  }
  
  ArchRProj
  
}


## Create fragment files -----------------------------------------------------------
.getGroupFragsFromProj <- function(ArchRProj = NULL,
                                   groupBy = NULL,
                                   outDir = file.path("Shiny", "fragments")) {
  dir.create(outDir, showWarnings = FALSE)
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters.
  clusters <- names(cellGroups)
  
  
  for (cluster in clusters) {
    cat("Making fragment file for cluster:", cluster, "\n")
    # get GRanges with all fragments for that cluster
    cellNames = cellGroups[[cluster]]
    fragments <-
      getFragmentsFromProject(ArchRProj = ArchRProj, cellNames = cellNames)
    fragments <- unlist(fragments, use.names = FALSE)
    # filter Fragments
    fragments <-
      GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
    saveRDS(fragments, file.path(outDir, paste0(cluster, "_cvg.rds")))
  }
}

addSeqLengths <- function (gr, genome) {
  gr <- ArchR:::.validGRanges(gr)
  genome <- validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <-
    seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

.getClusterCoverage <- function(ArchRProj = NULL,
                                tileSize = 100,
                                scaleFactor = 1,
                                groupBy = "Clusters",
                                outDir = file.path("Shiny", "coverage")) {
  fragfiles = list.files(path = file.path("Shiny", "fragments"),
                         full.names = TRUE)
  dir.create(outDir, showWarnings = FALSE)
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters.
  clusters <- names(cellGroups)
  
  chrRegions <- getChromSizes(ArchRProj)
  genome <- getGenome(ArchRProj)
  
  for (file in fragfiles) {
    fragments <- readRDS(file)
    #fragmentsToInsertions()
    left <- GRanges(seqnames = seqnames(fragments),
                    ranges = IRanges(start(fragments), width = 1))
    right <- GRanges(seqnames = seqnames(fragments),
                     ranges = IRanges(end(fragments), width = 1))
    # call sort() after sortSeqlevels() to sort also the ranges in addition
    # to the chromosomes.
    insertions <- c(left, right) %>% sortSeqlevels() %>%
      sort()
    
    cluster <- file %>% basename() %>% gsub("_.*", "", .)
    #binnedCoverage
    # message("Creating bins for cluster ",clusters[clusteridx], "...")
    bins <-
      unlist(slidingWindows(chrRegions, width = tileSize, step = tileSize))
    # message("Counting overlaps for cluster ",clusters[clusteridx], "...")
    bins$reads <-
      countOverlaps(
        bins,
        insertions,
        maxgap = -1L,
        minoverlap = 0L,
        type = "any"
      )
    addSeqLengths(bins, genome)
    # message("Creating binned coverage for cluster ",clusters[clusteridx], "...")
    #each value is multiplied by that weight.
    # TODO add scaleFactor
    # allCells as.vector(ArchRProj@cellColData$Sample, mode="any")
    clusterReadsInTSS <-
      ArchRProj@cellColData$ReadsInTSS[cells %in% cellGroups$cluster]
    # scaleFactor <- 5e+06 / sum(clusterReadsInTSS)
    binnedCoverage <-
      coverage(bins, weight = bins$reads * scaleFactor)
    saveRDS(binnedCoverage, file.path(outDir, paste0(cluster, "_cvg.rds")))
  }
  
}


#############################################################

# ArchRProj=myLoadArchRProject("./Shiny/inputData/")


# Load all hidden ArchR functions ------------------------------------------------
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
  tryCatch({
    eval(parse(text = paste0(fn[i], "<-", fn[i])))
  }, error = function(x) {
  })
}

# UMAP Visualization ------------------------------------------------------------





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
  
  # ArchRProj <- myLoadArchRProject(path = paste0("./",outputDir,"/inputData/"))
  # ArchRProj <- addImputeWeights(ArchRProj = ArchRProj)

  # TODO: Check that all columns exist in cellColData 
  if(is.null(groupBy)){
    stop("groupBy must be provided")
  } 
  else if(groupBy %ni% colnames(getCellColData(ArchRProj))){
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
  } 
  else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
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
  ## main umaps -----------------------------------------------------------------

  dir.create(file.path(getOutputDirectory(ArchRProj), outputDir, subOutputDir),showWarnings = FALSE)

  # need arrowFiles to getFeatures so need to save genes as RDS
  # TODO change hardcoding of these paths: Should be file.path(getOutputDirectory(ArchRProj), outputDir, subOutputDir, "features.rds"
  if(!file.exists(file.path(getOutputDirectory(ArchRProj), outputDir, subOutputDir, "gene_names.rds"))){
    gene_names <- getFeatures(ArchRProj = ArchRProj)
    saveRDS(gene_names, paste0("./", outputDir, "/", subOutputDir,"/gene_names.rds"))
  }else{
    message("gene_names already exists...")
    gene_names <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/gene_names.rds"))
    
  }
  
  if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/umaps.rds"))){  
    umaps <- list()
    umapNames <- colnames(ArchRProjShiny@cellColData)
    
    for(x in 1:length(umapNames)){
      tryCatch(
        umap <- plotEmbedding(
          ArchRProj = ArchRProjShiny,
          baseSize=12,
          colorBy = "cellColData",
          name = umapNames[x],
          embedding = embedding,
          rastr = FALSE,
          size=0.5,
        )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12),
                                                          legend.title = element_text(size = 12),legend.text = element_text(size = 6)),
        
        umaps[[umapNames[[x]]]] <- umap,
        error = function(e){
        print(e)
        })
    }
    
    saveRDS(umaps, paste0("./", outputDir, "/", subOutputDir,"/umaps.rds"))
    
    
    }else{
      message("umaps already exists...")
      umaps <- readRDS(paste0("./", outputDir, "/", subOutputDir,"/umaps.rds"))
    }
    
    # cluster_umap <- plotEmbedding(
    #   ArchRProj = ArchRProjShiny,
    #   baseSize=12,
    #   colorBy = "cellColData",
    #   name = "TSSEnrichment",
    #   embedding = embedding,
    #   rastr = FALSE,
    #   size=0.5,
    # )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12),
    #                                                   legend.title = element_text(size = 12),legend.text = element_text(size = 6))
    # umaps[["Clusters"]] <- cluster_umap
    # 
    # sample_umap <- plotEmbedding(
    #   ArchRProj = ArchRProj,
    #   baseSize=12,
    #   colorBy = "cellColData",
    #   name = "Sample",
    #   embedding = "UMAP",
    #   rastr = FALSE,
    #   size=0.5
    # )+ ggtitle("Colored by original identity")+theme(text=element_text(size=12),
    #                                                  legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    # umaps[["Sample"]] <- sample_umap
    # 
    # constrained_umap <- plotEmbedding(
    #   ArchRProj = ArchRProjShiny,
    #   colorBy = "cellColData",
    #   name = "predictedGroup_Co",
    #   rastr = FALSE,
    #   baseSize=12,
    #   size=0.5
    # )+ggtitle("UMAP: constrained integration")+theme(text=element_text(size=12),
    #                                                  legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    # umaps[["Constrained"]] <- constrained_umap
    # 
    # unconstrained_umap <- plotEmbedding(
    #   ArchRProj = ArchRProjShiny,
    #   embedding = "UMAP",
    #   colorBy = "cellColData",
    #   name = "predictedGroup_Un",
    #   baseSize=12,
    #   rastr = FALSE,
    #   size=0.5
    # )+ggtitle("UMAP: unconstrained integration")+theme(text=element_text(size=12),
    #                                                    legend.title = element_text(size = 12),legend.text = element_text(size = 6))
    # # saveRDS(unconstrained_umap, "./UMAPs/unconstrained_umap.rds")
    # umaps[["Unconstrained"]] <- unconstrained_umap
    # 
    # constrained_remapped_umap <- plotEmbedding(
    #   ArchRProj = ArchRProjShiny,
    #   colorBy = "cellColData",
    #   name = "Clusters2",
    #   rastr = FALSE,
    # )+ggtitle("UMAP: Constrained remapped clusters")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    # # saveRDS(constrained_remapped_umap, "./UMAPs/constrained_remapped_umap.rds")
    # umaps[["Constrained remap"]] <- constrained_remapped_umap
    # 
   
  
  ## colorMats without Impute Weights----------------------------------------------------------------
  
  #TODO check with matrices are available
  allMatrices <- getAvailableMatrices(ArchRProj)

  #Get gene and motif names and save as RDS
  if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/gene_names_GSM.rds"))){
    # TODO check if ArchRProj has GSM
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
  # TODO check if ArchRProj has GIM
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
     # TODO check if ArchRProj has MM
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
    
    #colorMatMM has 1740 rows because in name = getFeatures() returns the 870 z: + the 870 deviations:
    colorMatMM <- Matrix(.getMatrixValues(
      ArchRProj = ArchRProj,
      #name = getFeatures(ArchRProj, "MotifMatrix")
      name =  paste0("deviations:", motif_names), #used deviations:
      matrixName = "MotifMatrix",
      log2Norm = FALSE,
      threads = threads
    ), sparse = TRUE)
    matrices$"MotifMatrix" <- colorMatMM
    
    #TODO modify this so it only has the matrices we are actually supporting
    matrices$allColorBy=c("colData", "cellColData", .availableArrays(head(getArrowFiles(ArchRProj), 2)))
    # shouldn't save rds because it's too hefty for ShinyApps
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

    # Create an HDF5 containing the nativeRaster vectors for the main matrices
    if (!file.exists(file.path(outputDir, subOutputDir, "mainEmbeds.h5"))) {
      
      mainEmbed(ArchRProj = ArchRProj,
                outDirEmbed = file.path(outputDir, subOutputDir),
                names = as.list(colnames(ArchRProjShiny@cellColData))
                )
    } else{
      message("H5 for main embeds already exists...")
    }
    
    
    if(!file.exists(paste0("./", outputDir, "/", subOutputDir,"/plotBlank72.h5"))){
      
      shinyRasterUMAPs(
        ArchRProj = NULL,
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
          ArchRProj <- loadArchRProject('",getwd(),"') and 
          run shiny::runApp('", outputDir, "') from parent directory")
    #  runApp("myappdir")
  }  

}


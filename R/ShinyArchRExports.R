# Functions for exporting an ArchR-based Shiny app -----------------------------------------------------------
#'
#' Export a Shiny App based on an ArchRProj
#' 
#' Generate all files required for an autonomous Shiny app to display browser tracks and embeddings.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param outputDir The name (not the path!) of the directory for the Shiny App files. This will become a sub-directory of the ArchRProj output directory
#' given by `getOutputDirectory(ArchRProj)`.
#' @param subOutputDir Where data to upload to Shinyapps.io will be stored. Defaults to `inputData`.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for generating BigWig-style sequencing tracks. 
#' Only one cell grouping is allowed.
#' @param cellColEmbeddings A character vector of columns in `cellColData` to plot as part of the Shiny app. No default is provided so this must be set.
#' For ex. `c("Sample","Clusters","TSSEnrichment","nFrags")`.
#' @param embedding The name of the embedding from `ArchRProj` to be used for plotting embeddings in the Shiny app.
#' @param matsToUse A character vector containing the matrices that you want to include in the Shiny app. This should be used to limit
#' which matrices are included in the app. Matrices listed here must exist in your project (see `getAvailableMatrices()`).
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param force A boolean value that indicates whether to overwrite any relevant files during the `exportShinyArchR()` process.
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#' 
#' proj <- getTestProject(version = 2)
#' proj@geneAnnotation$genes <- proj@geneAnnotation$genes[which(proj@geneAnnotation$genes$symbol %in% c("CD14","CD3D","MS4A1","CD74"))]
#' proj <- addGeneScoreMatrix(input = proj, force = TRUE)
#' proj <- addImputeWeights(proj)
#'
#' exportShinyArchR(ArchRProj = proj,
#'                  mainDir = "Shiny",
#'                  subOutDir = "inputData",
#'                  savedArchRProjFile = "Save-ArchR-Project.rds",
#'                  groupBy = "Clusters",
#'                  cellColEmbeddings = "Clusters",
#'                  embedding = "UMAP",
#'                  matsToUse = "GeneScoreMatrix",
#'                  tileSize = 100,
#'                  force = FALSE,
#'                  threads = getArchRThreads(),
#'                  logFile = createLogFile("exportShinyArchR"))
#' 
#' @export
exportShinyArchR <- function(
  ArchRProj = NULL,
  mainDir = "Shiny",
  subOutDir = "inputData",
  savedArchRProjFile = "Save-ArchR-Project.rds",
  groupBy = "Clusters",
  cellColEmbeddings = "Clusters",
  embedding = "UMAP",
  matsToUse = NULL,
  tileSize = 100,
  force = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("exportShinyArchR")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = mainDir, name = "mainDir", valid = c("character"))
  .validInput(input = subOutDir, name = "subOutDir", valid = c("character"))
  .validInput(input = savedArchRProjFile, name = "savedArchRProjFile", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = cellColEmbeddings, name = "cellColEmbeddings", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = matsToUse, name = "matsToUse", valid = c("character","null"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')
  
  if(length(groupBy) > 1){
    stop("Only one value is allowed for groupBy.")
  }
  
  if(!all(cellColEmbeddings %in% colnames(ArchRProj@cellColData))){
    stop("Not all entries in cellColEmbeddings exist in the cellColData of your ArchRProj. Please check provided inputs.")
  }
  
  # Check that the embedding exists in ArchRProj@embeddings
  if(embedding %ni% names(ArchRProj@embeddings)){
    stop("embedding doesn't exist in ArchRProj@embeddings")
  }
  
  #check that groupBy column exists and doesn't have NA values
  if (groupBy %ni% colnames(ArchRProj@cellColData)) {
    stop("groupBy is not an entry in cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("Some entries in the column indicated by groupBy have NA values. Please subset your project using subsetArchRProject() to only contain cells with values for groupBy.")
  }
  
  supportedMatrices <- c("GeneScoreMatrix", "GeneIntegrationMatrix", "MotifMatrix") #only these matrices are currently supported for ShinyArchR
  #subset matrices for use in Shiny app
  allMatrices <- getAvailableMatrices(ArchRProj)
  if(!is.null(matsToUse)){
    if(!all(matsToUse %in% allMatrices)){
      stop("Not all matrices defined in matsToUse exist in your ArchRProject. See getAvailableMatrices().")
    } else {
      allMatrices <- allMatrices[which(allMatrices %in% matsToUse)]
    }
  }

  # get directories paths
  projDir <- getOutputDirectory(ArchRProj)
  mainOutputDir <- file.path(projDir, mainDir)
  subOutputDir <- file.path(projDir, mainDir, subOutDir)
  
  
  # Make directory for Shiny App  and download the app, global, server, and ui files if they dont already exist
  dir.create(mainOutputDir, showWarnings = FALSE)
  
  ## Check the links for the files
  filesUrl <- data.frame(
    fileUrl = c(
      "https://files.corces.gladstone.org/Users/rcorces/ArchR/Shiny/1.0.3/app.R",
      "https://files.corces.gladstone.org/Users/rcorces/ArchR/Shiny/1.0.3/global.R",
      "https://files.corces.gladstone.org/Users/rcorces/ArchR/Shiny/1.0.3/server.R",
      "https://files.corces.gladstone.org/Users/rcorces/ArchR/Shiny/1.0.3/ui.R"
    ),
    md5sum = c(
      "fe63ffcd28d04001997fe1efb900ac42",
      "df9a4773d7dd4b3954446618e29aa197",
      "8aa79954bfc191c0189d7ac657cb2b61",
      "95a811550e8b8577ead13c3a010c9939"
    ),
    stringsAsFactors = FALSE
  )
  
  dl <- .downloadFiles(filesUrl = filesUrl, pathDownload = mainOutputDir, threads = threads)
  
  dir.create(subOutputDir, showWarnings = FALSE)
  
  # Create a copy of the ArchRProj 
  ArchRProjShiny <- ArchRProj
  # Add metadata to ArchRProjShiny
  ArchRProjShiny@projectMetadata[["groupBy"]] <- groupBy
  ArchRProjShiny@projectMetadata[["tileSize"]] <- tileSize

  #copy the RDS corresponding to the ArchRProject to a new directory for use in the Shiny app
  file.copy(file.path(getOutputDirectory(ArchRProjShiny), savedArchRProjFile), file.path(mainOutputDir), recursive=FALSE)

  # Create fragment files - should be saved within a dir called ShinyFragments within the ArchRProjShiny output directory
  fragDir <- file.path(mainOutputDir, "ShinyFragments", groupBy)
  fragFiles <- list.files(path = file.path(fragDir, pattern = "\\_frags.rds$"))
  dir.create(file.path(mainOutputDir, "ShinyFragments"), showWarnings = FALSE)

  groups <- unique(ArchRProj@cellColData[,groupBy])

  #check for the existence of each expected fragment file and create if not found
  fragOut <- .safelapply(seq_along(groups), function(x){
    groupsx <- groups[x]
    if(!file.exists(file.path(fragDir,paste0(groupsx,"_frags.rds"))) | force){
      .exportGroupFragmentsRDS(ArchRProj = ArchRProjShiny, groupBy = groupBy, outDir = fragDir, threads = threads)
    } else {
      message(paste0("Fragment file for ", groupsx," already exist. Skipping fragment file generation..."))
    }
    return(NULL)
  }, threads = threads)
  
  # Create coverage objects - should be saved within a dir called ShinyCoverage within the mainOutputDir
  covDir <- file.path(mainOutputDir, "ShinyCoverage", groupBy)
  covFiles <- list.files(path = covDir, pattern = "\\_cvg.rds$")
  dir.create(file.path(mainOutputDir, "ShinyCoverage"), showWarnings = FALSE)
  
  covOut <- .safelapply(seq_along(groups), function(x){
    groupsx <- groups[x]
    if(!file.exists(file.path(covDir,paste0(groupsx,"_cvg.rds"))) | force){
      .exportClusterCoverageRDS(ArchRProj = ArchRProjShiny, tileSize = tileSize, groupBy = groupBy, outDir = covDir, fragDir = fragDir, threads = threads)
    } else {
      message(paste0("Coverage file for ", groupsx," already exist. Skipping coverage file generation..."))
    }
    return(NULL)
  }, threads = threads)
  
  #Create embedding plots for columns in cellColData
  message("Generating raster embedding images for cellColData entries...")
  
  # mainEmbeds will create an HDF5 file containing the nativeRaster vectors for data stored in cellColData
  if (!file.exists(file.path(subOutputDir, "mainEmbeds.h5"))) {
    .mainEmbeds(
      ArchRProj = ArchRProjShiny,
      outDirEmbed = file.path(subOutputDir),
      colorBy = "cellColData",
      cellColEmbeddings = cellColEmbeddings,
      embedding = embedding,
      logFile = logFile
    )
  } else{
    message("H5 for main embeddings already exists...")
  }
  
  #Create embedding plots for matrices
  message("Generating raster embedding images for matrix data...")
  
  # matrixEmbeds will create an HDF5 file containing he nativeRaster vectors for data stored in matrices
  if(!file.exists(file.path(subOutputDir, "plotBlank72.h5"))){
    .matrixEmbeds(
      ArchRProj = ArchRProj,
      outDirEmbed = file.path(subOutputDir),
      colorBy = intersect(supportedMatrices, allMatrices),
      embedding = embedding,
      threads = threads,
      verbose = TRUE,
      logFile = logFile
    )
    
  }else{
    message("H5 file already exists...")
  }  
  
  ## delete unnecessary files -----------------------------------------------------------------
  unlink(file.path(projDir, "ShinyFragments"), recursive = TRUE) 
  
  ## ready to launch ---------------------------------------------------------------
  message("App is created!", '\n',
          "Please run the following code chunk to launch the app:",'\n\n',
          
              "ArchRProj <- loadArchRProject('", projDir,"')\n", 
              "mainDir = ", "'", mainDir, "'" ,'\n',
              "subOutDir = ", "'",subOutDir,"'",'\n',
              "savedArchRProjFile = ", "'",savedArchRProjFile,"'",'\n',
              "groupBy = ", "'",groupBy,"'",'\n',
              "cellColEmbeddings = ", "c(",paste(shQuote(cellColEmbeddings, type = "cmd"), collapse=", "),")",'\n',
              "embedding = ", "'",embedding,"'",'\n', 
              "availableMatrices = ", "c(",paste(shQuote(allMatrices, type = "cmd"), collapse=", "),")",'\n',
              "shiny::runApp('", mainOutputDir, "')"
          
          )

}

#' Create an HDF5 file, mainEmbeds.h5, containing the nativeRaster vectors for the 5 main embeddings. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outDirEmbed Where the HDF5 and the jpgs will be saved.
#' @param colorBy `cellColData` ("cellColData") only.
#' @param cellColEmbeddings A character vector of columns in `cellColData` to plot as part of the Shiny app. No default is provided so this must be set.
#' For ex. `c("Sample","Clusters","TSSEnrichment","nFrags")`.
#' @param embedding The embedding to use. Default is "UMAP".
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
.mainEmbeds <- function(
  ArchRProj = NULL,
  outDirEmbed = NULL,
  colorBy = "cellColData", 
  cellColEmbeddings = NULL,
  embedding = "UMAP",
  threads = getArchRThreads(),
  logFile = createLogFile("mainEmbeds")
){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outDirEmbed, name = "outDirEmbed", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = cellColEmbeddings, name = "cellColEmbeddings", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("numeric"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "mainEmbeds Input-Parameters", logFile = logFile)
  
  if(!file.exists(file.path(outDirEmbed, "embeds.rds"))){  
    
    embeds <- .safelapply(1:length(cellColEmbeddings), function(x){ # 
      
      tryCatch({
        named_embed <- plotEmbedding(
          ArchRProj = ArchRProj,
          baseSize = 12,
          colorBy = colorBy,
          name = cellColEmbeddings[x],
          embedding = embedding,
          rastr = FALSE,
          size = 0.5,
        ) + ggtitle(paste0("Colored by ", cellColEmbeddings[x])) + 
        theme(
          text = element_text(size=12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 6)
          )
      }, error = function(x){
        print(x)
      })
      return(named_embed)
    }, threads = threads)
    
    names(embeds) <- cellColEmbeddings
    saveRDS(embeds, file.path(outDirEmbed, "embeds.rds"))
    
  } else {
    message("Main embeddings already exist. Skipping generation and reading in embeds.rds file...")
    embeds <- readRDS(file.path(outDirEmbed, "embeds.rds"))
  }
  
  h5closeAll()
  points <- H5Fcreate(name = file.path(outDirEmbed, "mainEmbeds.h5"))
  
  embed_legend <- list()
  embed_color <- list()
  
  
  for(i in 1:length(embeds)){
    
    embed_plot <- embeds[i]
    
    embed_plot[[1]]$labels$title <- NULL
    embed_plot_blank <- embed_plot[[1]] + theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.title = element_blank()) +
      theme(legend.position = "none") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        title=element_blank()
      )
    
    #save plot without axes etc as a jpg
    ggsave(filename = file.path(outDirEmbed, paste0(names(embeds)[i],"_blank72.jpg")),
           plot = embed_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
    
    #read back in that jpg because we need vector in native format
    blank_jpg72 <- jpeg::readJPEG(source = file.path(outDirEmbed, paste0(names(embeds)[[i]],"_blank72.jpg")), native = TRUE)
    
    # save the native raster vectors 
    h5createDataset(file = points, dataset = names(embeds)[i], dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = as.vector(blank_jpg72), h5loc = points, name= names(embeds)[i])
    
    # save legend and color scale
    embed_legend[[i]] <- levels(embed_plot[[1]]$data$color)
    names(embed_legend)[[i]] <- names(embed_plot)
    
    embed_color[[i]] <- unique(ggplot_build(embed_plot[[1]])$data[[1]][,"colour"])
    names(embed_color)[[i]] <- names(embed_plot)
    
  }
  
  saveRDS(embed_color, file.path(outDirEmbed, "embed_color.rds"))
  saveRDS(embed_legend, file.path(outDirEmbed, "embed_legend_names.rds"))
}

#' Create an HDF5 containing the nativeRaster vectors for all features for all feature matrices. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outDirEmbed Where the HDF5 and the jpgs will be saved.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData") or by
#' a data matrix in the corresponding ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix"). 
#' @param embedding The embedding to use. Default is "UMAP".
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
.matrixEmbeds <- function(
  ArchRProj = NULL,
  outDirEmbed = NULL,
  colorBy = NULL,
  embedding = "UMAP",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("matrixEmbeds")
){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outDirEmbed, name = "outDirEmbed", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("numeric"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "matrixEmbeds Input-Parameters", logFile = logFile)
  
  if (file.exists(file.path(outDirEmbed, "plotBlank72.h5"))){
    file.remove(file.path(outDirEmbed, "plotBlank72.h5"))
  }

  # save the scale
  embeds_min_max_list = list()
  # save the palette
  embeds_pal_list = list()
  
  for(mat in colorBy){
      
    dir.create(paste0(outDirEmbed, "/",mat, "/embeds"), showWarnings = FALSE)
    
    featureNames <- getFeatures(ArchRProj = ArchRProj, useMatrix = mat)
    featureNames <- featureNames[which(!is.na(featureNames))]

    message(paste0("Creating plots for ", mat,"..."))

    if(!is.null(featureNames)){
        
        featurePlots <- plotEmbedding(
          ArchRProj = ArchRProj,
          colorBy = mat,
          name = featureNames,
          quantCut = c(0.01, 0.95),
          imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
          plotAs = "points",
          rastr = TRUE,
          threads = threads
        )
        
        embeds_points <- .safelapply(seq_along(featurePlots), function(x){
          featurePlotx <- featurePlots[x][[1]]
          if(!is.null(featurePlotx)){
            
            featurePlotx_blank <- featurePlotx + theme(axis.title.x = element_blank()) +
              theme(axis.title.y = element_blank()) +
              theme(axis.title = element_blank()) +
              theme(legend.position = "none") +
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                title=element_blank()
              )
            
            #save plot without axes etc as a jpg.
            ggsave(filename = file.path(outDirEmbed, mat, "embeds", paste0(featureNames[x],"_blank72.jpg")),
                   plot = featurePlotx_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
            
            #read back in that jpg because we need vector in native format
            blank_jpg72 <- jpeg::readJPEG(source = file.path(outDirEmbed, mat, "embeds", paste0(featureNames[x],"_blank72.jpg")),
                                          native = TRUE)
            
            g <- ggplot_build(featurePlotx)
            
            res = list(list(plot=as.vector(blank_jpg72), min = round(min(featurePlotx$data$color),1),
                            max = round(max(featurePlotx$data$color),1), pal = unique(g$data[[1]][,"colour"])))
            
            return(res)
          }
        }, threads = threads)
        
        names(embeds_points) <- featureNames
        embeds_points = embeds_points[!unlist(lapply(embeds_points, is.null))]
        
        embeds_min_max <- data.frame(matrix(NA, 2, length(embeds_points)))
        colnames(embeds_min_max) <- names(embeds_points)[which(!unlist(lapply(embeds_points, is.null)))]
        rownames(embeds_min_max) <- c("min","max")
        
        h5closeAll()
        points =  H5Fcreate(name = file.path(outDirEmbed, paste0(mat,"_plotBlank72.h5")))
        H5Gcreate(points, mat)
        
        for(i in 1:length(embeds_points)){
          h5createDataset(file = points, dataset = paste0(mat,"/",featureNames[i]), dims = c(46656,1), storage.mode = "integer")
          h5writeDataset(obj = embeds_points[[i]][[1]]$plot, h5loc = points, name=paste0(mat,"/",featureNames[i]))
          embeds_min_max[1,i] = embeds_points[[i]][[1]]$min
          embeds_min_max[2,i] = embeds_points[[i]][[1]]$max
          
        }
        
        embeds_min_max_list[[mat]] =  embeds_min_max
        embeds_pal_list[[mat]] = embeds_points[[length(embeds_points)]][[1]]$pal

      }else{
      
      stop("Matrix ", mat,"has no features!")
    }  
    
  }

for(i in 1:length(embeds_pal_list)){
 
  cols = embeds_pal_list[[i]]
  rgb <- col2rgb(cols)
  lab <- convertColor(t(rgb), 'sRGB', 'Lab')
  embeds_pal_list[[i]] <- cols[order(lab[, 'L'])]
  
}
  
scale <- embeds_min_max_list
pal <- embeds_pal_list

saveRDS(scale, file.path(outDirEmbed, "scale.rds"))
saveRDS(pal, file.path(outDirEmbed, "pal.rds"))

}


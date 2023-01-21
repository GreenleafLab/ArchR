# Functions for exporting an ArchR-based Shiny app -----------------------------------------------------------
#'
#' Export a Shiny App based on an ArchRProj
#' 
#' Generate all files required for an autonomous Shiny app to display browser tracks and embeddings.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param outputDir The name (not the path!) of the directory for the Shiny App files. This will become a sub-directory of the ArchRProj output directory
#' given by `getOutputDirectory(ArchRProj)`.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for generating BigWig-style sequencing tracks. 
#' Only one cell grouping is allowed.
#' @param cellColEmbeddings A character vector of columns in `cellColData` to plot as part of the Shiny app. No default is provided so this must be set.
#' For ex. `c("Sample","Clusters","TSSEnrichment","nFrags")`.
#' @param embedding The name of the embedding from `ArchRProj` to be used for plotting embeddings in the Shiny app.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param force A boolean value that indicates whether to overwrite any relevant files during the `exportShinyArchR()` process.
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportShinyArchR <- function(
  ArchRProj = NULL,
  outputDir = "Shiny",
  groupBy = "Clusters",
  cellColEmbeddings = NULL,
  embedding = "UMAP",
  tileSize = 100,
  force = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("exportShinyArchR")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outputDir, name = "outputDir", valid = c("character"))
  .validInput(input = outputDir, name = "subOutputDir", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = cellColEmbeddings, name = "groupBy", valid = c("character", "null"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')
  
  if(length(groupBy) > 1){
    stop("Only one value is allowed for groupBy".)
  }
  
  if(is.null(cellColEmbeddings)){
    stop("The cellColEmbeddings parameter must be defined! Please see function input definitions.")
  } else if(!all(cellColEmbeddings %in% colnames(ArchRProj@cellColData)){
    stop("Not all entries in cellColEmbeddings exist in the cellColData of your ArchRProj. Please check provided inputs.")
  }
  
  # Check that the embedding exists in ArchRProj@embeddings
  if(embedding %ni% names(ArchRProj@embeddings)){
    stop("embedding doesn't exist in ArchRProj@embeddings")
  }else{
    print(paste0("embedding:", embedding))
  }
  
  #check that groupBy column exists and doesnt have NA values
  if (groupBy %ni% colnames(ArchRProj@cellColData)) {
    stop("groupBy is not part of cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("Some entries in the column indicated by groupBy have NA values. Please subset your project using subsetArchRProject() to only contain cells with values for groupBy")
  }
  
  subOutputDir <- "inputData" #hardcoded
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
  ArchRProjShiny@projectMetadata[["groupBy"]] <- groupBy
  ArchRProjShiny@projectMetadata[["tileSize"]] <- tileSize
  units <- tryCatch({
   .h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
  },error=function(e){
   "values"
  })
  ArchRProjShiny@projectMetadata[["units"]] <- units
  ArchRProjShiny <- saveArchRProject(ArchRProj = ArchRProjShiny, outputDirectory = 
  file.path(mainDir, outputDir, "Save-ArchRProjShiny"), dropCells = TRUE, overwrite = TRUE, load = TRUE)
  
  projDir <- getOutputDirectory(ArchRProj = ArchRProjShiny)
  
  # Create fragment files - should be saved within a dir called ShinyFragments within the ArchRProjShiny output directory
  fragDir <- file.path(projDir, "ShinyFragments", groupBy)
  fragFiles <- list.files(path = file.path(fragDir, pattern = "\\_frags.rds$"))
  #this is still a slightly dangerous comparison, better would be to compare for explicitly the file names that are expected
  if(length(fragFiles) == length(unique(ArchRProjShiny@cellColData[,groupBy]))){
    if(force){
      .getGroupFragsFromProj(ArchRProj = ArchRProjShiny, groupBy = groupBy, outDir = fragDir)
    } else{
      message("Fragment files already exist. Skipping fragment file generation...")
    }    
  }else{ 
    .getGroupFragsFromProj(ArchRProj = ArchRProjShiny, groupBy = groupBy, outDir = fragDir)
  }
  
  # Create coverage objects - should be saved within a dir called ShinyCoverage within the ArchRProjShiny output directory
  covDir <- file.path(projDir, "ShinyCoverage", groupBy)
  covFiles <- list.files(path = covDir, pattern = "\\_cvg.rds$")
  #this is still a slightly dangerous comparison, better would be to compare for explicitly the file names that are expected
  if(length(covFiles) == length(unique(ArchRProjShiny@cellColData[,groupBy]))){
    if(force){
      .getClusterCoverage(ArchRProj = ArchRProjShiny, tileSize = tileSize, groupBy = groupBy, outDir = covDir)
    } else{
      message("Coverage files already exist. Skipping fragment file generation...")
    }
  }else{
    .getClusterCoverage(ArchRProj = ArchRProjShiny, tileSize = tileSize, groupBy = groupBy, outDir = covDir)
  }

  # Create directory to save input data to Shinyapps.io (everything that will be preprocessed)
  dir.create(file.path(mainDir, outputDir, subOutputDir), showWarnings = TRUE)
 
  allMatrices <- getAvailableMatrices(ArchRProjShiny)
  matrices <- list()
  imputeMatrices <- list()
  imputeWeights <- getImputeWeights(ArchRProj = ArchRProjShiny)
  df <- getEmbedding(ArchRProjShiny, embedding = embedding, returnDF = TRUE)
  
  if(!file.exists(file.path(outputDir, subOutputDir, "matrices.rds")) && !file.exists(file.path(outputDir, subOutputDir, "imputeMatrices.rds"))){
    for(matName in allMatrices){
      featuresNames <- getFeatures(ArchRProj = ArchRProj, useMatrix = matName)
      saveRDS(featuresNames, file.path(outputDir, subOutputDir, matName, "_names.rds"))
    
      if(!is.null(featuresNames)){
        
         mat = Matrix(.getMatrixValues(
          ArchRProj = ArchRProjShiny, 
          name = featuresNames,
          matrixName = matName,
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
   
  # mainEmbeds will create an HDF5 file containing the nativeRaster vectors for data stored in cellColData
  if (!file.exists(file.path(mainDir, outputDir, subOutputDir, "mainEmbeds.h5"))) {
    .mainEmbeds(ArchRProj = ArchRProjShiny,
              outDirEmbed = file.path(mainDir, outputDir, subOutputDir),
              colorBy = "cellColData",
              cellColEmbeddings = cellColEmbeddings,
              embeddingDF = df,
              matrices =  matrices,
              imputeMatrices = imputeMatrices,
              logFile = createLogFile("mainEmbeds")
            )          
  } else{
    message("H5 for main embeddings already exists...")
  }

  # matrixEmbeds will create an HDF5 file containing he nativeRaster vectors for data stored in matrices
  supportedMatrices <- c("GeneScoreMatrix", "GeneIntegrationMatrix", "MotifMatrix") #only these matrices are currently supported for ShinyArchR
  if(!file.exists(file.path(outputDir, subOutputDir, "plotBlank72.h5"))){

    .matrixEmbeds(
      ArchRProj = ArchRProj,
      outDirEmbed = file.path(mainDir, outputDir, subOutputDir),
      colorBy = intersect(supportedMatrices, allMatrices),
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
  unlink(file.path(projDir, "ShinyFragments"), recursive = TRUE) 

  ## ready to launch ---------------------------------------------------------------
  message("App created! To launch, 
            ArchRProj <- loadArchRProject('", projDir,"') and 
            run shiny::runApp('", outputDir, "') from parent directory")
  #  runApp("myappdir")

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
#' @param matrices List of stored matrices to use for plotEmbedding so that it runs faster. 
#' @param imputeMatrices List of stored imputed matrices to use for plotEmbedding so that it runs faster. 
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
.mainEmbeds <- function(
  ArchRProj = NULL,
  outDirEmbed = NULL,
  colorBy = "cellColData", 
  cellColEmbeddings = NULL,
  embedding = "UMAP",
  matrices = NULL,
  imputeMatrices = NULL,
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

     # check all names exist in ArchRProj
      if(cellColEmbeddings %ni% colnames(ArchRProj@cellColData)){
        stop("All columns should be present in cellColData")
      }

    embeds <- .safelapply(1:length(cellColEmbeddings), function(x){
       name <- cellColEmbeddings[[x]]

       tryCatch({
         named_embed <- plotEmbedding(
           ArchRProj = ArchRProj,
           baseSize = 12,
           colorBy = colorBy,
           name = name,
           allNames = names,
           embedding = embedding,
           embeddingDF = df,
           rastr = FALSE,
           size = 0.5,
           # imputeWeights = NULL, # unsure if inputWeights needed for cellColData
           Shiny = TRUE
         )+ggtitle(paste0("Colored by ", name))+theme(text = element_text(size=12), 
                                                      legend.title = element_text(size = 12),legend.text = element_text(size = 6))
       }, error = function(x){
         print(x)
       })
        return(named_embed) 
    })

    names(embeds) <- names
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
    blank_jpg72 <- readJPEG(source = file.path(outDirEmbed, paste0(names(embeds)[[i]],"_blank72.jpg")), native = TRUE)

    # save the native raster vectors 
    h5createDataset(file = points, dataset = names(embeds)[i], dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = as.vector(blank_jpg72), h5loc = points, name= names[i])

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
#' @param matrices List of stored matrices to use for plotEmbedding so that it runs faster. 
#' @param imputeMatrices List of stored imputed matrices to use for plotEmbedding so that it runs faster. 
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
.matrixEmbeds <- function(
  ArchRProj = NULL,
  outDirEmbed = NULL,
  colorBy = c("GeneScoreMatrix", "GeneIntegrationMatrix", "MotifMatrix"),
  embedding = "UMAP",
  matrices = NULL,
  imputeMatrices = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("matrixEmbeds")
){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outDirEmbed, name = "outDirEmbed", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = matrices, name = "matrices", valid = c("list"))
  .validInput(input = imputeMatrices, name = "imputeMatrices", valid = c("list"))
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

  allMatrices <- getAvailableMatrices(ArchRProj)

  for(mat in colorBy){
      if(mat %ni% shinyMatrices){
        stop(mat,"not in ArchRProj")
      }

      if(file.exists(paste0(outDirEmbed, "/", mat, ".rds"))){
        featureNames <- readRDS(file.path(outputDir, subOutputDir, mat, "_names.rds"))
      }else{ 
        

        if(!is.null(featureNames)){

        embeds_points <- .safelapply(1:length(featureNames), function(x){ 

            print(paste0("Creating plots for ", mat,": ",x,": ",round((x/length(featureNames))*100,3), "%"))

          if(!is.na(matrices[[mat]][x])){

              gene_plot <- plotEmbedding(
                ArchRProj = ArchRProj,
                colorBy = mat,
                name = featureNames[x],
                embedding = embedding,
                quantCut = c(0.01, 0.95),
                imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
                plotAs = "points",
                matrices = mat,
                embeddingDF = df,
                imputeMatrices = imputeMatrices,
                rastr = TRUE
              )
          }else{
            gene_plot = NULL
          }

          if(!is.null(gene_plot)){

            gene_plot_blank <- gene_plot + theme(axis.title.x = element_blank()) +
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
                ggsave(filename = file.path(outDirEmbed, paste0(shinymatrices,"_embeds"), paste0(featureNames[x],"_blank72.jpg")),
                       plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)

                #read back in that jpg because we need vector in native format
                blank_jpg72 <- readJPEG(source = file.path(outDirEmbed, paste0(shinymatrices,"_embeds"),
                                                           paste0(featureNames[x],"_blank72.jpg")), native = TRUE)

            g <- ggplot_build(gene_plot)

            res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot$data$color),1),
                            max = round(max(gene_plot$data$color),1), pal = unique(g$data[[1]][,"colour"])))

            return(res)
          }


        }, threads = threads) 

          names(embeds_points) <- featureNames

        embeds_points = embeds_points[!unlist(lapply(embeds_points, is.null))]

        embeds_min_max <- data.frame(matrix(NA, 2, length(embeds_points)))
        colnames(embeds_min_max) <- names(embeds_points)[which(!unlist(lapply(embeds_points, is.null)))]
        rownames(embeds_min_max) <- c("min","max")

          h5closeAll()
          points =  H5Fcreate(name = file.path(outDirEmbed, paste0(shinymatrices,"_plotBlank72.h5")))
          h5createGroup(file.path(outDirEmbed, paste0(shinymatrices,"_plotBlank72.h5")), shinymatrices)

          for(i in 1:length(embeds_points)){

            print(paste0("Getting H5 files for embeds_points: ",i,": ",round((i/length(embeds_points))*100,3), "%"))
            h5createDataset(file = points, dataset = paste0(shinymatrices,"/",featureNames[i]), dims = c(46656,1), storage.mode = "integer")
            h5writeDataset(obj = embeds_points[[i]][[1]]$plot, h5loc = points, name=paste0(shinymatrices,"/",featureNames[i]))
            embeds_min_max[1,i] = embeds_points[[i]][[1]]$min
            embeds_min_max[2,i] = embeds_points[[i]][[1]]$max

          }

          embeds_min_max_list[[shinymatrices]] =  embeds_min_max
          embeds_pal_list[[shinymatrices]] = embeds_points[[length(embeds_points)]][[1]]$pal


        }else{

          message(matName,".rds file is NULL")

        }

        embeds_min_max_list[[matrix]] =  embeds_min_max
        embeds_pal_list[[matrix]] = embeds_points[[length(embeds_points)]][[1]]$pal

      }else{

        message(matName,".rds file does not exist")
      }
    }else{
      stop(matrixName,".rds file does not exist. This file should have been created previously be exportShinyArchR.")
    }


  }

  scale <- embeds_min_max_list
  pal <- embeds_pal_list

  saveRDS(scale, file.path(outDirEmbed, "scale.rds"))
  saveRDS(pal, file.path(outDirEmbed, "pal.rds"))

}

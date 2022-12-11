# mainEmbed function -----------------------------------------------------------
#' 
#' Create an HDF5, mainEmbeds.h5, containing the nativeRaster vectors for the 5 main embeddings. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outDirEmbed Where the HDF5 and the jpgs will be saved.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData") or by
#' a data matrix in the corresponding ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param names A list of the names of the columns in `cellColData` or the featureName/rowname of the data matrix to be used for plotting. 
#' For example if colorBy is "cellColData" then `names` refers to a column names in the `cellcoldata` (see `getCellcoldata()`). If `colorBy`
#' is "GeneScoreMatrix" then `name` refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param embedding The embedding to use. Default is "UMAP".
#' @param Shiny A boolean value that tells the function is calling for Shiny or not.
#' @param matrices A list that contains color matrices for genes.
#' @param imputeMatricesList A list that contains color matrices for genes after imputation.   
#' @param threads The number of threads to use for parallel execution.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
mainEmbed <- function(
  ArchRProj = NULL,
  outDirEmbed = NULL,
  colorBy = "cellColData", 
  names = NULL,
  embedding = "UMAP",
  Shiny = FALSE,
  matrices = matrices,
  imputeMatricesList = imputeMatricesList,
  threads = getArchRThreads(),
  logFile = createLogFile("mainEmbeds")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outDirEmbed, name = "outDirEmbed", valid = c("character"))
  .validInput(input = colorBy, name = "colorBy", valid = c("character"))
  .validInput(input = names, name = "names", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character"))
  .validInput(input = Shiny, name = "Shiny", valid = c("boolean"))
  .validInput(input = matrices, name = "matrices", valid = c("list"))
  .validInput(input = imputeMatricesList, name = "imputeMatricesList", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("numeric"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)

# check to see if the matrix exists using getAvailableMatrices()
# Check if colorBy is cellColData or Matrix (e.g. GSM, GIM, or MM)
# Check if embedding exists in ArchRProj@embeddings
# Check all names exist
  
  if(!file.exists(file.path(outDirEmbed, "embeds.rds"))){  
    
     # check all names exist in ArchRProj
      ccd <- getCellColData(ArchRProj)
      discreteCols <- lapply(seq_len(ncol(ccd)), function(x){
        .isDiscrete(ccd[, x])
      }) %>% unlist %>% {colnames(ccd)[.]}
      if("Clusters" %in% discreteCols){
        selectCols <- "Clusters"
      }else{
        selectCols <- "Sample"
      }

    embeds <- .safelapply(1:length(names), function(x){
       name <- names[[x]]
       print(name)
       
       tryCatch({
         named_embed <- plotEmbedding(
           ArchRProj = ArchRProj,
           baseSize = 12,
           colorBy = colorBy,
           name = name,
           embedding = embedding,
           rastr = FALSE,
           size = 0.5,
           matrices = matrices,
           imputeMatricesList = imputeMatricesList,
           Shiny = ShinyArchR
         )+ggtitle(paste0("Colored by ", name))+theme(text = element_text(size=12), 
                                                      legend.title = element_text(size = 12),legend.text = element_text(size = 6))
       }, error = function(x) {
       })
       
        
        return(named_embed)
    })

    saveRDS(embeds, file.path(outDirEmbed, "embeds.rds"))
 
  } else {
    message("embeddings already exist...")
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
    
    #save plot without axes etc as a jpg.
    ggsave(filename = file.path(outDirEmbed, paste0(names(embeds)[i],"_blank72.jpg")),
           plot = embed_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
    
    #read back in that jpg because we need vector in native format
    blank_jpg72 <- readJPEG(source = file.path(outDirEmbed, paste0(names(embeds)[[i]],"_blank72.jpg")), native = TRUE)
    
    h5createDataset(file = points, dataset = names(embeds)[i], dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = as.vector(blank_jpg72), h5loc = points, name= names[i])
    
    embed_legend[[i]] <- levels(embed_plot[[1]]$data$color)
    names(embed_legend)[[i]] <- names(embed_plot)
    
    
    embed_color[[i]] <- unique(ggplot_build(embed_plot[[1]])$data[[1]][,"colour"])
    names(embed_color)[[i]] <- names(embed_plot)
    
  }
  
  saveRDS(embed_color, file.path(outDirEmbed, "embeddings.rds"))
  saveRDS(embed_legend, file.path(outDirEmbed, "embed_legend_names.rds"))
}

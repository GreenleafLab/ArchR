# mainUmaps function -----------------------------------------------------------
#' 
#' Create an HDF5, mainUMAPs.h5, containing the nativeRaster vectors for the 5 main UMAPS. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outputDirUmaps Where the HDF5 and the jpgs will be saved.  
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
mainUMAPs <- function(
  ArchRProj = NULL,
  outputDirUmaps = "Shiny/inputData",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("mainUMAPs")
){
  
  if(!file.exists(file.path(outputDirUmaps, "umaps.rds"))){  
    umaps <- list()
    
    cluster_umap <- plotEmbedding(
      ArchRProj = ArchRProj,
      baseSize=12,
      colorBy = "cellColData",
      name = "Clusters",
      embedding = "UMAP",
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
    umaps[["Unconstrained"]] <- unconstrained_umap
    
    constrained_remapped_umap <- plotEmbedding(
      ArchRProj = ArchRProjShiny,
      colorBy = "cellColData",
      name = "Clusters2",
      rastr = FALSE,
    )+ggtitle("UMAP: Constrained remapped clusters")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))
    umaps[["Constrained remap"]] <- constrained_remapped_umap
    
    saveRDS(umaps, file.path(outputDirUmaps, "umaps.rds"))
  } else {
    message("umaps already exists...")
    umaps <- readRDS(file.path(outputDirUmaps, "umaps.rds"))
  }
  
  h5closeAll()
  
  points <- H5Fcreate(name = file.path(outputDirUmaps, "mainUMAPs.h5"))
  umap_legend <- list()
  umap_color <- list()
  for(i in 1:length(umaps)){
    
    umap_plot <- umaps[i]
    
    umap_plot[[1]]$labels$title <- NULL
    umap_plot_blank <- umap_plot[[1]] + theme(axis.title.x = element_blank()) +
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
    ggsave(filename = file.path(outputDirUmaps, paste0(names(umaps)[i],"_blank72.jpg")),
           plot = umap_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
    
    #read back in that jpg because we need vector in native format
    blank_jpg72 <- readJPEG(source = file.path(outputDirUmaps, paste0(names(umaps)[[i]],"_blank72.jpg")), native = TRUE)
    
    h5createDataset(file = points, dataset = names(umaps)[i], dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = as.vector(blank_jpg72), h5loc = points, name= names(umaps)[i])
    
    umap_legend[[i]] <- levels(umap_plot[[1]]$data$color)
    names(umap_legend)[[i]] <- names(umap_plot)
    
    
    umap_color[[i]] <- unique(ggplot_build(umap_plot[[1]])$data[[1]][,"colour"])
    names(umap_color)[[i]] <- names(umap_plot)
    
  }
  
  saveRDS(umap_color, file.path(outputDirUmaps, "color_umaps.rds"))
  saveRDS(umap_legend, file.path(outputDirUmaps, "umap_legend_names.rds"))
}

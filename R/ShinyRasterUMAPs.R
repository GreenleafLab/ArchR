# shinyRasterUmaps function -----------------------------------------------------------
#' 
#' 
#' Create an HDF5 containing the nativeRaster vectors for all features for all feature matrices. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outputDirUmaps Where the HDF5 and the jpgs will be saved.  
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
shinyRasterUMAPs <- function(
  ArchRProj = NULL,
  outputDirUmaps = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("ShinyRasterUMAPs")
){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outputDirUmaps, name = "outputDirUmaps", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("numeric"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  
  if (file.exists(file.path(outputDirUmaps, "plotBlank72.h5"))){
    
    file.remove(file.path(outputDirUmaps, "plotBlank72.h5"))
    
  }
  
  umaps_min_max_list = list()
  umaps_pal_list = list()
  
  shinyMatrices <- getAvailableMatrices(ArchRProj)
  
  for(shinymatrices in shinyMatrices){
    
      print(shinymatrices)
      matrixName = paste0(shinymatrices,"_names")
    
      if(file.exists(paste0(outputDirUmaps, "/", matrixName, ".rds"))){
        
        geneMatrixNames <- readRDS(paste0(outputDirUmaps, "/", matrixName, ".rds"))
        
        if(!is.null(geneMatrixNames)){
          
        umaps_points <- .safelapply(1:10, function(x){
            
            print(paste0("Creating plots for ",shinymatrices,": ",x,": ",round((x/length(geneMatrixNames))*100,3), "%"))
            
            gene_plot <- plotEmbedding(
              ArchRProj = ArchRProj,
              colorBy = shinymatrices,
              name = geneMatrixNames[x],
              embedding = "UMAP",
              quantCut = c(0.01, 0.95),
              imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
              plotAs = "points",
              matrices = matrices,
              imputeMatricesList = imputeMatricesList,
              rastr = TRUE
            )
            
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
              ggsave(filename = file.path(outputDirUmaps, paste0(shinymatrices,"_umaps"), paste0(geneMatrixNames[x],"_blank72.jpg")),
                     plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
              
              #read back in that jpg because we need vector in native format
              blank_jpg72 <- readJPEG(source = file.path(outputDirUmaps, paste0(shinymatrices,"_umaps"),
                                                         paste0(geneMatrixNames[x],"_blank72.jpg")), native = TRUE)
              
              g <- ggplot_build(gene_plot)
              
              res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot$data$color),1),
                              max = round(max(gene_plot$data$color),1), pal = unique(g$data[[1]][,"colour"])))
              
              return(res)
            }
          }, threads = threads) 
          names(umaps_points) <- geneMatrixNames[1:10]
          
        }else{
          
          message(matrixName,".rds file is NULL")
          
        }
      
      }else{
        
        message(matrixName,".rds file does not exist")
      }
      
      umaps_points = umaps_points[!unlist(lapply(umaps_points, is.null))]
      
      umaps_min_max <- data.frame(matrix(NA, 2, length(umaps_points)))
      colnames(umaps_min_max) <- names(umaps_points)[which(!unlist(lapply(umaps_points, is.null)))]
      rownames(umaps_min_max) <- c("min","max")
      
      h5closeAll()
      points <- H5Fcreate(name = file.path(outputDirUmaps,"plotBlank72.h5"))
      h5createGroup(file.path(outputDirUmaps, "plotBlank72.h5"), shinymatrices)
      
      for(i in 1:length(umaps_points)){
        
        print(paste0("Getting H5 files for umaps_points: ",i,": ",round((i/length(umaps_points))*100,3), "%"))
        
        h5createDataset(file = points, dataset = paste0(shinymatrices,"/",geneMatrixNames[i]), dims = c(46656,1), storage.mode = "integer")
        h5writeDataset(obj = umaps_points[[i]][[1]]$plot, h5loc = points, name=paste0(shinymatrices,"/",geneMatrixNames[i]))
        
        umaps_min_max[1,i] = umaps_points[[i]][[1]]$min
        umaps_min_max[2,i] = umaps_points[[i]][[1]]$max
        
      }
      
      umaps_min_max_list[[shinymatrices]] =  umaps_min_max
      umaps_pal_list[[shinymatrices]] = umaps_points[[1]][[1]]$pal
    }
  
  scale <- umaps_min_max_list
  pal <- umaps_pal_list
  
  saveRDS(scale, file.path(outputDirUmaps, "scale.rds"))
  saveRDS(pal, file.path(outputDirUmaps, "pal.rds"))
  
  # if(exists("umaps_points")){ rm(umaps_points) }
  # if(exists("GIM_umaps_points")){ rm(GIM_umaps_points) }
  # if(exists("MM_umaps_points")){ rm(MM_umaps_points) }
  
}


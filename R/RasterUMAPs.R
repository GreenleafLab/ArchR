# rasterUmaps function -----------------------------------------------------------
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
#' @export
rasterUMAPs <- function(
  ArchRProj = NULL,
  outputDirUmaps = "Shiny/inputData",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("rasterUMAPs")
){
  
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  
  if (file.exists(file.path(outputDirUmaps, "plotBlank72.h5"))){
    
    file.remove(file.path(outputDirUmaps, "plotBlank72.h5"))
    
  }
  
  h5closeAll()
  points <- H5Fcreate(name = file.path(outputDirUmaps,"plotBlank72.h5"))
  h5createGroup(file.path(outputDirUmaps, "plotBlank72.h5"),"GSM")
  h5createGroup(file.path(outputDirUmaps, "plotBlank72.h5"),"GIM")
  h5createGroup(file.path(outputDirUmaps, "plotBlank72.h5"),"MM")
  

  if(!exists("GSM_umaps_points")){ 
  
    GSM_umaps_points <- ArchR:::.safelapply(1:length(gene_names_GSM), function(x){
      
      print(paste0("Creating plots for GSM_umaps_points: ",x,": ",round((x/length(gene_names_GSM))*100,3), "%"))
      
      gene_plot <- plotEmbeddingShiny(
        ArchRProj = ArchRProj,
        colorBy = "GeneScoreMatrix",
        name = gene_names_GSM[x],
        embedding = "UMAP",
        embeddingDF = df,
        quantCut = c(0.01, 0.95),
        imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
        plotAs = "points",
        rastr = TRUE
      )
      
      if(!is.null(gene_plot)){
        
        gene_plot_blank <- gene_plot[[1]] + theme(axis.title.x = element_blank()) +
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
        ggsave(filename = file.path(outputDirUmaps, "GSM_umaps", paste0(gene_names_GSM[x],"_blank72.jpg")),
               plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
        
        #read back in that jpg because we need vector in native format
        blank_jpg72 <- readJPEG(source = file.path(outputDirUmaps, "GSM_umaps",
                                                   paste0(gene_names_GSM[x],"_blank72.jpg")), native = TRUE)
        
        res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot[[1]]$data$color),1),
                   max = round(max(gene_plot[[1]]$data$color),1), pal = gene_plot[[2]]))
        
        return(res)
      }
    }, threads = threads) 
    names(GSM_umaps_points) <- gene_names_GSM
  }else{
    message("GSM_umaps_points already exists. Skipping the loop...")
  }
  
  if(!exists("GIM_umaps_points")){ 
    GIM_umaps_points <- ArchR:::.safelapply(1:length(gene_names_GIM), function(x){ 
      
      print(paste0("Creating plots for GIM_umaps_points: ",x,": ",round((x/length(gene_names_GIM))*100,3), "%"))
      
      gene_plot <- plotEmbeddingShiny(
        ArchRProj = ArchRProj,
        colorBy = "GeneIntegrationMatrix",
        name = gene_names_GIM[x],
        embedding = "UMAP",
        embeddingDF = df,
        quantCut = c(0.01, 0.95),
        imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
        plotAs = "points",
        rastr = TRUE
      )
      
      if(!is.null(gene_plot)){
        gene_plot_blank <- gene_plot[[1]] + theme(axis.title.x = element_blank()) +
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
        ggsave(filename = file.path(outputDirUmaps, "GIM_umaps", paste0(gene_names_GIM[x],"_blank72.jpg")),
               plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
        
        #read back in that jpg because we need vector in native format
        blank_jpg72 <- readJPEG(source = file.path(outputDirUmaps, "GIM_umaps", 
                                                   paste0(gene_names_GIM[x],"_blank72.jpg")), native = TRUE)
        
        
        res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot[[1]]$data$color),1), 
                   max = round(max(gene_plot[[1]]$data$color),1), pal = gene_plot[[2]]))
        return(res)
      }
      
    }, threads = threads)
    names(GIM_umaps_points) <- gene_names_GIM
  }else{
    message("GIM_umaps_points already exists. Skipping the loop...")
  }
  
  if(!exists("MM_umaps_points")){ 

    MM_umaps_points <- ArchR:::.safelapply(1:length(motif_names), function(x){
      
      print(paste0("Creating plots for MM_umaps_points: ",x,": ",round((x/length(motif_names))*100,3), "%"))
      
      gene_plot <- plotEmbeddingShiny(
        ArchRProj = ArchRProj,
        colorBy = "MotifMatrix",
        name = motif_names[x],
        embedding = "UMAP",
        embeddingDF = df,
        quantCut = c(0.01, 0.95),
        imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
        plotAs = "points",
        rastr = TRUE
      )
      
      if(!is.null(gene_plot)){
        gene_plot_blank <- gene_plot[[1]] + theme(axis.title.x = element_blank()) +
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
        ggsave(filename = file.path(outputDirUmaps, "MM_umaps", paste0(motif_names[x],"_blank72.jpg")),
               plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
        
        
        #read back in that jpg because we need vector in native format
        blank_jpg72 <- readJPEG(source = file.path(outputDirUmaps, "MM_umaps", 
                                                   paste0(motif_names[x],"_blank72.jpg")), native = TRUE)
        
        res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot[[1]]$data$color),1), 
                   max = round(max(gene_plot[[1]]$data$color),1), pal = gene_plot[[2]]))
        
        names(res) = motif_names[x]
        return(res)
      }
    }, threads = threads)
    names(MM_umaps_points) <- motif_names
  }else{
    message("MM_umaps_points already exists. Skipping the loop...")
  }
  
  GSM_umaps_points = GSM_umaps_points[!unlist(lapply(GSM_umaps_points, is.null))]
  GIM_umaps_points = GIM_umaps_points[!unlist(lapply(GIM_umaps_points, is.null))]
  MM_umaps_points = MM_umaps_points[!unlist(lapply(MM_umaps_points, is.null))]
  
  GSM_min_max <- data.frame(matrix(NA, 2, length(GSM_umaps_points)))
  colnames(GSM_min_max) <- names(GSM_umaps_points)[which(!unlist(lapply(GSM_umaps_points, is.null)))]
  rownames(GSM_min_max) <- c("min","max")
  
  GIM_min_max <- data.frame(matrix(NA, 2, length(GIM_umaps_points)))
  colnames(GIM_min_max) <- names(GIM_umaps_points)[which(!unlist(lapply(GIM_umaps_points, is.null)))]
  rownames(GIM_min_max) <- c("min","max")
  
  MM_min_max <- data.frame(matrix(NA, 2, length(MM_umaps_points)))
  colnames(MM_min_max) <- names(MM_umaps_points)[which(!unlist(lapply(MM_umaps_points, is.null)))]
  rownames(MM_min_max) <- c("min","max")
  
  for(i in 1:length(GSM_umaps_points)){
    
    print(paste0("Getting H5 files for GSM_umaps_points: ",i,": ",round((i/length(GSM_umaps_points))*100,3), "%"))
    
    h5createDataset(file = points, dataset = paste0("GSM/", gene_names_GSM[i]), dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = GSM_umaps_points[[i]][[1]]$plot, h5loc = points, name=paste0("GSM/", gene_names_GSM[i]))
    
    GSM_min_max[1,i] = GSM_umaps_points[[i]][[1]]$min
    GSM_min_max[2,i] = GSM_umaps_points[[i]][[1]]$max
    
  }
  
  for(i in 1:length(GIM_umaps_points)){
    
    print(paste0("Getting H5 files for GIM_umaps_points: ",i,": ",round((i/length(GIM_umaps_points))*100,3), "%"))
    
    h5createDataset(file = points, dataset = paste0("GIM/", gene_names_GIM[i]), dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = GIM_umaps_points[[i]][[1]]$plot, h5loc = points, name=paste0("GIM/", gene_names_GIM[i]))
    
    GIM_min_max[1,i] = GIM_umaps_points[[i]][[1]]$min
    GIM_min_max[2,i] = GIM_umaps_points[[i]][[1]]$max
    
  }
  
  for(i in 1:length(MM_umaps_points)){
    
    print(paste0("Getting H5 files for MM_umaps_points: ",i,": ",round((i/length(MM_umaps_points))*100,3), "%"))
    h5createDataset(file = points, dataset = paste0("MM/", motif_names[i]), dims = c(46656,1), storage.mode = "integer")
    h5writeDataset(obj = MM_umaps_points[[i]][[1]]$plot, h5loc = points, name=paste0("MM/", motif_names[i]))
    MM_min_max[1,i] = MM_umaps_points[[i]][[1]]$min
    MM_min_max[2,i] = MM_umaps_points[[i]][[1]]$max
    
  }
  
  scale <- list(gsm = GSM_min_max, gim = GIM_min_max, mm = MM_min_max)
  pal <- list(gsm = GSM_umaps_points[[1]][[1]]$pal, gim = GIM_umaps_points[[1]][[1]]$pal, mm = MM_umaps_points[[1]][[1]]$pal)
  
  saveRDS(scale, file.path(outputDirUmaps, "scale.rds"))
  saveRDS(pal, file.path(outputDirUmaps, "pal.rds"))
  
  if(exists("GSM_umaps_points")){ rm(GSM_umaps_points) }
  if(exists("GIM_umaps_points")){ rm(GIM_umaps_points) }
  if(exists("MM_umaps_points")){ rm(MM_umaps_points) }
  
}


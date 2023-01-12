# matrixEmbeds function -----------------------------------------------------------
#' 
#' 
#' Create an HDF5 containing the nativeRaster vectors for all features for all feature matrices. 
#' This function will be called by exportShinyArchR()
#' 
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param outputDirEmbeds Where the HDF5 and the jpgs will be saved.  
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#'
matrixEmbeds <- function(
  ArchRProj = NULL,
  outputDirEmbeds = NULL,
  embedding = "UMAP",
  matrices = NULL,
  imputeMatricesList = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("matrixEmbeds")
){
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = outputDirEmbeds, name = "outputDirEmbeds", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("numeric"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  
  if (file.exists(file.path(outputDirEmbeds, "plotBlank72.h5"))){
    
    file.remove(file.path(outputDirEmbeds, "plotBlank72.h5"))
    
  }
  
  embeds_min_max_list = list()
  embeds_pal_list = list()
  
  shinyMatrices <- getAvailableMatrices(ArchRProj)
  
  for(matrix in shinyMatrices){
    
      matrixName = paste0(matrix,"_names")
    
      if(file.exists(paste0(outputDirEmbeds, "/", matrixName, ".rds"))){
        
        geneMatrixNames <- readRDS(paste0(outputDirEmbeds, "/", matrixName, ".rds"))
        
        if(!is.null(geneMatrixNames)){
          
        embeds_points <- .safelapply(1:length(geneMatrixNames), function(x){ 
            
            print(paste0("Creating plots for ", matrix,": ",x,": ",round((x/length(geneMatrixNames))*100,3), "%"))
            
          if(!is.na(matrices[[matrix]][x])){
          
              gene_plot <- plotEmbedding(
                ArchRProj = ArchRProj,
                colorBy = matrix,
                name = geneMatrixNames[x],
                embedding = embedding,
                quantCut = c(0.01, 0.95),
                imputeWeights = getImputeWeights(ArchRProj = ArchRProj),
                plotAs = "points",
                matrices = matrices,
                imputeMatricesList = imputeMatricesList,
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
                ggsave(filename = file.path(outputDirEmbeds, paste0(shinymatrices,"_embeds"), paste0(geneMatrixNames[x],"_blank72.jpg")),
                       plot = gene_plot_blank, device = "jpg", width = 3, height = 3, units = "in", dpi = 72)
                
                #read back in that jpg because we need vector in native format
                blank_jpg72 <- readJPEG(source = file.path(outputDirEmbeds, paste0(shinymatrices,"_embeds"),
                                                           paste0(geneMatrixNames[x],"_blank72.jpg")), native = TRUE)
                
                g <- ggplot_build(gene_plot)
                
                res = list(list(plot=as.vector(blank_jpg72), min = round(min(gene_plot$data$color),1),
                                max = round(max(gene_plot$data$color),1), pal = unique(g$data[[1]][,"colour"])))
                
                return(res)
              }
              
            
          }, threads = threads) 
          
          names(embeds_points) <- geneMatrixNames
          
          embeds_points = embeds_points[!unlist(lapply(embeds_points, is.null))]
          
          embeds_min_max <- data.frame(matrix(NA, 2, length(embeds_points)))
          colnames(embeds_min_max) <- names(embeds_points)[which(!unlist(lapply(embeds_points, is.null)))]
          rownames(embeds_min_max) <- c("min","max")
          
          h5closeAll()
          points =  H5Fcreate(name = file.path(outputDirEmbeds, paste0(shinymatrices,"_plotBlank72.h5")))
          h5createGroup(file.path(outputDirEmbeds, paste0(shinymatrices,"_plotBlank72.h5")), shinymatrices)
          
          for(i in 1:length(embeds_points)){
            
            print(paste0("Getting H5 files for embeds_points: ",i,": ",round((i/length(embeds_points))*100,3), "%"))
            h5createDataset(file = points, dataset = paste0(shinymatrices,"/",geneMatrixNames[i]), dims = c(46656,1), storage.mode = "integer")
            h5writeDataset(obj = embeds_points[[i]][[1]]$plot, h5loc = points, name=paste0(shinymatrices,"/",geneMatrixNames[i]))
            embeds_min_max[1,i] = embeds_points[[i]][[1]]$min
            embeds_min_max[2,i] = embeds_points[[i]][[1]]$max
            
          }
          
          embeds_min_max_list[[shinymatrices]] =  embeds_min_max
          embeds_pal_list[[shinymatrices]] = embeds_points[[length(embeds_points)]][[1]]$pal
          
          
        }else{
          
          message(matrixName,".rds file is NULL")
          
        }
        
        
      
      }else{
        
        message(matrixName,".rds file does not exist")
      }
      
    
    }
  
  scale <- embeds_min_max_list
  pal <- embeds_pal_list
  
  saveRDS(scale, file.path(outputDirEmbeds, "scale.rds"))
  saveRDS(pal, file.path(outputDirEmbeds, "pal.rds"))
  
}

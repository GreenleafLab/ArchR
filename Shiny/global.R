# Setting up ----------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(farver)
library(rhdf5)
library(plotfunctions)
library(jpeg)
library(ArchR)
library(htmltools)

#As best I can tell, these are not used explicitly.
#Some are dependencies of ArchR already. Others not.
# library(sparseMatrixStats)
# library(BiocManager)
# library(ComplexHeatmap)
# library(shinybusy)
# library(patchwork)
# library(hexbin)
# library(magick)
# library(shinycssloaders)
# library(ggpubr)
# library(raster)

############# NEW ADDITIONS (start) ###############################

# Calling ArchRProj
ArchRProj <- ArchR::loadArchRProject(path = ".", shiny = TRUE)
ArchRProj <- ArchR::addImputeWeights(ArchRProj = ArchRProj)
mainDir <- 'Shiny'
subOutDir <- 'inputData'
groupBy <- 'Clusters'
cellColEmbeddings <- 'Clusters'
embedding <- 'UMAP'
availableMatrices <- c("GeneScoreMatrix", "MotifMatrix", "PeakMatrix", "TileMatrix")
ShinyArchR <- TRUE
sampleLabels <- 'Clusters'



############# NEW ADDITIONS (end) ###############################

# EMBED Visualization ------------------------------------------------------------

# create a list of dropdown options for EMBED tab
EMBEDs_dropdown = base::colnames(ArchRProj@cellColData)[base::colnames(ArchRProj@cellColData) %in% groupBy]
matrices_dropdown = base::names(base::readRDS(base::file.path(subOutDir, "scale.rds")))

for(i in 1:base::length(matrices_dropdown)){
  
  if(base::file.exists(base::paste0(subOutDir, "/", base::paste0(matrices_dropdown[i], "/", matrices_dropdown[i],"_names"), ".rds"))){
    
    base::assign(base::paste0(matrices_dropdown[i], "_dropdown"), base::readRDS(base::paste0(subOutDir, "/", base::paste0(matrices_dropdown[i], "/", matrices_dropdown[i],"_names"), ".rds")))
    
  }
  
}

embed_legend = base::readRDS(base::paste0(subOutDir, "/embed_legend_names.rds"))
color_embeddings = base::readRDS(base::paste0(subOutDir, "/embed_color.rds"))

# define a function to get the EMBED for a feature/gene
getEMBEDplotWithCol<-function(gene,EMBEDList,scaffoldName,matrixType)
{
  gene_plot=EMBEDList[[gene]]
  
  p_template1=base::readRDS(base::paste0(subOutDir, "/" ,scaffoldName,".rds"))
  
  p_template1$scales$scales <- gene_plot$scale
  
  title=base::paste("EMBED of IterativeLSI colored by\n",matrixType," : ",sep="")
  
  p_template1$labels$title <- base::paste0(title, gene)
  
  return(p_template1)
}


# define a function to get the filename for a gene and then call get EMBED function
getEMBED<-function(gene,fileIndexer,folderName,scaffoldName,matrixType)
{
  # getFilename
  for(file in base::names(fileIndexer))
  {
    if(gene %in% fileIndexer[[file]])
    {
      EMBEDs_data_subset=base::readRDS(base::paste(base::paste0(subOutDir, "/" ,folderName),file,sep="/"))
      
      return(getEMBEDplotWithCol(gene,EMBEDs_data_subset,scaffoldName,matrixType))
    }
  }
}

# PlotBrowser ------------------------------------------------------------------

# create a list of dropdown options for plotbroswer tab
gene_names=base::readRDS(base::paste0(subOutDir, "/GeneScoreMatrix/GeneScoreMatrix_names.rds"))



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
library(ComplexHeatmap)
library(ArchR)


############# NEW ADDITIONS (start) ###############################

# Adjusting ArchR functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
  tryCatch({
    eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
  }, error = function(x) {
  })
}

source("AllClasses.R")
source("ArchRBrowser.R")
source("GgplotUtils.R")


# Calling ArchRProj
ArchRProj=loadArchRProject(path = ".", shiny = TRUE)
ArchRProj <- addImputeWeights(ArchRProj = ArchRProj)
mainDir = 'Shiny'
subOutDir = 'inputData'
groupBy = 'Clusters'
cellColEmbeddings = 'Clusters'
embedding = 'UMAP'
availableMatrices = c("GeneScoreMatrix", "MotifMatrix", "PeakMatrix", "TileMatrix")
ShinyArchR = TRUE
sampleLabels = 'Clusters'



############# NEW ADDITIONS (end) ###############################

# EMBED Visualization ------------------------------------------------------------

# create a list of dropdown options for EMBED tab
EMBEDs_dropdown=colnames(ArchRProj@cellColData)[colnames(ArchRProj@cellColData) %in% groupBy]
matrices_dropdown = names(readRDS(file.path(subOutDir, "scale.rds")))

for(i in 1:length(matrices_dropdown)){
  
  if(file.exists(paste0(subOutDir, "/", paste0(matrices_dropdown[i], "/", matrices_dropdown[i],"_names"), ".rds"))){
    
    assign(paste0(matrices_dropdown[i], "_dropdown"), readRDS(paste0(subOutDir, "/", paste0(matrices_dropdown[i], "/", matrices_dropdown[i],"_names"), ".rds")))
    
  }
  
}

embed_legend = readRDS(paste0(subOutDir, "/embed_legend_names.rds"))
color_embeddings = readRDS(paste0(subOutDir, "/embed_color.rds"))

# define a function to get the EMBED for a feature/gene
getEMBEDplotWithCol<-function(gene,EMBEDList,scaffoldName,matrixType)
{
  gene_plot=EMBEDList[[gene]]
  
  p_template1=readRDS(paste0(subOutDir, "/" ,scaffoldName,".rds"))
  
  p_template1$scales$scales <- gene_plot$scale
  
  title=paste("EMBED of IterativeLSI colored by\n",matrixType," : ",sep="")
  
  p_template1$labels$title <- paste0(title, gene)
  
  return(p_template1)
}


# define a function to get the filename for a gene and then call get EMBED function
getEMBED<-function(gene,fileIndexer,folderName,scaffoldName,matrixType)
{
  # getFilename
  for(file in names(fileIndexer))
  {
    if(gene %in% fileIndexer[[file]])
    {
      EMBEDs_data_subset=readRDS(paste(paste0(subOutDir, "/" ,folderName),file,sep="/"))
      
      return(getEMBEDplotWithCol(gene,EMBEDs_data_subset,scaffoldName,matrixType))
    }
  }
}

# PlotBrowser ------------------------------------------------------------------

# create a list of dropdown options for plotbroswer tab
gene_names=readRDS(paste0(subOutDir, "/GeneScoreMatrix/GeneScoreMatrix_names.rds"))



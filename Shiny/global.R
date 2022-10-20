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
library(ArchR)


# specify desired number of threads 
addArchRThreads(threads = 1) 
# specify genome version. Default hg19 set
addArchRGenome("hg19")
set.seed(1)

# Load all hidden ArchR functions ------------------------------------------------
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
  tryCatch({
    eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
  }, error = function(x) {
  })
}

# UMAP Visualization ------------------------------------------------------------

# create a list of dropdown options for umap tab
Umaps_dropdown=c("Clusters","Sample","Unconstrained","Constrained","Constrained remap")
MM_dropdown=readRDS("./inputData/motif_names.rds")
GSM_dropdown=readRDS("./inputData/gene_names_GSM.rds")
GIM_dropdown=readRDS("./inputData/gene_names_GIM.rds")
umap_legend_names = readRDS("./inputData/umap_legend_names.rds")
color_umaps=readRDS("./inputData/color_umaps.rds")


# define a function to get the umap for a gene
getUMAPplotWithCol<-function(gene,umapList,scaffoldName,matrixType)
{
  gene_plot=umapList[[gene]]
  
  p_template1=readRDS(paste("./inputData/",scaffoldName,".rds",sep=""))
  p_template1$scales$scales <- gene_plot$scale
  
  title=paste("UMAP of IterativeLSI colored by\n",matrixType," : ",sep="")
  
  p_template1$labels$title <- paste0(title, gene)
  
  return(p_template1)
}


# define a function to get the filename for a gene and then call get umap function
getUmap<-function(gene,fileIndexer,folderName,scaffoldName,matrixType)
{
  # getFilename
  for(file in names(fileIndexer))
  {
    if(gene %in% fileIndexer[[file]])
    {
      Umaps_data_subset=readRDS(paste(paste0("./inputData/",folderName),file,sep="/"))
      return(getUMAPplotWithCol(gene,Umaps_data_subset,scaffoldName,matrixType))
    }
  }
}

# PlotBrowser ------------------------------------------------------------------

# create a list of dropdown options for plotbroswer tab
gene_names=readRDS("./inputData/gene_names.rds")
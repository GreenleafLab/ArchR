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
library(AnnotationDbi)
library(BSgenome)
library(Biobase)
library(BiocGenerics)
library(BiocParallel)
library(Biostrings)
library(CNEr)
library(ComplexHeatmap)
library(ArchR)

#' # specify whether you use a local machine or the shiny app
#' ShinyArchR = TRUE
#' 
#' # specify desired number of threads
#' addArchRThreads(threads = 1)
#' # specify genome version. Default hg19 set
#' addArchRGenome("hg19")
#' set.seed(1)
#' 
#' ArchRProj=loadArchRProject(path = "Save-ArchRProjShiny/")
#' ArchRProj <- addImputeWeights(ArchRProj = ArchRProj)
#' 
#' ############################################################
#' 
#' # myLoadArchRProject -----------------------------------
#' #' Load Previous ArchRProject into R
#' #'
#' #' This function will load a previously saved ArchRProject and re-normalize paths for usage.
#' #'
#' #' @param path A character path to an `ArchRProject` directory that was previously saved using `saveArchRProject()`.
#' #' @param force A boolean value indicating whether missing optional `ArchRProject` components (i.e. peak annotations /
#' #' background peaks) should be ignored when re-normalizing file paths. If set to `FALSE` loading of the `ArchRProject`
#' #' will fail unless all components can be found.
#' #' @param showLogo A boolean value indicating whether to show the ascii ArchR logo after successful creation of an `ArchRProject`.
#' #' @export
#' myLoadArchRProject <- function(path = "./",
#'                                force = FALSE,
#'                                showLogo = TRUE) {
#'   .validInput(input = path,
#'                       name = "path",
#'                       valid = "character")
#'   .validInput(input = force,
#'                       name = "force",
#'                       valid = "boolean")
#'   .validInput(input = showLogo,
#'                       name = "showLogo",
#'                       valid = "boolean")
#' 
#'   path2Proj <- file.path(path, "Save-ArchR-Project.rds")
#' 
#'   if (!file.exists(path2Proj)) {
#'     stop("Could not find previously saved ArchRProject in the path specified!")
#'   }
#' 
#'   ArchRProj <- recoverArchRProject(readRDS(path2Proj))
#'   outputDir <- getOutputDirectory(ArchRProj)
#'   outputDirNew <- normalizePath(path)
#' 
#' 
#'   ArchRProj@projectMetadata$outputDirectory <- outputDirNew
#' 
#'   message("Successfully loaded ArchRProject!")
#'   if (showLogo) {
#'     .ArchRLogo(ascii = "Logo")
#'   }
#' 
#'   ArchRProj
#' 
#' }
#' 
#' 
#' ## Create fragment files -----------------------------------------------------------
#' .getGroupFragsFromProj <- function(ArchRProj = NULL,
#'                                    groupBy = NULL,
#'                                    outDir = file.path("Shiny", "fragments")) {
#'   dir.create(outDir, showWarnings = FALSE)
#' 
#'   # find barcodes of cells in that groupBy.
#'   groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
#'   cells <- ArchRProj$cellNames
#'   cellGroups <- split(cells, groups)
#' 
#'   # outputs unique cell groups/clusters.
#'   clusters <- names(cellGroups)
#' 
#' 
#'   for (cluster in clusters) {
#'     cat("Making fragment file for cluster:", cluster, "\n")
#'     # get GRanges with all fragments for that cluster
#'     cellNames = cellGroups[[cluster]]
#'     fragments <-
#'       getFragmentsFromProject(ArchRProj = ArchRProj, cellNames = cellNames)
#'     fragments <- unlist(fragments, use.names = FALSE)
#'     # filter Fragments
#'     fragments <-
#'       GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
#'     saveRDS(fragments, file.path(outDir, paste0(cluster, "_cvg.rds")))
#'   }
#' }
#' 
#' 
#' .getClusterCoverage <- function(ArchRProj = NULL,
#'                                 tileSize = 100,
#'                                 scaleFactor = 1,
#'                                 groupBy = "Clusters",
#'                                 outDir = file.path("Shiny", "coverage")) {
#'   fragfiles = list.files(path = file.path("Shiny", "fragments"),
#'                          full.names = TRUE)
#'   dir.create(outDir, showWarnings = FALSE)
#' 
#'   # find barcodes of cells in that groupBy.
#'   groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
#'   cells <- ArchRProj$cellNames
#'   cellGroups <- split(cells, groups)
#' 
#'   # outputs unique cell groups/clusters.
#'   clusters <- names(cellGroups)
#' 
#'   chrRegions <- getChromSizes(ArchRProj)
#'   genome <- getGenome(ArchRProj)
#' 
#'   for (file in fragfiles) {
#'     fragments <- readRDS(file)
#'     #fragmentsToInsertions()
#'     left <- GRanges(seqnames = seqnames(fragments),
#'                     ranges = IRanges(start(fragments), width = 1))
#'     right <- GRanges(seqnames = seqnames(fragments),
#'                      ranges = IRanges(end(fragments), width = 1))
#'     # call sort() after sortSeqlevels() to sort also the ranges in addition
#'     # to the chromosomes.
#'     insertions <- c(left, right) %>% sortSeqlevels() %>%
#'       sort()
#' 
#'     cluster <- file %>% basename() %>% gsub("_.*", "", .)
#'     #binnedCoverage
#'     # message("Creating bins for cluster ",clusters[clusteridx], "...")
#'     bins <-
#'       unlist(slidingWindows(chrRegions, width = tileSize, step = tileSize))
#'     # message("Counting overlaps for cluster ",clusters[clusteridx], "...")
#'     bins$reads <-
#'       countOverlaps(
#'         bins,
#'         insertions,
#'         maxgap = -1L,
#'         minoverlap = 0L,
#'         type = "any"
#'       )
#'     addSeqLengths(bins, genome)
#'     # message("Creating binned coverage for cluster ",clusters[clusteridx], "...")
#'     #each value is multiplied by that weight.
#'     # TODO add scaleFactor
#'     # allCells as.vector(ArchRProj@cellColData$Sample, mode="any")
#'     clusterReadsInTSS <-
#'       ArchRProj@cellColData$ReadsInTSS[cells %in% cellGroups$cluster]
#'     # scaleFactor <- 5e+06 / sum(clusterReadsInTSS)
#'     binnedCoverage <-
#'       coverage(bins, weight = bins$reads * scaleFactor)
#'     saveRDS(binnedCoverage, file.path(outDir, paste0(cluster, "_cvg.rds")))
#'   }
#' 
#' }
#' 
#' 
#' #############################################################
#' 
#' ArchRProj=loadArchRProject("~/Documents/upwork/Paulina Paiz/Shiny_28_11_2022/Save-ProjHeme5/")
#' 
#' 
#' # Load all hidden ArchR functions ------------------------------------------------
#' fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
#' for (i in seq_along(fn)) {
#'   tryCatch({
#'     eval(parse(text = paste0(fn[i], "<-", fn[i])))
#'   }, error = function(x) {
#'   })
#' }

# EMBED Visualization ------------------------------------------------------------

# create a list of dropdown options for EMBED tab
EMBEDs_dropdown=colnames(ArchRProj@cellColData)[colnames(ArchRProj@cellColData) %in% groupBy]
matrices_dropdown = names(readRDS(paste0("./", subOutputDir, "/scale.rds")))

for(i in 1:length(matrices_dropdown)){
  
  if(file.exists(paste0(outputDir, "/", subOutputDir, "/", paste0(matrices_dropdown[i],"_names"), ".rds"))){
    
    assign(paste0(matrices_dropdown[i], "_dropdown"), readRDS(paste0(outputDir, "/", subOutputDir, "/", paste0(matrices_dropdown[i],"_names"), ".rds")))
    
  }
  
}

# if("MotifMatrix" %in% matrices_dropdown){
#   Feature_dropdown = readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/motif_names.rds"))
# }
# 
# if("GeneScoreMatrix" %in% matrices_dropdown){
#   GSM_dropdown = readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/gene_names_GSM.rds"))
# }
# 
# if("GeneIntegrationMatrix" %in% matrices_dropdown){
#   GIM_dropdown = readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/gene_names_GIM.rds"))
# }
embed_legend = readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/embed_legend_names.rds"))
color_embeddings = readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/embeddings.rds"))


# define a function to get the EMBED for a gene
getEMBEDplotWithCol<-function(gene,EMBEDList,scaffoldName,matrixType)
{
  gene_plot=EMBEDList[[gene]]
  
  p_template1=readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/" ,scaffoldName,".rds"))
  
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
      EMBEDs_data_subset=readRDS(paste(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/" ,folderName),file,sep="/"))
      
      return(getEMBEDplotWithCol(gene,EMBEDs_data_subset,scaffoldName,matrixType))
    }
  }
}

# PlotBrowser ------------------------------------------------------------------

# create a list of dropdown options for plotbroswer tab
gene_names=readRDS(paste0(getOutputDirectory(ArchRProj),"/",outputDir, "/", subOutputDir, "/features.rds"))



#' Export a Shiny App based on ArchRProj
#' 
#' Generate all files required for an autonomous shiny app
#' This function will open an interactive shiny session in style of a browser track. It allows for normalization of the signal which
#' enables direct comparison across samples. Note that the genes displayed in this browser are derived from your `geneAnnotation`
#' (i.e. the `BSgenome` object you used) so they may not match other online genome browsers that use different gene annotations.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument
#' can be used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param browserTheme A `shinytheme` from shinythemes for viewing the ArchR Browser. If not installed this will be NULL.
#' To install try devtools::install_github("rstudio/shinythemes").
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportShinyArchR <- function(
  ArchRProj = NULL,
  outputDir = "Shiny",
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  # minCells = 25,
  # baseSize = 10,
  borderWidth = 0.5,
  tickWidth = 0.5,
  facetbaseSize = 12,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  browserTheme = "cosmo",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("exportShinyArchR")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  # .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  # .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  # .validInput(input = minCells, name = "minCells", valid = c("integer"))
  # .validInput(input = baseSize, name = "baseSize", valid = c("integer"))
  # .validInput(input = borderWidth, name = "borderWidth", valid = c("numeric"))
  # .validInput(input = tickWidth, name = "tickWidth", valid = c("numeric"))
  # .validInput(input = facetbaseSize, name = "facetbaseSize", valid = c("numeric"))
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = browserTheme, name = "browserTheme", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "ArchRBrowser Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')

# Make directory for Shiny App 
  if(!dir.exists(outputDir)) dir.create(outputDir)
  if(length(dir(outputDir,  all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
    stop("Please specify a new or empty directory")
  }

# Create fragment files 
.getGroupFragsFromProj(ArchRProj = ArchRProj, groupBy = groupBy)

# Create coverage objects
.getClusterCoverage(ArchRProj = ArchRProj, tileSize = tileSize, groupBy = groupBy)

# Create a copy of the ArchRProj object 
# Add metadata to ArchRProj   
  cells = ArchRProj$cellNames
  data = replicate(length(cells), tileSize)
  ArchRProj <- addCellColData(ArchRProj = ArchRProj, data = data, 
                              cells = cells, name = "tileSize")
  
#' Generate code files required for shiny app (one dataset)
#'
#' Generate code files required for shiny app containing only one dataset. In 
#' particular, two R scripts will be generated, namely \code{server.R} and 
#' \code{ui.R}.  

  ## Clone Github repo with ui.R and server.R ----------------------------------------

  ## ready to launch ---------------------------------------------------------------
  message("App created! To launch, run shiny::runApp('", outputDir, "')")
#  runApp("myappdir")
  
# Bulk Tracks Methods -----------------------------------------------------------

## Create fragment files -----------------------------------------------------------

keepFilteredChromosomes <- function (x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode = "coarse") {
  if (standard) {
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  seqNames <- seqlevels(x)
  chrRemove <- c()
  if (underscore) {
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  chrRemove <- c(chrRemove, which(seqNames %in% remove))
  if (length(chrRemove) > 0) {
    chrKeep <- seqNames[-chrRemove]
  }
  else {
    chrKeep <- seqNames
  }
  seqlevels(x, pruning.mode = pruning.mode) <- chrKeep
  return(x)
}

.getGroupFragsFromProj <- function(
  ArchRProj = NULL,
  groupBy = NULL,
  outDir = "Shiny/fragments"
){
  
  dir.create(outDir, showWarnings = FALSE) 
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters. 
  clusters <- names(cellGroups)
  
  for (cluster in clusters){
    cat("Making fragment file for cluster:", cluster,"\n")
    # get GRanges with all fragments for that cluster 
    cellNames = cellGroups[[cluster]]
    fragments <- getFragmentsFromProject(ArchRProj = ArchRProj, cellNames = cellNames)
    fragments <- unlist(fragments, use.names = FALSE) 
    # filter Fragments
    ATACFragments <- keepFilteredChromosomes(fragments)
    saveRDS(ATACFragments, paste0(outdir,"/fragments/",cluster, "_fragments.rds"))
  }
}

## Create coverage objects -----------------------------------------------------------

addSeqLengths <- function (gr, genome){
  gr <- validGRanges(gr)
  genome <- validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <- seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.getClusterCoverage <- function(
  ArchRProj = NULL,
  tileSize = 100, 
  groupBy = "Clusters",
  outDir = "Shiny/coverage"
  # geneAnnotation = getGeneAnnotation(ArchRProj),
){
  fragfiles = list.files(path = "./Shiny/fragments", full.names = TRUE)
  dir.create(outDir, showWarnings = FALSE) 
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters. 
  clusters <- names(cellGroups)
  
  chrRegions <- getChromSizes(ArchRProj)
  genome <- getGenome(ArchRProj)
  
  clusteridx=0
  for(file in fragfiles){
    clusteridx = clusteridx + 1
    ATACFragments <- readRDS(file) 
    
    #fragmentsToInsertions()
    left <- GRanges(seqnames = seqnames(ATACFragments), 
                    ranges = IRanges(start(ATACFragments), width = 1))
    right <- GRanges(seqnames = seqnames(ATACFragments), 
                     ranges = IRanges(end(ATACFragments), width = 1))
    # call sort() after sortSeqlevels() to sort also the ranges in addition 
    # to the chromosomes. 
    insertions <- c(left, right) %>% sortSeqlevels() %>% 
      sort()
    
    #binnedCoverage
    message("creating bins for cluster ",clusters[clusteridx], "...")
    bins <- unlist(slidingWindows(chrRegions, width = tileSize, step = tileSize))
    message("counting overlaps for cluster ",clusters[clusteridx], "...")
    bins$reads <- countOverlaps(bins, insertions, maxgap = -1L, minoverlap = 0L, type = "any")
    addSeqLengths(bins, genome)
    message("creating binned coverage for cluster ",clusters[clusteridx], "...")
    #each value is multiplied by that weight.
    binnedCoverage <- coverage(bins, weight = bins$reads)
    saveRDS(binnedCoverage, paste0(outdir,"/coverage/",clusters[clusteridx], "_cvg.rds"))
  }
}

  
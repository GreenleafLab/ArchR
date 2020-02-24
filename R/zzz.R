########################################################
# On Load / Attach
########################################################
ArchR_Defaults <- list(
  ArchR.threads = 1,
  ArchR.genome = NA
)

ArchR_Depends <- list(
  c(
    "ggplot2",
    "SummarizedExperiment",
    "data.table",
    "Matrix",
    "rhdf5",
    "magrittr",
    "S4Vectors",
    "BiocGenerics",
    "GenomicRanges"
  )
)

.onLoad <- function(libname, pkgname){
  packageStartupMessage("Loading Packages...")
  for(i in seq_along(ArchR_Depends)){
    suppressPackageStartupMessages(requireNamespace(ArchR_Depends[i], quietly = TRUE))
  }
  invisible()
}

.onAttach <- function(libname, pkgname){
  if(!interactive()) return()
  v <- packageVersion("ArchR")
  ArchR:::.ArchRLogo()
  packageStartupMessage("ArchR : Version ", v, "\nFor more information see our website : https://greenleaflab.github.io/ArchR_Website/\nIf you encounter a bug please report : https://github.com/GreenleafLab/ArchR/issues")
  op <- options()
  toset <- !(names(ArchR_Defaults) %in% names(op))
  if (any(toset)) options(ArchR_Defaults[toset])
  if(!.isWholenumber(options()[["ArchR.threads"]])){
    addArchRThreads()
  }else if(options()[["ArchR.threads"]] == 1){
    addArchRThreads()
  }
  invisible()
}

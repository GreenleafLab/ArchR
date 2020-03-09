########################################################
# On Load / Attach
########################################################
ArchR_Defaults <- list(
  ArchR.threads = 1,
  ArchR.genome = NA
)

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

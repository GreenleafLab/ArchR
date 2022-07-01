#' @include GgplotUtils.R

########################################################
# On Load / Attach
########################################################
ArchRDefaults <- list(
  ArchR.threads = 1,
  ArchR.logging = TRUE,
  ArchR.genome = NA,
  ArchR.chrPrefix = TRUE,
  ArchR.debugging = FALSE,
  ArchR.verbose = TRUE
)

ArchRDependency <- c(
  "grid",
  "gridExtra",
  "gtools",
  "gtable",
  "ggplot2",
  "magrittr",
  "plyr",
  "stringr",
  "data.table",
  "matrixStats",
  "S4Vectors",
  "GenomicRanges",
  "BiocGenerics",
  "Matrix",
  "Rcpp",
  "SummarizedExperiment",
  "rhdf5"
)

.onAttach <- function(libname, pkgname){
  
  #Logo
  .ArchRLogo()
  
  #Package Startup
  v <- packageVersion("ArchR")
  packageStartupMessage("ArchR : Version ", v, "\nFor more information see our website : www.ArchRProject.com\nIf you encounter a bug please report : https://github.com/GreenleafLab/ArchR/issues")
  
  #Load Packages
  packageStartupMessage("Loading Required Packages...")
  pkgs <- ArchRDependency
  for(i in seq_along(pkgs)){
    packageStartupMessage("\tLoading Package : ", pkgs[i], " v", packageVersion(pkgs[i]))
    tryCatch({
      suppressPackageStartupMessages(require(pkgs[i], character.only=TRUE))
    }, error = function(e){
      packageStartupMessage("\tFailed To Load Package : ", pkgs[i], " v", packageVersion(pkgs[i]))
    })
  }

  if(!interactive()) return()

  #Set Default Options
  op <- options()
  toset <- !(names(ArchRDefaults) %in% names(op))
  
  if (any(toset)) options(ArchRDefaults[toset])
  
  if(!.isWholenumber(options()[["ArchR.threads"]])){
    addArchRThreads()
  }else if(options()[["ArchR.threads"]] == 1){
    addArchRThreads()
  }
  
  if(!.checkCairo()){
    packageStartupMessage("WARNING : Cairo check shows Cairo is not functional.\n          ggplot2 rasterization will not be available without Cario.\n          This may cause issues editing plots with many thousands of points from single cells.")
  }
  
  if(.checkJupyter()){
    packageStartupMessage("Detected Jupyer Notebook session. Disabling Log Messages!\n\tIf this is undesired use `addArchRVerbose(TRUE)`")
    addArchRVerbose(verbose = FALSE)
  }
  
  invisible()

}

#Check Jupyer Status
.checkJupyter <- function(){
  tryCatch({
    sysID <- Sys.getenv("JPY_PARENT_PID")
    if(!is.character(sysID)){
      return(FALSE)
    }
    if(sysID == ""){
      FALSE
    }else{
      TRUE
    }
  },error= function(e){
    FALSE
  })
}

#' Install extra packages used in ArchR that are not installed by default
#' 
#' This function will install extra packages used in ArchR that are not installed by default.
#' 
#' @param force If you want to force a reinstall of these pacakges.
#' @export
installExtraPackages <- function(force = FALSE){

  n <- 0
  f <- 0
  fp <- c()

  # Seurat
  if(!requireNamespace("Seurat", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing Seurat..")
      install.packages('Seurat')
      if(!requireNamespace("Seurat", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "Seurat")
      message("Error With Seurat Installation")
      message("Try install.packages('Seurat')")
    })
  }

  # Cairo
  if(!requireNamespace("Cairo", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing Cairo..")
      install.packages('Cairo')
      if(!requireNamespace("Cairo", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "Cairo")
      message("Error With Cairo Installation")
      message("Try install.packages('Cairo')")
    })
  }

  # ggrastr plots
  # if(!requireNamespace("ggrastr", quietly = TRUE) || force){
  #   o <- tryCatch({
  #     message("Installing ggrastr..")
  #     devtools::install_github('VPetukhov/ggrastr')
  #     if(!requireNamespace("ggrastr", quietly = TRUE)){
  #       stop()
  #     }
  #     n <- n + 1
  #   }, error = function(e){
  #     message("Error With ggrastr Installation")
  #     message("Try devtools::install_github('VPetukhov/ggrastr')")
  #   })
  # }

  # harmony is a package that can correct batch effects
  if(!requireNamespace("harmony", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing harmony..")
      devtools::install_github('immunogenomics/harmony', repos = BiocManager::repositories())
      if(!requireNamespace("harmony", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "harmony")
      message("Error With harmony Installation")
      message("Try devtools::install_github('immunogenomics/harmony', repos = BiocManager::repositories())")
    })
  }

  # presto is a package that has efficient tools for wilcoxon tests on sparseMatrices
  if(!requireNamespace("presto", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing presto..")
      devtools::install_github('immunogenomics/presto', repos = BiocManager::repositories())
      if(!requireNamespace("presto", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "presto")
      message("Error With presto Installation")
      message("Try devtools::install_github('immunogenomics/presto', repos = BiocManager::repositories())")
    })
  }

  # shiny
  if(!requireNamespace("shiny", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing shiny..")
      install.packages('shiny')
      if(!requireNamespace("shiny", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "shiny")
      message("Error With shiny Installation")
      message("Try install.packages('shiny')")
    })
  }

  # rhandsontable
  if(!requireNamespace("rhandsontable", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing rhandsontable..")
      install.packages('rhandsontable')
      if(!requireNamespace("rhandsontable", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "rhandsontable")
      message("Error With rhandsontable Installation")
      message("Try install.packages('rhandsontable')")
    })
  }

  # shinythemes
  if(!requireNamespace("shinythemes", quietly = TRUE) || force){
    o <- tryCatch({
      message("Installing shinythemes..")
      install.packages('shinythemes')
      if(!requireNamespace("shinythemes", quietly = TRUE)){
        stop()
      }
      n <- n + 1
    }, error = function(e){
      f <- f + 1
      fp <- c(fp, "rhandsontable")
      message("Error With shinythemes Installation")
      message("Try install.packages('shinythemes')")
    })
  }

  message("Successfully installed ", n, " packages. ", f, " packages failed installation.")

  if(f > 0){
    message("Failed packages : ", paste0(fp, collapse=", "))
  }

  return(invisible(0))

}

##########################################################################################
# Chr Prefix
##########################################################################################

#' Add a globally-applied requirement for chromosome prefix
#' 
#' This function will set the default requirement of chromosomes to have a "chr" prefix.
#' 
#' @param chrPrefix A boolean describing the requirement of chromosomes to have a "chr" prefix.
#' @export
addArchRChrPrefix <- function(chrPrefix = TRUE){
  
  .validInput(input = chrPrefix, name = "chrPrefix", valid = "boolean")

  if(chrPrefix){
    message("ArchR is now requiring chromosome prefix = 'chr'")
  }else{
    message("ArchR is now disabling the requirement of chromosome prefix = 'chr'")
  }

  options(ArchR.chrPrefix = chrPrefix)

}

#' Get a globally-applied requirement for chromosome prefix
#' 
#' This function will get the default requirement of chromosomes to have a "chr" prefix.
#' 
#' @export
getArchRChrPrefix <- function(){
  
  .ArchRChrPrefix <- options()[["ArchR.chrPrefix"]]
  
  if(!is.null(.ArchRChrPrefix)){
    if(is.logical(.ArchRChrPrefix)){
      .ArchRChrPrefix
    }else{
      TRUE
    }
  }else{
    TRUE
  }

}

##########################################################################################
# Parallel Information
##########################################################################################

#' Add a globally-applied number of threads to use for parallel computing.
#' 
#' This function will set the number of threads to be used for parallel computing across all ArchR functions.
#' 
#' @param threads The default number of threads to be used for parallel execution across all ArchR functions.
#' This value is stored as a global environment variable, not part of the `ArchRProject`.
#' This can be overwritten on a per-function basis using the given function's `threads` parameter.
#' @param force If you request more than the total number of CPUs minus 2, ArchR will set `threads` to `(nCPU - 2)`.
#' To bypass this, setting `force = TRUE` will use the number provided to `threads`.
#' @export
addArchRThreads <- function(threads = floor(parallel::detectCores()/ 2), force = FALSE){
  
  .validInput(input = threads, name = "threads", valid = "integer")
  .validInput(input = force, name = "force", valid = "boolean")

  if(tolower(.Platform$OS.type) == "windows"){
    message("Detected windows OS, setting threads to 1.")
    threads <- 1
  }

  if(threads >= parallel::detectCores() - 1){
    if(force){
      message("Input threads is equal to or greater than ncores minus 1 (",parallel::detectCores()-1,")\nOverriding since force = TRUE.")
    }else{
      message("Input threads is equal to or greater than ncores minus 1 (",parallel::detectCores()-1,")\nSetting cores to ncores minus 2. Set force = TRUE to set above this number!")
      threads <- parallel::detectCores()-2
    }
  }
  if(threads > 1){
    RNGkind("L'Ecuyer-CMRG")
  }
  
  message("Setting default number of Parallel threads to ", threads, ".")
  options(ArchR.threads = as.integer(round(threads)))

}

#' Get globally-applied number of threads to use for parallel computing.
#' 
#' This function will get the number of threads to be used for parallel execution across all ArchR functions.
#' 
#' @export
getArchRThreads <- function(){
  .ArchRThreads <- options()[["ArchR.threads"]]
  if(!is.null(.ArchRThreads)){
    if(!.isWholenumber(.ArchRThreads)){
      message("option(ArchR.threads) : ", .ArchRThreads, " is not an integer. \nDid you mistakenly set this to a value without addArchRThreads? Reseting to default!")
      addArchRThreads()
      options()[["ArchR.threads"]]
    }else{
      .ArchRThreads
    }
  }else{
    1
  }
}

##########################################################################################
# Create Gene/Genome Annotation
##########################################################################################

#' Add a globally defined genome to all ArchR functions.
#' 
#' This function will set the genome across all ArchR functions.
#' 
#' @param genome A string indicating the default genome to be used for all ArchR functions.
#' Currently supported values include "hg19","hg38","mm9", and "mm10".
#' This value is stored as a global environment variable, not part of the `ArchRProject`.
#' This can be overwritten on a per-function basis using the given function's `geneAnnotation`and `genomeAnnotation` parameter.
#' For something other than one of the currently supported, see `createGeneAnnnotation()` and `createGenomeAnnnotation()`.
#' @param install  A boolean value indicating whether the `BSgenome` object associated with the provided `genome` should be
#' automatically installed if it is not currently installed. This is useful for helping reduce user download requirements.
#' @export
addArchRGenome <- function(genome = NULL, install = TRUE){
  
  .validInput(input = genome, name = "genome", valid = "character")
  .validInput(input = install, name = "install", valid = c("boolean"))

  supportedGenomes <- c("hg19","hg38","mm9","mm10","hg19test")
  
  if(tolower(genome) %ni% supportedGenomes){
    
    message("Genome : ", genome, " is not currently supported by ArchR.")
    message("Currently supported genomes : ", paste0(supportedGenomes, collapse = ","))
    message("To continue try building a custom geneAnnotation with createGeneAnnnotation(),\nand genomeAnnotation with createGenomeAnnotation()!")
 
  }else{

    #Check if BSgenome exists!
    if(tolower(genome)=="hg19"){
      if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
        if(install){
          message("BSgenome for hg19 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        }else{
          stop("BSgenome for hg19 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
        }
      }
    }else if(tolower(genome)=="hg19test"){
      if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
        if(install){
          message("BSgenome for hg19 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        }else{
          stop("BSgenome for hg19 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
        }
      }
    }else if(tolower(genome)=="hg38"){
      if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
        if(install){
          message("BSgenome for hg38 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
        }else{
          stop("BSgenome for hg38 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
        }
      }
    }else if(tolower(genome)=="mm9"){
      if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE)){
        if(install){
          message("BSgenome for mm9 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
          BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
        }else{
          stop("BSgenome for mm9 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
        }
      }
    }else if(tolower(genome)=="mm10"){
      if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
        if(install){
          message("BSgenome for mm10 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
          BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
        }else{
          stop("BSgenome for mm10 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
        }
      }
    }

    genome <- paste(toupper(substr(genome, 1, 1)), substr(genome, 2, nchar(genome)), sep="")
    
    message("Setting default genome to ", genome, ".")
    #assign(".ArchRGenome", genome, envir = .GlobalEnv)
    options(ArchR.genome = genome)

  }

  invisible(0)

}

#' Get the globally defined genome, the geneAnnotation, or genomeAnnotation objects associated with the globally defined genome.
#' 
#' This function will retrieve the genome that is currently in use by ArchR. Alternatively, this function can return either the `geneAnnotation`
#' or the `genomeAnnotation` associated with the globally defined genome if desired.
#' 
#' @param geneAnnotation A boolean value indicating whether the `geneAnnotation` associated with the ArchRGenome should be returned
#' instead of the globally defined genome. The `geneAnnotation` is used downstream to calculate things like TSS Enrichment Scores.
#'This function is not meant to be run with both `geneAnnotation` and `genomeAnnotation` set to `TRUE` (it is an either/or return value).
#' @param genomeAnnotation A boolean value indicating whether the `genomeAnnotation` associated with the ArchRGenome should be returned
#' instead of the globally defined genome. The `genomeAnnotation` is used downstream to determine things like chromosome sizes and nucleotide content.
#' This function is not meant to be run with both `geneAnnotation` and `genomeAnnotation` set to `TRUE` (it is an either/or return value).
#' @export
getArchRGenome <- function(
  geneAnnotation=FALSE, 
  genomeAnnotation=FALSE
  ){

  .validInput(input = geneAnnotation, name = "geneAnnotation", valid = "boolean")
  .validInput(input = genomeAnnotation, name = "genomeAnnotation", valid = "boolean")

  supportedGenomes <- c("hg19","hg38","mm9","mm10","hg19test")
  .ArchRGenome <- options()[["ArchR.genome"]]

  if(!is.null(.ArchRGenome)){

    ag <- .ArchRGenome
    
    if(!is.character(ag)){
    
      return(NULL)
    
    }else{
      
      if(tolower(ag) %in% supportedGenomes){
        
        genome <- paste(toupper(substr(ag, 1, 1)), substr(ag, 2, nchar(ag)), sep="")

        if(geneAnnotation & genomeAnnotation){
          stop("Please request either geneAnnotation or genomeAnnotation, not both!")
        }

        if(geneAnnotation){

          message("Using GeneAnnotation set by addArchRGenome(",ag,")!")

          geneAnno <- paste0("geneAnno", genome)
          eval(parse(text=paste0("data(geneAnno",genome,")")))
          return(eval(parse(text=geneAnno)))
        
        }else if(genomeAnnotation){

          message("Using GeneAnnotation set by addArchRGenome(",ag,")!")

          genomeAnno <- paste0("genomeAnno", genome)
          eval(parse(text=paste0("data(genomeAnno",genome,")")))
          return(eval(parse(text=genomeAnno)))

        }else{

          return(genome)

        }

      }else{
        
        stop("option(ArchR.genome) : ", ag, " is not currently supported by ArchR. \nDid you mistakenly set this to a value without addArchRGenome?")
      
      }
    }

  }else{
    
    return(NULL)

  }

}

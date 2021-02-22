####################################
# Log Tools
####################################

#' Set ArchR Logging
#' 
#' This function will set ArchR logging
#'
#' @param useLogs A boolean describing whether to use logging with ArchR.
#' @export
addArchRLogging <- function(useLogs = TRUE){
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  message("Setting ArchRLogging = ", useLogs)
  options(ArchR.logging = useLogs)
  return(invisible(0))
}

#' Get ArchR Logging
#' 
#' This function will get ArchR logging
#'
#' @export
getArchRLogging <- function(){
  ArchRLogging <- options()[["ArchR.logging"]]
  if(!is.logical(ArchRLogging)){
    options(ArchR.logging = TRUE)
    return(TRUE)
  }
  ArchRLogging
}

#' Set ArchR Debugging
#' 
#' This function will set ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @param debug A boolean describing whether to use logging with ArchR.
#' @export
addArchRDebugging <- function(debug = FALSE){
  .validInput(input = debug, name = "debug", valid = "boolean")
  message("Setting ArchRDebugging = ", debug)
  options(ArchR.logging = debug)
  return(invisible(0))
}

#' Get ArchR Debugging
#' 
#' This function will get ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @export
getArchRDebugging <- function(){
  ArchRDebugging <- options()[["ArchR.debugging"]]
  if(!is.logical(ArchRDebugging)){
    options(ArchR.debugging = FALSE)
    return(FALSE)
  }
  ArchRDebugging
}

#' Set ArchR Verbosity for Log Messaging
#' 
#' This function will set ArchR logging verbosity.
#'
#' @param verbose A boolean describing whether to printMessages in addition to logging with ArchR.
#' @export
addArchRVerbose <- function(verbose = TRUE){
  .validInput(input = verbose, name = "verbose", valid = "boolean")
  message("Setting addArchRVerbose = ", verbose)
  options(ArchR.verbose = verbose)
  return(invisible(0))
}

#' Set ArchR Verbosity for Log Messaging
#' 
#' This function will get ArchR logging verbosity.
#'
#' @export
getArchRVerbose <- function(){
  ArchRVerbose <- options()[["ArchR.verbose"]]
  if(!is.logical(ArchRVerbose)){
    options(ArchR.verbose = TRUE)
    return(TRUE)
  }
  ArchRVerbose
}

#' Create a Log File for ArchR
#' 
#' This function will create a log file for ArchR functions. If ArchRLogging is not TRUE
#' this function will return NULL.
#'
#' @param name A character string to add a more descriptive name in log file.
#' @param logDir The path to a directory where log files should be written.
#' @export
createLogFile <- function(
    name = NULL,
    logDir = "ArchRLogs",
    useLogs = getArchRLogging()
  ){
  
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = logDir, name = "logDir", valid = "character")
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")

  if(!useLogs){
    return(NULL)
  }
  dir.create(logDir, showWarnings = FALSE)
  if(is.null(name)){
    logFile <- .tempfile(pattern = "ArchR", fileext = ".log", tmpdir = logDir)
  }else{
    logFile <- .tempfile(pattern = paste0("ArchR-", name), fileext = ".log", tmpdir = logDir)
  }
  logFile
}

.messageDiffTime <- function(...){ #Deprecated
  .logDiffTime(...)
}

.logDiffTime <- function(
  main = "",
  t1 = NULL,
  verbose = TRUE,
  addHeader = FALSE,
  t2 = Sys.time(),
  units = "mins",
  header = "###########",
  tail = "elapsed.",
  precision = 3,
  logFile = NULL,
  useLogs = getArchRLogging()
  ){

  if(verbose){
    
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),precision))
      if(addHeader){
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
      }else{
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
      }
      if(getArchRVerbose()) message(msg)
    }, error = function(x){
      if(getArchRVerbose()) message("Time Error : ", x)
    })

  }

  if(!useLogs){
    return(invisible(0))
  }

  if(!is.null(logFile)){
    if(file.exists(logFile)){
      logStamp <- tryCatch({
        dt <- abs(round(difftime(t2, t1, units = units),precision))
        if(addHeader){
          msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
        }else{
          msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
        }
        cat(paste0(msg,"\n"), file = logFile, append = TRUE)
      }, error = function(x){
        0
      })
    }
  }
  
  return(invisible(0))

}

.startLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
  ){

  if(!useLogs){
    return(invisible(0))
  }

  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(file.exists(logFile)){
    return(invisible(0))
  }

  .getRam <- function(OS = .Platform$OS.type){
  if(grepl("linux", OS, ignore.case = TRUE)){
    ram <- paste0("Linux : ", as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)))
  }else if(grepl("unix", OS, ignore.case = TRUE)){
    ram <- system("/usr/sbin/system_profiler SPHardwareDataType", intern = TRUE)
    ram <- paste0("MAC : ", gsub("Memory:","",gsub(" ","", grep("Memory", ram, value = TRUE))))
  }else{
    ram <- NA
  }
  }

  if(getArchRVerbose()) message("ArchR logging to : ", logFile, 
    "\nIf there is an issue, please report to github with logFile!")
  
  #Begin With
  cat(.ArchRLogo(ascii = "Package", messageLogo = FALSE), file = logFile, append = FALSE) 
  cat("\nLogging With ArchR!\n\n", file = logFile, append = TRUE) 
  cat(paste0("Start Time : ",Sys.time(),"\n\n"), file = logFile, append = TRUE)

  #ArchR Info
  cat("------- ArchR Info\n\n", file = logFile, append = TRUE)
  cat(paste0("ArchRThreads = ", getArchRThreads()), file = logFile, append = TRUE)
  tryCatch({
    if(!is.null(getArchRGenome())){
      cat(paste0("\nArchRGenome = ", getArchRGenome()), file = logFile, append = TRUE)
    }
  }, error = function(x){
  })
  cat("\n\n", file = logFile, append = TRUE)

  #Add Info
  cat("------- System Info\n\n", file = logFile, append = TRUE)
  cat(paste0("Computer OS = ", .Platform$OS.type), file = logFile, append = TRUE)
  tryCatch({
    cat(paste0("\nTotal Cores = ", detectCores()), file = logFile, append = TRUE)
  }, error = function(x){
  })
  # tryCatch({
  #     cat(paste0("\nTotal RAM = ", .getRam()), file = logFile, append = TRUE)
  # }, error = function(x){
  # })
  cat("\n\n", file = logFile, append = TRUE)

  #Session Info
  cat("------- Session Info\n\n", file = logFile, append = TRUE)
  utils::capture.output(sessionInfo(), file = logFile, append = TRUE)
  cat("\n\n------- Log Info\n\n", file = logFile, append = TRUE)

  return(invisible(0))

}

.logMessage <- function(
  ..., 
  logFile = NULL,
  verbose = FALSE,   
  useLogs = getArchRLogging()
  ){

  if(getArchRVerbose()){
    msg <- utils::capture.output(message(...), type = "message")
    msg <- paste0(msg, collapse = "\n")
  }else{
    msg <- "SuppressedMessaged due to getArchRVerbose() is FALSE!"
  }

  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }

  if(verbose){
    message(sprintf("%s : %s", Sys.time(), msg))
  }

  if(!useLogs){
    return(invisible(0))
  }

  if(is.null(logFile)){
    return(invisible(0))
  }
    
  cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)

  return(invisible(0))

}

.logHeader <- function(
  name = NULL, 
  logFile = NULL,   
  useLogs = getArchRLogging()
  ){
  
  if(!useLogs){
    return(invisible(0))
  }

  if(is.null(logFile)){
    return(invisible(0))
  }

  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
    
  header <- "###########"
  cat(sprintf("\n%s\n%s : %s\n%s\n\n", header, Sys.time(), name, header), file = logFile, append = TRUE)

  return(invisible(0))
}

.logStop <- function(
  ..., 
  logFile = NULL,
  useLogs = getArchRLogging()
  ){

  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
    
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }

  if(useLogs){
    if(!is.null(logFile)){
      cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
    }
  }

  stop(sprintf("%s\n", msg), call. = FALSE)

  return(invisible(0))

}

.logError <- function(
  e = NULL,
  fn = NULL,
  info = NULL, 
  errorList = NULL,
  logFile = NULL,   
  useLogs = getArchRLogging(),
  throwError = TRUE,
  debug = getArchRDebugging()
  ){

  header <- "************************************************************"

  if(is.null(logFile)){
    useLogs <- FALSE
  }

  if(useLogs){
    #To Log File
    cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile), file = logFile, append = TRUE)

    utils::capture.output(print(e), file = logFile, append = TRUE)

    if(!is.null(errorList)){
      tryCatch({
        #.safeSaveRDS(errorList, "Save-Error.rds")
        .logThis(errorList, name = "errorList", logFile)
      }, error = function(e){
        cat("Error recording errorList", file = logFile, append = TRUE)
      })
    }

    cat(sprintf("\n%s\n\n", header), file = logFile, append = TRUE)
  }
  
  #To Console
  cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile))
  
  if(debug){
    if(!is.null(errorList)){
      debugFile <- paste0(gsub("\\.log", "", logFile), "-debug.rds")
      cat(sprintf("\n%s : ArchRDebugging is set to TRUE, DebugFile = %s\n", Sys.time(), debugFile))
      .safeSaveRDS(errorList, debugFile)
    }
  }

  print(e)

  cat(sprintf("\n%s\n\n", header))

  if(throwError) stop("Exiting See Error Above")

  return(invisible(0))

}


.logThis <- function(
    x = NULL, 
    name = NULL, 
    logFile = NULL, 
    useLogs = getArchRLogging()
  ){

  if(!useLogs){
    return(invisible(0))
  }

  if(is.null(logFile)){
      return(invisible(0))
  }

  if(!file.exists(logFile)){
    stop("logFile does not exist! Something may have deleted this file! Exiting...")
  }
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  cat(paste0("\n", Sys.time(), " : ", name, ", Class = ", class(x), "\n"), file = logFile, append = TRUE)

  if(missing(x)){
    cat("Data is Missing\n\n", file = logFile, append = TRUE)
    return(invisible(0))
  }

  if(is.matrix(x)){
    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)

  }else if(is.data.frame(x)){

    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)

  }else if(is(x, "dgCMatrix") | is(x, "dgeMatrix")){

    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    cat(paste0(name, ": NonZeroEntries = ", length(x@x), ", EntryRange = [ ", paste0(range(x@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)

  }else if(is(x, "GRanges")){

    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))

  }else if(is(x, "SummarizedExperiment")){

    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))

  }else if(is(x, "DataFrame")){

    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))

  }else if(is(x, "ArchRProj")){

    suppressMessages(utils::capture.output(print(proj), file = logFile, append = TRUE))

  }else if(is(x, "SimpleList") | is(x, "list")){

    for(i in seq_along(x)){

      y <- x[[i]]

      if(missing(y)){
        next
      }

      if(is.matrix(y)){

        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)

      }else if(is.data.frame(y)){

        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)

      }else if(is(y, "dgCMatrix")){

        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": NonZeroEntries = ", length(y@x), ", EntryRange = [ ", paste0(range(y@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)

      }else if(is(y, "SimpleList") | is(y, "list")){

        for(j in seq_along(y)){

          z <- y[[j]]

          if(missing(z)){
            next
          }

          if(is.matrix(z)){

            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)

          }else if(is.data.frame(z)){

            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)

          }else if(is(z, "dgCMatrix")){

            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": NonZeroEntries = ", length(z@x), ", EntrzRange = [ ", paste0(range(z@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)

          }else if(is(y, "SimpleList") | is(y, "list")){

            #Only print 2x nested lists

          }else{

            tryCatch({
              cat("\n", file = logFile, append = TRUE)
              cat(paste0(paste0(name,"$", names(y[j])), ": length = ", length(z), "\n"), file = logFile, append = TRUE)
              suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
              cat("\n", file = logFile, append = TRUE)
            }, error = function(q){
            })

          }

        }

      }else{

        tryCatch({
          cat("\n", file = logFile, append = TRUE)
          cat(paste0(paste0(name,"$", names(x[i])), ": length = ", length(y), "\n"), file = logFile, append = TRUE)
          suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
          cat("\n", file = logFile, append = TRUE)
        }, error = function(q){
        })

      }

    }

  }else{

    tryCatch({
      cat("\n", file = logFile, append = TRUE)
      cat(paste0(name, ": length = ", length(x), "\n"), file = logFile, append = TRUE)
      suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
      cat("\n", file = logFile, append = TRUE)
    }, error = function(q){
    })

  }

  cat("\n", file = logFile, append = TRUE)
  return(invisible(0))

}

.endLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
  ){
  
  if(!useLogs){
    return(invisible(0))
  }

  if(is.null(logFile)){
    return(invisible(0))
  }

  rL <- readLines(logFile)
  o <- tryCatch({
    t1 <- gsub("Start Time : ","", grep("Start Time", rL, ignore.case = TRUE, value = TRUE))
    mn <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "mins"))
    hr <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "hours"))
    cat("\n------- Completed\n\n", file = logFile, append = TRUE)
    cat(paste0("End Time : ",Sys.time(),"\n"), file = logFile, append = TRUE)
    cat(paste0("Elapsed Time Minutes = ", mn), file = logFile, append = TRUE)
    cat(paste0("\nElapsed Time Hours = ", hr), file = logFile, append = TRUE)
    cat("\n\n-------\n\n\n\n", file = logFile, append = TRUE)
    if(getArchRVerbose()) message("ArchR logging successful to : ", logFile)
  }, error = function(x){
  })

  # tryCatch({
  #   R.utils::gzip(logFile, paste0(logFile, ".gz"))
  #   message("ArchR logging successful to : ", paste0(logFile, ".gz"))
  # }, error = function(x){
  # })
  
  return(invisible(0))

}






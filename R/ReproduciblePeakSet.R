####################################################################
# Peak Set Creation Methods
####################################################################

#' Add a Reproducible Peak Set to an ArchRProject
#' 
#' This function will get insertions from coverage files, call peaks,
#' and merge peaks to get a "Union Reproducible Peak Set"
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for peak calling.
#' @param reproducibility A string that indicates how peak reproducibility should be handled. This string is dynamic and can be a function of `n` where `n` is the number of samples being assessed. For example, `reproducibility = "2"` means at least 2 samples must have a peak call at this locus and `reproducibility = "(n+1)/2"` means that the majority of samples must have a peak call at this locus.
#' @param peaksPerCell The limit of number of peaks that can be identified per cell (this is useful for controlling how many peaks can be called from low cell groups).
#' @param maxPeaks A numeric threshold for the maximum peaks to retain per group in `groupBy` in the union reproducible peak set.
#' @param minCells The minimum number of unique cells that was used to create the coverage files on which peaks are called. This is important to allow for exclusion of pseudo-bulk replicates derived from very low cell numbers.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from peak calling.
#' @param pathToMacs2 The full path to the MACS2 executable.
#' @param genomeSize The genome size to be used for MACS2 peak calling (see MACS2 documentation).
#' @param shift The number of basepairs to shift each Tn5 insertion. When combined with `extsize` this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).
#' @param extsize The number of basepairs to extend the MACS2 fragment after `shift` has been applied. When combined with `extsize` this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).
#' @param method The method to use for significance testing in MACS2. Options are "p" for p-value and "q" for q-value. When combined with `cutOff` this gives the method and significance threshold for peak calling (see MACS2 documentation).
#' @param cutOff The numeric significance cutoff for the testing method indicated by `method` (see MACS2 documentation).
#' @param additionalParams A string of additional parameters to pass to MACS2 (see MACS2 documentation).
#' @param extendSummits The number of basepairs to extend peak summits (in both directions) to obtain final fixed-width peaks. For example, `extendSummits = 250` will create 501-bp fixed-width peaks from the 1-bp summits.
#' @param promoterDist The maximum distance in basepairs from the peak summit to the nearest transcriptional start site to allow for a peak to be annotated as a "promoter" peak.
#' @param genomeAnno The genomeAnnotation (see createGenomeAnnotation) is used for downstream analyses for genome information such as nucleotide information (GC info) or chromosome sizes.
#' @param geneAnno The geneAnnotation (see createGeneAnnotation) is used for peak labeling as promoter etc.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the reproducible peak set to be overwritten if it already exist in the given `ArchRProject` peakSet.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... additional args
#' @export
addReproduciblePeakSet <- function(
	ArchRProj = NULL,
	groupBy = "Clusters",
	reproducibility = "2",
	peaksPerCell = 500,
	maxPeaks = 150000,
	minCells = 25,
	excludeChr = c("chrM","chrY"),
	pathToMacs2 = findMacs2(),
	genomeSize = NULL, 
	shift = -75, 
	extsize = 150, 
	method = "q",
	cutOff = 0.1, 
	additionalParams = "--nomodel --nolambda",
	extendSummits = 250,
	promoterDist = 500,
	genomeAnno = getGenomeAnnotation(ArchRProj),
	geneAnno = getGeneAnnotation(ArchRProj),
	threads = getArchRThreads(),
	parallelParam = "mclapply",
	force = FALSE,
	verboseHeader = TRUE,
	verboseAll = FALSE,
	...
	){

	.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
	.validInput(input = groupBy, name = "groupBy", valid = c("character"))
	.validInput(input = reproducibility, name = "reproducibility", valid = c("character"))
	.validInput(input = peaksPerCell, name = "peaksPerCell", valid = c("integer"))
	.validInput(input = maxPeaks, name = "maxPeaks", valid = c("integer"))
	.validInput(input = minCells, name = "minCells", valid = c("integer"))
	.validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
	.validInput(input = pathToMacs2, name = "pathToMacs2", valid = c("character"))
	.validInput(input = genomeSize, name = "genomeSize", valid = c("character", "numeric", "null"))
	.validInput(input = shift, name = "shift", valid = c("integer"))
	.validInput(input = extsize, name = "extsize", valid = c("integer"))
	.validInput(input = method, name = "method", valid = c("character"))
	.validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
	.validInput(input = additionalParams, name = "additionalParams", valid = c("character"))
	.validInput(input = extendSummits, name = "extendSummits", valid = c("integer"))
	.validInput(input = promoterDist, name = "promoterDist", valid = c("integer"))
	geneAnno <- .validGeneAnnotation(geneAnno)
	genomeAnno <- .validGeneAnnotation(genomeAnno)
	.validInput(input = threads, name = "threads", valid = c("integer"))
	.validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
	.validInput(input = force, name = "force", valid = c("boolean"))
	.validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
	.validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))

	tstart <- Sys.time()
	utility <- .checkPath(pathToMacs2)

	coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy, minCells = minCells)
	coverageParams <- .getCoverageParams(ArchRProj = ArchRProj, groupBy = groupBy)

	#####################################################
	# Peak Calling Summary
	#####################################################
	tableGroups <- table(getCellColData(ArchRProj, groupBy, drop = TRUE))
	groupSummary <- lapply(seq_along(coverageParams$cellGroups), function(y){
		x <- coverageParams$cellGroups[[y]]
		uniq <- unique(unlist(x))
		n <- lapply(x, length) %>% unlist %>% sum
		nmin <- lapply(x, length) %>% unlist %>% min
		nmax <- lapply(x, length) %>% unlist %>% max
		data.frame(
		  Group=names(coverageParams$cellGroups)[y], 
		  nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
		  nCellsUsed=length(uniq), 
		  nReplicates=length(x), 
		  nMin=nmin, 
		  nMax=nmax, 
		  maxPeaks = min(maxPeaks, length(uniq) * peaksPerCell)
		)
	}) %>% Reduce("rbind",.)

	.messageDiffTime("Peak Calling Parameters!", tstart)
	#printSummary <- groupSummary
	#rownames(printSummary) <- NULL
	print(groupSummary)

	#####################################################
	# Create Output Directory
	#####################################################
	outDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls")
	outSubDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "ReplicateCalls")
	outBedDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "InsertionBeds")
	dir.create(outDir, showWarnings = FALSE)
	dir.create(outSubDir, showWarnings = FALSE)
	dir.create(outBedDir, showWarnings = FALSE)

	#####################################################
	# Genome Size Presets
	#####################################################
	if(is.null(genomeSize)){
		if(grepl("hg19|hg38", getGenome(ArchRProj), ignore.case = TRUE)){
			genomeSize <- 2.7e9
		}else if(grepl("mm9|mm10", getGenome(ArchRProj), ignore.case = TRUE)){
			genomeSize <- 1.87e9
		}
	}

	#####################################################
	# Arguments for Peak Calling
	#####################################################
	coverageFiles <- coverageMetadata$File
	names(coverageFiles) <- coverageMetadata$Name
	args <- list()
	args$X <- seq_len(nrow(coverageMetadata))
	args$FUN <- .callSummitsOnCoverages
	args$coverageFiles <- coverageFiles
	args$outFiles <- file.path(outSubDir, paste0(coverageMetadata$Name,"-summits.rds"))
	args$bedDir <- outBedDir
	args$excludeChr <- excludeChr
	args$peakParams <- list(
			pathToMacs2 = pathToMacs2,
			genomeSize = genomeSize, 
			shift = shift, 
			extsize = extsize, 
			cutOff = cutOff, 
			method = method,
			additionalParams = additionalParams
		)
	args$parallelParam <- parallelParam
	args$threads <- threads
	args$registryDir <- file.path(outDir, "batchRegistry")

	#####################################################
	# Batch Call Peaks
	#####################################################
	.messageDiffTime("Batching Peak Calls!", tstart)
		
	#back lapply
	outSummmits <- unlist(.batchlapply(args))
	
	#Summarize Output
	outSummitList <- split(outSummmits, coverageMetadata$Group)
	summitNamesList <- split(coverageMetadata$Name, coverageMetadata$Group)

	#####################################################
	# BSgenome for Add Nucleotide Frequencies!
	#####################################################
	.requirePackage(genomeAnno$genome)
	.requirePackage("Biostrings")
	BSgenome <- eval(parse(text = genomeAnno$genome))
	BSgenome <- .validBSgenome(BSgenome)

	#####################################################
	# Identify Reproducible Peaks!
	#####################################################
	.messageDiffTime("Identifying Reproducible Peaks!", tstart)
	groupPeaks <- .safelapply(seq_along(outSummitList), function(i){
		.messageDiffTime(sprintf("Creating Reproducible Peaks for Group %s of %s", i, length(outSummitList)), tstart)
		peaks <- suppressMessages(.identifyReproduciblePeaks(
			summitFiles = outSummitList[[i]],
			summitNames = summitNamesList[[i]],
			reproducibility = reproducibility,
	    	extendSummits = extendSummits,
	    	blacklist = genomeAnno$blacklist
		))
		peaks <- sort(sortSeqlevels(peaks))
		peaks <- subsetByOverlaps(peaks, genomeAnno$chromSizes, type = "within")
		peaks <- .fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnno = geneAnno, promoterDist = promoterDist)
		peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
		peaks <- peaks[order(peaks$groupScoreQuantile, decreasing = TRUE)]
		peaks <- head(peaks, groupSummary[names(outSummitList)[i],"maxPeaks"])
		mcols(peaks)$N <- NULL #Remove N Column
		saveRDS(peaks, file.path(outDir, paste0(names(outSummitList)[i], "-reproduciblePeaks.gr.rds")))
		return(peaks)
	}, threads = threads)
	names(groupPeaks) <- names(outSummitList)

	#Construct Union Peak Set
	.messageDiffTime("Creating Union Peak Set!", tstart)
	unionPeaks <- Reduce("c",groupPeaks)
	unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)

	#Summarize Output
	peakDF <- lapply(seq_along(groupPeaks), function(x){
		data.frame(Group = names(groupPeaks)[x], table(groupPeaks[[x]]$peakType))
	}) %>% Reduce("rbind", .)
	peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
	peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
	peakDF$Freq <- peakDF$Freq / 1000
	metadata(unionPeaks)$PeakCallSummary <- peakDF

	#Add Peak Set
	ArchRProj <- addPeakSet(ArchRProj, unionPeaks, force = TRUE)

	plotPDF(.plotPeakCallSummary(ArchRProj), name = "Peak-Call-Summary", width = 6, height = 4, ArchRProj = ArchRProj)

	.messageDiffTime(sprintf("Finished Creating Union Peak Set (%s)!", length(unionPeaks)), tstart)

    suppressWarnings(sink())

	return(ArchRProj)

}

#' @export
.plotPeakCallSummary <- function(ArchRProj, pal = NULL){

  peakDF <- metadata(ArchRProj@peakSet)$PeakCallSummary
  if(is.null(peakDF)){
    stop("Error no Peak Call Summary available are you sure these peaks were called with CreateReproduciblePeakSet?")
  }
  if(is.null(pal)){
    pal <- paletteDiscrete(values=peakDF$Var1)
  }

  lengthMax <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% max

  colnames(peakDF)[colnames(peakDF)=="Var1"] <- "PeakAnno"

  p <- ggplot(peakDF, aes(x=Group, y=Freq, fill=PeakAnno)) + 
    geom_bar(stat = "identity") + 
    theme_ArchR(xText90 = TRUE) +
    ylab("Number of Peaks (x10^3)") +
    xlab("") + theme(legend.position = "bottom", legend.key = element_rect(size = 2), legend.box.background = element_rect(color = NA)) +
    scale_fill_manual(values=pal) +
    scale_y_continuous(
      breaks = seq(0, lengthMax * 2,50), 
      limits = c(0, lengthMax * 1.1), 
      expand = c(0,0)
      )

  attr(p, "ratioYX") <- 0.5

  return(p)

}

#####################
# Utility Functions
#####################

.fastAnnoPeaks <- function(peaks, BSgenome, geneAnno, promoterDist = 1000){

	#Validate
	peaks <- .validGRanges(peaks)
	peakSummits <- resize(peaks,1,"center")
	geneAnno$genes <- .validGRanges(geneAnno$genes)
	geneAnno$exons <- .validGRanges(geneAnno$exons)
	geneAnno$TSS <- .validGRanges(geneAnno$TSS)
	BSgenome <- .validBSgenome(BSgenome)

	#First Lets Get Distance to Nearest Gene Start
	distPeaks <- distanceToNearest(peakSummits, resize(geneAnno$genes, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
	mcols(peaks)$nearestGene <- mcols(geneAnno$genes)$symbol[subjectHits(distPeaks)]
	og <- overlapsAny(peakSummits, geneAnno$genes, ignore.strand = TRUE)
	oe <- overlapsAny(peakSummits, geneAnno$exons, ignore.strand = TRUE)
	type <- rep("Distal", length(peaks))
	type[which(og & oe)] <- "Exon"
	type[which(og & !oe)] <- "Intron"
	type[mcols(peaks)$distToGeneStart < promoterDist] <- "Promoter"
	mcols(peaks)$peakType <- type

	#First Lets Get Distance to Nearest TSS's
	distTSS <- distanceToNearest(peakSummits, resize(geneAnno$TSS, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToTSS <- mcols(distTSS)$distance
	if("symbol" %in% colnames(mcols(geneAnno$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnno$TSS)$symbol[subjectHits(distPeaks)]
	}else if("tx_name" %in% colnames(mcols(geneAnno$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnno$TSS)$tx_name[subjectHits(distPeaks)]
	}

	#Get NucleoTide Content
	nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  	mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  	mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  	peaks

}

.identifyReproduciblePeaks <- function(
	summitFiles = NULL,
	summitNames = NULL,
    reproducibility = 0.51,
    extendSummits = 250,
    blacklist
  ){

  start <- Sys.time()
  summits <- lapply(seq_along(summitFiles), function(x){
  	grx <- readRDS(summitFiles[x])
  	grx <- subsetByOverlaps(grx, blacklist, invert = TRUE) #Not Overlapping Blacklist!
  	grx$GroupReplicate <- paste0(summitNames[x])
  	grx
  })
  summits <- Reduce("c", as(summits, "GRangesList"))

  extendedSummits <- resize(summits, extendSummits * 2 + 1, "center")
  extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
    nonES <- nonOverlappingGR(x, by = "score", decreasing = TRUE)
    nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score),3)
    nonES
  })
  extendedSummits <- Reduce("c", as(extendedSummits, "GRangesList"))

  nonOverlapES <- nonOverlappingGR(extendedSummits, by = "replicateScoreQuantile", decreasing = TRUE)

  overlapMat <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
    overlapsAny(nonOverlapES, x)
  }) %>% Reduce("cbind", .)

  if(length(summitFiles) > 1){
    nonOverlapES$Reproducibility <- rowSums(overlapMat)
    nonOverlapES$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
    n <- length(summitFiles)
    minRep <- eval(parse(text=reproducibility))
    if(!is.numeric(minRep)){
    	stop("Error reproducibility not numeric when evaluated!")
    }
  	idxPass <- which(nonOverlapES$Reproducibility >= minRep)
  	nonOverlapPassES <- nonOverlapES[idxPass]
  }else{
    nonOverlapES$Reproducibility <- rep(NA, length(nonOverlapES))
    nonOverlapPassES <- nonOverlapES
  }

  nonOverlapPassES$groupScoreQuantile <- round(.getQuantiles(nonOverlapPassES$replicateScoreQuantile),3)
  mcols(nonOverlapPassES) <- mcols(nonOverlapPassES)[,c("score","replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "GroupReplicate")]

  return(nonOverlapPassES)

}

.callSummitsOnCoverages <- function(i, coverageFiles, outFiles, peakParams, bedDir, excludeChr, tstart, ...){
	
	.messageDiffTime(sprintf("Group %s of %s, Calling Peaks with MACS2!", i, length(coverageFiles)), tstart)

	################
	# Create Bed File from Coverage File
	################
	bedFile <- file.path(bedDir, paste0(names(coverageFiles)[i], ".insertions.bed"))
	o <- ArchR:::.writeCoverageToBed(coverageFiles[i], bedFile, excludeChr = excludeChr)
	peakParams$bedFile <- bedFile
	
	################
	# MACS2 Peak-Calling Leave Room For Other Options?
	################
	summits <- do.call(.callSummitsMACS2, peakParams)
	rmf <- file.remove(bedFile)
	
	################
	# Save output
	################
	saveRDS(summits, outFiles[i])

	.messageDiffTime(sprintf("Group %s of %s, Finished Calling Peaks with MACS2!", i, length(coverageFiles)), tstart)

	outFiles[i]

}

.callSummitsMACS2 <- function(
	bedFile = NULL,
	pathToMacs2 = "macs2",
	genomeSize = 2.7e9, 
	shift = -75, 
	extsize = 150, 
	cutOff = 0.05, 
	method = "q",
	additionalParams = "--nomodel --nolambda",
	...
	){

	stopifnot(tolower(method) %in% c("p","q"))
	stopifnot(!is.null(genomeSize))
	utility <- ArchR:::.checkPath(pathToMacs2)

	#Output Files
	bedName <- gsub("\\.insertions.bed", "", bedFile)
	summitsFile <- paste0(bedName, "_summits.bed")
	narrowPeaksFile <- paste0(bedName, "_peaks.narrowPeak")
	xlsFile <- paste0(bedName, "_peaks.xls")

	#Create MACS2 Command
	cmd <- sprintf("callpeak -g %s --name %s --treatment %s --outdir %s --format BED --call-summits --keep-dup all %s", 
		genomeSize, basename(bedName), bedFile, dirname(bedName), additionalParams)

	if(!is.null(shift) & !is.null(extsize)){
		cmd <- sprintf("%s --shift %s --extsize %s", cmd , shift, extsize)
	}

	if(tolower(method) == "p"){
		cmd <- sprintf("%s -p %s", cmd , cutOff)
	}else{
		cmd <- sprintf("%s -q %s", cmd , cutOff)
	}

	run <- system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)

	#Read Summits!
	out <- fread(summitsFile, select = c(1,2,3,5))
	out <- GRanges(out$V1, IRanges(out$V2 + 1, out$V3), score = out$V5)

	#Remove Files
	r2 <- suppressWarnings(file.remove(summitsFile, narrowPeaksFile, xlsFile))

	return(out)

}

#' Find the installed location of the MACS2 executable
#' 
#' This function attempts to find the path to the MACS2 executable by serting the path and python's pip.
#'
#' @export
findMacs2 <- function(){
  
  message("Searching For MACS2...")

  #Check if in path
  if(.suppressAll(.checkPath("macs2", throwError = FALSE))){
  	message("Found with $path!")
    return("macs2")
  }
  #Try seeing if its pip installed
  search2 <- tryCatch({system2("pip", "show macs2", stdout = TRUE, stderr = NULL)}, error = function(x){"ERROR"})
  search3 <- tryCatch({system2("pip3", "show macs2", stdout = TRUE, stderr = NULL)}, error = function(x){"ERROR"})
  
  message("Not Found in $Path")

  if(search2[1] != "ERROR"){
	  path2Install <- gsub("Location: ","",search2[grep("Location", search2, ignore.case=TRUE)])
	  path2Bin <- gsub("lib/python/site-packages", "bin/macs2",path2Install)
	  if(.suppressAll(.checkPath(path2Bin, throwError = error))){
	  	message("Found with pip!")
	    return(path2Bin)
	  }else{
  		message("Not Found with pip")
	  }
  }

  if(search3[1] != "ERROR"){
	  path2Install <- gsub("Location: ","",search3[grep("Location", search3, ignore.case=TRUE)])
	  path2Bin <- gsub("lib/python/site-packages", "bin/macs2",path2Install)
	  if(.suppressAll(.checkPath(path2Bin, throwError = error))){
	  	message("Found with pip3!")
	    return(path2Bin)
	  }else{
  		message("Not Found with pip3")
	  }
  }
  
  cat(out)

}




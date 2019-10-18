#' Add Reproducible Peak Set to ArchR Project
#' 
#' This function will get insertions from coverage files call peaks
#' and merge to get a Union Reproducible Peak Set
#'
#' @param ArchRProj ArchRProject
#' @param groupBy use groupings for peak calling matching group coverage files
#' @param reproducibility reproducibility for peak calling (string that is a function of n)
#' @param peaksPerCell number of peaks that can be identified per cell on average
#' @param excludeChr exclude chromosomes from peak calling
#' @param pathToMacs2 path to macs2 executable (see Macs2)
#' @param genomeSize genome size for peak calling (see Macs2)
#' @param shift shift of Tn5 insertions (<- |   ) (see Macs2)
#' @param extsize extension of Tn5 insertions (|<- ->|) (see Macs2)
#' @param method significance method for Macs2 (see Macs2)
#' @param cutOff significance cutoff for Macs2 (see Macs2)
#' @param extendSummits extend peak summits for final fixed-width peaks
#' @param promoterDist promoter distance from TSS for annotating peaks
#' @param genomeAnno genome annotation for ArchRProject
#' @param geneAnno gene annotation for ArchRProject
#' @param additionalParams additional parameters to pass to Macs2 (see Macs2)
#' @param threads number of threads for parallel execution
#' @param parallelParam parallel parameters for batch style execution
#' @param force force creating peakSet if existed
#' @param verboseHeader verbose sections
#' @param verboseAll verbose sections and subsections
#' @param ... additional args
#' @export
addReproduciblePeakSet <- function(
	ArchRProj = NULL,
	groupBy = "Clusters",
	reproducibility = "2",
	peaksPerCell = 500,
	maxPeaks = 150000,
	excludeChr = c("chrM","chrY"),
	pathToMacs2 = "macs2",
	genomeSize = NULL, 
	shift = -75, 
	extsize = 150, 
	method = "q",
	cutOff = 0.05, 
	extendSummits = 250,
	promoterDist = 500,
	genomeAnno = getGenomeAnnotation(ArchRProj),
	geneAnno = getGeneAnnotation(ArchRProj),
	additionalParams = "--nomodel --nolambda",
	threads = 1,
	parallelParam = "mclapply",
	force = FALSE,
	verboseHeader = TRUE,
	verboseAll = FALSE,
	...
	){

	tstart <- Sys.time()
	utility <- ArchR:::.checkPath(pathToMacs2)

	coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy)
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
		  nCells=tableGroups[y], 
		  nCellsUsed=length(uniq), 
		  nReplicates=length(x), 
		  nMin=nmin, 
		  nMax=nmax, 
		  maxPeaks = min(maxPeaks, nmin * peaksPerCell)
		)
	}) %>% Reduce("rbind",.)

	.messageDiffTime("Peak Calling Parameters!", tstart)
	printSummary <- groupSummary
	rownames(printSummary) <- NULL
	print(printSummary)

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
	unionPeaks <- nonOverlappingGRanges(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)

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

	pdf(file.path(outDir, "PeakCallSummary.pdf"), width = 6, height = 4, onefile=FALSE)
	print(plotPeakCallSummary(ArchRProj))
	dev.off()

	.messageDiffTime("Finished Creating Union Peak Set!", tstart)

	return(ArchRProj)

}

#' @export
plotPeakCallSummary <- function(ArchRProj, pal = NULL){

  peakDF <- metadata(ArchRProj@peakSet)$PeakCallSummary
  if(is.null(peakDF)){
    stop("Error no Peak Call Summary available are you sure these peaks were called with CreateReproduciblePeakSet?")
  }
  if(is.null(pal)){
    pal <- paletteDiscrete(values=peakDF$Var1)
  }

  lengthMax <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% max

  p <- ggplot(peakDF, aes(x=Group, y=Freq, fill=Var1)) + 
    geom_bar(stat = "identity") + 
    theme_ArchR(rotate_x_axis_text_90 = TRUE) +
    ylab("Number of Peaks (x10^3)") +
    xlab("") +
    scale_fill_manual(values=pal) +
    scale_y_continuous(
      breaks = seq(0, lengthMax * 2,50), 
      limits = c(0, lengthMax * 1.1), 
      expand = c(0,0)
      )

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
  }) %>% Reduce("c", .)

  extendedSummits <- resize(summits, extendSummits * 2 + 1, "center")
  extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
    nonES <- nonOverlappingGRanges(x, by = "score", decreasing = TRUE)
    nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score),3)
    nonES
  }) %>% Reduce("c", .)
  nonOverlapES <- nonOverlappingGRanges(extendedSummits, by = "replicateScoreQuantile", decreasing = TRUE)

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
	o <- .writeCoverageToBed(coverageFiles[i], bedFile, excludeChr = excludeChr)
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


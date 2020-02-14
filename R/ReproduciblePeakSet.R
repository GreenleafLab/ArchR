####################################################################
# Peak Set Creation Methods
####################################################################

#' Add a Reproducible Peak Set to an ArchRProject
#' 
#' This function will get insertions from coverage files, call peaks,
#' and merge peaks to get a "Union Reproducible Peak Set".
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param groupBy The name of the column in `cellColData` to use for grouping cells together for peak calling.
#' @param peakMethod The name of peak calling method to be used. Options include "Macs2" for using macs2 callpeak or "Tiles" for using a TileMatrix.
#' @param reproducibility A string that indicates how peak reproducibility should be handled. This string is dynamic and can be a
#' function of `n` where `n` is the number of samples being assessed. For example, `reproducibility = "2"` means at least 2 samples
#' must have a peak call at this locus and `reproducibility = "(n+1)/2"` means that the majority of samples must have a peak call at this locus.
#' @param peaksPerCell The upper limit of the number of peaks that can be identified per cell-grouping in `groupBy`. This is useful
#' for controlling how many peaks can be called from cell groups with low cell numbers.
#' @param maxPeaks A numeric threshold for the maximum peaks to retain per group from `groupBy` in the union reproducible peak set.
#' @param minCells The minimum allowable number of unique cells that was used to create the coverage files on which peaks are called.
#' This is important to allow for exclusion of pseudo-bulk replicates derived from very low cell numbers.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from peak calling.
#' @param pathToMacs2 The full path to the MACS2 executable.
#' @param genomeSize The genome size to be used for MACS2 peak calling (see MACS2 documentation).
#' @param shift The number of basepairs to shift each Tn5 insertion. When combined with `extsize` this allows you to create proper fragments,
#' centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).
#' @param extsize The number of basepairs to extend the MACS2 fragment after `shift` has been applied. When combined with `extsize` this
#' allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).
#' @param method The method to use for significance testing in MACS2. Options are "p" for p-value and "q" for q-value. When combined with
#' `cutOff` this gives the method and significance threshold for peak calling (see MACS2 documentation).
#' @param cutOff The numeric significance cutOff for the testing method indicated by `method` (see MACS2 documentation).
#' @param additionalParams A string of additional parameters to pass to MACS2 (see MACS2 documentation).
#' @param extendSummits The number of basepairs to extend peak summits (in both directions) to obtain final fixed-width peaks. For example,
#' `extendSummits = 250` will create 501-bp fixed-width peaks from the 1-bp summits.
#' @param promoterDist The maximum distance in basepairs from the peak summit to the nearest transcriptional start site to allow for a peak
#' to be annotated as a "promoter" peak.
#' @param genomeAnnotation The genomeAnnotation (see `createGenomeAnnotation()`) to be used for generating peak metadata such as nucleotide
#' information (GC content) or chromosome sizes.
#' @param geneAnnotation The geneAnnotation (see `createGeneAnnotation()`) to be used for labeling peaks as "promoter", "exonic", etc.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the reproducible peak set to be overwritten if it already exist in the given `ArchRProject` peakSet.
#' @param verboseHeader A boolean value that determines whether standard output includes verbose sections.
#' @param verboseAll A boolean value that determines whether standard output includes verbose subsections.
#' @param ... Additional parameters to be pass to `addGroupCoverages()` to get sample-guided pseudobulk cell-groupings. Only used for TileMatrix-based
#' peak calling (not for MACS2). See `addGroupCoverages()` for more info.
#' @export
addReproduciblePeakSet <- function(
	ArchRProj = NULL,
	groupBy = "Clusters",
	peakMethod = "Macs2",
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
	genomeAnnotation = getGenomeAnnotation(ArchRProj),
	geneAnnotation = getGeneAnnotation(ArchRProj),
	threads = getArchRThreads(),
	parallelParam = NULL,
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

	if(tolower(peakMethod) == "macs2"){
		.validInput(input = pathToMacs2, name = "pathToMacs2", valid = c("character"))
		.validInput(input = genomeSize, name = "genomeSize", valid = c("character", "numeric", "null"))
		.validInput(input = shift, name = "shift", valid = c("integer"))
		.validInput(input = extsize, name = "extsize", valid = c("integer"))
		.validInput(input = method, name = "method", valid = c("character"))
		.validInput(input = additionalParams, name = "additionalParams", valid = c("character"))
		.validInput(input = extendSummits, name = "extendSummits", valid = c("integer"))		
	}else if(tolower(peakMethod) == "tiles"){

	}else{
		stop("peakMethod not recognized! Supported peakMethods are Macs2 or Tiles!")
	}

	.validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
	.validInput(input = promoterDist, name = "promoterDist", valid = c("integer"))
	geneAnnotation <- .validGeneAnnotation(geneAnnotation)
	genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
	.validInput(input = threads, name = "threads", valid = c("integer"))
	.validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
	.validInput(input = force, name = "force", valid = c("boolean"))
	.validInput(input = verboseHeader, name = "verboseHeader", valid = c("boolean"))
	.validInput(input = verboseAll, name = "verboseAll", valid = c("boolean"))

	tstart <- Sys.time()

	if(tolower(peakMethod) == "macs2"){

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
		.requirePackage(genomeAnnotation$genome)
		.requirePackage("Biostrings",source="bioc")
		BSgenome <- eval(parse(text = genomeAnnotation$genome))
		BSgenome <- validBSgenome(BSgenome)

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
		    	blacklist = genomeAnnotation$blacklist
			))
			peaks <- sort(sortSeqlevels(peaks))
			peaks <- subsetByOverlaps(peaks, genomeAnnotation$chromSizes, type = "within")
			peaks <- .fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterDist = promoterDist)
			peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
			peaks <- peaks[order(peaks$groupScoreQuantile, decreasing = TRUE)]
			peaks <- head(peaks, groupSummary[names(outSummitList)[i],"maxPeaks"])
			mcols(peaks)$N <- NULL #Remove N Column
			saveRDS(peaks, file.path(outDir, paste0(names(outSummitList)[i], "-reproduciblePeaks.gr.rds")))
			return(peaks)
		}, threads = threads) %>% GRangesList()
		names(groupPeaks) <- names(outSummitList)

		#Construct Union Peak Set
		.messageDiffTime("Creating Union Peak Set!", tstart)
		unionPeaks <- unlist(groupPeaks)
		unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)

		#Summarize Output
		peakDF <- lapply(seq_along(groupPeaks), function(x){
			data.frame(Group = names(groupPeaks)[x], table(groupPeaks[[x]]$peakType))
		}) %>% Reduce("rbind", .)
		peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
		peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
		peakDF$Freq <- peakDF$Freq / 1000
		metadata(unionPeaks)$PeakCallSummary <- peakDF

	}else if(tolower(peakMethod) == "tiles"){

		useMatrix <- "TileMatrix"
		
		cellGroups <- addGroupCoverages(
			ArchRProj = ArchRProj, 
			groupBy = groupBy,
			returnGroups = TRUE,
			minCells = minCells,
			...
		)[[1]]$cellGroups

		#####################################################
		# Peak Calling Summary
		#####################################################
		tableGroups <- table(getCellColData(ArchRProj, groupBy, drop = TRUE))
		groupSummary <- lapply(seq_along(cellGroups), function(y){
			x <- cellGroups[[y]]
			uniq <- unique(unlist(x))
			n <- lapply(x, length) %>% unlist %>% sum
			nmin <- lapply(x, length) %>% unlist %>% min
			nmax <- lapply(x, length) %>% unlist %>% max
			data.frame(
			  Group=names(cellGroups)[y], 
			  nCells=tableGroups[names(cellGroups)[y]], 
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
		# Peak Calling from TileMatrix
		#####################################################

		#MatrixFiles
		ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]
		chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)
	    featureDF <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix)
	    
	    #Determine Resolution
	    d1 <- featureDF$start[2] - featureDF$start[1]
	    d2 <- featureDF$start[3] - featureDF$start[2]
	    if(d1 != d2){
	    	stop("Something is wrong with TileMatrix, Could not determine a resolution!")
	    }else{
	    	res <- d1
	    }

		#Compute Row Sums Across All Samples
		.messageDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = verboseAll, verbose = verboseHeader)
		totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
		nTiles <- nrow(totalAcc)
		gc()

		#Pre-Filter 0s
		topFeatures <- totalAcc[which(totalAcc$rowSums != 0), ]

		#Group Matrix
		#Consider reading in group-wise if this is getting too large?
		.messageDiffTime("Computing Pseudo-Grouped Tile Matrix", tstart, addHeader = verboseAll, verbose = verboseHeader)
	    groupMat <- .getGroupMatrix(
	      ArrowFiles = ArrowFiles, 
	      featureDF = topFeatures,
	      useMatrix = useMatrix, 
	      threads = threads,
	      groupList = unlist(cellGroups),
	      useIndex = FALSE,
	      asSparse = TRUE,
	      verbose = verboseAll
	    )
		.messageDiffTime(sprintf("Created Pseudo-Grouped Tile Matrix (%s GB)", round(object.size(groupMat) / 10^9, 3)), tstart, addHeader = verboseAll, verbose = verboseHeader)

		#Expectation, we lowered the expectation (by 25%) to better match 
		#the results between macs2 and tile peak calling methods.
	    exp <- 0.75 * Matrix::colSums(groupMat) / nTiles

		#####################################################
		# Peak Calling Per Group
		#####################################################

	    groupPeaks <- .safelapply(seq_along(cellGroups), function(i){

	    	.messageDiffTime(sprintf("Computing Peaks from Tile Matrix Group (%s of %s)", i, length(cellGroups)), tstart, addHeader = verboseAll, verbose = verboseHeader)

	    	gx <- grep(paste0(names(cellGroups)[i],"."), colnames(groupMat))
		    
		    pMat <- lapply(seq_along(gx), function(x){
			    pval <- ppois(q = groupMat[,gx[x]]-1, lambda = exp[gx[x]], lower.tail = FALSE, log=FALSE)
			    if(tolower(method) == "q"){
			    	pval <- p.adjust(pval, method = "fdr", n = nTiles)
			    }else if(tolower(method)=="p"){
			    }else{
			    	stop("method should be either p for p-value or q for adjusted p-value!")
			    }
			    as(as.matrix(-log10(pmax(pval, 10^-250))), "dgCMatrix")
		    }) %>% Reduce("cbind", .)

		    n <- ncol(pMat)
		    passPeaks <- Matrix::rowSums(pMat >= -log10(cutOff)) >= eval(parse(text=reproducibility))
		    mlogp <- Matrix::rowSums(Matrix::t(Matrix::t(pMat) / Matrix::colSums(pMat)) * 10^6) / ncol(pMat)

		    rm(pMat)
		    gc()

		    passPeaks <- passPeaks[order(mlogp, decreasing=TRUE)]
		    mlogp <- mlogp[order(mlogp, decreasing=TRUE)]

		    nMax <- groupSummary[names(cellGroups)[i], "maxPeaks"]
		    passPeaks <- head(passPeaks, nMax)
		    mlogp <- head(mlogp, nMax)

		    mlogp <- mlogp[which(passPeaks)]
		    passPeaks <- passPeaks[which(passPeaks)]

		    if(length(passPeaks) > 0){
			    DataFrame(
			    	Group = Rle(names(cellGroups)[i]), 
			    	peaks = names(passPeaks), 
			    	mlog10p = mlogp, 
			    	normmlogp = 10^6 * mlogp / sum(mlogp)
			    )
		    }else{
		    	NULL
		    }

	    }, threads = threads) %>% Reduce("rbind", .)

	    groupPeaks <- groupPeaks[order(groupPeaks$normmlogp, decreasing=TRUE), ]

		#####################################################
		# BSgenome for Add Nucleotide Frequencies!
		#####################################################
		.requirePackage(genomeAnnotation$genome)
		.requirePackage("Biostrings",source="bioc")
		BSgenome <- eval(parse(text = genomeAnnotation$genome))
		BSgenome <- validBSgenome(BSgenome)
		outDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls")

		#####################################################
		# Create Group Peaks
		#####################################################
		.messageDiffTime("Creating Group Peak Sets with Annotations!", tstart)

		peakDF <- .safelapply(seq_along(cellGroups), function(i){
			
			groupPeaksi <- groupPeaks[groupPeaks$Group == names(cellGroups)[i], ]

			if(nrow(groupPeaksi) == 0){
				return(NULL)
			}
			
			groupPeaksi <- cbind(
				topFeatures[as.integer(gsub("f", "", groupPeaksi$peaks)),],
				groupPeaksi
			)
			groupPeaksi <- groupPeaksi[order(groupPeaksi$seqnames, groupPeaksi$idx), ]
			groupPeaksi$peaks <- NULL
			groupPeaksi$normmlogp <- NULL
			groupPeaksi$mlog10p <- round(groupPeaksi$mlog10p, 3)
			groupPeaksGRi <- GRanges(groupPeaksi$seqnames, 
				IRanges(
					start = pmax((groupPeaksi$idx - 1) * res, 1),
					end = (groupPeaksi$idx) * res - 1
				)
			)
			mcols(groupPeaksGRi) <- groupPeaksi[,c("mlog10p", "Group")]

			groupPeaksGRi <- subsetByOverlaps(groupPeaksGRi, genomeAnnotation$chromSizes, type = "within")
			groupPeaksGRi <- subsetByOverlaps(groupPeaksGRi, genomeAnnotation$blacklist, invert = TRUE)
			groupPeaksGRi <- .fastAnnoPeaks(
				groupPeaksGRi, 
				BSgenome = BSgenome, 
				geneAnnotation = geneAnnotation, 
				promoterDist = promoterDist
			)
			groupPeaksGRi <- groupPeaksGRi[which(mcols(groupPeaksGRi)$N < 0.001)] #Remove N Containing Peaks
			mcols(groupPeaksGRi)$N <- NULL #Remove N Column

			#Table Peak Types
			tabPT <- data.frame(Group = names(cellGroups)[i], table(groupPeaksGRi$peakType))

			#Save
			saveRDS(groupPeaksGRi, file.path(outDir, paste0(names(cellGroups)[i], "-reproduciblePeaks.gr.rds")))

			#Remove
			rm(groupPeaksGRi)
			gc()

			tabPT

		}, threads = threads) %>% Reduce("rbind", .)

		#####################################################
		# Call Union Peaks
		#####################################################
		.messageDiffTime("Creating Union Peak Set with Annotations!", tstart)
	    unionPeaks <- groupPeaks[!duplicated(groupPeaks$peaks), ]
	    rm(groupPeaks)
	    gc()

	    unionPeaks <- cbind(
	    	topFeatures[as.integer(gsub("f", "", unionPeaks$peaks)),],
	    	unionPeaks
	    )
	    unionPeaks <- unionPeaks[order(unionPeaks$seqnames, unionPeaks$idx), ]
	    unionPeaks$peaks <- NULL
	    unionPeaks$normmlogp <- NULL
	    unionPeaks$mlog10p <- round(unionPeaks$mlog10p, 3)

	    unionPeaksGR <- GRanges(unionPeaks$seqnames, 
	    	IRanges(
	    		start = pmax((unionPeaks$idx - 1) * res, 1),
	    		end = (unionPeaks$idx) * res - 1
	    	)
	    )
	    mcols(unionPeaksGR) <- unionPeaks[,c("mlog10p", "Group")]

		unionPeaksGR <- subsetByOverlaps(unionPeaksGR, genomeAnnotation$chromSizes, type = "within")
		unionPeaksGR <- subsetByOverlaps(unionPeaksGR, genomeAnnotation$blacklist, invert = TRUE)
		unionPeaksGR <- .fastAnnoPeaks(
			unionPeaksGR, 
			BSgenome = BSgenome, 
			geneAnnotation = geneAnnotation, 
			promoterDist = promoterDist
		)
		unionPeaksGR <- unionPeaksGR[which(mcols(unionPeaksGR)$N < 0.001)] #Remove N Containing Peaks
		mcols(unionPeaksGR)$N <- NULL #Remove N Column

		#Set To unionPeaks
		unionPeaks <- unionPeaksGR
		rm(unionPeaksGR)
		gc()

		#Summarize Output
		peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
		peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
		peakDF$Freq <- peakDF$Freq / 1000
		metadata(unionPeaks)$PeakCallSummary <- peakDF


	}else{
		stop("method not recognized! Supported methods are Macs2 or Tiles!")
	}	

	#Add Peak Set
	ArchRProj <- addPeakSet(ArchRProj, unionPeaks, force = TRUE)

	plotPDF(.plotPeakCallSummary(ArchRProj), name = "Peak-Call-Summary", width = 6, height = 4, ArchRProj = ArchRProj)

	.messageDiffTime(sprintf("Finished Creating Union Peak Set (%s)!", length(unionPeaks)), tstart)

    suppressWarnings(sink())

	return(ArchRProj)

}

.plotPeakCallSummary <- function(ArchRProj = NULL, pal = NULL){

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

.fastAnnoPeaks <- function(peaks = NULL, BSgenome = NULL, geneAnnotation = NULL, promoterDist = 1000){

	#Validate
	peaks <- .validGRanges(peaks)
	peakSummits <- resize(peaks,1,"center")
	geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
	geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
	geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
	BSgenome <- validBSgenome(BSgenome)

	#First Lets Get Distance to Nearest Gene Start
	distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
	mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
	og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
	oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
	type <- rep("Distal", length(peaks))
	type[which(og & oe)] <- "Exon"
	type[which(og & !oe)] <- "Intron"
	type[mcols(peaks)$distToGeneStart < promoterDist] <- "Promoter"
	mcols(peaks)$peakType <- type

	#First Lets Get Distance to Nearest TSS's
	distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToTSS <- mcols(distTSS)$distance
	if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
	}else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
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

.callSummitsOnCoverages <- function(
	i = NULL,
	coverageFiles = NULL,
	outFiles = NULL,
	peakParams = NULL,
	bedDir = NULL,
	excludeChr = NULL,
	subThreads = 1,
	tstart = NULL
	){
	
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

	#.messageDiffTime(sprintf("Group %s of %s, Finished Calling Peaks with MACS2!", i, length(coverageFiles)), tstart)

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
	additionalParams = "--nomodel --nolambda"
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
  
  message("Searching For MACS2..")

  #Check if in path
  if(.suppressAll(.checkPath("macs2", throwError = FALSE))){
  	message("Found with $path!")
    return("macs2")
  }

  message("Not Found in $PATH")

  #Try seeing if its pip installed
  search2 <- tryCatch({system2("pip", "show macs2", stdout = TRUE, stderr = NULL)}, error = function(x){"ERROR"})
  search3 <- tryCatch({system2("pip3", "show macs2", stdout = TRUE, stderr = NULL)}, error = function(x){"ERROR"})
  
  if(length(search2) > 0){
	  if(search2[1] != "ERROR"){
		  path2Install <- gsub("Location: ","",search2[grep("Location", search2, ignore.case=TRUE)])
		  path2Bin <- gsub("lib/python/site-packages", "bin/macs2",path2Install)
		  if(.suppressAll(.checkPath(path2Bin, throwError = error))){
		  	message("Found with pip!")
		    return(path2Bin)
		  }
	  }
  }

  message("Not Found with pip")

  if(length(search3) > 0){
	  if(search3[1] != "ERROR"){
		  path2Install <- gsub("Location: ","",search3[grep("Location", search3, ignore.case=TRUE)])
		  path2Bin <- gsub("lib/python/site-packages", "bin/macs2",path2Install)
		  if(.suppressAll(.checkPath(path2Bin, throwError = error))){
		  	message("Found with pip3!")
		    return(path2Bin)
		  }
	  }
  }

  message("Not Found with pip3")

  stop("Could Not Find Macs2! Please install w/ pip, add to your $PATH, or just supply the macs2 path directly and avoid this function!")

}




# Tests for arrow 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

#Context
context("test_arrow")

#Fragments
fragments <- file.path(system.file("testdata", package="ArchR"), "PBSmall.tsv.gz")

#Check
test_that("Test Fragments Exist...", {
	expect_equal(file.exists(fragments), TRUE)
})

################################################
# Testing Create Arrow File
################################################

addArchRGenome("hg19test2")

test_that("Test Genome is Correct Exist...", {
	check1 <- all(paste0(seqnames(getGenes())) %in% c("chr5", "chr11"))
	check2 <- all(paste0(seqnames(getChromSizes())) %in% c("chr5", "chr11"))
	expect_equal(check1, TRUE)
	expect_equal(check2, TRUE)
})

#Create Arrow Files
arrowFiles <- createArrowFiles(
	inputFiles = fragments,
	sampleNames = "PBSmall",
	minFrags = 100,
	nChunk = 1,
	TileMatParams=list(tileSize=10000),
	force = TRUE
)

#Get Contents
test_that("Checking Arrow Contents...", {

	expect_equal(
		ArchR:::.validArrow(arrowFiles),
		arrowFiles
	)

	expect_equal(
		ArchR:::.availableArrays(arrowFiles),
		c("GeneScoreMatrix", "TileMatrix")
	)

	expect_equal(
		nrow(ArchR:::.getFeatureDF(arrowFiles, "TileMatrix")),
		31593
	)

	expect_equal(
		nrow(ArchR:::.getFeatureDF(arrowFiles, "GeneScoreMatrix")),
		2454
	)

	expect_equal(
		ArchR:::.availableSeqnames(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		ArchR:::.availableChr(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		as.vector(table(substr(ArchR:::.availableCells(arrowFiles), 9, 9))[c("B", "M", "T")]),
		c(45, 33, 49)
	)

	expect_equal(
		paste0(ArchR:::.sampleName(arrowFiles)),
		"PBSmall"
	)

})

################################################
# Testing Dropping Matrices
################################################

arrowFiles <- ArchR:::.dropGroupsFromArrow(
	ArrowFile = arrowFiles,
	dropGroups = c("GeneScoreMatrix", "TileMatrix")
)

#Get Contents
test_that("Checking Arrow Contents After Drop...", {

	expect_equal(
		ArchR:::.validArrow(arrowFiles),
		arrowFiles
	)

	expect_equal(
		paste0(ArchR:::.availableArrays(arrowFiles)),
		c("")[-1]
	)

	expect_equal(
		ArchR:::.availableSeqnames(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		ArchR:::.availableChr(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		as.vector(table(substr(ArchR:::.availableCells(arrowFiles), 9, 9))[c("B", "M", "T")]),
		c(45, 33, 49)
	)

	expect_equal(
		paste0(ArchR:::.sampleName(arrowFiles)),
		"PBSmall"
	)

})

################################################
# Adding Tile Matrix
################################################

arrowFiles <- addTileMatrix(
	input = arrowFiles,
	tileSize = 25000,
	chromSizes = getChromSizes(),
	force = TRUE
)

#Get Contents
test_that("Checking Arrow Contents after addTileMatrix", {

	expect_equal(
		ArchR:::.validArrow(arrowFiles),
		arrowFiles
	)

	expect_equal(
		paste0(ArchR:::.availableArrays(arrowFiles)),
		"TileMatrix"
	)

	expect_equal(
		nrow(ArchR:::.getFeatureDF(arrowFiles, "TileMatrix")),
		12638
	)

	expect_equal(
		ArchR:::.availableSeqnames(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		ArchR:::.availableChr(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		as.vector(table(substr(ArchR:::.availableCells(arrowFiles), 9, 9))[c("B", "M", "T")]),
		c(45, 33, 49)
	)

	expect_equal(
		paste0(ArchR:::.sampleName(arrowFiles)),
		"PBSmall"
	)

})


################################################
# Adding Gene Score Matrix
################################################

arrowFiles <- addGeneScoreMatrix(
	input = arrowFiles,
	genes = getGenes()
)

#Get Contents
test_that("Checking Arrow Contents after addGeneScoreMatrix...", {

	expect_equal(
		ArchR:::.validArrow(arrowFiles),
		arrowFiles
	)

	expect_equal(
		ArchR:::.availableArrays(arrowFiles),
		c("GeneScoreMatrix", "TileMatrix")
	)

	expect_equal(
		nrow(ArchR:::.getFeatureDF(arrowFiles, "TileMatrix")),
		12638
	)

	expect_equal(
		nrow(ArchR:::.getFeatureDF(arrowFiles, "GeneScoreMatrix")),
		2454
	)

	expect_equal(
		ArchR:::.availableSeqnames(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		ArchR:::.availableChr(arrowFiles),
		c("chr11", "chr5")
	)

	expect_equal(
		as.vector(table(substr(ArchR:::.availableCells(arrowFiles), 9, 9))[c("B", "M", "T")]),
		c(45, 33, 49)
	)

	expect_equal(
		paste0(ArchR:::.sampleName(arrowFiles)),
		"PBSmall"
	)

})

################################################
# Final Checks
################################################

test_that("Final Checks...", {

	expect_equal(
		length(getFragmentsFromArrow(arrowFiles)),
		26443
	)

	expect_equal(
		getMatrixFromArrow(arrowFiles, "GeneScoreMatrix") %>% {c(nrow(.), ncol(.))},
		c(2454, 127)
	)

	expect_equal(
		getMatrixFromArrow(arrowFiles, "TileMatrix", binarize = TRUE) %>% {c(nrow(.), ncol(.))},
		c(12638, 127)
	)

})

################################################
# Clear
################################################

files <- list.files()
files <- files[!grepl("\\.R", files)]
for(i in seq_along(files)){
	if(dir.exists(files[i])){
		unlink(files[i], recursive=TRUE)
	}else if(file.exists(files[i])){
		file.remove(files[i])
	}
}





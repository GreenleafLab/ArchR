# Tests for ArchR 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

#Context
context("test_workflow")

#Genome
addArchRGenome("hg19test2")

#Arrows
arrows <- file.path(system.file("testdata", package="ArchR"), "PBSmall.arrow")

#Project
proj <- ArchRProject(arrows, outputDirectory = "Test")

#LSI
proj <- addIterativeLSI(proj, dimsToUse = 1:5, varFeatures=1000, iterations = 2, force=TRUE)

#Clusters
proj <- addClusters(proj, force=TRUE, dimsToUse = 1:5)

#Cell type
proj$CellType <- substr(getCellNames(proj), 9, 9)

#Check Clusters
cM <- confusionMatrix(proj$CellType, proj$Clusters)

test_that("Test Cluster Purity...", {
	expect_equal(all(apply(cM / Matrix::rowSums(cM), 1, max) > 0.9), TRUE)
})

#UMAP
proj <- addUMAP(proj, nNeighbors = 40, dimsToUse = 1:5, minDist=0.1, force=TRUE)

#Plot UMAP
p1 <- plotEmbedding(proj, name = "CellType", size = 3)
p2 <- plotEmbedding(proj, name = "Clusters", size = 3)
plotPDF(p1, p2, name = "UMAP", addDOC = FALSE, ArchRProj = proj)
 
#Group Coverages
proj <- addGroupCoverages(proj)

#Custom Peak Calling
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    peakMethod = "tiles",
    minCells = 20,
    cutOff = 0.1
)

#Macs2 Peak Calling
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2,
    minCells = 20,
    cutOff = 0.1
)

#Add Peak Matrix
proj <- addPeakMatrix(proj)

#Motif Annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbpTest", name = "Motif", force=TRUE)

#Motif Deviations
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

#Check Matrices
se <- getMatrixFromProject(proj, "MotifMatrix")

pval_CEBP <- t.test(
	assays(se)$z["CEBPB_1", se$CellType=="T"],
	assays(se)$z["CEBPB_1", se$CellType=="M"]
)$p.value

pval_EOMES <- t.test(
	assays(se)$z["EOMES_6", se$CellType=="T"],
	assays(se)$z["EOMES_6", se$CellType=="M"]
)$p.value

pval_PAX <- t.test(
	assays(se)$z["PAX5_5", se$CellType=="M"],
	assays(se)$z["PAX5_5", se$CellType=="B"]
)$p.value

test_that("Test Motif Deviations...", {
	expect_equal(pval_CEBP < 0.01, TRUE)
	expect_equal(pval_EOMES < 0.01, TRUE)
	expect_equal(pval_PAX < 0.01, TRUE)
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


# Tests for functions in src directory 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
library(chromVAR)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

#Context
context("test_chromVAR")

#Test Project
proj <- getTestProject()

#Matrix
se <- getMatrixFromProject(proj, "PeakMatrix")
names(assays(se)) <- "counts"

#Bias
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
rowData(se)$bias <- round(rowData(se)$bias, 4)

#Motif Matches
matches <- getMatches(proj)

#Check
stopifnot(all(paste0(rowRanges(se))==paste0(rowRanges(matches))))

#Background Peaks
se2 <- SummarizedExperiment(assays=SimpleList(counts=cbind(Matrix::rowSums(assay(se)), 1)))
rowRanges(se2) <- rowRanges(se)
bgdPeaks <- getBackgroundPeaks(se2)

#Computing deviations
dev <- computeDeviations(
  object = se, 
  annotations = matches,
  background_peaks = bgdPeaks
)

#Background Peaks ArchR
bgdPeaks2 <- getBgdPeaks(proj, force = TRUE)

#Computing deviations
dev2 <- computeDeviations(
  object = se, 
  annotations = matches,
  background_peaks = assay(bgdPeaks2)
)

#Compute Deviations ArchR
proj <- ArchR::addDeviationsMatrix(proj, force=TRUE, bgdPeaks = bgdPeaks2)

#Compare
dev3 <- getMatrixFromProject(proj, "MotifMatrix")

#Check Deviations
corDev_12 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assay(dev2)[x,], assay(dev)[x,])
}) %>% unlist

corZ_12 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assays(dev2)[[2]][x,], assays(dev)[[2]][x,])
}) %>% unlist

corDev_13 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assay(dev3)[x,], assay(dev)[x,])
}) %>% unlist

corZ_13 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assays(dev3)[[2]][x,], assays(dev)[[2]][x,])
}) %>% unlist

corDev_23 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assay(dev2)[x,], assay(dev3)[x,])
}) %>% unlist

corZ_23 <- lapply(seq_len(nrow(dev)), function(x){
  cor(assays(dev2)[[2]][x,], assays(dev3)[[2]][x,])
}) %>% unlist

#Test
test_that("Testing ArchR ChromVAR", {
  expect_equal(all(corDev_12 > 0.95), TRUE)
  expect_equal(all(corZ_12 > 0.95), TRUE)
  expect_equal(all(corDev_13 > 0.95), TRUE)
  expect_equal(all(corZ_13 > 0.95), TRUE)
  expect_equal(all(corDev_23 > 0.99), TRUE)
  expect_equal(all(corZ_23 > 0.99), TRUE)
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

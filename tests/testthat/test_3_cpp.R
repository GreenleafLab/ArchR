# Tests for functions in src directory 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

#Context
context("test_cpp")

#Default Matrices
data("m1")
m2 <- m1[rev(1:10), rev(1:10)]

################################################
# Testing Row Correlation
################################################

#Correlations
c1 <- ArchR:::rowCorCpp(1:10, 1:10, m1, m2)
c2 <- lapply(1:10, function(x){
	cor(m1[x, ], m2[x, ])
}) %>% unlist

#Test
test_that("Row-wise Correlation is working...", {
	expect_equal(c1, c2)
})

################################################
# Testing KNN Utils
################################################

#KNN
knnObj <- ArchR:::.computeKNN(m1, m2, k = 5)

#Check Knn Overlap Cpp
overlapCpp <- ArchR:::determineOverlapCpp(knnObj, 3)

#Check Knn Overlap R
overlapR <- lapply(seq_len(nrow(knnObj)), function(x){
	o <- lapply(seq_len(x-1), function(y){
		sum(knnObj[x, ] %in% knnObj[y, ])
	}) %>% unlist
	if(any(o > 3)){
		-1
	}else{
		0
	}
}) %>% unlist

#Test
test_that("KNN Utils is working...", {
	expect_equal(overlapCpp, overlapR)
})

################################################
# Testing General Utils
################################################

#tabulate2dCpp
tab2d <- as.matrix(ArchR:::tabulate2dCpp(
	x = c(0,0,2,2,3),
	xmin = 0,
	xmax = 3,
	y = c(1,1,2,3,3),
	ymin = 1,
	ymax = 3
))

tabSm <- as.matrix(Matrix::sparseMatrix(
	i = c(1,1,2,3,3),
	j = c(0,0,2,2,3) + 1,
	x = c(1,1,1,1,1)
))

#Test
test_that("Tabulate Utils is working...", {
	expect_equal(max(tab2d-tabSm) <= testthat_tolerance(), TRUE)
})

#rowSparseVariances
sm1 <- as(m1, "dgCMatrix")

#Variances
var1 <- ArchR:::computeSparseRowVariances(
	sm1@i + 1, sm1@x, Matrix::rowMeans(sm1), ncol(sm1))

var2 <- apply(m1, 1, var)

test_that("Variance Utils is working...", {
	expect_equal(max(var1-var2) <= testthat_tolerance(), TRUE)
})





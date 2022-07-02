# Tests for ArchR 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

clear_test_dir <- function(){
	files <- list.files()
	files <- files[!grepl("\\.R", files)]
	for(i in seq_along(files)){
		if(dir.exists(files[i])){
			unlink(files[i], recursive=TRUE)
		}else if(file.exists(files[i])){
			file.remove(files[i])
		}
	}
}

################################################
# Test Rd
################################################

fn <- list.files(system.file("man", package="ArchR"), full.names=TRUE)
fn <- grep("\\.Rd", fn, value=TRUE)
fn <- c(
	fn[grepl("createArrow", fn)],
	fn[!grepl("createArrow", fn)]
)
for(i in seq_along(fn)){
	test_example(path = fn[i], basename(fn[i]))
	clear_test_dir()
}
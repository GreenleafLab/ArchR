<p align="center"><a href ="https://www.archrproject.com"><img src="Figures/ArchR_Logo_Integrated.png" alt="" width="350"></a></p>
<hr>

[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

### ArchR has new features available for scATAC-seq Analysis

**Paired scATAC-seq and scRNA-seq Analysis**

ArchR now supports paired scATAC-seq and scRNA-seq Analysis! <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;See updates with importFeatureMatrix, addGeneExpressionMatrix, addIterativeLSI, addCombinedDims <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For a brief tutorial of these features : https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html

**Trajectory Analysis**

ArchR now directly supports both monocle3 and Slingshot based trajectory analysis! <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;See updates with getMonocleTrajectories, addMonocleTrajectory, addSlingShotTrajectories <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For a brief tutorial of these features : https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Trajectory.html

Additionally ArchR now enables export of a peak matrix that is compatible with STREAM!<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;See updates with exportPeakMatrixForSTREAM <br />

### ArchR is currently in Beta and will be in active development through the peer review process.

ArchR is a full-featured R package for processing and analyzing single-cell ATAC-seq data. ArchR provides the most extensive suite of scATAC-seq analysis tools of any software available. Additionally, ArchR excels in both speed and resource usage, making it possible to analyze 1 million cells in 8 hours on a MacBook Pro laptop.

### For installation instructions and full documentation, visit www.ArchRProject.com.

<hr>

![](Figures/ArchR_Workflow_Horizontal.png)

# Quick Installation of ArchR
For a full walk through of installation and frequently related issues please visit www.ArchRProject.com.

**First, install devtools (for installing GitHub packages) if it isn't already installed:**
``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

**Then, install BiocManager (for installing bioconductor packages) if it isn't already installed:**
``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

**Then, install ArchR:**
``` r
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```

**Lastly, install all of the ArchR dependencies that aren't installed by default:**
``` r
library(ArchR)
ArchR::installExtraPackages()
```
If any of these steps fails, you should identify the offending package and troubleshoot that individual installation before proceeding. Additionally, please see the ArchR website (www.ArchRProject.com) where we have installation troubleshooting tips.

# Pre-compiled ArchR environment
We provide two methods in which a user can manage R dependencies.  

### Using renv to manage dependencies
The first is by using renv to manage a project's dependencies. To utilize this, make sure that the renv package is installed and loaded.  Before you are ready to use `renv`, you must ensure that you are working on the same R version that we used for the provided renv environment. 
The R versions we currently support are:
```
- R 4.4
- R 4.1
```
Secondly, make sure that the renv package is installed and loaded.
```
install.packages("renv")
library(renv)
```
Set a working directory for your project
```
dir.create(path = "./<project_name>", showWarnings = FALSE)
setwd("./<project_name>")
```
Then, lets download the lock file for the current master branch of ArchR.
For R 4.4:
```
download.file(url = "https://pub-9ae435458ecc412abbbc9420a502ec38.r2.dev/renv.lock", destfile = "./renv.lock")
```
For R 4.1:
```
download.file(url = "https://pub-9ae435458ecc412abbbc9420a502ec38.r2.dev/renv_4_1.lock", destfile = "./renv.lock")
```

Now, we can initiate our renv project environment, utilizing the renv.lock to bootstrap a new renv environment.
```
renv::init()
```

### Using Docker to manage dependencies
We also provide Docker images, built off of `rocker/rstudio`, that already have ArchR and all dependencies pre-loaded.

The latest version can be found at:
```
greenleaflab/archr:latest
```
and other versions, including images built with differing R versions, can be found at:
```
https://hub.docker.com/r/greenleaflab/archr/tags
```

To utilize these images, the user can first install Docker as mentioned in their [documentation](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository)

Following, create a container using the following command:
```
docker image pull greenleaflab/archr:latest
docker run -it --rm -v <your_workspace>:/workspace -p <your_port_of_interest>:8787
```
This will spin up a container that has Rstudio turned on by default. Rstudio can be accessed through:
```
localhost:<your_port_of_interest>
```
If you would like an interactive bash console instead, the following command can instead be called:
```
docker run -it --rm -v <your_workspace>:/workspace -p <your_port_of_interest>:8787 bash
```

# Issues using ArchR?

ArchR is currently in __beta__. We expect there to be bumps in the road. If you think you have found a bug, please first install the latest version of ArchR via
``` r
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```
If this does not fix your problem, please [report an issue on Github](https://github.com/GreenleafLab/ArchR/issues) with the __Bug Report__ form.

If you have questions about ArchR usage, please refer to [the searchable full user's manual](https://www.archrproject.com/bookdown/index.html), [the FAQ section](https://www.archrproject.com/articles/Articles/faq.html), and the [publication](https://greenleaf.stanford.edu/assets/pdf/). If you think the documentation on this website or in the function annotations is unclear, please [submit an issue on Github](https://github.com/GreenleafLab/ArchR/issues) with the __Documentation Request__ form. If there is a feature that you think is missing from ArchR _and you have already searched the user's manual_, [submit an issue on Github](https://github.com/GreenleafLab/ArchR/issues) with the __Feature Request__ form. If none of these options help, [send us an email](mailto:archr.devs@gmail.com). We will do our best to respond to questions that are not otherwise answered in the documentation.



<p align="center"><a href ="https://www.archrproject.com"><img src="Figures/ArchR_Logo_Integrated.png" alt="" width="350"></a></p>
<hr>

[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

### ArchR is currently in Beta and will be in active development through the peer review process.

ArchR is a full-featured R package for processing and analyzing single-cell ATAC-seq data. ArchR provides the most extensive suite of scATAC-seq analysis tools of any software available. Additionally, ArchR excels in both speed and resource usage, making it possible to analyze 1 million cells in 8 hours on a MacBook Pro laptop.

### All documentation is available at www.ArchRProject.com.

<hr>

![](Figures/ArchR_Workflow_Horizontal.png)

# Quick Installation of ArchR
For a full walk through of installation and frequently related issues please visit www.ArchRProject.com.

**First, install devtools (for installing GitHub packages) if it isn't already installed:**
```{r}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

**Then, install BiocManager (for installing bioconductor packages) if it isn't already installed:**
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

**Then, install ArchR:**
```{r}
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```

**Lastly, install all of the ArchR dependencies that arent installed by default:**
```{r}
library(ArchR)
ArchR::installExtraPackages()
```
If any of these steps fails, you should identify the offending package and troubleshoot that individual installation before proceeding. Additionally, please see the ArchR website (www.ArchRProject.com) where we have installation troubleshooting tips.

# Issues using ArchR?
If you find a bug, please report it as an issue on [Github](https://github.com/GreenleafLab/ArchR/issues). If you think the documentation on this website or in the function annotations is unclear, please submit this as an issue on [Github](https://github.com/GreenleafLab/ArchR/issues) with the __documentation__ tag. If you have questions about ArchR usage, please refer to the [the searchable full user's manual](https://www.archrproject.com/bookdown/index.html), [the FAQ section](https://www.archrproject.com/articles/faq.html), and the publication. If none of these options help, [send us an email](mailto:archr.devs@gmail.com). We will do our best to respond to questions that are not otherwise answered in the documentation.



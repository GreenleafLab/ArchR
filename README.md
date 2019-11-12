# ArchR ![](man/figures/ArchR_dartLogo_small.jpg)

ArchR is a full-featured R package for processing and analyzing single-cell ATAC-seq data. Its strength is speed and resource usage, making it possible to analyze 1 million cells in QQQ hours on a macbook pro laptop. It provides facile tools to do just about anything you would want to do with single-cell ATAC-seq data. For a more detailed description of the software, please see the [publication](https://greenleaf.stanford.edu/assets/pdf/) ([pdf](http://greenleaf.stanford.edu/assets/pdf/), [supplement](http://greenleaf.stanford.edu/assets/pdf/)) or the vignettes on this page.

# Installation of ArchR

ArchR installation currently requires devtools. The following commands will use the Bioconductor BiocManager to install required dependences:

```{r}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

devtools::install_github("GreenleafLab/ArchR",
	auth_token = token, #Need a token to download (see personalized access tokens)
	repos = BiocManager::repositories()
)
```

### Additional Packages that are used from github
To complete installation, you also must maually install the following packages using these devtools commands:

```{r}
# ggrastr is a package for plotting ggplots with rastr'd points which is super helpful for large UMAPs etc
# # You need to have Cairo installed for ggrastr to work
# # On Mac OSx you need to have XQuartz (https://www.xquartz.org/)
devtools::install_github('VPetukhov/ggrastr')

# harmony is a package that can correct batch effects
devtools::install_github("immunogenomics/harmony")

# presto is a package that has efficient tools for wilcoxon tests on sparseMatrices
devtools::install_github('immunogenomics/presto')
```

# ArchR Workflow
![ArchR Workflow](https://web.stanford.edu/~mcorces/ArchR_Workflow_v1.PNG)

# Documentation
Please see the navigation bar at the top of this page for links to background information about ArchR as well as detailed vignettes for each of the major ArchR analytical components.




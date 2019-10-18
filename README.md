# Installation of ArchR

```r

if (!requireNamespace("devtools", quietly = TRUE))
	install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

devtools::install_github("jgranja24/ArchR",
	auth_token = token, #please email me for a token
	repos = BiocManager::repositories()
)

```

# Additional Packages that are used from github


```r

# ggrastr is a package for plotting ggplots with rastr'd points which
# is super helpful for large UMAPs etc
# 
# You need to have Cairo for ggrastr
#
# On Mac OSx you need to have XQuartz (https://www.xquartz.org/)
#
devtools::install_github('VPetukhov/ggrastr')

# harmony is a package that can correct batch effects
devtools::install_github("immunogenomics/harmony")

# presto is a package that has efficient tools for wilcoxon tests on sparseMatrices
devtools::install_github('immunogenomics/presto')


```


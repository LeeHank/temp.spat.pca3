
<!-- README.md is generated from README.Rmd. Please edit that file -->
temp.spat.pca3
==============

The goal of temp.spat.pca3 is to implement the functions mentioned in the thesis: Regularized PCA on dynamic spatial pattern.

Installation
------------

You can install the released version of temp.spat.pca3 from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("temp.spat.pca3")
```

You can also install the development version of temp.spat.pca3 by:

``` r
#install.packages("devtools")
library(devtools)
install_github("LeeHank/temp.spat.pca3", build_vignettes = TRUE)
```

Tutorials
---------

After downloading the package, one can use following code to run tutorial:

``` r
vignette("simulation_example", package="temp.spat.pca3")
vignette("real_data_example", package="temp.spat.pca3")
```

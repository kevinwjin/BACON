# BACON: Bayesian Clustering of n-gons via a Double Dirichlet Mixture Model

## Introduction
`BACON` is an R package for landmark-based Bayesian clustering of closed polygonal chains, relying on intrinsic shape features (i.e. relative interior angles and relative side lengths). The algorithm accepts a matrix, list, or dataframe containing the coordinates of closed polygonal chains, extracts the aforementioned inherent shape features, and clusters them by approximating posterior distributions using a Gibbs sampler.

## Installation
Install BACON by one of the following methods, depending on your preference:

```R
# Install from The Comprehensive R Archive Network (CRAN)
install.packages("BACON")

# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools:install_github("kevinwjin/BACON")
```
## Usage
*(under development)*

## Arguments
*(under development)*

## Contents
* `code/simulate_shapes.R` - Data simulation
* `code/shape_simu_analysis.R` - MCMC algorithm execution
* `code/data_analysis.R` - Testing alternative clustering methods
* `code/functions/shape_generation.R` - Functions for data simulation
* `code/functions/shape_mcmc.cpp` - MCMC algorithm
* `data/` - Real-world datasets

## Prerequisites
* `sf` - Spatial point detection for interior angle calculation
* `MCMCpack` - For the `rdirichlet()` function
* `mcclust` - MCMC clustering sample processing
* `mclust` - GMM implementation for comparison purposes
* `Rcpp` - Faster MCMC approximation
* `RcppArmadillo` - Fast matrix operations (requires GNU Fortran library)
* `RcppDist` - Call statistical distributions from within C++
* `RcppEigen` - Faster matrix operations

## Examples
*(under development)*

## Contributors
* [Kevin Jin](https://www.linkedin.com/in/kevin-w-jin/)
* [Huimin Li](https://www.linkedin.com/in/huimin-li-19789248)
* [Stephen McKeown](https://personal.utdallas.edu/~sxm190098/)
* [Qiwei Li](https://sites.google.com/site/liqiwei2000/)

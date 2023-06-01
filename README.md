# BACON: Bayesian Clustering of n-gons via a Double Dirichlet Mixture Model

## Introduction
`BACON` is an R package for landmark-based Bayesian clustering of closed polygonal chains, relying on intrinsic shape features (i.e. proportions of interior angles and side lengths). The algorithm accepts the coordinates of closed polygonal chains, extracts the aforementioned inherent shape features, and clusters them using a Gibbs sampler.

## Installation (for future use)
Install the package by one of the following methods:

```R
# Install from The Comprehensive R Archive Network (CRAN)
install.packages("BACON")

# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools:install_github("kevinwjin/BACON")
```

## Usage
### Simulate shape data
To simulate shape data for clustering, begin by adding the relevant functions to
your R environment. The following file contains all functions necessary for shape 
data simulation:

```R
source("~/BACON/code/data_simulation/functions.R")
```

Next, generate some simulated shape data. Our data simulation code generates 
the Cartesian coordinates of the vertices of a shape. As an example, 
we will simulate 1000 20-gons belonging to 10 evenly-spaced 
clusters. Each cluster will have 100 20-gons within it:

```R
dataset <- simulate_shapes(x = 1000, # 1000 shapes total
                          z = 10, # 10 clusters
                          n = 100, # 100 shapes per cluster
                          k = 20, # Each shape is a 20-gon
                          jitter_factor = 0.01) # Each cluster has 0.01 jitter
```

As a sanity check, you may visualize the simulated shapes by plotting the 
coordinates of each shape in the dataset:

```R
for (i in seq_along(dataset)) {
  for (shape in dataset[[i]]) {
    plot(shape, type = "l")
  }
}
```

We cannot cluster the raw coordinate data, as it is sensitive to geometric
transformations; therefore, we must convert the raw coordinate data to 
normalized compositional data by extracting the intrinsic
transformation-invariant shape features (interior 
angles and side lengths, both normalized to 1) that our model will be able to 
cluster. To do this, call `get_interior_angles()` and `get_side_lengths()`, and
store the resultant interior angle and side length proportions into new matrices:

```R
angles <- matrix(nrow = x, ncol = k, byrow = TRUE)
side_lengths <- matrix(nrow = x, ncol = k, byrow = TRUE)

counter <- 1
for (i in seq_along(dataset)) {
  for (j in dataset[[i]]) {
    angles[counter, ] <- get_interior_angles(j)
    side_lengths[counter, ] <- get_side_lengths(j)
    counter <- counter + 1
  }
}

# Clean up variables
rm(i, j)
```

Finally, we finish the simulated dataset by generating the ground truth 
containing the cluster labels:

```R
ground_truth <- rep(1:z, each = n)
```

A template for the above procedure is provided in `data_simulation_demo.Rmd`.

### Cluster shape data
To use BACON, begin by sourcing the following file into your R environment.

```R
# Source the clustering function into R environment
setwd("~/Documents/Repositories/BACON/code/clustering/")
source("bacon.R")
```
Call the clustering function `bacon()`, passing it the following arguments:

* `side_lengths`: n x k matrix containing n samples of k-gon side length proportions (required)
* `angles`: n x k matrix containing n samples of k-gon angle proportions (required)
* `K`: Number of clusters *a priori* (required)
* `weight_L`: Numerical weight of the contribution of the side length proportions to the mixture model (tuning parameter [0, 1] differing between datasets; default is 1)
* `weight_A`: Numerical weight of the contribution of the angle proportions to the mixture model (tuning parameter [0, 1] differing between datasets; default is 1)
* `estimate.s`: Whether to estimate the `s` parameter for shape registration (default is `TRUE`)
* `estimate.r`: Whether to estimate the `r` parameter for shape registration (default is `TRUE`)
* `iter`: Number of MCMC iterations (default is 2000)
* `burn`: Number of MCMC iterations for burn-in (default is `iter`/2, or 1000)

```R
# Execute the clustering function and save results
res <- bacon(side_lengths, angles, K, weight_L, weight_A, estimate.s, estimate.r, iter, burn)
```

Finally, evaluate the clustering accuracy by calculating the ARI between the estimated clusters and the ground truth. We will do this by calling the the `adjustedRandIndex()` function from the `mclust` package.

```R
# Return the ARI of BACON-derived clusters and the ground truth
mclust::adjustedRandIndex(res$cluster, ground_truth)
```

## Directory contents
* `code/clustering/` - Code for shape clustering
* `code/data_simulation/` - Code for simulating data for model testing
* `code/model_testing/` - Code for model testing on data
* `data/` - Simulated and real datasets
* `figures/` - Figures from model testing

## Prerequisites
* `Rcpp` - C++ implementation of the Markov chain Monte Carlo (MCMC) algorithm
* `RcppArmadillo` - Fast matrix operations (requires GNU Fortran library)
* `RcppEigen` - Fast matrix operations
* `RcppDist` - Call probability distributions from within C++
* `MCMCpack` - Required by `rdirichlet()` function
* `mcclust` - MCMC clustering sample processing
* `sf` - Spatial point detection for interior angle calculation
* `mclust` - Implementation of Gaussian mixture model for comparison of clustering methods

## Collaborators
* [Kevin Jin](https://kevinwjin.com/)
* [Huimin Li](https://www.linkedin.com/in/huimin-li-19789248)
* [Bryn Brakefield](https://github.com/brakefieb)
* [Stephen McKeown](https://personal.utdallas.edu/~sxm190098/)
* [Qiwei Li](https://sites.google.com/site/liqiwei2000/)

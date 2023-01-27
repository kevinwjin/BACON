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
To simulate shape data for clustering, begin by adding `shape_simulation.R` to 
your R environment. This file contains all functions necessary for shape 
data simulation:

```R
source("~/BACON/code/shape_simulation.R")
```

Next, generate some simulated shape data. We will use the `generate()` function,
which simulates the Cartesian coordinates of the vertices of a shape.
As an example, we will simulate 1000 20-gons belonging to 10 clusters. Each
cluster will have 100 20-gons.

We define a function that will create a list of sublists where each sublist is
a cluster and each element within each sublist is the a (k + 1) x n matrix 
containing the Cartesian coordinates of a shape. Each row will have some jitter 
(the degree, or factor, of jitter can be specified) applied to it to 
distinguish the individual shapes within a cluster:

```R
simulate_shapes <- function(x, z, n, k, jitter_factor) {
  # Entire dataset, a list of length z containing z sublists/clusters
  dataset <- list()
  # Each sublist is a cluster containing n cluster members
  for (cluster in 1:z) {
    # Generate random k-gon
    shape <- generate(k = k)
    # Replicate k-gon n times
    shapes <- replicate(n, shape, simplify = FALSE)
    # Apply jitter to each k-gon
    for (i in seq_along(shapes)) {
      shapes[[i]] <- 
        jitter(shapes[[i]], random = c("vertices"), factor = jitter_factor)
    }
    # Add to main shape dataset
    dataset[[length(dataset) + 1]] <- shapes
  }
  return(dataset)
}
```

Now we call the function and generate 1000 20-gons, storing
the simulated shape coordinate data into a list:

```R
dataset <- simulate_shapes(x = 1000, # 1000 shapes total
                          z = 10, # 10 clusters
                          n = 100, # 100 shapes per cluster (could differ per cluster)
                          k = 20, # Each shape is a 20-gon
                          jitter_factor = 0.01) # Each cluster has 0.01 jitter
```

For a sanity check, you may visualize the simulated 20-gons by plotting the 
coordinates of each shape:

```R
for (i in seq_along(dataset)) {
  for (shape in dataset[[i]]) {
    plot(shape, type = "l")
  }
}
```

We cannot cluster the raw coordinate data, as it is sensitive to geometric
transformations; therefore, we must convert the 
raw coordinate data to compositional data by extracting the intrinsic
transformation-invariant shape features (interior 
angles and side lengths, both normalized to 1) that our model will be able to 
cluster. To do this, call `get_interior_angles()` and `get_side_lengths()`, and
store the resultant interior angle and side length datasets into new matrices:

```R
angles <- list()
side_lengths <- list()

for (i in seq_along(dataset)) { # For all clusters
  a <- list()
  l <- list()
  i <- 1 # Loop counter
  for (shape in dataset[[i]]) { # Within one cluster
    a[[i]] <- get_interior_angles(shape) # Extract interior angle proportions
    l[[i]] <- get_side_lengths(shape) # Extract side length proportions
    i <- i + 1
  }
  angles[[length(angles) + 1]] <- a # Add to main angle proportion dataset
  side_lengths[[length(side_lengths) + 1]] <- l # Add to main side length dataset
}
```

Finally, we finish the simulated data by generating a matrix containing the
cluster labels:

```R
clusters <- c()
for (i in 1:z) {
  clusters <- append(clusters, rep.int(i, times = (n / z)))
}
```

A template for the above procedure of simulating 1000 20-gons is
provided in `shape_simulation_example.R`.

### Cluster shape data
*(under development)*

## Directory contents
* `code/clustering/bacon.R` - Clustering example
* `code/clustering/bacon.R` - Clustering execution function
* `code/clustering/BACONmcmc.cpp` - MCMC clustering algorithm
* `code/data_simulation/shape_simulation.R` - Data simulation functions
* `code/data_simulation/shape_simulation_example.R` - Data simulation template
* `data/demo.RData` - Demo dataset

## Prerequisites
* `sf` - Spatial point detection for interior angle calculation
* `MCMCpack` - For the `rdirichlet()` function
* `mcclust` - MCMC clustering sample processing
* `mclust` - GMM implementation for comparison purposes
* `Rcpp` - Fast MCMC operations
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

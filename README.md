# BACON: Bayesian Clustering of n-gons via a Double Dirichlet Mixture Model

## Introduction
`BACON` is an R package for landmark-based Bayesian clustering of closed polygonal chains, relying on intrinsic shape features (i.e. relative interior angles and relative side lengths). The algorithm accepts a matrix, list, or dataframe containing the coordinates of closed polygonal chains, extracts the aforementioned inherent shape features, and clusters them using a Gibbs sampler.

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
To simulate shape data for clustering, begin by adding `shape_generation.R` to 
your R environment. This file contains all functions necessary for shape 
data simulation:

```R
source("~/BACON/code/shape_generation.R")
```

Next, generate some simulated shape data. We will use the `generate()` function,
which simulates the Cartesian coordinates of the vertices of a shape.
As an example, we will simulate two hundred 20-gons belonging to 10 clusters. 
Let's define a function that will create a 200 x 20 matrix where each row 
is a separate 20-gon and each column contains the Cartesian coordinates of a 
vertex in a given 20-gon. Each row will have some jitter (the degree, or factor, 
of jitter can be specified) applied to it to distinguish the individual shapes
and clusters:

```R
n <- 200 # Total number of shapes
k <- 20 # Number of vertices
z <- 10 # Pre-specified number of clusters

simulate_shapes <- function(f) {
  # Generate random shape with k vertices
  shape <- generate(k = k)
  
  # Replicate shape n / z times
  shapes <- do.call(rbind, replicate(n / z, shape, simplify = FALSE))
  
  # Count by k + 1 because last vertex repeats due to closedness
  for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  
    # Apply jitter to each chain
    shapes[i:(i + k), ] <- jitter(shapes[i:(i + k), ], random = c("vertices"),
                                    factor = 0.01)
  }
  return(shapes)
}
```

Now we run the function and generate two hundred 20-gons, storing
the simulated shape coordinate data into one three-dimensional array:

```R
# Simulate z clusters of shapes into 3-D array [(k + 1) * z, xy, z]
# To access first shape: array[, , 1]
shapes <- sapply(1:z, simulate_shapes, simplify = "array")
```

If you wish, you may visualize the simulated 20-gons by plotting the 
coordinates of each shape:

```R
# Plot each shape
for (i in seq_len(z)) {
  for (j in seq(1, (n / z * (k + 1)), by = (k + 1))) {
    plot(shapes[j:(j + k), , i], type = "l")
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
# Extract relative interior angles of the shapes
angles <- matrix(nrow = n, ncol = k, byrow = TRUE)
row <- 1 # Row counter
for (i in seq_len(z)) { # Extract all clusters
  for (j in seq(1, n / z * (k + 1), by = (k + 1))) { # Extract all shapes per cluster
    print(j)
    angles[row, ] <- get_interior_angles(shapes[j:(j + k), , i])
    row <- row + 1
  }
}

# Extract relative side lengths of the shapes
side_lengths <- matrix(nrow = n, ncol = k, byrow = TRUE)
row <- 1 # Row counter
for (i in seq_len(z)) { # Extract all clusters
  for (j in seq(1, n / z * (k + 1), by = (k + 1))) { # Extract all shapes per cluster
    print(j)
    side_lengths[row, ] <- get_side_lengths(shapes[j:(j + k), , i])
    row <- row + 1
  }
}
```

Finally, we finish the simulated data by appending a column containing the
cluster labels to each new dataset:

```R
clusters <- c()
for (i in 1:z) {
  clusters <- append(clusters, rep.int(i, times = (n / z)))
}
angles <- cbind(angles, clusters) # Completed angle dataset
side_lengths <- cbind(side_lengths, clusters) # Completed side length dataset
```

A template for the above procedure of simulating two hundred 20-gons is
provided in `simulate_shapes.R`.

### Cluster shape data
*(under development)*

## Directory contents
* `code/functions/shape_generation.R` - Data simulation functions
* `code/simulate_shapes.R` - Data simulation template
* `code/functions/shape_mcmc.cpp` - MCMC clustering algorithm
* `code/shape_simu_analysis.R` - Clustering template
* `code/data_analysis.R` - Cluster validation and method comparison
* `data/` - Real-world datasets

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

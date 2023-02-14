# Shape data simulation for BACON
# Author: Kevin Jin

#### Load data simulation functions ####
source("~/Documents/Repositories/BACON/code/data_simulation/shape_simulation.R")

#### Generate simulated data ####
## Step 1: Simulate z clusters of a total of x k-gons, each cluster having
## n k-gons each with modifiable jitter

# x: Total number of shapes in the dataset
# z: Number of clusters in the shape dataset
# n: Number of shapes per cluster
# k: Number of vertices in each shape
# jitter_factor: Amount of jitter to apply to each shape within a cluster
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
dataset <- simulate_shapes(x = 1000, # 1000 shapes total
                           z = 10, # 10 clusters
                           n = 100, # 100 shapes per cluster (could differ per cluster)
                           k = 20, # Each shape is a 20-gon
                           jitter_factor = 0.01) # Each cluster has 0.01 jitter

# Plot each shape as a sanity check
for (i in seq_along(dataset)) {
  for (shape in dataset[[i]]) {
    plot(shape, type = "l")
  }
}

## Step 2: Extract the shapes' normalized interior angles and side lengths
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

## Done! Your simulated data is stored in the lists 'angles' and 
## 'side_lengths'. Within each list, each sublist represents its own cluster.
## Within each cluster, each element represents a shape.

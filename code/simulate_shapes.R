# Shape data simulation for BACON
# Author: Kevin Jin

#### Load data simulation functions ####
source("~/Documents/Repositories/BACON/code/shape_generation.R")

#### Generate simulated data ####
## Step 1: Simulate 10 shapes 20 times each with jitter (200 total)
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
                                    factor = 0.01) # Jitter factor is sensitive
  }
  return(shapes)
}

# Simulate z clusters of shapes into 3-D array [(k + 1) * z, xy, z]
# To access first shape: array[, , 1]
shapes <- sapply(1:z, simulate_shapes, simplify = "array")

# Plot each shape
for (i in seq_len(z)) {
  for (j in seq(1, (n / z * (k + 1)), by = (k + 1))) {
    plot(shapes[j:(j + k), , i], type = "l")
  }
}

## Step 2: Extract the angle vectors of the shapes
angles <- matrix(nrow = n, ncol = k, byrow = TRUE)
row <- 1 # Row counter
for (i in seq_len(z)) { # Extract all clusters
  for (j in seq(1, n / z * (k + 1), by = (k + 1))) { # Extract all shapes per cluster
    print(j)
    angles[row, ] <- get_interior_angles(shapes[j:(j + k), , i])
    row <- row + 1
  }
}

## Step 3: Extract side lengths of the shapes
side_lengths <- matrix(nrow = n, ncol = k, byrow = TRUE)
row <- 1 # Row counter
for (i in seq_len(z)) { # Extract all clusters
  for (j in seq(1, n / z * (k + 1), by = (k + 1))) { # Extract all shapes per cluster
    print(j)
    side_lengths[row, ] <- get_side_lengths(shapes[j:(j + k), , i])
    row <- row + 1
  }
}

## Step 3: Data analysis
# Create and append cluster labels to data
clusters <- c()
for (i in 1:z) {
  clusters <- append(clusters, rep.int(i, times = (n / z)))
}
angles <- cbind(angles, clusters)
side_lengths <- cbind(side_lengths, clusters)

# Optional step: PCA dimensional reduction
require(ggfortify)
require(compositions)
angles <- cbind(angles, clusters)
angles <- as.data.frame(angles)
autoplot(stats::prcomp(clr(angles[, 1:k])), 
         data = angles,
         colour = factor(clusters),
         main = 'PCA of Simulated Shapes - Relative Interior Angles')
#biplot(stats::prcomp(angles[, 1:k]))
#plot(stats::prcomp(angles[, 1:k]), type = 'l')

side_lengths <- cbind(side_lengths, clusters)
side_lengths <- as.data.frame(side_lengths)
autoplot(stats::prcomp(clr(side_lengths[, 1:k])), 
         data = side_lengths, 
         colour = factor(clusters),
         main = 'PCA of Simulated Shapes - Relative Side Lengths')
#biplot(stats::prcomp(side_lengths[, 1:k]))
#plot(stats::prcomp(side_lengths[, 1:k], type = 'l'))

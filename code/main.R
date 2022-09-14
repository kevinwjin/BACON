# Execution script (data generation, model testing)
# Author: Kevin Jin

#### Preparation ####
set.seed(13579)
source("~/Documents/Programming/Repositories/bsclust/code/shape_generation.R")

#### Test shape generation ####
k <- 5
shape <- generate(k = k)
plot(shape, type = "l")
text(shape[1:nrow(shape), ], labels = 1:nrow(shape))
sum_interior_angles(shape)
sum(get_interior_angles(shape))
get_interior_angles(shape)
get_side_lengths(shape)
#polygon(jitter(shape, random = c("angles"), factor = 0.01))
polygon(translate(shape, x = -0.5, y = 1))
polygon(dilate(shape, factor = 0.5))
polygon(reflect(shape, direction = c("vertical")))
polygon(rotate(shape, angle = 180, clockwise = FALSE))

#### Generate data: Two different shapes, 30 each with slight variations ####
## Step 1: Replicate two shapes 5 times each and add jitter
n <- 60 # Total number of shapes
k <- 3 # Number of vertices
z <- 2 # Pre-specified number of clusters

# Generate first shape
shape_1 <- generate(k = k)
# Replicate n / 2 times (half of the total)
shapes_1 <- do.call(rbind, replicate(n / 2, shape_1, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / 2 * (k + 1)),
              by = (k + 1))) {
  shapes_1[i:(i + k), ] <- jitter(shapes_1[i:(i + k), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
# Generate second shape (z = 2 clusters)
shape_2 <- generate(k = k)
# Replicate n / 2 times (half of the total)
shapes_2 <- do.call(rbind, replicate(n / 2, shape_2, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / 2 * (k + 1)), 
              by = (k + 1))) {
  shapes_2[i:(i + k), ] <- jitter(shapes_2[i:(i + k), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
# Combine into 1 matrix
shapes <- rbind(shapes_1, shapes_2)

# Plot shapes - draw one shape per plot
for (i in seq(1, (n / 2 * (k + 1)), 
              by = (k + 1))) {
  plot(shapes[i:(i + k), ], type = "l")
}

## Step 2: Extract the angle vectors of the shapes (data to be clustered)
p <- k # Number of categories (length of angle vector)
x <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (angle vectors)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  x[row, ] <- get_interior_angles(shapes[i:(i + k), ])
  row <- row + 1
}

## Step 3: Extract side lengths of the shapes (data to be clustered)
s <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (side length vectors)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  s[row, ] <- get_side_lengths(shapes[i:(i + k), ])
  row <- row + 1
}

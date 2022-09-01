# Execution script (data generation, model testing)
# Author: Kevin Jin

## To-do:
## 1. Combine triangle generation within each cluster into one for loop

# Preparation
set.seed(13579)
source("~/Documents/bsclust/code/shape_generation.R")

# Test Rcpp functions
library(Rcpp)
sourceCpp("code/model.cpp")

# Generate test shapes
shape <- generate(k = 3)
plot(shape, xlim = c(-2, 2), ylim = c(-2, 2))
polygon(shape)

polygon(jitter(shape, random = c("vertices"), factor = 0.01)) # Factor should be small
#polygon(jitter(shape, random = c("angles"), factor = 0.01))
polygon(dilate(shape, factor = 0.5))
polygon(rotate(shape, angle = 90, clockwise = TRUE))
polygon(flip(shape, direction = c("horizontal")))
polygon(flip(shape, direction = c("vertical")))
polygon(rotate(shape, angle = 90, clockwise = FALSE))
polygon(translate(shape, x = -0.5, y = 1))

# Multinomial angle vector-based triangle generation
dmultinom(c(29, 61, 90), size = (sides - 2) * 180, prob = probs) # Density

sides <- 3
probs <- c(1 / 3, 1 / 3, 1 / 3)
triangles <- rmultinom(n = 100, # Number of polygons
                       size = (sides - 2) * 180, # Total internal angle
                       prob = probs) # Probability for K angles
t(triangles) # Generated vectors

# Generate 30 equilateral and right triangles each with slight variations
# Method 1: Generate 60 different triangles
n <- 60 # Number of samples
k <- 3 # Number of vertices
shapes <- matrix(nrow = n * k, ncol = 2, byrow = TRUE) # Generated shapes
for (i in seq(1, n * k, by = k)) { 
  shapes[i:(i + (k - 1)), ] <- generate(k = k)
}

# Method 2: Replicate two triangles 30 times each and add jitter
n <- 60 # Number of samples
k <- 3 # Number of vertices
z <- 2 # Pre-specified number of clusters

shape_1 <- generate(k = k) # Generate a triangle
shapes_1 <- do.call(rbind, replicate(n / 2, shape_1, simplify = FALSE)) # Replicate n / 2 times
for (i in seq(1, (n * k) / 2, by = k)) { # Apply jitter
  shapes_1[i:(i + (k - 1)), ] <- jitter(shapes_1[i:(i + (k - 1)), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
shape_2 <- generate(k = k) # Generate another triangle
shapes_2 <- do.call(rbind, replicate(n / 2, shape_2, simplify = FALSE)) # Replicate n / 2 times
for (i in seq(1, (n * k) / 2, by = k)) { # Apply jitter
  shapes_2[i:(i + (k - 1)), ] <- jitter(shapes_2[i:(i + (k - 1)), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
shapes <- rbind(shapes_1, shapes_2) # Combine into 1 matrix

# Draw one triangle per plot
for (i in seq(1, n * k, by = k)) {
  plot(NULL, xlim = c(-2, 2), ylim = c(-2, 2))
  polygon(shapes[i:(i + (k - 1)), ])
}

# Extract the angle vectors of the triangles (x, the data to be clustered)
p <- 3 # Number of categories (length of angle vector)
x <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (angle vectors)

x_row <- 1 # Independent row counter for X
for (i in seq(1, n * k, by = k)) {
  x[x_row, ] <- get_internal_angles(shapes[i:(i + (k - 1)), ])
  x_row <- x_row + 1
}
round(x) # Round to integers since we are treating the angle vector as a discrete variable

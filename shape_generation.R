# Random polygon generation
# Author: Kevin Jin

# Preparation
set.seed(13579)
source("~/Documents/bsclust/functions.R")

# Create empty canvas
plot.new()
#plot(1, 1, col = "white", xlab = "X", ylab = "Y") 

# Draw simple polygons
a <- data.frame(c(0, 0, 0.2, 0.2), c(0, 0.7, 0.7, 0))
polygon(a)

b <- data.frame(c(0.3, 0.4, 0.5), c(0.05, 0.4, 0.05))
polygon(b, lwd = 4) # Line weight

c <- data.frame(c(0, 0.1, 0.5, 0.5, 0.1), c(0.4, 0.5, 0.5, 0.3, 0.3))
polygon(c, col = "green") # Fill color

# Generate random data frames
for (shape in 1:3) {
  d <- matrix(rbeta(n = 2 * 5, # n = a dimensions * b vertices
                    shape1 = 1, 
                    shape2 = 5), ncol = 2) # Matrices are faster than data frames
  polygon(d)
}
 
# Generate regular polygons
plot.new()
e <- generate_polygon()
f <- generate_polygon(x_center = 1, y_center = 1)
g <- generate_polygon(N = 5, x_center = 1, y_center = 1)
plot(e)
plot(f)
plot(g)
polygon(g)

# Vary side length
for (r in sample(c(1:10))) {
  h <- generate_polygon(r = r, N = 3)
  plot(h)
  polygon(h)
}

# Vary rotation
for (theta in sample(seq(from = 30, to = 360, by = 30))) {
  h <- generate_polygon(N = 3, theta = theta)
  plot(h)
  polygon(h)
}

# Vary inner angles
for (n in sample(c(1:10))) {
  h <- generate_polygon(n = n, N = 3)
  polygon(h)
}

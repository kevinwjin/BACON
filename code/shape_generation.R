# Random polygon generation
# Author: Kevin Jin

# Preparation
set.seed(13579)
source("~/Documents/bsclust/functions.R")

# Test shape generation
shape <- generate(k = 3)

plot(shape, xlim = c(-2, 2), ylim = c(-2, 2))

hpts <- chull(shape)
hpts <- c(hpts, hpts[1]) # Close convex hull
lines(shape[hpts, ]) # Trace convex hull

lines(shape) # Trace polygon (alternate, less good)
trace(shape) # Trace polygon (alternate, less good)
polygon(shape) # Trace polygon (gold standard)

polygon(dilate(shape, factor = 0.5))
polygon(rotate(shape, angle = 90, clockwise = TRUE))
polygon(flip(shape, direction = c("horizontal")))
polygon(flip(shape, direction = c("vertical")))
polygon(rotate(shape, angle = 90, clockwise = FALSE))
polygon(translate(shape, x = -0.5, y = 1))

polygon(jitter(shape, random = c("vertices"), factor = 0.1))
polygon(jitter(shape, random = c("vertices"), factor = 0.2))
polygon(jitter(shape, random = c("vertices"), factor = 0.3))
polygon(jitter(shape, random = c("vertices"), factor = 0.4))

polygon(jitter(shape, random = c("angles")))

# Probabilistic triangle generation via multinomially distributed angle vectors
sides <- 3
probs <- c(1/3, 1/3, 1/3)
triangles <- rmultinom(n = 100, # Number of polygons
                       size = (sides - 2) * 180, # Total internal angle
                       prob = probs) # Probability for K angles
t(triangles)

dmultinom(c(29, 61, 90), size = (sides - 2) * 180, prob = probs)

# Test Rcpp functions
library(Rcpp)
sourceCpp("code/bsclust.cpp")

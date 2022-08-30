# Testing and execution script
# Author: Kevin Jin

# Preparation
set.seed(13579)
source("~/Documents/bsclust/code/shape_generation.R")

# Test shape generation
shape <- generate(k = 3)

plot(shape, xlim = c(-2, 2), ylim = c(-2, 2))
polygon(shape)

polygon(jitter(shape, random = c("vertices"), factor = 0.01)) # Should be small
#polygon(jitter(shape, random = c("angles")))
polygon(dilate(shape, factor = 0.5))
polygon(rotate(shape, angle = 90, clockwise = TRUE))
polygon(flip(shape, direction = c("horizontal")))
polygon(rotate(shape, angle = 90, clockwise = FALSE))
polygon(translate(shape, x = -0.5, y = 1))

# Probabilistic triangle generation from parameters (multinomial angle vectors)
sides <- 3
probs <- c(1 / 3, 1 / 3, 1 / 3)
triangles <- rmultinom(n = 100, # Number of polygons
                       size = (sides - 2) * 180, # Total internal angle
                       prob = probs) # Probability for K angles
t(triangles)

dmultinom(c(29, 61, 90), size = (sides - 2) * 180, prob = probs)

# Test Rcpp functions
library(Rcpp)
sourceCpp("code/bsclust.cpp")

# Random polygon generation
# Author: Kevin Jin

# Preparation
set.seed(13579)
source("~/Documents/bsclust/functions.R")

# Create empty canvas and test functions
plot.new()

pentagon <- generate(k = 5)

plot(pentagon, xlim = c(-2, 2), ylim = c(-2, 2))
polygon(pentagon)

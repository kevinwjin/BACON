# bsclust: Bayesian landmark-based shape clustering
# Authors: Kevin Jin, Qiwei Li

## Read data
## Preparation
set.seed(13579)
source("~/Documents/bsclust/code/shape_generation.R")

# Generate 30 equilateral and 30 right triangles with slight variations
n <- 60 # Number of samples
p <- 3 # Number of categories
k <- 2 # Pre-specified number of clusters
x <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (angle samples)
omega <- matrix(nrow = n, ncol = p, byrow = TRUE) # Parameters (angle proportions)
z <- matrix(nrow = n, ncol = 1) # Parameters (Clustering vector)
z_prob <- matrix(nrow = n, ncol = k) # Parameters (Clustering probabilities)

## Formulate data likelihood
l <- 0

for (i in seq_len(n)) {
  for (j in seq_len(p)) {
    factorial(prod(x[n, ]))
  }
} 

for (i in seq_len(n)) {
  for (j in seq_len(p)) {
    factorial(prod(x[n, ]))
  }
} 

## Formulate priors
omega_0 <- 
z_0 <-
  
## Formulate posteriors

## MCMC approximation of joint posterior via Gibbs sampling
# Initialize parameters
t <- 10000 # MCMC iterations
z_initial <- 1
omega_initial <- 1

z_store <- rep(NA, t) # Create space to store the results
omega_store <- rep(NA, t)
z <- z_initial
omega <- omega_initial

# Posterior approximation via Gibbs sampling
for (iter in seq_len(t)) {
  z_n <- 
  omega_n <- 

  z_store[iter] <- z
  omega_store[iter] <- omega
}

# Monitor convergence


## Posterior inference
# Visualize joint and marginal posteriors

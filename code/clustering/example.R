#### Run BACON on Huimin's demo data. ####
## The following guide demonstrates clustering of exmaple data using BACON.
## BACON requires three inputs:
## 1. L: a m-by-n matrix of length proportions, where m is the number of n-gons 
## and n is the number of gons.
## 2. A: a m-by-n matrix of angle proportions.
## 3. K: the number of clusters in the dataset.

## Load required packages
setwd("~/Documents/Repositories/BACON/code/clustering")
source("bacon.R")

## Load demo data (100 triangles with 3 clusters)
setwd("~/Documents/Repositories/BACON/data/simulated/Data")
load("demo.RData")
head(L)
head(A)

## Run the model (Default iterations = 2000; burn in = 1000)
res <- bacon(L, A, K = 3, iter = 10000, burn = 5000) # ARI = 0.536 ~ 0.731
res <- bacon(L2, A2, K = 3) # ARI = 1

## BACON produces the following output variables:
## cluster - the estimated cluster assignment.
## s_map - the estimated starting vertex indicators. 
## r_map - the estimated reverse indicators r.
res$cluster
res$s_map
res$r_map

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", 
     ylab = "Number of samples with starting point as 0", xlab = "Iteration")

## Check clustering accuracy
mclust::adjustedRandIndex(res$cluster, z)

#### Run BACON on ADHD-200 data. ####
## Load required packages
setwd("~/Documents/Repositories/BACON/code/clustering")
source("bacon.R")

## Load ADHD-200 data (647 50-gons with 4 clusters)
setwd("~/Documents/Repositories/BACON/data/real/adhd_200")
load("ADHD-200_compositional.Rdata")
head(side_lengths)
head(angles)

## Reduce number of vertices to speed up computation
# Apply BayesLASA to reduce 50-gons to 10-gons - also for other datasets like MPEG-7 and ETH-80 (?)
# What is the minimum number of vertices needed to describe any object in this dataset? It may vary.

## Run the model
res <- bacon(side_lengths[, 1:50], 
             angles[, 1:50], 
             K = 4,
             weight_L = 1,
             weight_A = 1,
             estimate.s = TRUE,
             estimate.r = TRUE,
             iter = 10000, 
             burn = 5000) 

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", 
     ylab = "Number of samples with starting point as 0", xlab = "Iteration")

## Check clustering accuracy
mclust::adjustedRandIndex(res$cluster, phenotype)

# ARI: 0.027; User runtime: 95 minutes
# ARI: 0.022 at 10000/5000 iter; Elapsed runtime: 104 minutes

#### Run BACON on simulated 20-gon ("decagons") data. ####
## Load required packages
setwd("~/Documents/Repositories/BACON/code/clustering")
source("bacon.R")

## Load simulated shape data (100 20-gons with 10 clusters)
setwd("~/Documents/Repositories/BACON/data/simulated/Data")
load("decagons.Rdata")
head(side_lengths)
head(angles)

## Run the model
res <- bacon(side_lengths[, 1:20], 
             angles[, 1:20], 
             K = 10, 
             weight_A = 1, 
             weight_L = 1,
             estimate.s = TRUE, 
             estimate.r = TRUE,
             iter = 10000, 
             burn = 5000)

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", 
     ylab = "Number of samples with starting point as 0", xlab = "Iteration")

## Check clustering accuracy
mclust::adjustedRandIndex(res$cluster, angles[, 21]) 

# ARI: 0.46 at 10000/5000 iter

#### Run BACON on MPEG-7 data. ####
## Load required packages
setwd("~/Documents/Repositories/BACON/code/clustering")
source("bacon.R")

## Load MPEG-7 data
setwd("~/Documents/Repositories/BACON/data/real/MPEG-7")
load("MPEG-7_k=20.Rdata")
head(side_lengths)
head(angles)

## Run the model
res <- bacon(side_lengths, 
             angles, 
             K = 20,
             weight_L = 1,
             weight_A = 1,
             estimate.s = FALSE, 
             estimate.r = FALSE,
             iter = 10000, 
             burn = 5000)

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", 
     ylab = "Number of samples with starting point as 0", xlab = "Iteration")

## Check clustering accuracy
mclust::adjustedRandIndex(res$cluster, phenotype)

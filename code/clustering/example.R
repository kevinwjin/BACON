## Run BACON on demo data
## The following section will guide to run a exemplary data using BACON.

## Set directory
# setwd("/Users/huiminli/Desktop/example")

## Load required packages
source("bacon.R")

## Load demo data
## BACON requires three inputs:
## 1. L: a m-by-n matrix of length proportions, where m is the number of n-gons 
## and n is the number of gons.
## 2. A: a m-by-n matrix of angle proportions.
## 3. K: the number of clusters in the dataset.
load("demo.RData")
head(L)
head(A)

## Run the model
## We run BACON with default settings on the above demo data.
res = bacon(side_lengths[, 1:50], angles[, 1:50], K = 4)

## The following is the output of BACON. The important output includes:
## cluster: the estimated cluster assignment.
## s_map: the estimated starting vertex indicators. 
## r_map: the estimated reverse indicators r.
res$cluster
res$s_map
res$r_map

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", 
     ylab = "Number of samples with starting point as 0", xlab = "Iteration")

## ADHD-200 dataset (subsampled to 50 vertices)
save("res", file = "ADHD-200_clustering_k=50.Rdata")
 
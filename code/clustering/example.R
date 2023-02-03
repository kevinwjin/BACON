## Run BACON on demo data
## The following section will guide to run a exemplary data using BACON.

## Set directory
# setwd("/Users/huiminli/Desktop/example")

## Load required packages
source("bacon.R")

## Load demo data
## BACON requires three inputs:
## 1. L: a m-by-n matrix of length proportions, where m is the number of n-gons and n is the number of vertices.
## 2. A: a m-by-n matrix of angle proportions.
## 3. K: Number of clusters
load("demo.RData")
head(L)
head(A)

## Run the model
## We run bacon with its defaulting setting on the above demo data.
res = bacon(L, A, K = 10)

## Following is the output of bacon. The important output including:
## cluster is the estimated cluster assignment. s_map is the estimated starting vertex indicators. 
## r_map is the estimated reverse indicators r. 
res$cluster
res$s_map
res$r_map

## Check convergence
plot(rowSums(res$s_store == 0), type = "l", ylab = "Number of samples with starting point as 0", xlab = "Iteration")

  
 
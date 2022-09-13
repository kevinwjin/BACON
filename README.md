# BLAST: Bayesian Landmark-based Shape Clustering

## Introduction
`BLAST` is an R package for landmark-based Bayesian clustering of closed polygonal chains, relying on intrinsic shape features (i.e. interior angles abd relative side lengths). The algorithm accepts a matrix, list, or dataframe containing the coordinates of closed polygonal chains, extracts the aforementioned inherent shape features, and clusters them by approximating a posterior distribution using a Markov chain Monte Carlo method.

## Installation
Execute one of the following commands in R, depending on your preference:

```R
# Install from The Comprehensive R Archive Network (CRAN)
install.packages("BLAST")

# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools:install_github("kevinwjin/BLAST")
```
## Usage
*(under development)*

## Arguments
*(under development)*

## Contents
* `code/main.R` - Execution script
* `code/model.R` - Implementation of the model in R
* `code/model.cpp` - Faster implementation of the model in C++ *(recommended)*
* `code/shape_generation.R` - Functions for generating simulated data
* `data/` - Real-world datasets

## Prerequisites
* `Rcpp` - Faster MCMC approximation

## Examples
*(under development)*

## Contributors
* [Qiwei Li](https://profiles.utdallas.edu/qiwei.li)
* [Kevin Jin](https://www.linkedin.com/in/kevin-w-jin/)
* [Huimin Li](https://www.linkedin.com/in/huimin-li-19789248)

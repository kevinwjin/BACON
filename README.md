# bsclust: Bayesian landmark-based shape clustering

## Introduction
`bsclust` is an Rcpp package for Bayesian clustering of landmark-based closed polygonal chains, relying on intrinsic data (e.g. the interior angles of the chains). The algorithm accepts a matrix, list, or dataframe containing the coordinates of closed polygonal chains and clusters them probabilistically using Markov chain Monte Carlo approximation by way of Gibbs sampling.

## Usage
*(under development)*

## Contents
* `code/main.R` - Execution script
* `code/model.R` - Implementation of the model in R
* `code/model.cpp` - Faster implementation of the model in C++ *(recommended)*
* `code/shape_generation.R` - Functions for generating example data

## Prerequisites
* `Rcpp` - Faster simulation code

## Contributors
* [Qiwei Li](https://profiles.utdallas.edu/qiwei.li)
* [Kevin Jin](https://www.linkedin.com/in/kevin-w-jin/)

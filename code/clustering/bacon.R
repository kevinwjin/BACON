################################################################################
## BACON  : Bayesian Clustering of n-gons via a Double Dirichlet Mixture Model
## Authors: Kevin Jin, Huimin Li, Stephen McKeown, and Qiwei Li
##
## Observed data is involved:
## 1) L: a m-by-n matrix of length proportions, where m is the number of n-gons and n is the number of gons
## 2) A: a m-by-n matrix of angle proportions
##
## Parameters are to be estimated:
## 3) z: a vector with m elements indicating cluster allocation parameters
## 4) pi: a vector with K elements indicating underlying group proportion
## 5) s: a vector with m elements indicating starting vertex indicators
## 6) r: a vector with m elements indicating reverse indicators
## 7) Lambda: a K-by-n matrix of concentration parameter of L
## 8) Theta: a K-by-n matrix of concentration parameter of A
## 
## The main outputs including cluster assignment, estimation of s, and 
## estimation of r (posterior probability of inclusion (PPI))
################################################################################

Rcpp::sourceCpp("BACONmcmc.cpp")

## Label switching
switch_label <- function(Theta, z, K){
  iter <- dim(Theta)[3]
  
  for (i in 1:iter){
    z_temp <- z[i,] + 1
    Theta_temp <- Theta[1, ,i]
    order_temp <- order(Theta_temp)
    
    for (j in 1:K){    
      z[i, z_temp == order_temp[j]] <- j
    }
  }
  
  return(z)
}

## Function to get estimated clusters using pairwise probability matrix (PPM)
get.ppm <- function(z_store, burn, iter, K) {
  library(mcclust)
  n = ncol(z_store)
  ppm <- matrix(0, nrow = n, ncol = n);
  
  for (ii in (burn + 1):iter) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (z_store[ii, i] == z_store[ii, j]) {
          ppm[i, j] <- ppm[i, j] + 1;
          ppm[j, i] <- ppm[j, i] + 1;
        }
      }
    }
  }
  
  ppm <- ppm/(iter - burn);
  diag(ppm) <- rep(1, n);
  
  z_ppm <- minbinder(ppm, method = "comp", max.k = K)$cl
  return(list(z_ppm = z_ppm, ppm = ppm))
  
}

## Main function
## K is the number of clusters
bacon <- function(L, A, K, iter = 2000, burn = 1000) {
  ## Initials of clusters obtained by K-means
  set.seed(12345)
  z <- kmeans(cbind(L, A), centers = K)$cluster - 1

  ## Run model
  start_time = proc.time()
  res = BACONmcmc(L, A, K, z, iter, burn)
  end_time = proc.time()
  run_time = as.numeric((end_time - start_time)[3], "secs")
  print(paste0("run time is ", round(run_time, digits = 1), "s"))
  
  ## Label switching of z by Theta
  z_store <- res$z_store
  z_store <- switch_label(res$Theta_store, res$z_store, K)
  
  ## Estimation of clusters
  tt <- get.ppm(z_store, burn, iter, K)
  z_ppm <- tt$z_ppm
  ppm <- tt$ppm
  z_map <- as.vector(res$z_map)
  
  ## Estimation of starting vertex indicators s
  s_store <- res$s_store
  s_map <- res$s_map
  
  ## Estimation of reverse indicators r (PPI)
  PPI <- res$r_ppi
  
  return(list("z_store" = res$z_store, "cluster" = z_ppm, "s_store" = s_store, "s_map" = s_map, "PPI" = PPI, "time" = run_time))
}



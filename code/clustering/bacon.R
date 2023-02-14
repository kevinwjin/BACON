################################################################################
## BACON   : Bayesian Clustering of n-gons via a Double Dirichlet Mixture Model
## Authors : Kevin Jin, Huimin Li, Stephen McKeown, and Qiwei Li
## Modified: 2023-02-14
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
switch_label <- function(z_store,  pi_store, Lambda_store, Theta_store, K){
  iter <- dim(z_store)[1]
  
  for (i in 1:iter){
    z_temp <- z_store[i,] + 1
    pi_temp <- pi_store[i, ]
    Lambda_temp <- Lambda_store[, , i]
    Theta_temp <- Theta_store[, , i]
    order_temp <- order(Theta_store[, 1, i])
    
    pi_store[i, ] <- pi_temp[order_temp]
    Lambda_store[, , i] <- Lambda_temp[order_temp, ]
    Theta_store[, , i] <- Theta_temp[order_temp, ]
    
    for (j in 1:K){    
      z_store[i, z_temp == order_temp[j]] <- j
    }
  }
  
  return(list("z_store" = z_store, "Theta_store" = Theta_store, "Lambda_store" = Lambda_store, "pi_store" = pi_store))
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

################################################################################
## Main input arguments:
## 1) L: length proportions
## 2) A: angle proportions
## 3) K: number of clusters
## 4) weight: weight of length proportions, default is 1
## 5) estimate.s: logical, TRUE (default) or FALSE (s = 0)
## 6) estimate.r: logical, TRUE (default) or FALSE (r = 0)
################################################################################

bacon <- function(L, A, K, weight = 1, estimate.s = TRUE, estimate.r = TRUE, iter = 2000, burn = 1000) {
  ## Run model
  start_time = proc.time()
  res = BACONmcmc(L, A, K, weight, estimate_s=estimate.s, estimate_r=estimate.r, iter, burn)
  end_time = proc.time()
  run_time = as.numeric((end_time - start_time)[1:3], "secs")
  print(paste0(paste0(c("user", "system", "elapsed"), " time is "), round(run_time, digits = 3), "s"))
  
  ## Label switching of z, pi, Lambda, and Theta by Theta
  tt <- switch_label(res$z_store,  res$pi_store, res$Lambda_store, res$Theta_store, K)
  z_store <- tt$z_store
  Lambda_store <- tt$Lambda_store 
  Theta_store <- tt$Theta_store
  pi_store <- tt$pi_store
  
  ## Evaluate the underlying group proportion pi 
  pi_hat = colMeans(pi_store[(burn + 1):iter, ])
  
  ## Estimation of clusters
  temp <- get.ppm(z_store, burn, iter, K)
  z_ppm <- temp$z_ppm
  ppm <- temp$ppm
  z_map <- as.vector(res$z_map)
  
  ## Estimation of starting vertex indicators s
  s_store <- res$s_store
  s_map <- res$s_map
  
  ## Estimation of reverse indicators r (PPI)
  r_store <- res$r_store
  PPI <- res$r_ppi
  r_map <- res$r_map
  
  
  ## Evaluate the concentration parameters #######################################
  Lambda_store = lapply(1:iter, function(x) Lambda_store[ , , x])
  Theta_store = lapply(1:iter, function(x) Theta_store[ , , x])
  
  Lambda_hat = Reduce("+", Lambda_store[(burn + 1):iter]) / (iter - burn)
  Theta_hat = Reduce("+", Theta_store[(burn + 1):iter]) / (iter - burn)
  
  return(list("cluster" = z_ppm, "z_map" = z_map, "s_map" = s_map, "r_map" = r_map, "PPI" = PPI,
              "pi_hat" = pi_hat, "Lambda_hat" = Lambda_hat, "Theta_hat" = Theta_hat,
              "z_store" = z_store, "pi_store" = pi_store, "s_store" = s_store, 
              "r_store" = r_store, "Lambda_store" = Lambda_store, "Theta_store" = Theta_store,
              "accept_lambda" = res$accept_lambda, "accept_theta" = res$accept_theta,
              "time" = run_time))
}



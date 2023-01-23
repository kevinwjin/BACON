################################################################################
## Project : BACON: Bayesian clustering of n-gons
## Authors : Kevin Jin, Huimin Li, Stephen McKeown, and Qiwei Li
## Goal    : Generate simulated data
## Modified: 2023-01-07
################################################################################

rm(list= ls())

## Write function to generate  simulated data ##################################
generate_simu <- function(m, n, K, alpha, eta, omega, seed = NA){
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  message(paste0("Generating simulated dataset: m = ", m ,", n = ", n, ", K = ", K, ", replicate = ", replicate))  
  
  ## Underlying group proportion
  pi <- rdirichlet(1, alpha = alpha)
  if(min(pi) < 0.01){## in case on proportion is too small, generate K-1 groups
    pi <- rdirichlet(1, alpha = alpha)
    
  }
  
  ## Cluster parameters z
  z <- sample(1:K, m, prob = pi, replace = TRUE)

  ## Starting vertex indicators s
  s <- sample(0:(n-1), size = m, replace = TRUE, prob = c(0.5, rep(0.5/(n-1), n-1)))
  
  ## Reverse indicators r
  r <- rbern(m, prob = omega)
  
  ## Concentration parameter Lambda for length L
  # Lambda <- matrix(c(c(rep(1, n-1), 5), c(rep(1, n/2), rep(5, n/2)), c(rep(5, n-1), 1)), byrow = TRUE, nrow = K)
  # Lambda <- matrix(c(c(1,3,5,7), c(7,5,3,1), c(1,5,1,5)), byrow = TRUE, nrow = K)
  Lambda <- matrix(c(c(100, 100, 100), c(10000, 1, 1), c(1000, 2000, 1000*sqrt(3))), byrow = TRUE, nrow = K)
  
  ## Concentration parameter Theta for angle A
  # Theta <- matrix(c(rep(1, n), c(rep(1, n/2), rep(5, n/2)), c(rep(5, n/2), rep(1, n/2))), byrow = TRUE, nrow = K)
  Theta <- Lambda

  ## Rearranged length proportion L'
  L2 <- matrix(NA, nrow = m, ncol = n)
  for (j in 1:m) {
    L2[j, ] <- rdirichlet(1, Lambda[z[j], ])
  }
  
  ## Rearranged angle proportion A'
  A2 <- matrix(NA, nrow = m, ncol = n)
  for (j in 1:m) {
    A2[j, ] <- rdirichlet(1, Theta[z[j], ])
  }
  
  ## Length proportion L
  L <- matrix(NA, nrow = m, ncol = n)
  
  ## Angle proportion A
  A <- matrix(NA, nrow = m, ncol = n)
  
  for (j in 1:m) {
    if (r[j] == 0 & s[j] > 0){
      L[j,] = L2[j, ][(((n-s[j]) + c(0:(n-1))*(1 - 2*r[j])) %% n) + 1]
      A[j,] = A2[j, ][(((n-s[j]) + c(0:(n-1))*(1 - 2*r[j])) %% n) + 1]
      
    }else{
      L[j,] = L2[j, ][((s[j] + c(0:(n-1))*(1 - 2*r[j])) %% n) + 1]
      A[j,] = A2[j, ][((s[j] + c(0:(n-1))*(1 - 2*r[j])) %% n) + 1]
      
    }
  }

  return(list("L"= L, "A"= A, "L2"= L2, "A2"= A2, "z"= z, "pi"= pi, "s"= s, "r"= r, "Lambda"= Lambda, "Theta" = Theta))
  
}


## Hyperparameters settings 
link <- "/Users/Huimin Li/PhD/Project/Project 3/Data/simulated_data/data/"
m <- 100   # 100 n-gons
n <- 3    # 3-gons
K <- 3     # 3 clusters
alpha <- rep(1, K)
eta <- 0.5
omega <- 0.5
seed = 99

for (replicate in 1:30) {
  data <- generate_simu(m, n, K, alpha, eta, omega, seed = seed + replicate)
  L <- data$L
  A <- data$A
  L2 <- data$L2
  A2 <- data$A2
  z <- data$z
  pi <- data$pi
  s <- data$s 
  r <- data$r 
  Lambda <- data$Lambda 
  Theta <- data$Theta
  
  save(L, A, L2, A2, z, pi, s, r, Lambda, Theta, alpha, eta, omega,
       file = paste0(link, "m_", m ,"_n_", n, "_K_", K, "_replicate_", replicate, ".RData"))
}



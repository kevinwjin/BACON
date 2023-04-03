#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppEigen)]]

using namespace Rcpp;
using namespace arma;

static int num(IntegerVector x, int c);
static arma::uvec which(IntegerVector x, int c);
static arma::uvec organize_label(IntegerVector x);

// [[Rcpp::export]]
Rcpp::List BACONmcmc(arma::mat L, arma::mat A, int K, double weight_L, double weight_A, bool estimate_s, bool estimate_r, int iter, int burn) {
  // Read data information
  int m = L.n_rows;
  int n = L.n_cols;

  // Hyperparameters 
  NumericVector alpha(K);
  for (int k = 0; k < K; k++){
    alpha(k) = 1;
  }
  double eta = 0.5; // Truncated Poisson parameter - controls shift value probability
  double omega = 0.5;
  double a_lambda = 0.001;  
  double b_lambda = 0.001;  
  double a_theta = 0.001;  
  double b_theta = 0.001;  
  
  // Tunning parameters 
  double tau_lambda = 0.1;
  double tau_theta = 0.1;
  
  // Set temporary variables
  int i, ii, j, jj, k, t, tt, g, gg, num = 10;
  double map_z, map_z_temp, s_new, map_s, map_s_temp, map_r, map_r_temp, lambda_temp, theta_temp, lambda_new, theta_new, hastings;
  double try_lambda = 0, accept_lambda = 0, try_theta = 0, accept_theta = 0;
  IntegerVector z_map(m);
  IntegerVector s_map(m);
  IntegerVector r_map(m);
  NumericMatrix Lambda_new(K, n);
  NumericMatrix Theta_new(K, n);
  bool double_dirichlet = true;
  bool joint_update = true;
  int min_element = 2;
  
  // Create the spaces to store the results 
  IntegerMatrix z_store(iter, m);
  NumericMatrix pi_store(iter, K);
  IntegerMatrix s_store(iter, m);
  IntegerMatrix r_store(iter, m);
  IntegerVector r_sum(iter);
  NumericVector r_ppi(m);
  arma::cube Lambda_store(K, n, iter);
  arma::cube Theta_store(K, n, iter);
  NumericVector map_z_store(iter);
  NumericVector map_s_store(iter);
  NumericVector map_r_store(iter);
  
  // Initialization 
  IntegerVector clusters(K);
  NumericVector prop_init(K);
  for (k = 0; k < K; k++){
    clusters(k) = k;
    prop_init(k) = pow(K, -1);
  }
  
  IntegerVector z(m);
  for (j = 0; j < m; j++){
    z(j) = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prop_init));
  }
  
  // Initials of underlying group proportions pi
  NumericVector pi(K);
  for (k = 0; k < K; k++){
    pi(k) = pow(K, -1);
  }
  
  // Initials of starting point indicators s
  IntegerVector startings(n);
  for (i = 0; i < n; i++){
    startings(i) = i;
  }
  NumericVector prob_s(n);
  for (int i = 0; i < n; i++){
    prob_s(i) = pow(n, -1);
  }
  IntegerVector s(m);
  for (j = 0; j < m; j++){
    s(j) = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(startings, 1, TRUE, prob_s));
  }
  
  IntegerVector combination(2*n);
  for (i = 0; i < 2*n; i++){
    combination(i) = i;
  }
  
  // Initials of reverse indicators r
  IntegerVector r(m);
  
  // Initials of the concentration parameters
  NumericMatrix Lambda(K, n);
  NumericMatrix Theta(K, n);
  for (k = 0; k < K; k++){
    for (i = 0; i < n; i++){
      Lambda(k, i) = 1;
      Theta(k, i) = 1;
    }
  }
  
  // Rearranged length proportion and angle proportions L2 and A2
  NumericMatrix L2(m, n);
  NumericMatrix A2(m, n);
  
  for (j = 0; j < m; j++){
    for (i = 0; i < n; i++){
      for (ii = 0; ii < n; ii++){
        int temp = ((s(j) + i*(1-2*r(j))) + n) % n;  // in case negative value, -1%4=-1, -1%4=3 in R
        if (temp == ii){
          L2(j, i) = L(j, ii);
          A2(j, i) = A(j, ii);
        }
      }
    }
  }
  
  
  // Start MCMC algorithms 
  for (t = 0; t < iter; t++) {
    // Gibbs Sampler for updating allocation parameter z
    for (j = 0; j < m; j++){
      // Calculate posterior probability of z
      NumericVector loglklh(K);

      for (k = 0; k < K; k++){
        // Prior
        loglklh(k) = pi(k);

        // Data likelihood of A2 (add weight_A >= 0)
        loglklh(k) = loglklh(k) + weight_A*(lgamma(sum(Theta(k, _))));
        for (i = 0; i < n; i++) {
          loglklh(k) = loglklh(k) - weight_A*lgamma(Theta(k, i)) + weight_A*(Theta(k, i) - 1)*log(A2(j, i));
        }

        if (double_dirichlet){
          // Data likelihood of L2 (add weight_L >= 0)
          loglklh(k) = loglklh(k) + weight_L*(lgamma(sum(Lambda(k, _))));
          for (i = 0; i < n; i++) {
            loglklh(k) = loglklh(k) - weight_L*lgamma(Lambda(k, i)) + weight_L*(Lambda(k, i) - 1)*log(L2(j, i));
          }
          
        }
      }

      NumericVector prob(K);

      for (k = 0; k < K; k++) {
        double diff = 0;
        for (int tt = 0; tt < K; tt++) {
          if (tt != k) {
            diff = diff + exp(loglklh(tt) - loglklh(k));
          }
        }
        prob(k) = 1/(1 + diff);
      }

      // Generate new sample
      int z_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prob));

      // Make sure each group has at least multiple observations
      k = z(j);
      if (sum(z == k) == min_element && z_new != k) {
        z(j) = k;
      }else {
        z(j) = z_new;
      }
      
    } // End of updating cluster allocation parameter z

    
    // Gibbs Sampler for updating the underlying group proportion pi ###########
    NumericVector alpha_new(K);

    for (k = 0; k < K; k++){
      alpha_new(k) = alpha(k) + sum(z == k);
    }

    // Sampling from a Dirichlet distribution is equivalent to sampling from following gamma distribution
    for (k = 0; k < K; k++){
      pi(k) = rgamma(1, alpha_new(k), 1)(0);
    }
    double temp = sum(pi);
    for (k = 0; k < K; k++){
      pi(k) = pi(k)/temp;
    }
    
    // End of updating the underlying group proportion pi ######################

    
    // Gibbs Sampler for jointly updating s and r
    if (joint_update && estimate_s && estimate_r) {
      for (j = 0; j < m; j++) {
        // Calculate joint posterior probability of s and r
        NumericMatrix loglklh(n, 2);
        
        for (g = 0; g < 2; g++) {
          for (i = 0; i< n; i++){
            NumericVector L2_new(n);
            NumericVector A2_new(n);
            
            for (int tt = 0; tt < n; tt++){
              for (ii = 0; ii < n; ii++){
                int temp = ((i + tt*(1-2*g)) + n) % n;  
                if (temp == ii){
                  L2_new(tt) = L(j, ii);
                  A2_new(tt) = A(j, ii);
                }
              }
            }
            
            // Prior of s
            loglklh(i, g) = loglklh(i, g) + i*log(eta) - lgamma(i + 1);
            
            // Prior
            loglklh(i, g) = loglklh(i, g) + g*log(omega) + (1 - g)*log(1 - omega);
            
            // Data likelihood of A2
            for (ii = 0; ii< n; ii++) {
              loglklh(i, g) = loglklh(i, g) + weight_A*((Theta(z(j), ii) - 1)*log(A2_new(ii)));
            }
            
            
            if (double_dirichlet){
              // Data likelihood of L2
              for (ii = 0; ii< n; ii++) {
                loglklh(i, g) = loglklh(i, g) + weight_L*((Lambda(z(j), ii) - 1)*log(L2_new(ii)));
              }
            }
          }
        }
        
        // Calculate probability of each combination
        NumericVector prob(2*n);
        
        for (g = 0; g < 2; g++) {
          for (i = 0; i < n; i++){
            double diff = 0;
            
            for (int gg = 0; gg < 2; gg++) {
              for (ii = 0; ii < n; ii++){
                if (ii != i || gg != g) {
                  diff = diff + exp(loglklh(ii, gg) - loglklh(i, g));
                }
              }
            }
            prob(i + g*n) = 1/(1 + diff);
          }
        }
        
        tt = which_max(prob);
        int temp = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(combination, 1, TRUE, prob));
        
        if (temp < n){
          r(j) = 0;
        }else{
          r(j) = 1;
        }
        s_new = temp - r[j]*n;
        
        // Make sure each group has at least multiple observations
        i = s(j);
        if (sum(s == i) == min_element && s_new != i) {
          s(j) = i;
        }else {
          s(j) = s_new;
        }
        
      }
    } // End of jointly updating the s and r 
    
    
    // Gibbs Sampler for updating the starting vertex indicators s
    if (estimate_s){
      for (j = 0; j < m; j++) {
        // Calculate posterior probability of s
        NumericVector loglklh(n);
        
        for (i = 0; i < n; i++){
          NumericVector L2_new(n);
          NumericVector A2_new(n);
          
          for (int tt = 0; tt < n; tt++){
            for (ii = 0; ii < n; ii++){
              int temp = ((i + tt*(1-2*r(j))) + n) % n;
              if (temp == ii){
                L2_new(tt) = L(j, ii);
                A2_new(tt) = A(j, ii);
              }
            }
          }
          
          // Prior
          loglklh(i) = loglklh(i) + i*log(eta) - lgamma(i+1);
          
          // Data likelihood of A2
          for (ii = 0; ii < n; ii++) {
            loglklh(i) = loglklh(i) + weight_A*((Theta(z(j), ii) - 1)*log(A2_new(ii)));
          }
          
          if (double_dirichlet){
            // Data likelihood of L2
            for (ii = 0; ii < n; ii++) {
              loglklh(i) = loglklh(i) + weight_L*((Lambda(z(j), ii) - 1)*log(L2_new(ii)));
            }
          }
        }
        
        NumericVector prob(n);
        for (i = 0; i < n; i++) {
          double diff = 0;
          for (int tt = 0; tt < n; tt++) {
            if (tt != i) {
              diff = diff + exp(loglklh(tt) - loglklh(i));
            }
          }
          prob(i) = 1/(1 + diff);
        }
        
        // Generate new sample
        int s_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(startings, 1, TRUE, prob));
        
        // Make sure each group has at least multiple observations
        i = s(j);
        if (sum(s == i) == min_element && s_new != i) {
          s(j) = i;
        }else {
          s(j) = s_new;
        }
        
      }
      
    }else{
      for (j = 0; j < m; j++) {
        s(j) = 0;
      }
      
    }
    
    // End of updating the starting vertex indicators s


    // Gibbs Sampler for updating the reverse indicators r
    if (estimate_r){
      for (j = 0; j < m; j++) {
        // Calculate posterior probability of r
        NumericVector loglklh(2);
        
        for (g = 0; g < 2; g++) {
          NumericVector L2_new(n);
          NumericVector A2_new(n);
          
          for (int tt = 0; tt < n; tt++){
            for (ii = 0; ii < n; ii++){
              int temp = ((s(j) + tt*(1-2*g)) + n) % n;
              if (temp == ii){
                L2_new(tt) = L(j, ii);
                A2_new(tt) = A(j, ii);
              }
            }
          }
          
          // Prior
          loglklh(g) = g*log(omega) + (1 - g)*log(1 - omega);
          
          // Data likelihood of A2
          for (ii = 0; ii < n; ii++) {
            loglklh(g) = loglklh(g) + weight_A*((Theta(z(j), ii) - 1)*log(A2_new(ii)));
          }
          
          if (double_dirichlet){
            // Data likelihood of L2
            for (ii = 0; ii < n; ii++) {
              loglklh(g) = loglklh(g) + weight_L*((Lambda(z(j), ii) - 1)*log(L2_new(ii)));
            }
          }
        }
        
        // Posterior probability of r = 0 and r = 1
        NumericVector prob(2);
        
        prob(0) = 1/(1 + exp(loglklh(1) - loglklh(0)));
        prob(1) = 1/(1 + exp(loglklh(0) - loglklh(1)));
        r(j) = rbinom(1, 1, prob(1))(0);
      }
    }else{
      for (j = 0; j < m; j++) {
        r(j) = 0;
      }
    }
    // End of updating the reverse indicators r ################################
    
    
    // Rearranged length proportion and angle proportions L2 and A2
    for (j = 0; j < m; j++){
      for (i = 0; i < n; i++){
        for (ii = 0; ii < n; ii++){
          int temp = ((s(j) + i*(1-2*r(j))) + n) % n;  // in case negative value, -1%4=-1, -1%4=3 in R
          if (temp == ii){
            L2(j, i) = L(j, ii);
            A2(j, i) = A(j, ii);
          }
        }
      }
    }
    
    
    // RWMH to update the concentration parameter Lambda of L (length)
    if (double_dirichlet){
      for (k = 0; k < K; k++){
        for (i = 0; i < n; i++){
          lambda_temp = Lambda(k, i);
          Lambda_new(k, _) = Lambda(k, _);
          
          // lambda_new = exp(rnorm(1, log(lambda_temp), tau_lambda)(0));
          lambda_new = exp(r_truncnorm(log(lambda_temp), tau_lambda, log(1), log(1000)));
          Lambda_new(k, i) = lambda_new;
          
          hastings = 0;
          
          // Data (Length) likelihood ratio
          for (j = 0; j < m; j++){
            if (z(j) == k){
              hastings = hastings + weight_L*(lgamma(sum(Lambda_new(k, _))));
              hastings = hastings - weight_L*(lgamma(sum(Lambda(k, _))));
              
              for (ii = 0; ii < n; ii++) {
                hastings = hastings + weight_L*(- lgamma(Lambda_new(k, ii)) + (Lambda_new(k, ii) - 1)*log(L2(j, ii)));
                hastings = hastings - weight_L*(- lgamma(Lambda(k, ii)) + (Lambda(k, ii) - 1)*log(L2(j, ii)));
              }
            }
          }
          
          // Prior ratio
          hastings = hastings + (a_lambda - 1)*log(lambda_new) - b_lambda*lambda_new;
          hastings = hastings - ((a_lambda - 1)*log(lambda_temp) - b_lambda*lambda_temp);
          
          // Check if accent the nronosed new values
          if (hastings >= log(double(rand()%10001)/10000)) {
            Lambda(k, i) = lambda_new;
            if (t >= burn){
              accept_lambda = accept_lambda + 1;
            }
          }
          
          if (t >= burn){
            try_lambda = try_lambda + 1;
          }
        }
      } 
      
    } // End of updating the concentration parameter Lambda of length
    
    
    // RWMH to update the concentration parameter Theta of angle
    for (k = 0; k < K; k++){
      for (i = 0; i < n; i++){
        theta_temp = Theta(k, i);
        Theta_new(k, _) = Theta(k, _);
        
        // theta_new = exp(rnorm(1, log(theta_temp), tau_theta)(0));
        theta_new = exp(r_truncnorm(log(theta_temp), tau_theta, log(1), log(1000)));
        Theta_new(k, i) = theta_new;
        hastings = 0;
        
        // Data (Angle) likelihood ratio
        for (j = 0; j < m; j++){
          if (z(j) == k){
            hastings = hastings + weight_A*(lgamma(sum(Theta_new(k, _))));
            hastings = hastings - weight_A*(lgamma(sum(Theta(k, _))));
            
            for (ii = 0; ii < n; ii++) {
              hastings = hastings + weight_A*(- lgamma(Theta_new(k, ii)) + (Theta_new(k, ii) - 1)*log(A2(j, ii)));
              hastings = hastings - weight_A*(- lgamma(Theta(k, ii)) + (Theta(k, ii) - 1)*log(A2(j, ii)));
            }
          }
        }
        
        // Prior ratio
        hastings = hastings + (a_theta - 1)*log(theta_new) - b_theta*theta_new;
        hastings = hastings - ((a_theta - 1)*log(theta_temp) - b_theta*theta_temp);
        
        // Check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)) {
          Theta(k, i) = theta_new;
          if (t >= burn){
            accept_theta = accept_theta + 1;
          }
        }
        
        if (t >= burn){
          try_theta = try_theta + 1;
        }
      }
    } // End of updating the concentration parameter Theta of angle

    
    // Calculate the marginal posterior probability of z, s and r to obtain MAP 
    map_z = 0;
    map_s = 0;
    map_r = 0;
    
    for (j = 0; j < m; j++) {
      // Prior
      map_z = map_z + pi(z(j));
      map_s = map_s + (s(j))*log(eta) - lgamma(s(j)+1);
      map_r = map_r + r(j)*log(omega) + (1 - r(j))*log(1 - omega);
      
      // Data likelihood of L2
      map_z = map_z + lgamma(sum(Lambda(z(j), _)));
      map_s = map_s + lgamma(sum(Lambda(z(j), _)));
      map_r = map_r + lgamma(sum(Lambda(z(j), _)));
      
      for (i = 0; i < n; i++) {
        map_z = map_z - lgamma(Lambda(z(j), i)) + (Lambda(z(j), i) - 1)*log(L2(j, i));
        map_s = map_s - lgamma(Lambda(z(j), i)) + (Lambda(z(j), i) - 1)*log(L2(j, i));
        map_r = map_r - lgamma(Lambda(z(j), i)) + (Lambda(z(j), i) - 1)*log(L2(j, i));
        
      }
      
      if (double_dirichlet){
        // Data likelihood of A2
        map_z = map_z + lgamma(sum(Theta(z(j), _)));
        map_s = map_s + lgamma(sum(Theta(z(j), _)));
        map_r = map_r + lgamma(sum(Theta(z(j), _)));
        
        for (i = 0; i < n; i++) {
          map_z = map_z - lgamma(Theta(z(j), i)) + (Theta(z(j), i) - 1)*log(A2(j, i));
          map_s = map_s - lgamma(Theta(z(j), i)) + (Theta(z(j), i) - 1)*log(A2(j, i));
          map_r = map_r - lgamma(Theta(z(j), i)) + (Theta(z(j), i) - 1)*log(A2(j, i));
          
        }
      }
    }
    
    // Obtain estimate z_map, s_map, and r_map
    map_z_store(t) = map_z;
    map_s_store(t) = map_s;
    map_r_store(t) = map_r;
    
    if (t == burn) {
      z_map = z;
      map_z_temp = map_z;
      s_map = s;
      map_s_temp = map_s;
      r_map = r;
      map_r_temp = map_r;
      
    }else if (t > burn) {
      if (map_z > map_z_temp) {
        z_map = z;
        map_z_temp = map_z;
      }
      if (map_s > map_s_temp) {
        s_map = s;
        map_s_temp = map_s;
      }
      if (map_r > map_r_temp) {
        r_map = r;
        map_r_temp = map_r;
      }
    } // End of calculate MAP of z, s, and r
    

    // Store the results
    r_sum(t) = sum(r);
    
    for (j = 0; j < m; j++){
      z_store(t, j) = z(j);
      s_store(t, j) = s(j);
      r_store(t, j) = r(j);
      
      // Calculate PPI of r
      if (t >= burn){
        r_ppi(j) = r_ppi(j) + r(j);
      }
      
    }
    
    for (k = 0; k < K; k++){
      pi_store(t, k) = pi(k);
    }
    
    Lambda_store.slice(t) = as<arma::mat>(Lambda);
    Theta_store.slice(t) = as<arma::mat>(Theta);
    
    // Monitor the process
    if((t*100/(iter-1)) == num) {
      Rcout<<num<< "% has been done\n";
      num = num + 10;
    }

  } // End of iterations
  
  // Calculate PPI  of r
  for (j = 0; j < m; j++) {
    r_ppi(j) = r_ppi(j)/(iter - burn);
  }
  
  // Calculate acceptance rate
  accept_lambda = accept_lambda/try_lambda;
  accept_theta = accept_theta/try_theta;

  return Rcpp::List::create(Rcpp::Named("z_store") = z_store, Rcpp::Named("map_z_store") = map_z_store, Rcpp::Named("z_map") = z_map,
                            Rcpp::Named("pi_store") = pi_store, 
                            Rcpp::Named("s_store") = s_store, Rcpp::Named("map_s_store") = map_s_store, Rcpp::Named("s_map") = s_map,
                            Rcpp::Named("r_store") = r_store, Rcpp::Named("map_r_store") = map_r_store, Rcpp::Named("r_map") = r_map, Rcpp::Named("r_ppi") = r_ppi, Rcpp::Named("r_sum") = r_sum,
                            Rcpp::Named("Lambda_store") = Lambda_store, Rcpp::Named("Theta_store") = Theta_store,
                            Rcpp::Named("accept_lambda") = accept_lambda, Rcpp::Named("accept_theta") = accept_theta);
}


// Used functions ##############################################################
int num(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      count++;
    }
  }
  return count;
}

arma::uvec which(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  int m = num(x, c);
  arma::uvec index(m);
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      index(count) = i;
      count++;
    }
  }
  return index;
}

Rcpp::Environment base("package:base");
Function do_unique = base["unique"];

arma::uvec organize_label(IntegerVector x) {
  int n = x.size();
  arma::uvec x_new(n);
  int count = 0;
  int K = Rcpp::unique(x).size();
  IntegerVector cluster = do_unique(x);
  
  for (int k = 0; k < K; k++) {
    for (int i = 0; i < n; i++) {
      if (x(i) == cluster(k)) {
        x_new(i) = count;
      }
    }
    count = count + 1;
  }
  return x_new;
}


  
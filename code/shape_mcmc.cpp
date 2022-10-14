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
Rcpp::List mcmc(arma::mat X, arma::mat Y, IntegerVector z, String double_dirichlet, int min_element, int iter, int burn) {
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;

  // Cluster initials, the first cluster is indexed as 0
  z = z - 1;
  int K = Rcpp::unique(z).size();

  // Hyperparameters 
  double a_alpha = 1;  
  double b_alpha = 1;  
  double a_beta = 1;  
  double b_beta = 1;  
  
  // Tunning parameters 
  double tau_alpha = 1;
  double tau_beta = 1;
  
  // Set temporary variables
  int i, j, jj, k, t, num = 10;
  double map_z_temp, alpha_temp, beta_temp, alpha_new, beta_new, hastings;
  double try_alpha = 0, accept_alpha = 0, try_beta = 0, accept_beta = 0;
  IntegerVector z_map(n);
  IntegerVector clusters(K);
  
  for (k = 0; k < K; k++){
    clusters(k) = k;
  }
  
  NumericMatrix A_new(K, p);
  NumericMatrix B_new(K, p);
  
  // Create the spaces to store the results 
  IntegerMatrix z_store(iter, n);
  arma::cube A_store(K, p, iter);
  arma::cube B_store(K, p, iter);
  NumericVector map_z_store(iter);
  
  // Initialization 
  NumericMatrix A(K, p);
  NumericMatrix B(K, p);
  for (k = 0; k < K; k++){
    for (j = 0; j < p; j++){
      A(k, j) = 1;
      B(k, j) = 1;
    }
  }
  
  
  // Start MCMC algorithms 
  for (t = 0; t < iter; t++) {
    // Gibbs Sampler for updating allocation parameter z 
    for (i = 0; i < n; i++){
      // Calculate posterior probability of z
      NumericVector loglklh(K);
      
      for (k = 0; k < K; k++){
        // Data likelihood of X
        loglklh(k) = loglklh(k) + lgamma(sum(A(k, _)));
        for (j = 0; j < p; j++) {
          loglklh(k) = loglklh(k) - lgamma(A(k, j)) + (A(k, j) - 1)*log(X(i, j));
        }
        
        if (double_dirichlet == "TRUE"){
          // Data likelihood of Y
          loglklh(k) = loglklh(k) + lgamma(sum(B(k, _)));
          for (j = 0; j < p; j++) {
            loglklh(k) = loglklh(k) - lgamma(B(k, j)) + (B(k, j) - 1)*log(Y(i, j));
          }
        }
      }
      
      NumericVector prob(K);
      for (k = 0; k < K; k++) {
        double diff = 0;
        for (int m = 0; m < K; m++) {
          if (m != k) {
            diff = diff + exp(loglklh(m) - loglklh(k));
          }
        }
        prob(k) = 1/(1 + diff);
      }
      
      // Generate new sample
      int z_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prob));
      
      // Make sure each group has at least multiple observations
      k = z(i);
      if (sum(z == k) == min_element && z_new != k) {
        z(i) = k;
      }else {
        z(i) = z_new;
      }
      
      // To make the first observation is always in the group 1
      z = organize_label(z);
      
    } // End of updating cluster allocation parameter z
    
    
    // Calculate the marginal posterior probability to obtain MAP
    double map_z = 0;
    for (i = 0; i < n; i++) {
      z_store(t, i) = z(i);
        
      // Data likelihood of X
      map_z = map_z + lgamma(sum(A(z(i), _)));
      for (j = 0; j < p; j++) {
          map_z = map_z - lgamma(A(z(i), j)) + (A(z(i), j) - 1)*log(X(i, j));
      }
        
      if (double_dirichlet == "TRUE"){
      // Data likelihood of Y
        map_z = map_z + lgamma(sum(B(z(i), _)));
        for (j = 0; j < p; j++) {
          map_z = map_z - lgamma(B(z(i), j)) + (B(z(i), j) - 1)*log(Y(i, j));
        }
      }
      
    }
    
    
    // Obtain estimate z_map
    map_z_store(t) = map_z;
    if (t == burn) {
      z_map = z;
      map_z_temp = map_z;
      
    }else if (t > burn) {
      if (map_z > map_z_temp) {
        z_map = z;
        map_z_temp = map_z;
      }
    } // End of calculate MAP of z
    
    
    // RWMH to update concentration parameter A  
    for (k = 0; k < K; k++){
      for (j = 0; j < p; j++){
        alpha_temp = A(k, j);
        A_new(k, _) = A(k, _);
        
        alpha_new = exp(r_truncnorm(log(alpha_temp), tau_alpha, log(1), log(10)));
        A_new(k, j) = alpha_new;
        hastings = 0;
        
        // Data (X) likelihood ratio
        for (i = 0; i < n; i++){
          if (z(i) == k){
            hastings = hastings + lgamma(sum(A_new(k, _)));
            hastings = hastings - lgamma(sum(A(k, _)));
            
            for (jj = 0; jj < p; jj++) {
              hastings = hastings + (- lgamma(A_new(k, jj)) + (A_new(k, jj) - 1)*log(X(i, jj)));
              hastings = hastings - (- lgamma(A(k, jj)) + (A(k, jj) - 1)*log(X(i, jj)));
            }
          }
        }
        
        // Prior ratio
        hastings = hastings + (a_alpha - 1)*log(alpha_new) - b_alpha*alpha_new;
        hastings = hastings - ((a_alpha - 1)*log(alpha_temp) - b_alpha*alpha_temp);
        
        // Check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)) {
          A(k, j) = alpha_new;
          if (t >= burn){
            accept_alpha = accept_alpha + 1;
          }
        }
        
        if (t >= burn){
          try_alpha = try_alpha + 1;
        }
      }
    } // End of updating concentration parameter A
    
    
    // RWMH to update concentration parameter B 
    if (double_dirichlet == "TRUE"){
      for (k = 0; k < K; k++){
        for (j = 0; j < p; j++){
          beta_temp = B(k, j);
          B_new(k, _) = B(k, _);
          
          beta_new = exp(r_truncnorm(log(beta_temp), tau_beta, log(1), log(10)));
          B_new(k, j) = beta_new;
          
          hastings = 0;
          
          // Data (Y) likelihood ratio
          for (i = 0; i < n; i++){
            if (z(i) == k){
              hastings = hastings + lgamma(sum(B_new(k, _)));
              hastings = hastings - lgamma(sum(B(k, _)));
              
              for (jj = 0; jj < p; jj++) {
                hastings = hastings + (- lgamma(B_new(k, jj)) + (B_new(k, jj) - 1)*log(Y(i, jj)));
                hastings = hastings - (- lgamma(B(k, jj)) + (B(k, jj) - 1)*log(Y(i, jj)));
              }
            }
          }
          
          // Prior ratio
          hastings = hastings + (a_beta - 1)*log(beta_new) - b_beta*beta_new;
          hastings = hastings - ((a_beta - 1)*log(beta_temp) - b_beta*beta_temp);
          
          // Check if accept the proposed new values
          if (hastings >= log(double(rand()%10001)/10000)) {
            B(k, j) = beta_new;
            
            if (t > burn){
              accept_beta = accept_beta + 1;
            }
          }
          
          if (t > burn){
            try_beta = try_beta + 1;
          }
        }
      }
    } // End of updating concentration parameter B
    
    
    // Monitor the process
    if((t*100/(iter-1)) == num) {
      Rcout<<num<< "% has been done\n";
      num = num + 10;
    }
    
    // Store the results
    A_store.slice(t) = as<arma::mat>(A);
    B_store.slice(t) = as<arma::mat>(B);
        

  } // End of iterations
  
  accept_alpha = accept_alpha/try_alpha;
  accept_beta = accept_beta/try_beta;

  return Rcpp::List::create(Rcpp::Named("z_store") = z_store, Rcpp::Named("map_z_store") = map_z_store, Rcpp::Named("z_map") = z_map,
                                        Rcpp::Named("A_store") = A_store, Rcpp::Named("B_store") = B_store,
                                        Rcpp::Named("accept_alpha") = accept_alpha, Rcpp::Named("accept_beta") = accept_beta);
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
  int K = (Rcpp::unique(x)).size();
  IntegerVector cluster = do_unique(x);
  
  for (int j = 0; j < K; j++) {
    for (int i = 0; i < n; i++) {
      if (x(i) == cluster[j]) {
        x_new(i) = count;
      }
    }
    count = count + 1;
  }
  return x_new;
}


  
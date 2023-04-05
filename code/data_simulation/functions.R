# Functions for random polygonal chain generation and manipulation
# Authors: Kevin Jin, Bryn Brakefield

#' Generate random compositional data for BACON by drawing from a truncated 
#' Dirichlet distribution (TDD)
#' 
#' @author Bryn Brakefield
#'
#' @description Randomly generates a compositional angle or side length vector
#' for BACON by drawing from a TDD n times. Implemented in Rcpp, with a 
#' quantile function for a truncated beta distribution (qtbeta1) and a 
#' random generator of one sample from a TDD (rtdirichlet1).
#' 
#' @param n Number of samples to generate.
#' @param eta Vector of size k determining balance of the compositional data (
#' the shapes being generated).
#' @param a Vector of size k of minimum parameters for the TDD, usually all 0
#' @param b Vector of size k of maximum parameters for the TDD; no side length
#' or angle may exceed 0.5 due to geometric constraints
#'
#' @return A n x k(dimension of TDD) matrix containing the compositional data 
#' of one cluster of data for BACON.
#' 
require(Rcpp)
require(RcppArmadillo)
sourceCpp(code = '
  #include <RcppArmadillo.h>
  
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace arma;
  
  // [[Rcpp::export]]
  double qtbeta1(double p, double alpha, double beta, double a, double b) {
    double Fa = R::pbeta(a, alpha, beta, true, false);
    double Fb = R::pbeta(b, alpha, beta, true, false);
    double u = Fa + p*(Fb - Fa);
    double quan = R::qbeta(u, alpha, beta, true, false);
    return quan;
  }
  
  // [[Rcpp::export]]
  mat rtdirichlet1(vec eta, vec a, vec b) {
    int k = eta.size();
    vec u = randu(k);
    vec x = zeros(k);
    x[k - 2] = qtbeta1(u[k - 2], eta[k - 2], sum(eta) - eta[k - 2], std::max(a[k - 2], 1 - sum(b) + b[k - 2]), std::min(b[k - 2], 1 - sum(a) + a[k - 2]));
    for (int i = k - 3; i >= 0; i--) {
      double a0 = std::max(a[i]/(1 - sum(x(span(i + 1, k - 2)))), 1 - (sum(b) - sum(b(span(i, k - 2))))/(1 - sum(x(span(i + 1, k - 2)))));
      double b0 = std::min(b[i]/(1 - sum(x(span(i + 1, k - 2)))), 1 - (sum(a) - sum(a(span(i, k - 2))))/(1 - sum(x(span(i + 1, k - 2)))));
      x[i] = (1 - sum(x(span(i + 1, k - 2))))*qtbeta1(u[i], eta[i], sum(eta) - sum(eta(span(i, k - 2))), a0, b0);
    }
    x[k - 1] = 1 - sum(x);
    return x;
  }')
# Randomly generates n samples from TDD
rtdirichlet <- function(n, eta, a, b) {
  return(t(replicate(n, as.vector(rtdirichlet1(eta, a, b)))))
}

#' Generate a random closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description Randomly generates the vertices of a closed polygonal chain, 
#' sorting the vertices by angle sweep from the x-axis.
#'
#' @param k Number of vertices to generate.
#' @param min Minimum parameter for the uniform distribution.
#' @param max Maximum parameter for the uniform distribution.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
generate <- function(k = 3, min = 0, max = 1) {
  chain <- matrix(nrow = k, ncol = 2, byrow = FALSE)
  colnames(chain) <- c("x", "y")
  chain[, "x"] <- runif(k, min = min, max = max)
  chain[, "y"] <- runif(k, min = min, max = max)
  
  # Compute the centroid of the chain
  centroid <- rowSums(t(chain)) / nrow(chain)
  
  # Draw distance vectors between centroid and vertices and calculate angles
  sorting_angles <- matrix(nrow = nrow(chain), ncol = 1)
  for (n in seq_len(nrow(chain))) {
    dist <- chain[n, ] - centroid
    sorting_angles[n] <- atan2(dist[2], dist[1])
  }
  
  # Vertices will have anticlockwise order
  chain <- chain[order(sorting_angles, decreasing = FALSE), ]
  
  # Repeat first row at end to form closed chain
  chain <- rbind(chain, chain[1, ])
  
  return(chain)
}

#' Generate a dataset of closed polygonal chains
#' 
#' @author Kevin Jin
#'
#' @description Generates a dataset of x closed polygonal chains (k-gons), with 
#' z evenly spaced clusters, n k-gons per cluster, and a modifiable amount of
#' jitter.
#'
#' @param x Number of shapes total in the dataset.
#' @param z Number of clusters in the dataset.
#' @param n Number of shapes in each evenly-spaced cluster.
#' @param k Number of vertices in each shape.
#' @param jitter_factor Amount of jitter to apply to each shape within a cluster.
#'
#' @return A list of sublists (matrices) containing the dataset. Each main list 
#' represents a cluster, and each sublist (matrix) represents a shape within 
#' that cluster.
simulate_shapes <- function(x, z, n, k, jitter_factor) {
  dataset <- list()
  for (cluster in 1:z) {
    shape <- generate(k = k)
    shapes <- replicate(n, shape, simplify = FALSE)
    for (i in seq_along(shapes)) {
      shapes[[i]] <- 
        jitter(shapes[[i]], factor = jitter_factor)
    }
    dataset[[length(dataset) + 1]] <- shapes
  }
  return(dataset)
}

#' Check whether a polygonal chain is closed
#' 
#' @author Kevin Jin
#'
#' @description
#' Determines if the given chain of coordinates is a polygonal chain or not
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#' @param reconstruct A boolean indicating whether the calling function is
#' the mapping function, which requires a rounding in the closed chain check.
#'
#' @return A logical value indicating whether the given chain is a polygonal
#' chain.
is_closed <- function(chain, reconstruct) {
  # Order of vertices must remain the same and might be disturbed after jitter.
  if (reconstruct) {
    # Also, you must allow some room for error because the reconstruct() 
    # functionand side length/angle proportion extraction processes might 
    # introduce error.
    if (identical(chain[1, ], round(chain[nrow(chain), ], digits = 10))) {
      is_closed <- TRUE
    } else {
      is_closed <- FALSE
    }
  } else {
    if (identical(chain[1, ], chain[nrow(chain), ])) {
      is_closed <- TRUE
    } else {
      is_closed <- FALSE
    }
  }
  return(is_closed)
}

#' Reverse the orientation of a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Reverse the orientation, or matrix row order, of a closed polygonal chain.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return The same polygonal chain, except now in reverse orientation.
reverse_orientation <- function(chain) {
  if (is_closed(chain, reconstruct = FALSE)) {
    # Reverse clockwise orientation
    chain <- chain[nrow(chain):1, ]
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(chain)
}

#' Calculate the sum of the interior angles of a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Given a closed polygonal chain, return the sum of the interior angles.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A numeric containing the sum of the interior angles of the chain.
sum_interior_angles <- function(chain) {
  return((nrow(chain) - 3) * 180)
}

#' Calculate the angle between three points
#' 
#' @author Kevin Jin
#'
#' @description
#' Given a set of three points, calculate the angle between them using the 
#' Law of Cosines.
#'
#' @param points A 3 x 2 matrix containing the x-y coordinates of the vertices
#' of three points.
#'
#' @return A numeric containing the angle, in degrees, between the three points.
three_point_angle <- function(points) {
  pointA <- points[1, ]
  pointB <- points[2, ]
  pointC <- points[3, ]
  
  x1x2s <- (pointA[1] - pointB[1]) ^ 2
  x1x3s <- (pointA[1] - pointC[1]) ^ 2
  x2x3s <- (pointB[1] - pointC[1]) ^ 2
  
  y1y2s <- (pointA[2] - pointB[2]) ^ 2
  y1y3s <- (pointA[2] - pointC[2]) ^ 2
  y2y3s <- (pointB[2] - pointC[2]) ^ 2
  
  angle <- acos((x1x2s + y1y2s + x2x3s + y2y3s - x1x3s - y1y3s) /
                  (2 * sqrt(x1x2s + y1y2s) * sqrt(x2x3s + y2y3s)))
  
  return(angle * (180 / pi))
}

#' Check whether a vertex produces a reflex angle in a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Given the index of a vertex within a closed polygonal chain and the chain 
#' itself, close the chain without the vertex by drawing a third line and 
#' determine whether the given vertex is inside or outside of the resultant 
#' new polygonal chain.
#'
#' @param index A numeric containing the index of the vertex within the chain.
#' 
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A logical value indicating whether the given vertex is in the
#' polygonal chain.
is_reflex <- function(index, chain) {
  require(sf)
  
  # Get offending vertex from the index
  point <- chain[index, ]
  
  # Cast the offending vertex to sf point object
  point_sf <- st_as_sf(data.frame(t(point)), coords = c("x", "y"))
  
  # Create a new polygonal chain without the offending vertex
  new_chain <- chain[-index, ]
  new_chain <- rbind(new_chain, new_chain[1, ]) # Duplicate first vertex
  
  # Cast the new polygonal chain to sf polygon object
  chain_sf <- st_sfc(st_polygon(list(new_chain)))
  
  # Determine whether offending vertex is inside the new polygonal chain
  is_reflex <- st_within(point_sf, chain_sf, sparse = FALSE)[, 1]
  
  return(is_reflex)
}

#' Calculate the relative interior angles of a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' If the given chain of coordinates is a closed polygonal chain, 
#' return a vector of its relative interior angles.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A vector of length k containing the interior angles of the 
#' vertices of the polygonal chain.
get_interior_angles <- function(chain) {
  if (is_closed(chain, reconstruct = FALSE)) {
    # Number of vertices/rows is (k + 1) because the chain is closed
    n <- nrow(chain)
    
    # Loop over entire chain, checking each vertex for reflex angle
    reflex <- c()
    for (i in 1:n) {
      reflex <- c(reflex, is_reflex(i, chain))
    }
    
    # Create empty vector of length k to store interior angles
    angle <- rep(NA, n - 1)
    
    # Add penultimate vertex to first position for looping purposes
    a_chain <- rbind(chain[n - 1, ], chain)
    a_reflex <- append(reflex[n - 1], reflex)
    
    # Loop over entire chain, calculating the interior angles
    for (i in 2:n) {
      if (a_reflex[i]) { # Vertex produces a reflex angle; take the complement
        angle[i - 1] <- (360 - three_point_angle(a_chain[(i - 1):(i + 1), ]))
      } else { # Vertex does not produce a reflex angle
        angle[i - 1] <- three_point_angle(a_chain[(i - 1):(i + 1), ])
      }
    }
    
    # Normalize interior angles by the total to get relative interior angle
    angle <- compositional(angle, sum(angle))
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(angle)
}

#' Calculate the relative side lengths of a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' If the given chain of coordinates is a closed polygonal chain, 
#' return a vector of its side lengths relative to the perimeter.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A vector of length k containing the side lengths of the polygonal 
#' chain.
get_side_lengths <- function(chain) {
  if (is_closed(chain, reconstruct = FALSE)) {
    # Eliminate repeated row for now
    chain <- chain[-nrow(chain), ]
    
    # Store side lengths in matrix
    side_lengths <- matrix(nrow = nrow(chain), ncol = 1)
    
    # Calculate the perimeter of the whole chain via Euclidean distance
    perimeter <- 0
    for (i in seq_len(nrow(chain))) {
      j <- i + 1
      if (i == (nrow(chain))) {
        j <- 1 # Loop back to first vertex once end of chain is reached
      }
      side_length <- sqrt(sum((chain[i, ] - chain[j, ]) ^ 2))
      side_lengths[i, ] <- side_length
      perimeter <- perimeter + side_length
    }
    
    # Normalize side lengths by the perimeter to get relative length
    side_lengths <- compositional(side_lengths, perimeter)
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(c(side_lengths))
}

#' Normalize a numerical vector by its total to produce its compositional data
#' 
#' @author Kevin Jin
#'
#' @description
#' Given a numerical vector, return a vector of its relative compositional data.
#'
#' @param data A numeric vector of length k containing data.
#' @param sum A numeric containing the sum of the vector.
#'
#' @return A vector of length k containing compositional data.
compositional <- function(data, sum) {
  normalized <- data / sum
  return(normalized)
}

#' Reconstruct all possible closed unit polygonal chains given relative interior 
#' angles and relative side lengths
#' 
#' @author Bryn Brakefield, Kevin Jin
#'
#' @description
#' Given two numerical vectors of compositional data containing relative 
#' interior angles and relative side lengths, reconstruct all possible 
#' unit polygonal chains and return the closed ones.
#'
#' @param a A numeric vector of length n containing interior angles.
#' @param l A numeric vector of length n containing relative side lengths.
#'
#' @return A list or (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the unit polygonal chain.
reconstruct <- function(a, l) {
  # Get number of vertices
  n <- length(a)
  
  # Convert compositional angle proportions back to radians
  a <- a * (n - 2) * pi
  
  # If the shape is a triangle, then ignore the provided side length
  # proportions and manually recalculate using the law of sines
  if (n == 3) {
    l <- sin(a) / sum(sin(a))
    l <- c(l[3], l[-3]) # Change ordering to be correct
  }
  
  pc0 <- expand.grid(replicate(n - 1, 0:1, simplify = FALSE))
  pc <- vector(mode = "list", length = 2^(n - 1))
  for (i in 1:(2 ^ (n - 1))) {
    V <- matrix(0, 2, n + 1)
    V[1, 2] <- l[1]
    theta <- rep(0, n + 1)
    for (j in 1:(n - 1)) {
      index <- j + 1
      if (pc0[i, j] == 1) {
        theta[index] <- theta[index - 1] + (pi - a[index])
        V[, index + 1] <- c(l[index] * cos(theta[index]), 
                            l[index] * sin(theta[index])) + V[, index]
      }
      if (pc0[i, j] == 0) {
        theta[index] <- theta[index - 1] - (pi - a[index])
        V[, index + 1] <- c(l[index] * cos(theta[index]), 
                            l[index] * sin(theta[index])) + V[, index]
      }
    }
    pc[[i]] <- V
  }
  
  # Create new list of closed chains (check is rounded to 10 digits)
  pc_closed <- list()
  for (i in 1:length(pc)) {
    if (is_closed(t(pc[[i]]), reconstruct = TRUE)) {
      pc_closed[[length(pc_closed) + 1]] <- list(t(pc[[1]]))
    } else {
      next
    }
  }
  
  return(pc_closed)
}

#' Add jitter to a polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description Randomizes the vertices or angles of a polygonal chain. If
#' output is not a valid polygonal chain, re-jitter a maximum of 10 times 
#' before giving up.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#' @param factor Floating point number from 0 to 1 exclusive as a percentage
#' of the perimeter of the polygonal chain to randomize by.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
jitter <- function(chain, factor = 0.01) {
  loops <- 0 # Jitter loop counter
  max_loops <- 10 # Maximum jitter loops
  
  # Eliminate repeated row for now
  chain <- chain[-nrow(chain), ]
  
  repeat {
    # Calculate the perimeter of the whole chain via Euclidean distance
    perimeter <- 0
    for (i in seq_len(nrow(chain))) {
      j <- i + 1
      if (i == (nrow(chain))) {
        j <- 1 # Loop back to first vertex once end of chain is reached
      }
      perimeter <- perimeter + sqrt(sum((chain[i, ] - chain[j, ]) ^ 2))
    }
    # Choose a random point in a radius of uncertainty around the vertex
    for (n in seq_len(nrow(chain))) {
      radius <- perimeter * factor
      r <- radius * sqrt(runif(1))
      theta <- runif(1) * 2 * pi
      
      # Convert polar to Cartesian coordinates for easier plotting
      chain[n, 1] <- chain[n, 1] + r * cos(theta) # Newly randomized x
      chain[n, 2] <- chain[n, 2] + r * sin(theta) # Newly randomized y
    }
    
    # If jittered output is not a chain, then re-jitter at most max_loops times
    if (is_closed(chain, reconstruct = FALSE) == FALSE || loops == max_loops) {
      if (loops == max_loops) {
        warning("Could not create valid chain in 10 attempts; please lower the jitter factor.")
      }
      break
    } else {
      loops <- loops + 1
    }
  }
  # Add repeated row back
  chain <- rbind(chain, chain[1, ])
  return(chain)
}

#' Translate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Moves every point of a polygonal chain by the same distance in a given
#' direction.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param x Horizontal distance to translate the chain by.
#' @param y Vertical distance to translate the chain by.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the translated chain.
translate <- function(chain, x = 0, y = 0) {
  # Horizontal translation
  chain[, "x"] <- chain[, "x"] + x
  # Vertical translation
  chain[, "y"] <- chain[, "y"] + y
  return(chain)
}

#' Dilate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Scales the size of a polygonal chain to be greater or smaller, centered
#' about its centroid.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param factor Positive floating-point number to scale the chain by.
#' A factor > 1 dilates the chain, while a factor > 0 and < 1 shrinks the chain.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the dilated chain.
dilate <- function(chain, factor = 1) {
  if (factor <= 0) {
    stop("Dilation factor must be greater than 0.\n")
  } else {
    # Calculate the centroid vector
    centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
    chain <- (chain - centroid) * factor + centroid
  }
  return(chain)
}

#' Reflect a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Reflects a polygonal chain vertically or horizontally.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param direction String containing direction in which to reflect the chain.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the reflected chain.
#'
reflect <- function(chain, direction = c("horizontal", "vertical")) {
  centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
  if (direction == "horizontal") {
    # Reflect across the y-axis
    chain <- (chain - centroid) %*%
      matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) + centroid
    
    # For the sf function
    colnames(chain) <- c("x", "y")
    
  } else if (direction == "vertical") {
    # Reflect across the x-axis
    chain <- (chain - centroid) %*%
      matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) + centroid
    
    # For the sf function
    colnames(chain) <- c("x", "y")
    
  } else {
    stop("Invalid or no direction provided.\n")
  }
  return(chain)
}

#' Rotate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Rotates a polygonal chain by a specified angle about its centroid.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param angle Rotation angle in degrees.
#' @param clockwise Rotate the chain clockwise if true, counterclockwise if false.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the rotated chain.
#'
rotate <- function(chain, angle, clockwise = TRUE) {
  # Convert argument angle to radians, as R's trigonometric functions use radians
  angle <- angle * (pi / 180)
  
  # Eliminate repeated row for now
  chain <- chain[-nrow(chain), ]
  
  # Calculate centroid vector
  centroid <- matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))]
  if (clockwise) {
    # Calculate counterclockwise rotation matrix
    rotation <- matrix(c(cos(angle), -sin(angle),
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
    
  } else {
    # Calculate clockwise rotation matrix
    rotation <- matrix(c(cos(angle), sin(angle),
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
  }
  
  # Add repeated row back
  chain <- t(chain)
  chain <- rbind(chain, chain[1, ])
  
  # For the sf function
  colnames(chain) <- c("x", "y")
  
  # Transpose chain to original form
  return(chain)
}

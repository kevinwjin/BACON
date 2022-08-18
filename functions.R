# Functions for polygonal chain generation
# Author: Kevin Jin

#' Generate a random polygonal chain
#' 
#' @description Randomly generates the vertices of a polygonal chain by drawing 
#' from a uniform distribution, sorting vertices by angle sweep from the x-axis.
#' 
#' @param k Number of vertices to generate.
#' @param min Minimum parameter for the uniform distribution.
#' @param max Maximum parameter for the uniform distribution.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of the 
#' polygonal chain.
generate <- function(k = 3, min = 0, max = 1) {
  chain <- matrix(nrow = k, ncol = 2, byrow = FALSE)
  colnames(chain) <- c("x", "y")
  chain[, "x"] <- runif(k, min = min, max = max) # In increasing order
  chain[, "y"] <- runif(k, min = min, max = max)
  
  centroid <- rowSums(t(chain))/nrow(chain) # Compute centroid
  chain_p <- matrix(nrow = k, ncol = 2, byrow = FALSE) # Convert Cartesian to polar
  colnames(chain_p) <- c("r", "t")
  for (n in 1:nrow(chain)) {
    dist <- chain[n, ] - centroid # Draw vector between points and centroid
    r <- sqrt(dist[1] ^ 2 + dist[2] ^ 2)
    t <- atan2(dist[2], dist[1]) * (180 / pi)
    chain_p[n, "r"] <- r
    chain_p[n, "t"] <- t
  }
  chain <- chain[order(chain_p[, "t"]), ] # Sort vertices by increasing angle sweep from x-axis
  return(chain)
}

#' Randomize a polygonal chain
#' 
#' @description Randomizes the vertices or angles of a polygonal chain.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param random String containing polygon parameter to randomize.
#' @param factor Floating point number from 0 to 1 exclusive as a percentage
#' of the perimeter of the polygonal chain to randomize by.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of the 
#' polygonal chain.
jitter <- function(chain, random = c("vertices", "angles"), factor = 0.1) {
  if (random == "vertices") {
    for (n in 1:nrow(chain)) { # For each vertex...
      # Choose a random point in a radius of uncertainty around the vertex
      perimeter <- 1 # Placeholder; use Euclidean distance
      radius <- perimeter * factor
      r <- radius * sqrt(runif(1))
      theta <- runif(1) * 2 * pi
      
      # Convert to Cartesian
      chain[n, 1] <- chain[n, 1] + r * cos(theta) # Newly randomized x
      chain[n, 2] <- chain[n, 2] + r * sin(theta) # Newly randomized y
      
      # Check if the chain is valid after jittering (call separate function)
      # Order of the vertices must remain the same and might be disturbed after jitter
      # Re-jitter if the chain if invalid (set maximum # loops for while loop and warn to shrink factor if unsuccessful)
    }
  } else if (random == "angles") {
    stop("Not implemented yet.")
  } else {
    stop("Invalid or no randomization argument provided.")
  }
  return(chain)
}

#' Translate a polygonal chain
#' 
#' @description 
#' Moves every point of a polygonal chain by the same distance in a given
#' direction.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param x Horizontal distance to translate the chain by.
#' @param y Vertical distance to translate the chain by.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the translated chain.
translate <- function(chain, x = 0, y = 0) {
  chain[, "x"] <- chain[, "x"] + x # Horizontal
  chain[, "y"] <- chain[, "y"] + y # Vertical
  return(chain)
}

#' Dilate a polygonal chain
#' 
#' @description 
#' Scales the size of a polygonal chain to be greater or smaller around its 
#' centroid.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param factor Positive floating-point number to scale the chain by. 
#' A factor > 1 dilates the chain, while a factor > 0 and < 1 shrinks the chain.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the dilated chain.
dilate <- function(chain, factor = 1) {
  if (factor <= 0) {
    stop("Dilation factor must be greater than 0.\n")
  } else {
    centroid <- t(matrix(rowSums(t(chain))/nrow(chain))[, rep(1, each = nrow(chain))]) # Create centroid matrix
    chain <- (chain - centroid) * factor + centroid
  }
  return(chain)
}

#' Flip a polygonal chain
#' 
#' @description 
#' Flips the polygonal chain vertically or horizontally.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param direction String containing direction in which to flip the chain.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the flipped chain.
#' 
flip <- function(chain, direction = c("horizontal", "vertical")) {
  centroid <- t(matrix(rowSums(t(chain))/nrow(chain))[, rep(1, each = nrow(chain))])
  if (direction == "horizontal") {
    chain <- (chain - centroid) %*% 
      matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) + centroid # Flip across the y-axis
  } else if (direction == "vertical") {
    chain <- (chain - centroid) %*% 
      matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) + centroid # Flip across the x-axis
  } else {
    stop("Invalid or no direction provided.\n")
  }
  return(chain)
}

#' Rotate a polygonal chain
#' 
#' @description 
#' Rotates a polygonal chain by a specified angle about its centroid.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param angle Rotation angle in degrees.
#' @param clockwise Rotate the chain clockwise if true, counterclockwise if false.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the rotated chain.
#' 
rotate <- function(chain, angle, clockwise = TRUE) {
  angle <- angle * (pi / 180) # Convert angle to radians
  centroid <- matrix(rowSums(t(chain))/nrow(chain))[, rep(1, each = nrow(chain))]
  if (clockwise) {
    rotation <- matrix(c(cos(angle), -sin(angle), 
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- rotation %*% (t(chain) - centroid) + centroid
  } else {
    rotation <- matrix(c(cos(angle), sin(angle), 
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- (t(chain) - centroid) %*% rotation + centroid
  }
  return(t(chain)) # Transpose back to original form
}
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
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of the 
#' polygonal chain.
randomize <- function(chain, random = c("vertices, angles")) {
  if (random == "vertices") {
    chain <- jitter(chain, factor = 10)
  } else if (random == "angles") {
    stop("Angle variation not implemented yet.")
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
#' Scales the size of a polygonal chain to be greater or smaller.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param factor Positive or negative integer to scale the chain by.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the dilated chain.
dilate <- function(chain, factor = 1) {
  if (factor == 0) {
    stop("Cannot dilate chain by zero.\n")
  } else {
    dilation <- matrix(c(factor, 0, 0, factor), nrow = 2, ncol = 2, byrow = TRUE)
    chain <- chain %*% dilation
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
  if (direction == "horizontal") {
    chain <- chain %*% matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) # Flip across the y-axis
  } else if (direction == "vertical") {
    chain <- chain %*% matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) # Flip across the x-axis
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
  centroid <- matrix(rowSums(t(chain))/nrow(chain))[, rep(1, each = 3)] # Create centroid matrix
  if (clockwise) {
    rotation <- matrix(c(cos(angle), -sin(angle), 
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- rotation %*% (t(chain) - centroid) + centroid
  } else {
    rotation <- matrix(c(cos(angle), sin(angle), 
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- rotation %*% (t(chain) - centroid) + centroid
  }
  return(t(chain)) # Transpose back to original form
}
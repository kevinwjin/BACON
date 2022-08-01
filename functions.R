# Functions for polygonal chain generation
# Author: Kevin Jin

#' Generate a random polygonal chain
#' 
#' @description Randomly generates the vertices of a polygonal chain by drawing 
#' from a uniform distribution.
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
  chain[, "x"] <- sort(runif(k, min = min, max = max)) # In increasing order
  chain[, "y"] <- sort(runif(k, min = min, max = max))
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

#' Scale a polygonal chain
#' 
#' @description 
#' Scales the size of a polygonal chain to be greater or smaller around its
#' centroid.
#' 
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
#' @param factor Positive or negative integer to scale the chain by.
#' 
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the scaled chain.
scale <- function(chain, factor = 1) {
  if (factor < 0) {
    chain <- chain * (1 / -factor)
  } else if (factor > 0) {
    chain <- chain * factor
  } else {
    cat("Error: Cannot scale chain by zero.\n")
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
    chain <- triangle %*% matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) # Flip across the y-axis
  } else if (direction == "vertical") {
    chain <- triangle %*% matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) # Flip across the x-axis
  } else {
    cat("Invalid direction.\n")
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
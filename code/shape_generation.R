# Probabilistic polygonal chain generation and manipulation
# Author: Kevin Jin

## To-do:
## 1. Determine whether validate function is necessary, because it currently
## checks if the first and last vertices are the same. Note that this may
## necessitate rewriting all functions that check for nrows. Every chain,
## as long as they're not literally all in a line, they should be able to 
## assume a closed polygonal chain.
##
## 2. Finish interior angle calculation function.
##
## 3. Finish angle portion of jitter function.

#' Generate a random polygonal chain
#'
#' @description Uniformly generates the vertices of a polygonal chain, sorting
#' the vertices by angle sweep from the x-axis.
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
  chain[, "x"] <- runif(k, min = min, max = max)
  chain[, "y"] <- runif(k, min = min, max = max)

  # Compute the centroid of the chain
  centroid <- rowSums(t(chain)) / nrow(chain)
  # Convert Cartesian to polar coordinates for easier handling of angles
  chain_p <- matrix(nrow = k, ncol = 2, byrow = FALSE)
  colnames(chain_p) <- c("r", "t")
  for (n in seq_len(nrow(chain))) {
    # Draw the vector between the current point and the centroid
    dist <- centroid - chain[n, ]
    # Take the Euclidean distance between the points
    r <- sqrt(dist[1] ^ 2 + dist[2] ^ 2)
    # Calculate angle sweep from the x-axis
    t <- atan2(dist[2], dist[1]) * (180 / pi)
    chain_p[n, "r"] <- r
    chain_p[n, "t"] <- t
  }
  # Sort vertices by increasing angle sweep from x-axis
  chain <- chain[order(chain_p[, "t"]), ]
  
  return(chain)
}

#' Validate whether a polygonal chain is closed
#'
#' @description
#' Determines if the given chain of coordinates is a polygonal chain or not
#'
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#'
#' @return Boolean value indicating whether the given chain is a polygonal
#' chain.
validate <- function(chain) {
  # Order of vertices must remain the same and might be disturbed after jitter
  #if (chain[1, ] == shape[nrow(chain), ]) {
  if (identical(chain[1, ], chain[1, ])) {
    is_chain <- TRUE
  } else {
    is_chain <- FALSE
  }
  return(is_chain)
}

#' Calculate the internal angles of a closed polygonal chain
#'
#' @description
#' If the given chain of coordinates is a closed polygonal chain, 
#' return a vector of its internal angles.
#'
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#'
#' @return A vector of length k containing the internal angles of the 
#' vertices of the polygonal chain.
get_internal_angles <- function(chain) {
  if (validate(chain)) {
    require(pracma) # For calculating the vector cross product
    angles <- matrix(nrow = nrow(chain), ncol = 1) # Extracted angles
    
    for (i in seq_len(nrow(chain))) {
      j <- i + 1
      
      if (i == nrow(chain)) {
        j <- 1 # Loop back to first vertex once end of chain is reached
      }
      
      #angles[i] <- Pi + atan2(V[i] %*% V[i + 1], V[i] * V[i + 1])
      
      v1 <- chain[i, ] - chain[nrow(chain), ] # Incoming vector
      v2 <- chain[j, ] - chain[i, ] # Outgoing vector
      
      angles[i] <- pi + atan2(cross(c(v1, 0), c(v2, 0))[3], c(v1) %*% c(v2))
    }
    
    angles <- angles * (180 / pi) # Convert angles to degrees
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(angles)
}

#' Add jitter to a polygonal chain
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
jitter <- function(chain, random = c("vertices", "angles"), factor = 0.01) {
  if (random == "vertices") {
    loops <- 0 # Jitter loop counter
    max_loops <- 10 # Maximum jitter loops
    repeat {
      # Calculate the perimeter of the whole chain via Euclidean distance
      perimeter <- 0
      for (i in seq_len(nrow(chain))) {
        j <- i + 1
        if (i == nrow(chain)) {
          j <- 1 # Loop back to first vertex once end of chain is reached
        }
        perimeter <- perimeter + sqrt(sum((chain[i, ] - chain[j, ])^2))
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
      if (validate(chain) != FALSE || loops == max_loops) {
        if (loops == max_loops) {
          warning("Output is not a valid polygonal chain. Please lower the jitter factor.")
        }
        break
      } else {
        loops <- loops + 1
      }
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
  # Horizontal translation
  chain[, "x"] <- chain[, "x"] + x
  # Vertical translation
  chain[, "y"] <- chain[, "y"] + y
  return(chain)
}

#' Dilate a polygonal chain
#'
#' @description
#' Scales the size of a polygonal chain to be greater or smaller, centered
#' about its centroid.
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
    # Calculate the centroid vector
    centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
    chain <- (chain - centroid) * factor + centroid
  }
  return(chain)
}

#' Flip a polygonal chain
#'
#' @description
#' Inverts a polygonal chain vertically or horizontally.
#'
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param direction String containing direction in which to flip the chain.
#'
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the flipped chain.
#'
flip <- function(chain, direction = c("horizontal", "vertical")) {
  centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
  if (direction == "horizontal") {
    # Flip across the y-axis
    chain <- (chain - centroid) %*%
      matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) + centroid
  } else if (direction == "vertical") {
    # Flip across the x-axis
    chain <- (chain - centroid) %*%
      matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) + centroid
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
  # Convert argument to radians, as R's trigonometric functions use radians
  angle <- angle * (pi / 180)
  # Calculate centroid vector
  centroid <- matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))]
  if (clockwise) {
    rotation <- matrix(c(cos(angle), -sin(angle),
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- rotation %*% (t(chain) - centroid) + centroid
  } else {
    rotation <- matrix(c(cos(angle), sin(angle),
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    chain <- (t(chain) - centroid) %*% rotation + centroid
  }
  # Transpose chain to original form
  return(t(chain))
}

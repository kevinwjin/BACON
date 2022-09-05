# Probabilistic polygonal chain generation and manipulation
# Author: Kevin Jin

#' Generate a random closed polygonal chain
#'
#' @description Randomly generates the vertices of a closed polygonal chain, 
#' sorting the vertices by angle sweep from the x-axis.
#'
#' @param k Number of vertices to generate.
#' @param min Minimum parameter for the uniform distribution.
#' @param max Maximum parameter for the uniform distribution.
#'
#' @return A 2 x (k + 1) matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
generate <- function(k = 3, min = 0, max = 1) {
  chain <- matrix(nrow = k, ncol = 2, byrow = FALSE)
  colnames(chain) <- c("x", "y")
  chain[, "x"] <- runif(k, min = min, max = max)
  chain[, "y"] <- runif(k, min = min, max = max)

  # Compute the centroid of the chain
  centroid <- rowSums(t(chain)) / nrow(chain)
  
  # Sort vertices by decreasing angle
  sorting_angles <- matrix(nrow = nrow(chain), ncol = 1)
  for (n in seq_len(nrow(chain))) {
    dist <- chain[n, ] - centroid
    sorting_angles[n] <- atan2(dist[2], dist[1])
  }
  chain <- chain[order(sorting_angles, decreasing = TRUE), ]
  
  # Repeat first row at end to form closed chain
  chain <- rbind(chain, chain[1, ])
  
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
  if (identical(chain[1, ], chain[nrow(chain), ])) {
    is_chain <- TRUE
  } else {
    is_chain <- FALSE
  }
  return(is_chain)
}

#' Calculate the interior angles of a closed polygonal chain
#'
#' @description
#' If the given chain of coordinates is a closed polygonal chain, 
#' return a vector of its interior angles.
#'
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#'
#' @return A vector of length k containing the interior angles of the 
#' vertices of the polygonal chain.
get_interior_angles <- function(chain) {
  if (validate(chain)) {
    # Represent sides of the polygonal chain as vectors
    vectors <- matrix(nrow = nrow(chain), ncol = 2)
    for (i in seq_len(nrow(chain))) { # Loop over vertices clockwise
      j <- i + 1 # i = initial point; j = terminal point 
      if (i == nrow(chain)) {
        j <- 1 
      }
      vectors[i, ] <- chain[j, ] - chain[i, ]
    }
    
    # Extract the interior angles of the polygonal chain
    angles <- matrix(nrow = nrow(chain), ncol = 1)
    for (i in seq_len(nrow(chain))) {
      j <- i + 1 # i = incoming vector; j = outgoing vector
      if (i == nrow(chain)) { # Return to vertex 1 once end of chain is reached
        j <- 1 
      }
      # Calculate angle between incoming and outgoing vectors
      #angles[i] <- pi + atan2(v1[1] * v2[2] - v2[1] * v1[2],
      #                        v1[1] * v2[1] + v1[2] * v2[2])
      
      angles[i] <- pi + atan2(chain[i, 1] * chain[j, 2] - chain[j, 1] * chain[i, 2],
                              chain[i, 1] * chain[j, 1] + chain[i, 2] * chain[j, 2])
    }
    angles <- c(t(angles * (180 / pi))) # Convert to degrees and convert to vector
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(angles)
}

#' Add jitter to a polygonal chain
#'
#' @description Randomizes the vertices or angles of a polygonal chain. If
#' output is not a valid polygonal chain, re-jitter a maximum of 10 times 
#' before giving up.
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
      if (validate(chain) == FALSE || loops == max_loops) {
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

#' Reflect a polygonal chain
#'
#' @description
#' Reflects a polygonal chain vertically or horizontally.
#'
#' @param chain A 2 x k matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param direction String containing direction in which to reflect the chain.
#'
#' @return A 2 x k matrix containing the x-y coordinates of the vertices of
#' the reflected chain.
#'
reflect <- function(chain, direction = c("horizontal", "vertical")) {
  centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
  if (direction == "horizontal") {
    # Reflect across the y-axis
    chain <- (chain - centroid) %*%
      matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) + centroid
  } else if (direction == "vertical") {
    # REflect across the x-axis
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
  # Convert argument angle to radians, as R's trigonometric functions use radians
  angle <- angle * (pi / 180)
  
  # Eliminate repeated row for now
  chain <- chain[-nrow(chain), ]
  
  # Calculate centroid vector
  centroid <- matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))]
  if (clockwise) {
    # Calculate clockwise rotation matrix
    rotation <- matrix(c(cos(angle), -sin(angle),
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
  } else {
    # Calculate counter-clockwise rotation matrix
    rotation <- matrix(c(cos(angle), sin(angle),
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
  }

  # Add repeated row back
  chain <- t(chain)
  chain <- rbind(chain, chain[1, ])
  
  # Transpose chain to original form
  return(chain)
}

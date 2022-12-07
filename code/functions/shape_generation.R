# Functions for random polygonal chain generation and manipulation
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

#' Check whether a polygonal chain is closed
#'
#' @description
#' Determines if the given chain of coordinates is a polygonal chain or not
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A logical value indicating whether the given chain is a polygonal
#' chain.
is_closed <- function(chain) {
  # Order of vertices must remain the same and might be disturbed after jitter
  if (identical(chain[1, ], chain[nrow(chain), ])) {
    is_closed <- TRUE
  } else {
    is_closed <- FALSE
  }
  return(is_closed)
}

#' Reverse the orientation of a closed polygonal chain
#'
#' @description
#' Reverse the orientation, or matrix row order, of a closed polygonal chain.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return The same polygonal chain, except now in reverse orientation.
reverse_orientation <- function(chain) {
  if (is_closed(chain)) {
    # Reverse clockwise orientation
    chain <- chain[nrow(chain):1, ]
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(chain)
}

#' Calculate the sum of the interior angles of a closed polygonal chain
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
  if (is_closed(chain)) {
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
  if (is_closed(chain)) {
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

#' Reproduce the orientation of a a polygonal chain relative interior angles 
#' and relative side lengths
#'
#' @description
#' Given two numerical vectors of compositional data containing relative 
#' interior angles and relative side lengths, return a vector of length k
#' containing which directions to turn when redrawing the unit polygonal chain. 
#'
#' @param data A numeric vector of length k containing data.
#' @param sum A numeric containing the sum of the vector.
#'
#' @return A boolean vector of length k containing the directions to turn when
#' redrawing the polygonal chain.
turning_orientation <- function(angles, side_lengths) {
  v1 <- 0
  v2 <- side_lengths[1]
  m <- (v2y - v1y) / (v2x - v1x)
  b <- 0
  
  y <- v3y
  x <- v3x
  # The first direction will always be left (FALSE) due to anticlockwise
  # orientation.
  orientation[1] <- FALSE
  
  for (i in seq_len(angles)) {
    if (y > mx + b) {
      orientation[i] <- FALSE # Turn left
    } else {
      orientation[i] <- TRUE # Turn right
    }
  }
  return(orientation)
}

#' Reproduce the unit polygonal chain from relative interior angles and relative 
#' side lengths
#'
#' @description
#' Given two numerical vectors of compositional data containing relative 
#' interior angles and relative side lengths, redraw the unit polygonal chain, 
#' centered at the origin.
#'
#' @param angles A numeric vector of length n containing interior angles.
#' @param side_lengths A numeric vector of length n containing relative side lengths.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the unit polygonal chain.
unit_chain <- function(angles, side_lengths) {
  if (length(angles) == length(side_lengths)) {
    # x-y coordinates of the vertices of the unit chain, centered at the origin
    chain <- matrix(0, nrow = length(angles) + 1, ncol = 2, byrow = TRUE)
    colnames(chain) <- c("x", "y") # For the sf function
    
    # Calculate total interior angle for chain in degrees
    total_angle <- (length(angles) - 2) * 180
    
    # Try calculating sides one by one at first
    side <- side_lengths[1] # Draw side 1 vector (horizontal along the x-axis)
    chain[2, ] <- c(side, 0)
    
    angle <- total_angle * angles[2] * (pi / 180) # Calculate angle 1-2 in radians
    angle <- pi - angle # Sweep anticlockwise from the x-axis
    
    # Rotate side 1 vector anticlockwise and translate into side 2
    rotation <- matrix(c(cos(angle), -sin(angle),
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    side <- t(rotation %*% chain[2, ])
    side[1, 1] <- side[1, 1] + chain[2, 1] # Currently, side 1 = 2, which is not right... vertex 2's coordinates might be wrong 
    
    chain[3, ] <- side
    
    # Iterate over all sides
    for (side in seq_len(length(angles))) {
      angle <- total_angle * angles[side + 1]
    }
    
    
    
    # for (i in seq_len(length(angles))) {

      
      # Rotate the chain
      side <- side_lengths[2]
      chain[3, ] <- c(rotation %*% chain[2, ]) + c(side, 0)


      # Take second angle (the second vertex)
      angle <- sum * angles[3] # Calculate current interior angle
      angle <- angle * (pi / 180) # R's trigonometric functions use radians
      angle <- pi - angle
      rotation <- matrix(c(cos(angle), -sin(angle),
                           sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
      side <- side_lengths[3]
      chain[4, ] <- c(rotation %*% chain[3, ]) + c(side, 0)

      # Take third angle (the third vertex)
      angle <- sum * angles[4] # Calculate current interior angle
      angle <- angle * (pi / 180) # R's trigonometric functions use radians
      rotation <- matrix(c(cos(angle), -sin(angle),
                           sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
      side <- side_lengths[4]
      chain[5, ] <- c(rotation %*% chain[2, ]) + c(side, 0)

      # Take fourth angle (the fourth vertex)
      angle <- sum * angles[5] # Calculate current interior angle
      angle <- angle * (pi / 180) # R's trigonometric functions use radians
      rotation <- matrix(c(cos(angle), -sin(angle),
                           sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
      side <- side_lengths[5]
      chain[3, ] <- c(rotation %*% chain[2, ]) + c(side, 0)

    # }
    
    # Repeat first row at end to form closed chain
    # chain <- rbind(chain, chain[1, ])
    
  } else {
    stop("Angle and side length vectors differ in length.")
  }
  
  return(chain)
}

#' Add jitter to a polygonal chain
#'
#' @description Randomizes the vertices or angles of a polygonal chain. If
#' output is not a valid polygonal chain, re-jitter a maximum of 10 times 
#' before giving up.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#' @param random String containing polygon parameter to randomize.
#' @param factor Floating point number from 0 to 1 exclusive as a percentage
#' of the perimeter of the polygonal chain to randomize by.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices 
#' of the polygonal chain.
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
      if (is_closed(chain) == FALSE || loops == max_loops) {
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

#' Translate a closed polygonal chain
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

# Shape generation functions
# Author: Kevin Jin

generate_polygon <- function(r = 1, # radius
                             N = 3, # number of total vertices
                             theta = 0, # angle of rotation in radians
                             x_center = 0, # x-coordinate center
                             y_center = 0) { # y-coordinate center
  polygon <- matrix(nrow = N, ncol = 2)
  colnames(polygon) <- c("x", "y")
  for (n in 1:N) {
    polygon[n, "x"] <- r * cos(2 * pi * n / N + theta) + x_center 
    polygon[n, "y"] <- r * sin(2 * pi * n / N + theta) + y_center
  }
  return(polygon)
}
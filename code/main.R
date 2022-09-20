# Execution script
# Author: Kevin Jin

#### Preparation ####
source("~/Documents/Programming/Repositories/bsclust/code/shape_generation.R")

#### Test shape generation ####
k <- 50
chain <- generate(k = k)
plot(chain, type = "l")
text(chain, labels = 1:nrow(chain)) # Label vertices in order
get_interior_angles(chain)

#### Generate data: Two different shapes, 30 each with slight variations ####
## Step 1: Replicate two shapes 5 times each and add jitter
n <- 60 # Total number of shapes
k <- 5 # Number of vertices
z <- 2 # Pre-specified number of clusters

# Generate first shape
shape_1 <- generate(k = k)
# Replicate n / 2 times (half of the total)
shapes_1 <- do.call(rbind, replicate(n / 2, shape_1, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / 2 * (k + 1)),
              by = (k + 1))) {
  shapes_1[i:(i + k), ] <- jitter(shapes_1[i:(i + k), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
# Generate second shape (z = 2 clusters)
shape_2 <- generate(k = k)
# Replicate n / 2 times (half of the total)
shapes_2 <- do.call(rbind, replicate(n / 2, shape_2, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / 2 * (k + 1)), 
              by = (k + 1))) {
  shapes_2[i:(i + k), ] <- jitter(shapes_2[i:(i + k), ],
                                        random = c("vertices"),
                                        factor = 0.01)
}
# Combine into 1 matrix
shapes <- rbind(shapes_1, shapes_2)

# Plot shapes - draw one shape per plot
for (i in seq(1, (n / 2 * (k + 1)), 
              by = (k + 1))) {
  plot(shapes[i:(i + k), ], type = "l")
}

## Step 2: Extract the angle vectors of the shapes (data to be clustered)
p <- k # Number of categories (length of angle vector)
x <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (angle vectors)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  x[row, ] <- get_interior_angles(shapes[i:(i + k), ])
  row <- row + 1
}

## Step 3: Extract side lengths of the shapes (data to be clustered)
s <- matrix(nrow = n, ncol = p, byrow = TRUE) # Data (side length vectors)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  s[row, ] <- get_side_lengths(shapes[i:(i + k), ])
  row <- row + 1
}

#### Load ADHD-200 data ####
# Test data extraction
chain <- plg.chains[[1]] # Unlist polygonal chain
colnames(chain) <- c("x", "y") # Rename columns, otherwise sf will complain
k <- nrow(chain) - 1 # Number of vertices
plot(chain, type = "l")
text(chain, labels = 1:nrow(chain))
get_interior_angles(chain)
get_side_lengths(chain)

# Load all chains in the folder
setwd("~/Documents/Programming/Repositories/CAPoly/data/adhd_200/shapes/polygonal")
file_names <- list.files(pattern = "*.Rdata", full.names = TRUE)
chains <- lapply(file_names, function(x) {
  load(file = x)
  mget(ls()[ls()!= "filename"])
})

# Construct matrices for data
angles <- matrix(nrow = length(chains), ncol = 50, byrow = TRUE)
side_lengths <- matrix(nrow = length(chains), ncol = 50, byrow = TRUE)

# Extract data
for (i in 1:length(chains)) {
  chain <- chains[[i]]$plg.chains$`1` # Extract chain
  colnames(chain) <- c("x", "y")
  angles[i, ] <- get_interior_angles(chain)
  side_lengths[i, ] <- get_side_lengths(chain)
}

# Name data rows and columns
names <- sapply(1:length(chains), function(i) chains[[i]]$id[2])
rownames(angles) <- names
rownames(side_lengths) <- names
colnames(angles) <- 1:50
colnames(side_lengths) <- 1:50

#### Load MPEG-7 data ####
library(SAFARI) # Image processing
library(dplyr) # Data handling
library(parallel)

# Load all images in the folder
# Retrieve list of all images (Set directory to image folder
setwd("~/Documents/Programming/Repositories/SAFARI-cluster-analysis/data/MPEG-7/images")
file_list <- dir(pattern = "gif$") # Choose appropriate image extension

# Convert to polygonal chains
extract_chains <- function(img) {
  if (img == "Glas-12.gif") { # Special case of inverted image in MPEG-7
    this_img <- read.image(img, invert = TRUE)
  } else {
    this_img <- read.image(img)
  }
  img_segs <- binary.segmentation(this_img, id = c(""), filter = 150, k = 3,
                                  categories = c("geometric", "boundary",
                                                 "topological"))
  chains <- data.frame(img_segs$plg.chains$`1`)
  return(chains)
}

# Process all images in parallel
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, varlist = c("read.image", "binary.segmentation"))
chains <- parSapply(cl, file_list, FUN = extract_chains)
stopCluster(cl)

# Construct matrices for data
angles <- matrix(nrow = nrow(chains), 
                 ncol = 100, # Uncertain maximum number of vertices
                 byrow = TRUE)
side_lengths <- matrix(nrow = nrow(chains), 
                       ncol = 100, 
                       byrow = TRUE)

# Extract data
chains <- t(chains)
for (i in 1:length(chains)) {
  chain <- cbind(unlist(chains[i, "X"]), unlist(chains[i, "Y"]))
  chain <- chain[seq(1, dim(chain)[1], length.out = (50 + 1)), ] # Subsample
  colnames(chain) <- c("x", "y")
  angles[i, ] <- get_interior_angles(chain)
  side_lengths[i, ] <- get_side_lengths(chain)
}

# Name data rows and columns
rownames(angles) <- rownames(chains)
rownames(side_lengths) <- rownames(chains)
colnames(angles) <- 1:3000
colnames(side_lengths) <- 1:3000

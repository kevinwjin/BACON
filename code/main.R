# Execution script
# Author: Kevin Jin

#### Load data simulation functions ####
source("~/Documents/Programming/Repositories/bsclust/code/shape_generation.R")

#### Test data generation ####
k <- 20
chain <- generate(k = k)
plot(chain, type = "l")
text(chain, labels = 1:nrow(chain)) # Label vertices in order
get_interior_angles(chain)

#### Generate simulated data ####
## Step 1: Replicate two shapes 5 times each and add jitter
n <- 120 # Total number of shapes
k <- 10 # Number of vertices
z <- 4 # Pre-specified number of clusters

# Generate first shape
shape_1 <- generate(k = k)
# Replicate n / z times
shapes_1 <- do.call(rbind, replicate(n / z, shape_1, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_1[i:(i + k), ] <- jitter(shapes_1[i:(i + k), ], random = c("vertices"),
                                        factor = 0.01)
}

# Generate second shape
shape_2 <- generate(k = k)
# Replicate n / z times
shapes_2 <- do.call(rbind, replicate(n / z, shape_2, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_2[i:(i + k), ] <- jitter(shapes_2[i:(i + k), ], random = c("vertices"),
                                        factor = 0.01)
}

# Generate third shape
shape_3 <- generate(k = k)
# Replicate n / z times
shapes_3 <- do.call(rbind, replicate(n / z, shape_3, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_3[i:(i + k), ] <- jitter(shapes_3[i:(i + k), ], random = c("vertices"),
                                  factor = 0.01)
}

# Generate fourth shape
shape_4 <- generate(k = k)
# Replicate n / z times
shapes_4 <- do.call(rbind, replicate(n / z, shape_4, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_4[i:(i + k), ] <- jitter(shapes_4[i:(i + k), ], random = c("vertices"),
                                  factor = 0.01)
}

# Combine into 1 matrix
shapes <- rbind(shapes_1, shapes_2, shapes_3, shapes_4)

# Plot shapes - draw one shape per plot
for (i in seq(1, (n / 2 * (k + 1)), 
              by = (k + 1))) {
  plot(shapes[i:(i + k), ], type = "l")
}

## Step 2: Extract the angle vectors of the shapes (data to be clustered)
p <- k # Number of categories (length of angle vector)
angles <- matrix(nrow = n, ncol = p, byrow = TRUE)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  angles[row, ] <- get_interior_angles(shapes[i:(i + k), ])
  row <- row + 1
}

## Step 3: Extract side lengths of the shapes (data to be clustered)
side_lengths <- matrix(nrow = n, ncol = p, byrow = TRUE)
row <- 1 # Row counter
for (i in seq(1, n * (k + 1), by = (k + 1))) {
  side_lengths[row, ] <- get_side_lengths(shapes[i:(i + k), ])
  row <- row + 1
}

## Step 3: Data analysis
# Create cluster labels
clusters <- c()
for (i in 1:z) {
  clusters <- append(clusters, rep.int(i, times = 30))
}

# PCA dimensional reduction
require(ggfortify)
require(compositions)
angles <- cbind(angles, clusters)
angles <- as.data.frame(angles)
autoplot(stats::prcomp(clr(angles[, 1:k])), 
         data = angles,
         colour = factor(clusters),
         main = 'PCA of Simulated Shapes Relative Interior Angles')
#biplot(stats::prcomp(angles[, 1:k]))
#plot(stats::prcomp(angles[, 1:k]), type = 'l')

side_lengths <- cbind(side_lengths, clusters)
side_lengths <- as.data.frame(side_lengths)
autoplot(stats::prcomp(clr(side_lengths[, 1:k])), 
         data = angles, 
         colour = factor(clusters),
         main = 'PCA of Simulated Shapes Relative Side Lengths')
#biplot(stats::prcomp(side_lengths[, 1:k]))
#plot(stats::prcomp(side_lengths[, 1:k], type = 'l'))

## Test reflex angle counterexample
chain <- rbind(c(0, 0), c(5, 0), c(5, 5), c(2.5, 0.5), c(0, 5), c(0, 0))
get_interior_angles(chain) # Remove normalization first

#### Load MPEG-7 data ####
require(SAFARI) # Image processing
require(dplyr) # Data handling
require(parallel)

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


#### Load ADHD-200 data ####
# Load all chains in the folder
setwd("~/Documents/Programming/Repositories/CAPoly/data/adhd_200/shapes/polygonal")
file_names <- list.files(pattern = "*.Rdata", full.names = TRUE)
chains <- lapply(file_names, function(x) {
  load(file = x)
  mget(ls()[ls()!= "filename"])
})

# Test data extraction
chain <- chains[[1]]$plg.chains$`1` # Unlist polygonal chain
colnames(chain) <- c("x", "y") # Rename columns, otherwise sf will complain
plot(chain, type = "l")
text(chain, labels = 1:nrow(chain))

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


#### ADHD-200 data analysis ####
# Load clinical information (ground truth cluster labels)
# gender (discrete, binary): 0 = female; 1 = male
# age (continuous): rounded and discretized by factor()
# phenotype (discrete): 0 = typically developing children; 1 = ADHD-combined
#                       2 = ADHD-hyperactive/impulsive; 3 = ADHD-inattentive

gender <- clinical_info$gender[order(clinical_info$patient.id)]
age <- factor(round(clinical_info$age[order(clinical_info$patient.id)]))
phenotype <- clinical_info$phenotype[order(clinical_info$patient.id)]

# Centered & isometric log-ratio transformation
require(compositions)
angles.clr.pca <- prcomp(clr(angles), scale = TRUE)
autoplot(angles.clr.pca)

angles.ilr.pca <- prcomp(ilr(angles), scale = TRUE)
autoplot(angles.ilr.pca)

# PCA dimensional reduction on first 50 columns
require(ggfortify)
angles <- cbind(angles, phenotype, age, gender)
angles <- as.data.frame(angles)
autoplot(stats::prcomp(clr(angles[, 1:50])), 
         data = angles,
         colour = factor(phenotype),
         main = 'PCA of ADHD-200 Relative Interior Angles (Phenotype)')
autoplot(stats::prcomp(clr(angles[, 1:50])), 
         data = angles, 
         colour = 'age',
         main = 'PCA of ADHD-200 Relative Interior Angles (Age)')
autoplot(stats::prcomp(clr(angles[, 1:50])), 
         data = angles, 
         colour = factor(gender + 1),
         main = 'PCA of ADHD-200 Relative Interior Angles (Gender)')
biplot(stats::prcomp(angles[, 1:50]))
plot(stats::prcomp(angles[, 1:50]), type = 'l')

side_lengths <- cbind(side_lengths, phenotype, age, gender)
side_lengths <- as.data.frame(side_lengths)
autoplot(stats::prcomp(clr(side_lengths[, 1:50])), 
         data = angles, 
         colour = factor(phenotype),
         main = 'PCA of ADHD-200 Relative Side Lengths (Phenotype)')
autoplot(stats::prcomp(clr(side_lengths[, 1:50])), 
         data = angles, 
         colour = 'age',
         main = 'PCA of ADHD-200 Relative Side Lengths (Age)')
autoplot(stats::prcomp(clr(side_lengths[, 1:50])), 
         data = angles, 
         colour = factor(gender + 1),
         main = 'PCA of ADHD-200 Relative Side Lengths (Gender)')
biplot(stats::prcomp(side_lengths[, 1:50]))
plot(stats::prcomp(side_lengths[, 1:50], type = 'l'))


#### Bijective mapping for angle and side lengths

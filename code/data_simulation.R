# Data simulation and generation
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
n <- 200 # Total number of shapes
k <- 20 # Number of vertices
z <- 10 # Pre-specified number of clusters

# Generate first shape
shape_1 <- generate(k = k)
# Replicate n / z times
shapes_1 <- do.call(rbind, replicate(n / z, shape_1, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_1[i:(i + k), ] <- jitter(shapes_1[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate second shape
shape_2 <- generate(k = k)
# Replicate n / z times
shapes_2 <- do.call(rbind, replicate(n / z, shape_2, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_2[i:(i + k), ] <- jitter(shapes_2[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate third shape
shape_3 <- generate(k = k)
# Replicate n / z times
shapes_3 <- do.call(rbind, replicate(n / z, shape_3, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_3[i:(i + k), ] <- jitter(shapes_3[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate fourth shape
shape_4 <- generate(k = k)
# Replicate n / z times
shapes_4 <- do.call(rbind, replicate(n / z, shape_4, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_4[i:(i + k), ] <- jitter(shapes_4[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate fifth shape
shape_5 <- generate(k = k)
# Replicate n / z times
shapes_5 <- do.call(rbind, replicate(n / z, shape_5, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_5[i:(i + k), ] <- jitter(shapes_5[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate sixth shape
shape_6 <- generate(k = k)
# Replicate n / z times
shapes_6 <- do.call(rbind, replicate(n / z, shape_6, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_6[i:(i + k), ] <- jitter(shapes_6[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate seventh shape
shape_7 <- generate(k = k)
# Replicate n / z times
shapes_7 <- do.call(rbind, replicate(n / z, shape_7, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_7[i:(i + k), ] <- jitter(shapes_7[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate eighth shape
shape_8 <- generate(k = k)
# Replicate n / z times
shapes_8 <- do.call(rbind, replicate(n / z, shape_8, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_8[i:(i + k), ] <- jitter(shapes_8[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate ninth shape
shape_9 <- generate(k = k)
# Replicate n / z times
shapes_9 <- do.call(rbind, replicate(n / z, shape_9, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_9[i:(i + k), ] <- jitter(shapes_9[i:(i + k), ], random = c("vertices"),
                                  factor = 0.05)
}

# Generate tenth shape
shape_10 <- generate(k = k)
# Replicate n / z times
shapes_10 <- do.call(rbind, replicate(n / z, shape_10, simplify = FALSE))
# Apply jitter (each k + 1 rows is a shape)
for (i in seq(1, (n / z * (k + 1)), by = (k + 1))) {
  shapes_10[i:(i + k), ] <- jitter(shapes_10[i:(i + k), ], random = c("vertices"),
                                   factor = 0.05)
}


# Combine into 1 matrix
shapes <- rbind(shapes_1, shapes_2, shapes_3, shapes_4, shapes_5, shapes_6,
                shapes_7, shapes_8, shapes_9, shapes_10)

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
  clusters <- append(clusters, rep.int(i, times = (n / z)))
}
angles <- cbind(angles, clusters)
side_lengths <- cbind(side_lengths, clusters)

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
file_names <- dir(pattern = "gif$") # Choose appropriate image extension

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
chains <- parSapply(cl, file_names, FUN = extract_chains)
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


#### Convert MPEG-7 to polygonal chains ####
angles <- matrix(nrow = length(mpeg_7_twenty), ncol = 20, byrow = TRUE)
side_lengths <- matrix(nrow = length(mpeg_7_twenty), ncol = 20, byrow = TRUE)

for(i in 1:length(mpeg_7_twenty)) {
  if (i == 824) next # Temporarily skip, chain looks weird for some reason
  chain <- mpeg_7_twenty[[i]]
  colnames(chain) <- c("x", "y")
  angles[i, ] <- get_interior_angles(chain)
  side_lengths[i, ] <- get_side_lengths(chain)
}

#### Load ADHD-200 data ####
# Load all chains in the folder
setwd("~/Documents/Programming/Repositories/BACON/data/adhd_200/shapes/polygonal")
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


#### Bijection from data to coordinates ####

#### Create binary images from ADHD-200 chains ####
## Test conversion of 1 chain
adhd <- matrix(unlist(plg.chains), ncol = 2)
png(filename = "ADHD-200_1.png", width = 256, height = 256, units = "px")
par(bg = 'black', fg = 'white')
plot(adhd, type = 'l', xaxt = 'n', yaxt = 'n', ann = FALSE, frame.plot = FALSE,
     asp = 1)
polygon(adhd, col = "white")
dev.off()

library(EBImage)
img <- readImage("ADHD-200_1.png")
img <- channel(img, "gray")
img <- img > otsu(img, levels = 256)
display(img)
writeImage(img, files = "ADHD-200_1_binary.png")

library(SAFARI)
library(parallel)
img <- read.image(file = "ADHD-200_1_binary.png")
data <- binary.segmentation(img, 
                            id = "adhd", 
                            filter = 150, 
                            k = 3, 
                            categories = c("geometric", 
                                           "boundary", 
                                           "topological"))
data$desc
data$plg.chains

## Convert all chains
# Unlist entire list of chains and put it in one list
setwd("~/Documents/Programming/Repositories/BACON/data/adhd_200/shapes/polygonal")
file_names <- list.files(pattern = "*.Rdata", full.names = TRUE)
chains <- lapply(file_names, function(x) {
  load(file = x)
  mget(ls()[ls()!= "filename"])
})

# Loop through list of chains and convert all to binary images
setwd("~/Documents/Programming/Repositories/BACON/data/adhd_200/binary_images")

for (i in 1:(length(chains) - 1)) {
  chain <- chains[[i]]$plg.chains$`1` # Take 1 chain
  colnames(chain) <- c("x", "y")
  png(filename = paste0("ADHD-200_", i, ".png"), width = 256, height = 256, 
      units = "px")
  par(bg = 'black', fg = 'white')
  plot(chain, type = 'l', xaxt = 'n', yaxt = 'n', ann = FALSE, 
       frame.plot = FALSE, asp = 1)
  polygon(chain, col = "white")
  dev.off()
}

for (i in 1:(length(chains) - 1)) {
  require(EBImage)
  img <- readImage(paste0("ADHD-200_", i, ".png"))
  img <- channel(img, "gray")
  img <- img > otsu(img, levels = 256)
  writeImage(img, files = paste0("ADHD-200_", i, ".png"))
}

# Extract SAFARI features
require(SAFARI)
extract_features <- function(img) {
  this_img <- read.image(img)
  img_segs <- binary.segmentation(this_img, 
                                  id = "adhd", 
                                  filter = 150, 
                                  k = 3, 
                                  categories = c("geometric", 
                                                 "boundary", 
                                                 "topological"))
  features <- data.frame(img_segs[["desc"]][-c(1:2, 30)])
  features <- cbind(c(img), features)
  return(features)
}

require(parallel)
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, varlist = c("read.image", "binary.segmentation"))
file_names <- list.files(pattern = "*png", full.names = TRUE)
features <- parSapply(cl, file_names, FUN = extract_features)
stopCluster(cl)

require(dplyr)
features <- data.frame(t(features), row.names = 1)
features <- data.frame(lapply(features, unlist))
features_scaled <- features %>% mutate_all(scale)

# Extract ground truth
clinical_info <- read.csv("clinical_info.csv")
gender <- clinical_info$gender[order(clinical_info$patient.id)]
age <- factor(round(clinical_info$age[order(clinical_info$patient.id)]))
phenotype <- clinical_info$phenotype[order(clinical_info$patient.id)]

# k-means 
max_k <- 100
kmeans_mat <- matrix(nrow = max_k, ncol = length(file_names), byrow = TRUE)

# Generate k-means clustering from scaled features (number of holes is NaN, k-means doesn't work)
for (k in 1:max_k) {
  kmeans_mat[k, ] <- kmeans(
    x = features_scaled,
    centers = k
    # nstart = 20,
    # iter.max = 30
  )[["cluster"]]
}
kmeans_scaled <- as.data.frame(kmeans_mat)

# Generate k-means clustering from unscaled features
for (k in 1:max_k) {
  kmeans_mat[k, ] <- kmeans(
    x = features,
    centers = k
    # nstart = 20,
    # iter.max = 30
  )[["cluster"]] # Select cluster results
}
kmeans_unscaled <- as.data.frame(kmeans_mat)

# Calculate ARI values for k-means over k = 1:100
require(mclust)
kmeans_ari <- matrix(nrow = max_k, ncol = 2, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for k-means clustering of scaled features
  kmeans_ari[i, 1] <- adjustedRandIndex(
    unlist(kmeans_scaled[i, ]),
    phenotype
  )
  # ARI for k-means clustering of unscaled features
  kmeans_ari[i, 2] <- adjustedRandIndex(
    unlist(kmeans_unscaled[i, ]),
    phenotype
  )
}
kmeans_ari <- as.data.frame(kmeans_ari)
names(kmeans_ari) <- c("ARI_scaled", "ARI_unscaled")
k_values <- 1:max_k
kmeans_ari <- mutate(kmeans_ari, k_values = as.numeric(row.names(kmeans_ari)))

# Visualize k-means clustering accuracy
require(ggplot2)
ggplot(kmeans_ari, aes(x = k_values)) +
  geom_point(aes(
    y = ARI_scaled,
    color = 3
  )) +
  geom_smooth(aes(
    y = ARI_scaled,
    color = 3
  )) +
  # geom_point(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  # geom_line(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  geom_vline(
    xintercept = 70, # Ground truth
    color = "red"
  ) +
  labs(
    title = "k-means Clustering Accuracy (Data: Normalized ADHD-200 Shape Features)",
    x = "k-value",
    y = "Adjusted Rand Index"
    # color = "Features"
  ) +
  theme(legend.position = "none")
# scale_color_hue(labels = c("Scaled", "Unscaled")) +
# scale_x_continuous(breaks = seq(0, max_k, by = 10)) +
# scale_y_continuous(breaks = seq(0, 1, by = 0.1))

# Hierarchical
# Generate hierarchical clusters
hier_scaled_tree <- hclust(dist(features_scaled, method = "euclidean"),
                           method = "complete" # Complete linkage
)
hier_unscaled_tree <- hclust(dist(features, method = "euclidean"),
                             method = "complete"
)

max_k <- 100
hier_mat <- matrix(nrow = max_k, ncol = length(file_names), byrow = TRUE)

# Generate hierarchical clustering from scaled features
for (k in 1:max_k) {
  hier_mat[k, ] <- cutree(hier_scaled_tree, k = k)
}
hier_scaled <- as.data.frame(hier_mat)

# Generate hierarchical clustering from unscaled features
for (k in 1:max_k) {
  hier_mat[k, ] <- cutree(hier_unscaled_tree, k = k)
}
hier_unscaled <- as.data.frame(hier_mat)

# Calculate ARI values for hierarchical over k = 1:100 (Ground truth: k = 70)
hier_ari <- matrix(nrow = max_k, ncol = 2, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for hierarchical clustering of scaled features
  hier_ari[i, 1] <- adjustedRandIndex(
    unlist(hier_scaled[i, ]),
    phenotype
  )
  # ARI for hierarchical clustering of unscaled features
  hier_ari[i, 2] <- adjustedRandIndex(
    unlist(hier_unscaled[i, ]),
    phenotype
  )
}
hier_ari <- as.data.frame(hier_ari)
names(hier_ari) <- c("ARI_scaled", "ARI_unscaled")
k_values <- 1:max_k
hier_ari <- mutate(hier_ari, k_values = as.numeric(row.names(hier_ari)))

# Visualize hierarchical clustering accuracy
ggplot(hier_ari, aes(x = k_values)) +
  geom_point(aes(
    y = ARI_scaled,
    color = 3
  )) +
  geom_smooth(aes(
    y = ARI_scaled,
    color = 3
  )) +
  # geom_point(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  # geom_line(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  geom_vline(
    xintercept = 70, # Ground truth
    color = "red"
  ) +
  labs(
    title = "Hierachical Clustering Accuracy (Data: Normalized ADHD-200 Shape Features)",
    x = "k-value",
    y = "Adjusted Rand Index",
    # color = "Features"
  ) +
  theme(legend.position = "none")
# scale_color_hue(labels = c("Scaled", "Unscaled")) +
# scale_x_continuous(breaks = seq(0, max_k, by = 10)) +
# scale_y_continuous(breaks = seq(0, 1, by = 0.1))

# GMM
max_k <- 100
gmm_mat <- matrix(nrow = max_k, ncol = length(file_names), byrow = TRUE)

# Generate GMM clustering from scaled features
for (k in 1:max_k) {
  gmm_mat[k, ] <- Mclust(features_scaled, G = k)$classification
}
gmm_scaled <- as.data.frame(gmm_mat)

# Generate GMM clustering from unscaled features
for (k in 1:max_k) {
  gmm_mat[k, ] <- Mclust(features, G = k)$classification
}
gmm_unscaled <- as.data.frame(gmm_mat)

# Calculate ARI values for GMM over k = 1:100 (Ground truth: k = 70)
gmm_ari <- matrix(nrow = max_k, ncol = 2, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for GMM clustering of scaled features
  gmm_ari[i, 1] <- adjustedRandIndex(
    unlist(gmm_scaled[i, ]),
    phenotype
  )
  # ARI for GMM clustering of unscaled features
  gmm_ari[i, 2] <- adjustedRandIndex(
    unlist(gmm_unscaled[i, ]),
    phenotype
  )
}
gmm_ari <- as.data.frame(gmm_ari)
names(gmm_ari) <- c("ARI_scaled", "ARI_unscaled")
k_values <- 1:max_k
gmm_ari <- mutate(gmm_ari, k_values = as.numeric(row.names(gmm_ari)))

# Visualize GMM clustering accuracy
ggplot(gmm_ari, aes(x = k_values)) +
  geom_point(aes(
    y = ARI_scaled,
    color = 3
  )) +
  geom_smooth(aes(
    y = ARI_scaled,
    color = 3
  )) +
  # geom_point(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  # geom_line(aes(
  #   y = ARI_unscaled,
  #   color = "steelblue"
  # )) +
  geom_vline(
    xintercept = 70, # Ground truth
    color = "red"
  ) +
  labs(
    title = "GMM Clustering Accuracy (Data: Normalized ADHD-200 Shape Features)",
    x = "k-value",
    y = "Adjusted Rand Index",
    # color = "Features"
  ) +
  theme(legend.position = "none")
# scale_color_hue(labels = c("Scaled", "Unscaled")) +
# scale_x_continuous(breaks = seq(0, max_k, by = 10)) +
# scale_y_continuous(breaks = seq(0, 1, by = 0.1))

# Visualize accuracy of all clustering methods
accuracy <- right_join(kmeans_ari, gmm_ari, by = "k_values")
accuracy <- right_join(accuracy, hier_ari, by = "k_values")

names(accuracy) <- c(
  "kmeans_ari_scaled",
  "kmeans_ari_unscaled",
  "k_values",
  "gmm_ari_scaled",
  "gmm_ari_unscaled",
  "hier_ari_scaled",
  "hier_ari_unscaled"
)

ggplot(accuracy, aes(x = k_values)) +
  geom_smooth(aes(
    y = kmeans_ari_scaled,
    color = "darkred"
    # linetype = "solid"
  )) +
  geom_point(aes(
    y = kmeans_ari_scaled,
    color = "darkred"
  )) +
  # geom_line(aes(
  #   y = kmeans_ari_unscaled,
  #   color = "darkred",
  #   linetype = "twodash"
  # )) +
  geom_smooth(aes(
    y = hier_ari_scaled,
    color = "steelblue"
    # linetype = "solid"
  )) +
  geom_point(aes(
    y = hier_ari_scaled,
    color = "steelblue"
  )) +
  # geom_line(aes(
  #   y = hier_ari_unscaled,
  #   color = "steelblue",
  #   linetype = "twodash"
  # )) +
  geom_smooth(aes(
    y = gmm_ari_scaled,
    color = "seagreen"
    # linetype = "solid"
  )) +
  geom_point(aes(
    y = gmm_ari_scaled,
    color = "seagreen"
  )) +
  # geom_line(aes(
  #   y = gmm_ari_unscaled,
  #   color = "seagreen",
  #   linetype = "twodash"
  # )) +
  geom_vline(
    xintercept = 70,
    color = "red"
  ) +
  labs(
    title = "Accuracy of Several Clustering Methods (Data: Normalized ADHD-200 Shape Features)",
    x = "k-value",
    y = "Adjusted Rand Index"
  ) +
  # scale_x_continuous(breaks = seq(0, max_k, by = 10)) +
  # scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_hue(labels = c("k-means", "GMM", "Hierarchical"))
# scale_linetype(labels = c("Scaled", "Unscaled"))

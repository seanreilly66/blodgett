# R script for generating and plotting spatial patterns using spatstat

# Install and load spatstat (uncomment if not already installed)
# install.packages("spatstat")
library(spatstat)

# Set common parameters
win <- owin(c(0, 100), c(0, 100))
n_points <- 1830

# Set plot layout for patterns
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))

# 1. Complete Spatial Randomness (CSR)
csr_pattern <- runifpoint(n_points, win = win)
plot(csr_pattern, main = "Complete Spatial Randomness")

# 2. Evenly Spaced (Regular Grid)
grid_side <- ceiling(sqrt(n_points))
regular_coords <- expand.grid(x = seq(0, 100, length.out = grid_side),
                              y = seq(0, 100, length.out = grid_side))
regular_pattern <- as.ppp(regular_coords, W = win)
plot(regular_pattern, main = "Evenly Spaced")

# 3. Even S87tronger Clustering (Thomas Process)
set.seed(12345)
clustered_pattern <- rThomas(kappa = 0.001, scale = 6, mu = 250, win = win)
plot(clustered_pattern, main = "Clustering")

# Set plot layout for K functions
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))

# 1. CSR - K function (default)

plot(Kest(csr_pattern, correction = 'border'), main = "Complete Spatial Randomness", legend = FALSE)

# 2. Regular Grid - K function (default)
plot(Kest(regular_pattern, correction = 'border'), main = "Regular Grid", legend = FALSE)

# 3. Even Stronger Thomas - K function
plot(Kest(clustered_pattern, correction = 'border'), main = "Clustered", legend = FALSE)

# Reset plot layout
par(mfrow = c(1, 1))

# ==============================================================================
#
# Lidar tree segmentation
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 17 Sept 2024
# Last commit: 17 Sept 2024
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# ==============================================================================
#
# User inputs:
#
# ==============================================================================
#
# Package dependencies:
#
# ==============================================================================
#
# Known problems:
#
# ==============================================================================

library(lidR)
library(tidyverse)
library(sf)
library(glue)
library(doParallel)

# ================================= User inputs ================================

als_folder <- 'data/las/als_2022'
id <- 'bldgt_als.+hnorm_tin\\.las'

stem_map_file <- 'data/field_stem_maps/bldgt_field_stem_pts.shp'

# ==============================================================================
# ================================ Stem map eq =================================
# ==============================================================================

stem_map <- st_read(stem_map_file)

f_lm <- lm(crrad_m ~ h_m, data = stem_map)
summary(f_lm)

# ==============================================================================
# =============================== LMF Tree seg =================================
# ==============================================================================

als_files <- list.files(path = als_folder,
                        pattern = id,
                        full.names = T)

cl <- makeCluster(12)
registerDoParallel(cl)

foreach(
  file_c = als_files,
  .packages = c('lidR', 'tidyverse', 'sf', 'glue')
) %dopar% {
# ) %do% {
  
  # ===================================== ALS ====================================
  
  als_yr <- str_extract(file_c, '(?<=_als)[:digit:]+')
  
  comp <- str_extract(file_c, '(?<=_c)[:digit:]+')
  
  als_c <- readLAS(file_c)
  
  crs <- st_crs(als_c)
  
  # ==============================================================================
  # ========================= Testing tree segmentation ==========================
  # ==============================================================================

  kernel <- matrix(1, 3, 3)
  chm_c <- rasterize_canopy(als_c, 0.25, p2r(subcircle = 0.2), pkg = "terra") %>%
    terra::focal(w = kernel, fun = median, na.rm = TRUE)
  
  terra::writeRaster(
    chm_c, 
    glue('data/tree_seg/chm/als{als_yr}_c{comp}_chm.tif'),
    overwrite = T)
  
  # ------------------------------------------------------------------------------
  # -------------- LMF with lm function from height to radius ratio --------------
  # ------------------------------------------------------------------------------
  
  ws_f <- function(x) {x * f_lm$coefficients[2] + f_lm$coefficients[1]}
  ttop_c <- locate_trees(chm_c, lmf(ws_f))

  als_seg <- segment_trees(als_c, dalponte2016(chm_c, ttop_c, max_cr = 40))

  crown_poly_lmf <- crown_metrics(als_seg, func = .stdtreemetrics, geom = "convex")

  st_write(crown_poly_lmf,
           glue('data/tree_seg/lmf/als{als_yr}_c{comp}_lmf_chm.shp'),
           overwrite = T)
  
}

stopCluster(cl)

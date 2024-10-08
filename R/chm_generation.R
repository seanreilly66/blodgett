# ==============================================================================
#
# CHM generation
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 2 Oct 2024
# Last commit: 2 Oct 2024
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
library(glue)
library(doParallel)

# ================================= User inputs ================================

als_folder <- 'data/las'
start_pattern <- 'bldgt_als'

dtm_resolution <- 1
chm_resolution <- 0.5

subcircle_r <- 0.2

smooth <- TRUE

dtm_output <- 'data/dtm/{id}_dtm_tin_{res}.tif'
chm_output <- 'data/chm/{id}_chm_p2r_{res}.tif'

# ==============================================================================
# ============================== CHM generation ================================
# ==============================================================================

als_files <- list.files(
  als_folder,
  pattern = start_pattern,
  recursive = T,
  full.names = T
)

cl <- makeCluster(6)
registerDoParallel(cl)

output <- foreach(
  als_i = als_files,
  .packages = c('lidR', 'tidyverse', 'glue'),
  .export = c('start_pattern')
) %dopar% {
  
  # als_i = als_files [1]
  
  las <- als_i %>%
    readLAS() %>%
    filter_duplicates() %>%
    filter_poi(Classification != 7, Classification != 18) # Remove errors
  
  dtm <- rasterize_terrain(
    las,
    res = dtm_resolution,
    pkg = 'terra')
  
  chm <- (las - dtm) %>%
    rasterize_canopy(
      res = chm_resolution,
      algorithm = p2r(subcircle = subcircle_r),
      pkg = "terra"
    )
  
  if (smooth) {
    kernel <- matrix(1, 3, 3)
    chm <- terra::focal(chm,
                        w = kernel,
                        fun = median,
                        na.rm = TRUE)
    
  }
  
  id = str_remove(als_i, '.las$') |>
    str_remove(glue('[:graph:]+(?={start_pattern})'))
  
  chm_file <- glue(
    chm_output,
    res = chm_resolution |>
      as.character() |>
      str_remove('\\.')
  )
  
  dtm_file <- glue(
    dtm_output,
    res = dtm_resolution |>
      as.character() |>
      str_remove('\\.')
  )
  
  terra::writeRaster(chm, chm_file, overwrite = T)
  terra::writeRaster(dtm, dtm_file, overwrite = T)
  
  return(chm_file)

}

stopCluster(cl)

# ==============================================================================
# ==============================================================================
#
# Tree polygon metric extraction
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 18 Dec 2024
# Last commit: 18 Dec 2024
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

library(sf)
library(lidR)
library(tidyverse)
library(glue)
library(doParallel)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

crown_file <- 'data/tree_seg/bldgt_ffs_als22_merged_crowns.shp'

las_folder <- 'data/las/als_2022'
las_id <- 'bldgt_als22_.+_hnorm_tin\\.las$'

metric_func <- 'R/las_metric_function.R'

output <- 'data/poly_metrics/bldgt_ffs_als22_crown_metrics.shp'

# ==============================================================================
# ============================= Calculate metrics ============================== 
# ==============================================================================


crown_shp <- st_read(crown_file, quiet = T)

las_files <- list.files(
  las_folder, 
  pattern = las_id,
  full.names = T
) 

cl <- makeCluster(12)
registerDoParallel(cl)

metric <- foreach(
  las = las_files,
  .combine = 'rbind',
  .packages = c('tidyverse', 'sf', 'lidR')
) %dopar% {
  
  source(metric_func)
  
  compartment <- str_extract(las, '(?<=_c)[:digit:]+')
  
  shp <- crown_shp %>%
    filter(comp == compartment) %>%
    mutate(area = st_area(geometry))
  
  las <- readLAS(las)
  
  metric = lidR::plot_metrics(las, ~ las_cld_metrics(z = Z, r = ReturnNumber), shp)
  
  return(metric)
  
}

stopCluster(cl)

st_write(metric, output, append = F)

# ==============================================================================
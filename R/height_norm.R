# ==============================================================================
#
# ALS height normalization
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 12 Oct 2024
# Last commit: 12 Oct 2024
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

als_folder <- 'data/las'
id <- 'bldgt_als.+raw\\.las'

# ==============================================================================
# ================================= Data prep ==================================
# ==============================================================================

als_files <- list.files(als_folder, recursive = T, pattern = id, full.name = T) 

cl <- makeCluster(4)
registerDoParallel(cl)

foreach(
  als_i = als_files,
  .packages = c('lidR', 'tidyverse', 'sf', 'glue')
) %dopar% {
  
  las <- als_i %>%
    readLAS() %>%
    filter_duplicates() %>%
    filter_poi(Classification != 7, Classification != 18) %>% # Remove errors
    normalize_height(tin())
  
  file_out <- als_i %>% 
    str_replace('raw', 'hnorm_tin')
  
  writeLAS(las, file_out)
  
  return(NULL)
  
}

# ==============================================================================
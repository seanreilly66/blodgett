# ==============================================================================
#
# LAS catalog extract by polygon
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 30 July 2024
# Last commit: 30 July 2024
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# Generates las files for given polygon geometry (e.g. plot boundaries) from 
# catalog of las files (e.g. wall-to-wall als coverage)
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
library(sf)
library(tidyverse)
library(glue)

# ================================= User inputs ================================

# Files

las_folder <- 'data/las/als/raw/ctg'
file_ext <- '.laz$'
shp_file <- 'data/gis/BlodgettFFS_UnitBoundaries/FFS_COMPS_only_Alb83.shp'

# Output
file_skeleton <- 'data/las/als/bldgt_als_{id}.las'



# Convert laz to las

laz_files <- list.files(
  path = las_folder,
  pattern = '.laz$',
  full.names = T
)

# las_files <- list.files(
#   path = las_folder,
#   pattern = '.las$',
#   full.names = T
# ) 
# 
# prog <- laz_files[!(laz_files %in% str_replace(las_files, 'las$', 'laz'))] %>%
#   sample()



for (i in prog) {
  
  output_file <- str_replace(i, 'laz$', 'las')
  
  readLAS(i) %>%
    writeLAS(output_file)
  
}



ctg <- list.files(
    path = las_folder,
    pattern = '.las$',
    full.names = T
  ) %>%
  readLAScatalog()

boundary <- st_read(shp_file) %>%
  st_transform(st_crs(ctg))

opt_output_files(ctg) <- 'data/las/als/bldgt_als18_{id}'

x <- clip_roi(ctg, boundary)


#

als <- list.files(path = 'data/las/als',
                  pattern = '.las$',
                  full.names = T)

i = 1

x = readLAS(als[i]) %>%
  filter_duplicates() %>%
  filter_poi(Classification != 7, Classification != 18) %>%
  normalize_height(tin()) 

ttop <- locate_trees(x, lmf(ws = 5))

p <- plot(x, bg = "white", size = 4)
add_treetops3d(p, ttop)

# Look at variable window size by height



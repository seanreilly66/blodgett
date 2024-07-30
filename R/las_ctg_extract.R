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


ctg <- list.files(
  path = las_folder,
  pattern = file_ext,
  full.names = T
) %>%
  readLAScatalog()

boundary <- st_read(shp_file) %>%
  st_transform(st_crs(ctg))

x <- clip_roi(ctg, boundary)

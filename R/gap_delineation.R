# ==============================================================================
#
# Forest gap delineation
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 12 Dec 2024
# Last commit: 20 Jan 2025
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

library(tidyverse)
library(terra)
library(ForestGapR)
# library(viridis)
library(spatstat)
library(sf)
library(glue)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

chm_folder <- 'data/chm'
chm_id <- 'bldgt_als22_.+_chm_p2r_025\\.tif$'

study_shp <- 'data/gis/ffs_10m_buffer.shp'

gap_output <- 'bldgt_als22_c{comp}_t{treatment}_{type}_thresh{threshold}.{ext}'
gap_folder <- 'data/canopy_gaps'

threshold_heights <- c(1, 3, 5, 10)
size <- c(1, 10^10)


# ==============================================================================
# ============================= Gap identification =============================
# ==============================================================================

chm_files <- list.files(chm_folder, pattern = chm_id, full.names = T)

i = 1

for (chm_i in chm_files) {
  
  message('Processing: ', i, ' of ', length(chm_files))
  
  chm <- rast(chm_i)
  
  comp <- chm_i %>%
    str_extract('(?<=_c)[:digit:]+')
  
  treatment <- str_extract(chm_i, '(?<=_t)[:alpha:]+')
  
  analysis_window <- study_shp %>%
    st_read(quiet = T) %>%
    st_transform(st_crs(chm)) %>%
    filter(COMP == comp)
  
  chm <- mask(chm, analysis_window)
  
  # Generate gaps
  
  for (threshold in threshold_heights) {
    
    message('Threshold: ', threshold)
    
    gap_rast <- ForestGapR::getForestGaps(chm_layer = chm,
                              threshold = threshold,
                              size = size)
    
    gap_stats <- ForestGapR::GapStats(gap_layer = gap_rast, chm_layer = chm)
    
    gap_shp <- ForestGapR::GapSPDF(gap_rast) %>%
      merge(gap_stats, by = 'gap_id') %>%
      st_as_sf()
    
    shp_output <- glue(gap_folder, '/',
                       glue(gap_output, type = 'gapstatshp', ext = 'shp'))
    rast_output <- glue(gap_folder, '/',
                        glue(gap_output, type = 'gaprast', ext = 'tif'))
    
    writeRaster(gap_rast, 
                filename = rast_output,
                overwrite = T)
    
    st_write(gap_shp,
             dsn = shp_output,
             quiet = T)

  }
  
  i <- i + 1
  
}

# ==============================================================================
# ============================ Merged shape dataset ============================
# ==============================================================================

shp_files <- list.files(gap_folder, pattern = '.shp$', full.names = T)

# ==============================================================================
# ==============================================================================
# ==============================================================================
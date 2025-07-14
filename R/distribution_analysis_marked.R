# ==============================================================================
#
# Forest gap analysis
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 12 Oct 2024
# Last commit: 12 Dec 2024
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
# library(terra)
library(spatstat)
library(sf)
library(glue)
library(doParallel)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

canopy_shp_file <- 'data/poly_metrics/bldgt_ffs10m_als22_crown_metrics.shp'
window_shp_file <- 'data/gis/ffs_10m_buffer.shp'

t_i = 'f'

# ==============================================================================
# ============================= Centroid clustering ============================
# ==============================================================================

canopy_shp <- st_read(canopy_shp_file, quiet = T) %>%
  st_centroid() %>%
  group_by(comp) %>%
  mutate(quart = ntile(z_p95, 4) |>
           as.factor()) %>%
  filter(quart != 1) %>%
  mutate(quart = droplevels(quart))

window_shp <- st_read(window_shp_file, quiet = T) %>%
  st_transform(st_crs(canopy_shp))

canopy_shp <- canopy_shp %>%
  filter(lengths(st_within(geometry, window_shp)) > 0)

t_canopy <- canopy_shp %>%
  filter(trtmnt == t_i) 

t_window <- window_shp %>%
  filter(COMP %in% unique(t_canopy$comp)) 


t_df_q4 <- t_canopy %>%
  filter(quart == 4) %>%
  st_geometry %>%
  as.ppp(W = as.owin(t_window))

k_est_q4 <- envelope(t_df_q4, nsim = 25) 
allstat_q4 <- allstats(t_df_q4)




t_df <- t_canopy %>%
  st_geometry %>%
  as.ppp(W = as.owin(t_window))



t_quart <- t_df %mark% t_canopy$quart
k_quart <- alltypes(t_quart, Kcross, nsim = 25, envelope = TRUE)
ncor_quart <- nncorr(t_quart)
mcorr_quart <- markcorr(t_quart)
mconnect_quart <- alltypes(t_quart, markconnect)
mtable_quart <- marktable(t_quart, N = 1, collapse = T)



t_height <- t_df %mark% t_canopy$z_p95
ncor_height <- nncorr(t_height)
mcorr_height <- markcorr(t_height)
k_mark <- envelope(t_height, Kmark, nsim = 25) 

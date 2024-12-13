# ==============================================================================
#
# Polygon metric extraction
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

library(sf)
library(tidyverse)
library(glue)
library(doParallel)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

# File parameters

als_yr <- 22

tree_poly_folder <- 'data/tree_seg/lmf'
tree_poly_id = 'als{als_yr}.+treeseg_lmf\\.shp$'

chm_folder <- 'data/chm'
chm_id <- 'als{als_yr}_c{compartment}_t{treatment}_chm_p2r_025\\.tif$'

site_boundary <- 'data/gis/ffs_10m_buffer.shp'

# Polygon filter thresholds

h_thresh = 3
chm_cell_thresh = 2
n_point_thresh = 10

# Output

comp_output <- 'data/tree_seg/merged_crowns/bldgt_als{als_yr}_c{compartment}_t{treatment}_merged_crowns.shp'
combined_output <- 'data/tree_seg/bldgt_als{als_yr}_merged_crowns.shp'




las_folder <- 'data/las/als_2022'
las_id <- 'bldgt_als{als_yr}_c{compartment}_t{treatment}_hnorm_tin\\.las'




# ==============================================================================
# ====================== Extract metrics per merged group ======================
# ==============================================================================

source('R/las_metric_function.R')

las <- list.files(las_folder,
                  pattern = glue(las_id),
                  full.names = T) %>%
  readLAS(select = '')

poly_metrics <- polygon_metrics(las,
                                func = ~ las_cld_metrics(z = Z),
                                geometry = merged_poly)

poly_metrics <- poly_metrics %>%
  rowid_to_column(var = 'crown_id') %>%
  mutate(c = compartment, t = treatment, .before = 1) %>%
  mutate(area = st_area(geometry))




# merged_poly <- terra::extract(chm, merged_poly, method = 'simple') %>%
#   rename(merge_group = ID) %>%
#   nest(.by = merge_group, .key = 'chm') %>%
#   left_join(merged_poly, .)
# 
# theme_set(
#   theme(
#     text = element_text(family = 'serif', face = 'plain'),
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     line = element_line(linewidth = 1),
#     axis.line = element_line(),
#     panel.background = element_rect(color = 'white'),
#     legend.title = element_text(size = 16),
#     legend.text = element_text(size = 14),
#     legend.key = element_blank(),
#     legend.spacing = unit(0, "cm"),
#     legend.margin = margin(0, 5, 0, 5),
#     title = element_text(size = 12.8),
#     legend.position = 'bottom',
#     strip.background = element_blank(),
#     strip.text = element_text(size = 16,
#                               vjust = 1)
#   ))
# 
# ggplot(
#   data = poly_metrics,
#   mapping = aes(
#     x = area,
#     y = z_p95
#   )
# ) + 
#   geom_point() +
#   facet_wrap(~t, ncol = 2,  
#              strip.position = 'top', 
#              axes = 'all') +
#   labs(
#     x = 'Area (m2)',
#     y = 'p95 Height'
#   )
# 
# 
# 
# 
# col_pal = c(
#   'mf' = '#BB8B46',
#   'm' = '#356D34',
#   'f' = '#882255',
#   'c' = 'grey40'
# )
# 
# ggplot(
#   data = poly_metrics,
#   mapping = aes(
#     x = area,
#     y = z_p95,
#     color = t
#   )
# ) + 
#   geom_point(alpha = 0.75) +
#   facet_wrap(~t, ncol = 2,  
#              strip.position = 'top', 
#              axes = 'all') +
#   labs(
#     x = 'Area (m2)',
#     y = 'p95 Height'
#   )
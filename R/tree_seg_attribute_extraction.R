# ==============================================================================
#
# Segmented tree polygon and stem map overlap dataset generation
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
# From a given tree polygon shpfile, extracts chm and stem map attributes.
# Records all dbh and height values for stem map trees contained within each
# tree polygon. From CHM, records Z statistics per polygn.
#
# If supplied, extracts statistics from chm change raster
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

library(terra)
library(sf)
library(tidyverse)
library(lidR)
library(glue)
library(doParallel)
library(future)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

tree_poly_folder <- 'data/tree_seg/testing/lmf'
tree_poly_id = 'lmf\\.shp'

tree_z_thresh = 3

chm_folder <- 'data/chm'
chm_id <- 'als22_c{comp}.+p2r_025\\.tif$'

chm_diff_folder <- 'data/chm/diff'
chm_diff_id <- 'diff_c{comp}.+\\.tif$'

stem_map_folder <- 'data/field_stem_maps'
stem_map_id <- 'pts.shp$'

# las_folder <- 'data/las/als_2022'
# las_id <- 'bldgt_als22_c{comp}.+hnorm_tin\\.las'

# ==============================================================================
# =============================== Data wrangling ===============================
# ==============================================================================

tree_poly <- list.files(tree_poly_folder, pattern = tree_poly_id, full.names = T)

#testing
# poly_file <- tree_poly[1]

cl <- makeCluster(4)
registerDoParallel(cl)

df <- foreach(
  poly_file = tree_poly,
  .combine = 'rbind',
  .packages = c('lidR', 'tidyverse', 'sf', 'glue', 'terra')
) %dopar% {
  
  comp <- poly_file %>%
    str_extract('(?<=_c)[:digit:]+')
  
  sf <- st_read(poly_file, quiet = T) %>%
    add_column(comp = comp, .before = 1) %>%
    filter(Z > tree_z_thresh) %>%
    rowid_to_column(var = 'ID')
  
  # ----------------------------------------------------------------------------
  # ---------------------------------- CHM Z -----------------------------------
  # ----------------------------------------------------------------------------
  
  # Read in values
  
  chm <- list.files(chm_folder,
                    pattern = glue(chm_id),
                    full.names = T) %>%
    rast()
  
  names(chm) <- 'chm_z'
  
  chm_diff <- list.files(chm_diff_folder,
                         pattern = glue(chm_diff_id),
                         full.names = T) %>%
    rast()
  
  names(chm_diff) <- 'chm_diff'
  
  chm_stack = c(chm, chm_diff)
  
  # Extract values
  
  sf <- terra::extract(chm_stack, sf, method = 'simple') %>%
    nest(.by = ID, .key = 'chm') %>%
    left_join(sf, .)
  
  # Remove rasters
  rm(chm, chm_diff, chm_stack)
  gc()
  
  # ----------------------------------------------------------------------------
  # -------------------------------- Stem map ----------------------------------
  # ----------------------------------------------------------------------------
  
  # Read in values
  
  stem_map <- list.files(stem_map_folder, 
                         pattern = stem_map_id, 
                         full.names = T) %>%
    st_read(quiet = TRUE) %>%
    filter(comp == comp) %>%
    st_transform(st_crs(sf)) %>%
    rename(stem_map_treeID = treeID, stem_map_comp = comp)
  
  sf <- st_join(sf, stem_map) %>%
    as_tibble() %>%
    nest(stem_map = names(stem_map) |>
           str_subset('geometry', negate = T))
  
  
  
  # las <- list.files(
  #   las_folder,
  #   pattern = glue(las_id),
  #   full.names = T) %>%
  #   readLAS()
  
}


# ==============================================================================
# ================================== Analysis ==================================
# ==============================================================================


sm_func <- function(.df) {
  
  .df %>%
    summarize(
      sm_n = sum(!(is.na(stem_map_treeID))),
      sm_z_max = max(h_m, na.rm = F),
      sm_z_dif = diff(sort(h_m, decreasing = T))[1]
    )
  
}

df <- df %>%
  mutate(sm_stats = map(
    .x = stem_map,
    .f = ~ sm_func(.df = .x)
  )) %>%
  unnest(sm_stats)

# group_stats <- df %>%
#   group_by(
#     n_tree = sm_n > 0
#   ) %>%
#   summarize(
#     n_trees = n(),
#     mean_pts = mean(npoints),
#     min_pts = min(npoints),
#     max_pts = max(npoints),
#     mean_sm_trees = mean(sm_n),
#     min_sm_trees = min(sm_n),
#     max_sm_trees = max(sm_n),
#     mean_sm_z = mean(sm_z_dif, na.rm = T),
#     min_sm_z = min(sm_z_dif, na.rm = T),
#     max_sm_z = max(sm_z_dif, na.rm = T)
#   )
# 
# true_seg <- df %>%
#   filter(
#     sm_n > 0,
#     sm_z_dif < -0.1 * Z | is.na(sm_z_dif)
#   ) %>%
#   select(-chm, -stem_map)
# 
# st_write(true_seg, 'data/tree_seg/testing/true_seg.shp', append = FALSE)
# 
# ggplot(data = true_seg, 
#        mapping = aes(x = Z,
#                      y = npoints)) +
#   geom_point() +
#   geom_hline(yintercept = 100) +
#   theme_classic()


# ==============================================================================
# ================================== Analysis ==================================
# ==============================================================================

df = df %>%
  st_as_sf() %>%
  mutate(intersect = st_intersects(.)) %>%
  select(-ID) %>%
  rowid_to_column()

chm_sim_func <- function(i, x) {
  
  x = x %>%
    sf::st_drop_geometry() %>%
    as_tibble()
  
  intersect_id = x %>%
    filter(rowid == i) %>%
    pull(intersect) %>%
    simplify()
  
  if (length(intersect_id) < 2) {
    
    res <- tibble(
      rowid = i,
      chm_z_stat = NA,
      chm_diff_stat = NA,
      shared_stat = NA
    )
    
    return(res)
    
  }
  
  overlap <- x %>%
    filter(rowid %in% intersect_id) %>%
    unnest(chm)
  
  chm_z <- DescTools::DunnettTest(chm_z ~ rowid, data = overlap, control = i)
  
  chm_z = chm_z[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(pval < 0.05) %>%
    mutate(poly_id = str_extract(rowname, '^[:digit:]+(?=-)')) %>%
    pull(poly_id)
  
  chm_diff <- DescTools::DunnettTest(chm_diff ~ rowid, data = overlap, control = i)
  
  chm_diff = chm_diff[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(pval < 0.05) %>%
    mutate(poly_id = str_extract(rowname, '^[:digit:]+(?=-)')) %>%
    pull(poly_id)
  
  shared_stat <- chm_z[chm_z %in% chm_diff]
  
  res <- tibble(
    rowid = i,
    chm_z_stat = list(chm_z),
    chm_diff_stat = list(chm_diff),
    shared_stat = list(shared_stat)
  )
  
  return(res)
  
}



cl <- makeCluster(15)
registerDoParallel(cl)

sim_stat <- foreach(
  i = df$rowid,
  .combine = 'rbind',
  .packages = c('tidyverse')
) %dopar% {
  
  chm_sim_func(i, df)
  
}
  
stopCluster(cl)


# ==============================================================================
#
# Merged tree polygons
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

library(terra)
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

# ==============================================================================
# ======================== Similarity testing function =========================
# ==============================================================================

similarity_test <- function(row_id, data) {
  #testing
  # data = crown_poly
  # row_id = 27
  #
  
  intersect_id = data %>%
    filter(crown_id == row_id) %>%
    pull(intersect) %>%
    list_c()
  
  df = data %>%
    sf::st_drop_geometry() %>%
    as_tibble() %>%
    filter(crown_id %in% intersect_id) %>%
    mutate(across(crown_id, as_factor)) %>%
    unnest(chm) %>%
    mutate(crown_id = relevel(crown_id, ref = as.character(row_id)))
  
  if (length(intersect_id) == 1) {
    similarity_stats <- tibble(crown_id = row_id,
                               p_val = NA,
                               merge_crowns = NA)
    
  } else if (length(intersect_id) == 2) {
    t_test <- t.test(chm_z ~ crown_id, data = df)
    
    t_p = tibble(comparison_id = setdiff(intersect_id, row_id),
                 p = t_test$p.value)
    
    merge_crowns = t_p %>%
      filter(p > 0.05) %>%
      pull(comparison_id)
    
    if (length(merge_crowns) == 0) {
      merge_crowns = NA
      
    }
    
    similarity_stats <- tibble(
      crown_id = row_id,
      p_val = list(t_p),
      merge_crowns = list(merge_crowns)
    )
    
  } else {
    dunnett <- PMCMRplus::dunnettTest(chm_z ~ crown_id, data = df)
    
    dunnett_p <- dunnett$p.value %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      rename(comparison_id = 1, p = 2)
    
    merge_crowns <- dunnett_p %>%
      filter(p > 0.05) %>%
      pull(comparison_id)
    
    if (length(merge_crowns) == 0) {
      merge_crowns = NA
      
    }
    
    similarity_stats <- tibble(
      crown_id = row_id,
      p_val = list(dunnett_p),
      merge_crowns = list(merge_crowns)
    )
    
  }
  
  return(similarity_stats)
  
}


# ==============================================================================
# ============================ Merge crown polygons ============================
# ==============================================================================


tree_poly <- list.files(tree_poly_folder,
                        pattern = glue(tree_poly_id),
                        full.names = T)

i = 0

#testing
# poly_file <- tree_poly[7]

merged_poly <- foreach(poly_file = tree_poly, .combine = 'rbind') %do% {

  i = i + 1
  message(glue('Processing {i} of {length(tree_poly)}'))
  tictoc::tic()
  
  # ============================================================================
  # ============================== Data wrangling ==============================
  # ============================================================================

  # Extract attributes from file name
  
  compartment <- str_extract(poly_file, '(?<=_c)[:digit:]+')
  treatment <- str_extract(poly_file, '(?<=_t)[:alpha:]+')
  
  # Read in crown polygon data
  
  crown_poly <- st_read(poly_file, quiet = T) %>%
    add_column(c = compartment, .before = 1, t = treatment) %>%
    filter(Z >= h_thresh) %>%
    filter(npoints >= n_point_thresh) %>%
    select(-treeID)
  
  # Read in study site boundary and filter to only crowns contained within it
  
  site <- st_read(site_boundary, quiet = T) %>%
    st_transform(st_crs(crown_poly)) %>%
    filter(COMP == compartment)

  crown_poly <- crown_poly[st_within(crown_poly$geometry, site$geometry, sparse = F), ]
  
  crown_poly <- crown_poly %>%
    filter(!is.na(Z)) %>% # Compartment 40 creating empty geometries
    rowid_to_column(var = 'crown_id')
  
  # Read and extract CHM data
  
  chm <- list.files(chm_folder,
                    pattern = glue(chm_id),
                    full.names = T) %>%
    rast()
  
  names(chm) <- 'chm_z'
  
  crown_poly <- terra::extract(chm, crown_poly, method = 'simple') %>%
    rename(crown_id = ID) %>%
    nest(.by = crown_id, .key = 'chm') %>%
    left_join(crown_poly, ., by = join_by(crown_id))
  
  crown_poly <- crown_poly %>%
    mutate(n_chm = sapply(chm, nrow)) %>%
    filter(n_chm >= chm_cell_thresh) %>%
    select(-crown_id) %>%
    rowid_to_column(var = 'crown_id')
  
  
  # ============================================================================
  # ======================== ID similar crown polygons =========================
  # ============================================================================
  
  crown_poly <- crown_poly %>%
    st_as_sf() %>%
    mutate(intersect = st_intersects(.))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  similarity_stats <- foreach(
    row_id = crown_poly$crown_id,
    .combine = 'rbind',
    .packages = c('tidyverse', 'sf')
  ) %dopar% {
    
    similarity_test(row_id, data = crown_poly)
    
  }
  
  stopCluster(cl)
  
  crown_poly <- crown_poly %>%
    left_join(similarity_stats, by = join_by(crown_id))
  
  
  # ============================================================================
  # ============== Generate merge groups and combine geometries ================
  # ============================================================================
  
  crown_poly <- crown_poly %>%
    mutate(merge_group = crown_id)
  
  for (id in crown_poly$crown_id) {
    
    merge_group = crown_poly$merge_group[crown_poly$crown_id == id]
    
    merge_ids = crown_poly %>%
      filter(crown_id == id) %>%
      pull(merge_crowns) %>%
      list_c() %>%
      as.numeric()
    
    crown_poly$merge_group[crown_poly$crown_id %in% merge_ids] = merge_group
    
  }
  
  merged_poly <- crown_poly %>%
    group_by(merge_group) %>%
    summarize(geometry = st_union(geometry), .groups = "drop") %>%
    st_cast('MULTIPOLYGON') %>%
    add_column(
      comp = compartment,
      trtmnt = treatment,
      .before = 1
    )
  
  st_write(merged_poly, glue(comp_output), append = F)
  
  return(merged_poly)
  
  message(tictoc::toc())
  
}

# ==============================================================================

st_write(merged_poly, glue(combined_output))

# ==============================================================================
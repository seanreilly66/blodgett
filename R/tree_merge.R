# ==============================================================================
#
# Identify tree polygon regions
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

library(terra)
library(sf)
library(tidyverse)
library(glue)
library(doParallel)

# library(future)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

als_yr <- 22

tree_poly_folder <- 'data/tree_seg/lmf'
tree_poly_id = 'als{als_yr}.+treeseg_lmf\\.shp$'

chm_folder <- 'data/chm'
chm_id <- 'als{als_yr}_c{compartment}.+p2r_025\\.tif$'

# tree_z_thresh = 5
# las_folder <- 'data/las/als_2022'
# las_id <- 'bldgt_als22_c{comp}.+hnorm_tin\\.las'

# ==============================================================================
# =============================== Data wrangling ===============================
# ==============================================================================

tree_poly <- list.files(tree_poly_folder,
                        pattern = glue(tree_poly_id),
                        full.names = T)

#testing
poly_file <- tree_poly[1]

# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# df <- foreach(
#   poly_file = tree_poly,
#   .combine = 'rbind',
#   .packages = c('lidR', 'tidyverse', 'sf', 'glue', 'terra')
# ) %dopar% {


# Extract attributes from file name

compartment <- str_extract(poly_file, '(?<=_c)[:digit:]+')
treatment <- str_extract(poly_file, '(?<=_t)[:alpha:]+')

# Read in crown polygon data

crown_poly <- st_read(poly_file, quiet = T) %>%
  add_column(c = compartment, .before = 1,
             t = treatment) %>%
  # filter(Z > tree_z_thresh) %>%
  rowid_to_column(var = 'crown_id') %>%
  select(-treeID)

# Read in CHM data

chm <- list.files(chm_folder,
                  pattern = glue(chm_id),
                  full.names = T) %>%
  rast()

names(chm) <- 'chm_z'

# Extract chm values

crown_poly <- terra::extract(chm, crown_poly, method = 'simple') %>%
  rename(crown_id = ID) %>%
  nest(.by = crown_id, .key = 'chm') %>%
  left_join(crown_poly, .) %>%
  mutate(across(crown_id, as_factor))

# Remove chm and clean memory
rm(chm)
gc()


# ==============================================================================
# ========================= ID similar crown polygons ========================== 
# ==============================================================================

crown_poly <- crown_poly %>%
  st_as_sf() %>%
  mutate(intersect = st_intersects(.))

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
    unnest(chm) %>%
    mutate(crown_id = relevel(crown_id, ref = row_id))
  
  if (length(intersect_id) == 1) {
    
    similarity_stats <- tibble(
      crown_id = row_id,
      p_val = NA,
      merge_crowns = NA
    )
    
  } else if (length(intersect_id) == 2) {
    
    t_test <- t.test(chm_z ~ crown_id, data = df)
    
    t_p = tibble(
      comparison_id = setdiff(intersect_id, row_id),
      p = t_test$p.value
    )
    
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



cl <- makeCluster(detectCores())
registerDoParallel(cl)

similarity_stats <- foreach(
  row_id = as.numeric(crown_poly$crown_id),
  .combine = 'rbind',
  .packages = c('tidyverse', 'sf')
) %dopar% {
  
  similarity_test(row_id, data = crown_poly)
  
}

stopCluster(cl)

crown_poly <- crown_poly %>%
  left_join(similarity_stats |>
              mutate(across(crown_id, as.factor)))








x = slice(crown_poly, 1:100)

# Create a new column for the merge group
polygons <- x %>%
  mutate(
    merge_group = sapply(merge_crowns, function(ids) {
      if (is.na(ids)) {
        return(crown_id)  # Single id if no merge
      } else {
        return(unlist(ids))  # Merge with multiple ids
      }
    })
  )

# Flatten the merge_group column to a single ID for each group (just to make grouping easier)
polygons$merge_group <- sapply(polygons$merge_group, function(ids) {
  if (length(ids) > 1) {
    return(min(ids))  # Use the minimum id as the group identifier
  } else {
    return(ids)
  }
})

print(polygons)

# ==============================================================================
#
# Allometric equations for tree crown radii and height
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 19 Sept 2024
# Last commit: 24 Sept 2024
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# Computes tree height and crown radii using allometric equations from DBH
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

# ================================= User inputs ================================

stem_map_folder <- 'data/field_stem_maps/connie_stem_maps'

# ================================ Tree metrics ================================ 

h_eq <- function(eq_type, a, b, dbh) {
  
  # Equations from John Battles, personal communication
  # DBH in cm
  
  if (eq_type == 'linear') {
    
    h_m <- a + b * dbh
    
  } else if (eq_type == 'michment') {
    
    h_m <- (a * dbh)/(b * dbh)
    
  } else if (eq_type == 'power') {
    
    h_m <- a * dbh^b
    
  } else if (eq_type == 'siccama') {
    
    h_m <- 1.37 + a * (1 - exp(-b * dbh))
    
  } else {
    
    h_m = NA
    
  }
  
}

tree_calc <- function(species, dbh_cm, dbh_in) {
  
  # Gill et al 2000 crown radii
  
  if (species == 'DF') {
    
    crrad_m <- 1.6654 + 0.0355*dbh_cm
    h_m <- h_eq('siccama', a = 77.7282, b = 0.0088, dbh = dbh_cm)
    
    
  } else if (species == 'IC') {
    
    crrad_m <- 1.2960 + 0.0256*dbh_cm
    h_m <- h_eq('power', a = 0.7171, b = 0.8832, dbh = dbh_cm)
    

  # } else if (species == 'LP') {
  #   
  #   crown_r <- 0.5230 + 0.0440*dbh_cm
    
  } else if (species == 'F' | species == 'OC' | species == 'UKN') {
    
    crrad_m <- 1.1817 + 0.0265*dbh_cm
    h_m <- NA
    
  # } else if (species == 'OG') {
  #   
  #   crown_r <- 2.5067 + 0.0244*dbh_cm
    
  } else if (species == 'PP') {
    
    crrad_m <- 0.9488 + 0.0356*dbh_cm
    h_m <- h_eq('power', a = 0.5956, b = 0.9758, dbh = dbh_cm)
    
  # } else if (species == 'RF') {
  #   
  #   crown_r <- 1.0171 + 0.0259*dbh_cm
  #   
  # } else if (species == 'RW') {
  #   
  #   crown_r <- 1.7371 + 0.0311*dbh_cm
    
  } else if (species == 'SP') {
    
    crrad_m <- 0.9906 + 0.0398*dbh_cm
    h_m <- h_eq('siccama', a = 161.0388, b = 0.0033, dbh = dbh_cm)
    
    
  } else if (species == 'WF') {
    
    crrad_m <- 1.2256 + 0.0299*dbh_cm
    h_m <- h_eq('siccama', a = 109.7001, b = 0.0058, dbh = dbh_cm)
    
  # Bechtold 2004
    
  } else if (species == 'BO') {
    
    crdiam_ft <- 7.0284 + 1.0470*dbh_in 
    crrad_m <- crdiam_ft * 0.3048 / 2
    h_m <- h_eq('siccama', a = 35.2608, b = 0.0206, dbh = dbh_cm)
    
  } else if (species == 'TO') {
    
    crdiam_ft <- 6.7864 + 0.8443*dbh_in
    crrad_m <- crdiam_ft * 0.3048 / 2
    h_m <- h_eq('siccama', a = 51.3747, b = 0.0111, dbh = dbh_cm)
    
    
  } else {
    
    crrad_m <- NA
    h_m <- NA
    
  }
  
  return(tibble(crrad_m, h_m))
  
}



# ==========================-=== Data manipulation ============================= 

read_shp <- function(file) {
  
  shp <- st_read(file, quiet = TRUE)
  
  c <- str_extract(file, '(?<=_c)[:digit:]+')
  
  shp <- shp %>%
    add_column(
      comp = c
    )
  
}

stem_map <- list.files(path = stem_map_folder,
                       pattern = 'stem_pts.shp$',
                       full.names = TRUE) %>%
  lapply(read_shp) %>%
  bind_rows() %>%
  select(
    comp,
    treeID = TO,
    species = SPP,
    status = STATUS,
    dbh_cm = DBH_cm,
    dbh_in = DBH,
    x_m = X_m,
    y_m = Y_m,
    crrad_m,
    geometry
  ) 


stem_map <- stem_map %>%
  rowwise() %>%
  mutate(
    tree_calc(species = species, dbh_cm = dbh_cm, dbh_in = dbh_in) |>
      rename(crad_m_rep = crrad_m),
    .before = 'geometry'
  ) 

st_write(stem_map, 
         'data/field_stem_maps/bldgt_field_stem_pts.shp',
         append = FALSE) 

# ==============================================================================
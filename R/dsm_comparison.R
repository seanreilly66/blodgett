# ==============================================================================
#
# DSM time comparison
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 2 Oct 2024
# Last commit: 2 Oct 2024
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
library(glue)
library(doParallel)
library(terra)

# ================================= User inputs ================================

chm_folder <- 'data/chm'
dtm_folder <- 'data/dtm'
p2r_res <- 0.25

chm_output <- 'data/chm/diff'
dtm_output <- 'data/dtm/diff'

# ==============================================================================
# ================================= Data prep ================================== 
# ==============================================================================


chm_files <- list.files(
  chm_folder,
  pattern = glue('{str_remove(p2r_res, "[:punct:]")}.tif$'),
  full.names = T
)

dtm_files <- list.files(dtm_folder, pattern = '.tif$', full.names = T)

comps <- str_extract(chm_files, pattern = '(?<=_c)[:digit:]+') %>% unique

# ==============================================================================
# ========================== Difference calculation ============================ 
# ==============================================================================


cl <- makeCluster(12)
registerDoParallel(cl)

foreach(
  comp = comps,
  .packages = c('terra', 'tidyverse', 'glue')
) %do% {

  # comp = comps[2]
  
  chm22 <- chm_files %>%
    str_subset('als22') %>%
    str_subset(glue('_c{comp}_')) %>%
    rast()
  
  chm18 <- chm_files %>%
    str_subset('als18') %>%
    str_subset(glue('_c{comp}_')) %>%
    rast()
  
  chm18 <- extend(chm18, chm22)
  chm22 <- extend(chm22, chm18)
  chm_dif <- chm22 - chm18
  
  writeRaster(
    chm_dif,
    filename = chm_files |>
      str_subset('als22') |>
      str_subset(glue('_c{comp}_')) |>
      str_replace('als22', 'diff') |>
      str_replace(chm_folder, chm_output),
    overwrite = T
  )
  
  dtm22 <- dtm_files %>%
    str_subset('als22') %>%
    str_subset(glue('_c{comp}_')) %>%
    rast()
  
  dtm18 <- dtm_files %>%
    str_subset('als18') %>%
    str_subset(glue('_c{comp}_')) %>%
    rast()
  
  dtm18 <- extend(dtm18, dtm22)
  dtm22 <- extend(dtm22, dtm18)
  dtm_dif <- dtm22 - dtm18
  
  writeRaster(
    dtm_dif,
    filename = dtm_files |>
      str_subset('als22') |>
      str_subset(glue('_c{comp}_')) |>
      str_replace('als22', 'diff') |>
      str_replace(dtm_folder, dtm_output),
    overwrite = T
  )

}

stopCluster(cl)

# ==============================================================================
# ==============================================================================
#
# Lidar tree segmentation
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 17 Sept 2024
# Last commit: 17 Sept 2024
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

library(lidR)
library(tidyverse)
library(sf)
library(glue)
library(doParallel)

# ================================= User inputs ================================

als_folder <- 'data/las/als_2022'
id <- 'bldgt_als.+hnorm_tin\\.las'

stem_map_folder <- 'data/field_stem_maps'

test_compartment_list <- c(240, 350, 380, 400)

# ==============================================================================
# ================================= Data prep ==================================
# ==============================================================================

# Lazy code for looping through compartments
# if (exists('i')) {
# 
#   i = i + 1
# 
#   if (i > 4) {
# 
#     rm(i)
#     stop('Done')
# 
#   }
# 
# } else {
# 
#   i = 1
# 
# }
# 
# test_compartment <- test_compartment_list[i]

cl <- makeCluster(4)
registerDoParallel(cl)

foreach(
  test_compartment = test_compartment_list,
  .packages = c('lidR', 'tidyverse', 'sf', 'glue')
) %dopar% {
  
  # ===================================== ALS ====================================
  
  als_files <- list.files(path = als_folder,
                          pattern = id,
                          full.names = T)
  
  als_yr <- str_extract(als_folder, '[:digit:]{2}$')
  
  als_c <- als_files %>%
    str_subset(glue('c{test_compartment}')) %>%
    readLAS()
  
  crs <- st_crs(als_c)
  
  # ===================================== SHP ====================================
  
  stmap_shp_files <- list.files(
    path = stem_map_folder,
    pattern = '.shp$',
    full.names = T,
    recursive = T
  )
  
  stem_map <- stmap_shp_files %>%
    str_subset(glue('bldgt_field_stem_pts.shp')) %>%
    st_read() %>%
    st_transform(crs)
  
  plt_bndry <- stmap_shp_files %>%
    str_subset('bldgt_stem_map_plts.shp') %>%
    st_read() %>%
    filter(plot_name == str_subset(plot_name, glue('_{test_compartment}'))) %>%
    st_transform(crs) %>%
    st_zm() # Drop Z from plot boundary
  
  
  # ==============================================================================
  # ========================= Crown rad to height ratio ==========================
  # ==============================================================================
  
  # f_exp <- mosaic::fitModel(crrad_m ~ A + B * h_m ^ C,  data = stem_map)
  # summary(f_exp)

  f_lm <- lm(crrad_m ~ h_m, data = stem_map)
  summary(f_lm)
  
  # Plot of crown height vs radius relationship
  
  # theme_set(
  #   theme(
  #     text = element_text(family = 'serif', face = 'plain'),
  #     axis.title = element_text(size = 16),
  #     axis.text = element_text(size = 14),
  #     line = element_line(linewidth = 0.5),
  #     axis.line = element_line(),
  #     panel.background = element_rect(color = 'white'),
  #     legend.title = element_text(size = 16),
  #     legend.text = element_text(size = 14),
  #     legend.key = element_blank(),
  #     legend.spacing = unit(0, "cm"),
  #     legend.margin = margin(0, 5, 0, 5),
  #     title = element_text(size = 12.8)
  #   )
  # )
  #
  # ggplot(
  #   data = stem_map |>
  #     filter(status == 'L'),
  #   mapping = aes(
  #     x = h_m,
  #     y = crrad_m
  #   )
  # ) +
  #   geom_smooth(aes(color = species), se = F, linewidth = 1.5) +
  #   geom_abline(slope = f_lm$coefficients[2],
  #               intercept = f_lm$coefficients[1],
  #               linewidth = 1,
  #               linetype = 'dashed') +
  #   scale_color_brewer(palette = 'Pastel2') +
  #   labs(x = 'Tree height (m)',
  #     y = 'Crown radius (m)',
  #     color = 'Species'
  #   ) +
  #   lims(
  #     x = c(0,60),
  #     y = c(0,8)
  #   ) +
  #   guides(color = guide_legend(override.aes = list(linewidth = 6)))
  #
  # ggsave('figures/temp/tree_height_v_crown.png',
  #        width = 6.5,
  #        height = 6,
  #        units = 'in',
  #        dpi = 700,
  #        bg = 'white')
  
  # ==============================================================================
  # ========================= Testing tree segmentation ==========================
  # ==============================================================================
  
  als_p <- clip_roi(als_c, plt_bndry)
  
  kernel <- matrix(1, 3, 3)
  chm_p <- rasterize_canopy(als_p, 0.25, p2r(subcircle = 0.2), pkg = "terra") %>%
    terra::focal(w = kernel, fun = median, na.rm = TRUE)
  
  terra::writeRaster(
    chm_p, 
    glue('data/tree_seg/testing/chm/als{als_yr}_c{test_compartment}_chm_p.tif'),
    overwrite = T)
  
  
  # ------------------------------------------------------------------------------
  # -------------- LMF with lm function from height to radius ratio --------------
  # ------------------------------------------------------------------------------
  
  ws_f <- function(x) {x * f_lm$coefficients[2] + f_lm$coefficients[1]}
  ttop_p <- locate_trees(als_p, lmf(ws_f))

  als_seg <- segment_trees(als_p, dalponte2016(chm_p, ttop_p, max_cr = 40))

  crown_poly_lmf <- crown_metrics(als_seg, func = .stdtreemetrics, geom = "convex")

  st_write(crown_poly_lmf,
           glue('data/tree_seg/testing/lmf/als{als_yr}_c{test_compartment}_crown_poly_lmf.shp'),
           overwrite = T)
  
  # ------------------------------------------------------------------------------
  # --------------------------------- Watershed ----------------------------------
  # ------------------------------------------------------------------------------

  # als_seg <- segment_trees(als_p, watershed(chm_p, tol = 1))
  # 
  # crown_poly_wtrshd <- crown_metrics(als_seg, func = .stdtreemetrics, geom = "convex")
  # 
  # st_write(
  #   crown_poly_wtrshd,
  #   glue(
  #     'data/temp/tree_seg/c{test_compartment}_crown_poly_wtrshd_1tol.shp'
  #   )
  # )
  
  # ------------------------------------------------------------------------------
  # --------------------------------- Li ----------------------------------
  # ------------------------------------------------------------------------------

  # als_seg <- segment_trees(als_p, li2012(dt1 = 1, dt2 = 1.5))
  # 
  # crown_poly_wtrshd <- crown_metrics(als_seg, func = .stdtreemetrics, geom = "convex")
  # 
  # st_write(
  #   crown_poly_wtrshd,
  #   glue(
  #     'data/temp/tree_seg/als22_c{test_compartment}_crown_poly_li_halfparam.shp'
  #   )
  # )
  
}

stopCluster(cl)



# x = plot(crown_poly)
#
# plot(crwn_met)
#
# plot_dtm3d(chm_p)
#
# plot(crowns_p)
#
# plot(als_seg_p, color = 'treeID')
#
# p = plot(als_seg_p, color = 'treeID')
# add_treetops3d(p, ttop_p)
#
#
# terra::writeRaster(chm_p,
#                    glue('data/temp/tree_seg/c{test_compartment}_chm_p.tif'))
#
# # Segment with stem map locations
#
# stmap_ttop <- stmap_pt_c %>%
#   select(TO, geometry) %>%
#   rename(treeID = TO) %>%
#   st_zm()
#
# stmap_z = terra::extract(chm_p, stmap_ttop) %>%
#   rename(treeID = ID,
#          Z = focal_median)
#
# stmap_ttop <- stmap_ttop %>%
#   left_join(stmap_z)
#
# stmap_algo <- dalponte2016(chm, stmap_ttop, max_cr = 40)
# stmap_seg <- segment_trees(als_p, stmap_algo)
# stmap_crowns <- stmap_algo()
#
# plot(stmap_crowns)
#
# p = plot(stmap_seg, color = 'treeID')
# add_treetops3d(p, stmap_ttop)
#
#
#
#
# st_write(st_zm(ttop), glue('data/tree_seg/c{test_compartment}_ttop_test.shp'))
#
#
# p1 = plot(als_seg_p, bg = "white", size = 4)
# add_treetops3d(p1, ttop, col = 'black')
#
#
# x = plot(rbind(st_bbox(als_c), st_bbox(plt_bndry)))
# plot(st_bbox(plt_bndry), add = T)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# plot(filter_poi(mf180_seg, Z > 3), bg = "white", size = 4, color = "treeID")
#
#
# mf180_mtrc <- crown_metrics(mf180_seg, func = .stdtreemetrics, geom = "convex")
#
# plot(mf180_mtrc["Z"], main = "Height")
#
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
#     title = element_text(size = 12.8)
#   )
# )
#
# mf180_mtrc <- mf180_mtrc %>%
#   filter(npoints > 1000)
#
# ggplot(
#   data = as.data.frame(mf180_mtrc),
#   aes(x = Z)) +
#   geom_histogram(fill = '#77B68A', color = 'black', binwidth = 5, boundary = 0) +
#   labs(y = 'Count',
#        x = 'Height (m)') +
#   scale_x_continuous(breaks = seq(0,50,5))

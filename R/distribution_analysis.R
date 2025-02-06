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

# ==============================================================================
# ============================= Centroid clustering ============================
# ==============================================================================

canopy_shp <- st_read(canopy_shp_file, quiet = T) %>%
  st_centroid()

window_shp <- st_read(window_shp_file, quiet = T) %>%
  st_transform(st_crs(canopy_shp))

# k_est = list()

# c = c(240, 340, 190, 180)

cl <- makeCluster(2)
registerDoParallel(cl)

k_est <- foreach(
  c_i = unique(canopy_shp$comp),
  .combine = 'rbind',
  .packages = c('tidyverse', 'sf', 'spatstat')
) %dopar% {

  c_window <- window_shp %>%
    filter(COMP == c_i)
  
  c_canopy <- canopy_shp %>%
    filter(comp == c_i) 
  
  c_df <- c_canopy %>%
    st_geometry %>%
    as.ppp(W = as.owin(c_window))
  
  base_k <- envelope(c_df, nsim = 25) %>%
    as.data.frame() %>%
    add_column(
      mark = 'none',
      .before = 1
    )

  # marks(c_df) <- c_canopy$z_p95
  # 
  # zp95_k = envelope(c_df, nsim = 25, marks = TRUE) %>%
  #   add_column(
  #     mark = 'z_p95',
  #     .before = 1)
  # 
  # marks(c_df) <- c_canopy$area
  # marks(c_df) <- c_canopy$lddr_vn
  # 
  # full_k = envelope(c_df, nsim = 25, marks = TRUE) %>%
  #   add_column(
  #     mark = 'z_p95 area lddr_vn',
  #     .before = 1)
  
  base_k <- base_k %>%
    add_column(
      comp = c_i,
      trtmnt = unique(c_canopy$trtmnt)
    )

}

write_rds(k_est, 'canopy_p90_k.RData')

k_est <- read_rds('canopy_p90_k.RData')

x = lapply(k_est, as.data.frame)

x = map2(x, names(x), ~mutate(.x, trtmnt = .y)) %>%
  bind_rows() %>%
  mutate(trtmnt = case_match(
    trtmnt,
    'c' ~ 'Control',
    'f' ~ 'Fire',
    'm' ~ 'Mech',
    'mf' ~ 'Mech + Fire'
  ))

theme_set(
  theme(
    text = element_text(family = 'serif', face = 'plain'),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    line = element_line(linewidth = 1),
    axis.line = element_line(),
    panel.background = element_rect(color = 'white'),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.spacing = unit(0, "cm"),
    legend.margin = margin(0, 5, 0, 5),
    title = element_text(size = 12.8),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, hjust = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
    
  )
)

ggplot(
  data = x, 
  aes(x = r)) +
  geom_ribbon(
    aes(ymin = lo,
        ymax = hi),
    fill = "grey80") +  # Create the shaded region
  geom_line(
    aes(y = theo), 
    color = "firebrick",
    linetype = 'dashed') +  
  geom_line(
    aes(y = obs), 
    color = "black") +  
  facet_wrap(~trtmnt) +
  labs(title = 'Distribution of top 10% TAO by height',
       y = 'K(r)',
       x = 'radius (m)')

ggsave(
  filename = 'figures/kest_p90h.png',
  width = 8,
  height = 7,
  units = 'in',
  dpi = 700
)





# plot(density(cent_df))
# 
# 
# 
# 
# 
# 
# 
# # Raster cell clustering
# 
# 
# prop_samp = 0.01
# 
# pt_df <- as.data.frame(gap_rast, xy = T)
# 
# set.seed(111)
# 
# samp_pt_df <- pt_df %>%
#   group_by(gaps) %>%
#   slice_sample(prop = prop_samp) %>%
#   ungroup()
# 
# samp_pt_df <- samp_pt_df %>%
#   bind_rows(
#     pt_df |>
#       group_by(gaps) |>
#       filter(!gaps %in% samp_pt_df$gaps) |>
#       slice_sample(n = 1)  # Retaining at least one per group) %>%
#     ungroup() %>%
#       select(-gaps)
#     
#     pts <- as.ppp(samp_pt_df, W = as.owin(study_region))
#     
#     pts <- subset(pts, inside = T)
#     
#     plot(pts)
#     
#     pts_k = envelope(pts)
#     
#     plot(pts_k)
#     
#     
#     
#     
#     
#     
#     plot(pt_df)
#     
#     
#     
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(oldpar))
#     P <- spatstat.geom::as.ppp(sp::coordinates(gap_SPDF_layer), terra::ext(chm_layer)[])
#     K <- spatstat.explore::envelope(
#       P,
#       spatstat.explore::Kest,
#       nsim = 99,
#       verbose = F
#     )
#     L <- spatstat.explore::envelope(
#       P,
#       spatstat.explore::Lest,
#       nsim = 99,
#       verbose = F
#     )
#     graphics::par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))
#     graphics::plot(K)
#     graphics::plot(L)
#     CE <- spatstat.explore::clarkevans.test(P)
#     return(CE)
#     
#     
#     plot(gap_sp, col = 'darkred')
#     plot(ppp_df, add = T)
#     plot(gap_sp, col = 'grey50', add = T)
#     
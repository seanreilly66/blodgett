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
library(sf)
library(glue)
library(lme4)
library(lmerTest)
library(emmeans)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

# gap_output <- 'bldgt_als22_c{comp}_t{treatment}_{type}_thresh{threshold}.{ext}'
study_shp <- 'data/gis/ffs_10m_buffer.shp'
gap_folder <- 'data/canopy_gaps/shp'
primary_gap_sizes <- c(1,3,5,10)

figure_output <- 'figures/gap_dist'

plot_helpers <- 'R/plot_helpers.R'

# ==============================================================================
# ============================ Merged shape dataset ============================
# ==============================================================================

source(plot_helpers)

shp_read <- function(file) {
  
  c = str_extract(file, '(?<=_c)[:digit:]+')
  t = str_extract(file, '(?<=[:digit:]_t)[:alpha:]+')
  thresh = str_extract(file, '(?<=_thresh)[:digit:]+')
  
  shp <- st_read(file, quiet = T) %>%
    add_column(
      comp = c,
      trtmnt = t,
      gap_thresh = thresh,
      .before = 1
    )
  
  
}

gap_shp = gap_folder %>%
  list.files(pattern = '.shp$', full.names = T) %>%
  lapply(shp_read) %>%
  bind_rows() %>%
  mutate(
    trtmnt = case_match(
      trtmnt,
      'c' ~ 'Control',
      'mf' ~ 'Mech + Fire',
      'f' ~ 'Fire',
      'm' ~ 'Mech'
    ) |>
      fct_relevel(
        'Control', 
        'Fire', 
        'Mech', 
        'Mech + Fire'
      ),
    comp = as.integer(comp)
  ) %>%
  st_drop_geometry()

comp <- gap_shp %>%
  st_drop_geometry() %>%
  group_by(comp, trtmnt) %>%
  summarize(n = n())

comp_area = st_read(study_shp, quiet = T) %>%
  mutate(comp_area = st_area(geometry)) %>%
  left_join(comp, by = join_by(COMP == comp)) %>%
  rename(comp = COMP) %>%
  select(trtmnt, comp, comp_area)

trtmnt_area <- comp_area %>%
  group_by(trtmnt) %>%
  st_drop_geometry() %>%
  summarize(
    trtmnt_area = as.numeric(sum(comp_area))
  )

gap_shp <- gap_shp %>%
  left_join(comp_area, by = join_by(trtmnt, comp)) %>%
  left_join(trtmnt_area, by = join_by(trtmnt))

# ==============================================================================
# =========================== Gap area by threshold ============================
# ==============================================================================

# ----------------------------- Data manipulation ------------------------------ 

# Percent coverage by treatment 

trtmnt_gap_p <- gap_shp %>%
  st_drop_geometry() %>%
  group_by(trtmnt, gap_thresh, trtmnt_area) %>%
  summarize(
    area = sum(gap_area),
    p_area = sum(gap_area)/median(trtmnt_area)
  ) %>%
  mutate(gap_thresh = as.numeric(gap_thresh))

# Percent coverage by compartment

comp_gap_p <- gap_shp %>%
  st_drop_geometry() %>%
  group_by(trtmnt, comp, gap_thresh) %>%
  summarize(
    p_area = sum(gap_area)/median(as.numeric(comp_area))) %>%
  mutate(gap_thresh = as.numeric(gap_thresh),
         comp = as.factor(comp)) %>%
  ungroup()

# --------------------------------- Statistics --------------------------------- 

# Mixed effect model with gap_thresh and compartment as effects

mod_pct <- lmer(p_area ~ trtmnt * gap_thresh + (gap_thresh | comp),
                data = comp_gap_p)
mod_pct
anova(mod_pct, type = 3)

# Slope comparison

slope_trt <- emtrends(mod_pct, ~ trtmnt, var = "gap_thresh")
pairs(slope_trt, adjust = "tukey")


#  Predicted intercepts at height 1 m and 15 m

extent_trt <- emmeans(mod_pct, ~ trtmnt | gap_thresh, 
                       at = list(gap_thresh = c(1,15)))
pairs(extent_trt, adjust = "tukey")

# conditional and marginal R2
performance::r2_nakagawa(mod_pct)



# ------------------------------------ Plots ----------------------------------- 

# Percent coverage by treatment 

ggplot(
  data = trtmnt_gap_p,
  mapping = aes(
    x = gap_thresh,
    y = p_area,
    fill = trtmnt,
    color = trtmnt
  )
) + 
  geom_line(linewidth = 1,
            alpha = 0.8) +
  geom_point(size = 2) +
  labs(y = 'Proportion of study area',
       x = 'Height cutoff (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) +
  theme(legend.position = 'bottom') +
  lims(x = c(0,15))

ggsave(glue('{figure_output}/gaps_by_thresh_trtmnt.png'), 
       width = 6, height = 4, units = 'in', dpi = 700)


# Percent coverage by compartment

ggplot(
  data = comp_gap_p,
  mapping = aes(
    x = gap_thresh,
    y = p_area,
    group = interaction(trtmnt, comp),
    fill = trtmnt,
    color = trtmnt
  )
) + 
  geom_line(linewidth = 1,
            alpha = 0.4) +
  geom_point(size = 2) +
  labs(y = 'Proportion of compartment area',
       x = 'Height cutoff (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.65)) +
  theme(legend.position = 'bottom') +
  lims(x = c(0,15))

ggsave(glue('{figure_output}/gaps_by_thresh_comp.png'), 
       width = 6, height = 4, units = 'in', dpi = 700)



# ==============================================================================
# ================================= Gap height =================================
# ==============================================================================


# --------------------------------- Statistics --------------------------------- 

# # AOV - no comp effect
# 
# aov_h = aov(chm_mean ~ trtmnt, data = gap_shp |> filter(gap_thresh == 15))
# TukeyHSD(aov_h)
# 
# mean_groups_aov = tribble(
#   ~trtmnt, ~group,
#   'Control', 'a',
#   'Fire', 'b', 
#   'Mechanical', 'a',
#   'Mech + Fire', 'c'
# )

# Mixed effect model with compartment as random effect

h_mfx  <- lmer(chm_mean ~ trtmnt + (1 | comp),
                data = gap_shp |> filter(gap_thresh == 15))
anova(h_mfx, type = 3)
emmeans(h_mfx, pairwise ~ trtmnt)
performance::icc(h_mfx)

mean_groups_mixedfx = tribble(
  ~trtmnt, ~group,
  'Control', 'a',
  'Fire', 'a, b', 
  'Mech', 'a, b',
  'Mech + Fire', 'b'
)


# ------------------------------------ Plots ----------------------------------- 

ggplot(
  data = gap_shp |>
    filter(gap_thresh == 15),
  mapping = aes(
    x = trtmnt,
    y = chm_mean,
    fill = trtmnt,
    color = trtmnt
  )
) + 
  geom_jitter(alpha = 0.2,
              width = 0.25,
              shape = 20) +
  geom_boxplot(alpha = 0.5,
               color = 'black') +
  geom_text(
    data = mean_groups_mixedfx,
    mapping = aes(
      x = trtmnt,
      y = 15.5,
      label = group
    ),
    inherit.aes = FALSE,
    fontface = 'plain',
    family = 'serif',
    size = 18/ggplot2::.pt
  ) +
  labs(y = 'Gap height (m)',
       x = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  theme(legend.position = 'none')

ggsave(glue('{figure_output}/gap_height_15m_mixedfx.png'), 
       width = 6, height = 4, units = 'in', dpi = 700)


# ==============================================================================
# ============================= Gap heterogeneity ==============================
# ==============================================================================


# --------------------------------- Statistics --------------------------------- 

# # AOV - no comp effect
# 
# aov_gini = aov(chm_gini ~ trtmnt, data = gap_shp |> filter(gap_thresh == 15))
# TukeyHSD(aov_gini)
# 
# gini_groups_aov = tribble(
#   ~trtmnt, ~group,
#   'Control', 'a',
#   'Fire', 'b', 
#   'Mechanical', 'a',
#   'Mech + Fire', 'b'
# )

# Mixed effect model with compartment as random effect

gini_mfx  <- lmer(chm_gini ~ trtmnt + (1 | comp),
                data = gap_shp |> filter(gap_thresh == 15))
gini_mfx
anova(gini_mfx, type = 3)
emmeans(gini_mfx, pairwise ~ trtmnt)
performance::icc(gini_mfx)

gini_groups_mixedfx = tribble(
  ~trtmnt, ~group,
  'Control', 'a',
  'Fire', 'b', 
  'Mech', 'a, b',
  'Mech + Fire', 'b'
)


# ------------------------------------ Plots ----------------------------------- 

ggplot(
  data = gap_shp |> filter(gap_thresh == 15),
  mapping = aes(
    x = trtmnt,
    y = chm_gini,
    fill = trtmnt,
    color = trtmnt
  )
) + 
  geom_jitter(alpha = 0.2,
              width = 0.25,
              shape = 20) +
  geom_boxplot(alpha = 0.5,
               color = 'black') +
  geom_text(
    data = gini_groups_mixedfx,
    mapping = aes(
      x = trtmnt,
      y = 1.05,
      label = group
    ),
    inherit.aes = FALSE,
    fontface = 'plain',
    family = 'serif',
    size = 18/ggplot2::.pt
  ) +
  labs(y = 'Gap Gini coefficient',
       x = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  theme(legend.position = 'none')

ggsave(glue('{figure_output}/gap_gini_15m_mixedfx.png'), 
       width = 6, height = 4, units = 'in', dpi = 700)





# ==============================================================================
# ================================ Density plot ================================
# ==============================================================================

# primary_gap_shp <- gap_shp %>%
#   filter(as.numeric(gap_thresh) %in% primary_gap_sizes)  %>% 
#   mutate(
#     gap_thresh = as.factor(gap_thresh) |>
#       case_match(
#         '1' ~ 'a) Height cutoff: 1 m',
#         '3' ~ 'b) Height cutoff: 3 m',
#         '5' ~ 'c) Height cutoff: 5 m',
#         '10' ~ 'd) Height cutoff: 10 m'
#       ) |>
#       fct_relevel(
#         'a) Height cutoff: 1 m',
#         'b) Height cutoff: 3 m',
#         'c) Height cutoff: 5 m',
#         'd) Height cutoff: 10 m'
#       )
#   )
# 
# 
# ggplot(
#   data = primary_gap_shp,
#   mapping = aes(
#     x = gap_area,
#     fill = trtmnt,
#     color = trtmnt
#   )
# ) +
#   # geom_histogram(alpha = 0.3, linewidth = 1, position = 'identity') +
#   # geom_freqpoly(linewidth = 1, position = 'identity') +
#   # geom_freqpoly(aes(y = after_stat(ncount)), linewidth = 1) +
#   geom_density(linewidth = 1, alpha = 0.05) +
#   labs(y = 'Density',
#        x = expression(Gap~size~(m^2)),
#        fill = 'Treatment',
#        color = 'Treatment') +
#   scale_fill_manual(values = fill_pal) +
#   scale_color_manual(values = col_pal) +
#   theme(legend.position = 'bottom') +
#   # theme(legend.justification.top = "left",
#   #       legend.justification.left = "top",
#   #       legend.justification.bottom = "right",
#   #       legend.justification.inside = c(1, 1),
#   #       legend.location = 'plot',
#   #       legend.position = 'inside') +
#   lemon::facet_rep_wrap(~gap_thresh, ncol = 2) +
#   scale_x_log10()
# 
# ggsave(glue('{figure_output}/gap_size_density_curve.png'), 
#        width = 8, height = 6, units = 'in', dpi = 700)



# # ==============================================================================
# # ======================== Binned normalized area plot =========================
# # ==============================================================================
# 
# 
# log_bin_to_val = function(bin) {
#   
#   val_df = bin %>%
#     str_remove_all(pattern = '[^[:digit:],\\.]') %>%
#     str_split(',', simplify = T) %>%
#     as_tibble %>%
#     rename('log_bin_min' = V1, 'log_bin_max' = V2) %>%
#     mutate(across(everything(), as.numeric)) %>%
#     mutate(
#       log_bin_mid = mean(c_across(everything())),
#       bin_mid = 10^log_bin_mid)
#   
#   return(val_df)
#   
# } 
# 
# gap_binned = primary_gap_shp %>%
#   st_drop_geometry %>%
#   mutate(
#     log_area = log10(gap_area),
#     log_bin = cut_interval(log_area, n = 15)) %>%
#   group_by(log_bin, trtmnt, gap_thresh) %>%
#   summarize(
#     count = n(),
#     scaled_count = n()/median(trtmnt_area) * 1000,
#     area_sum = sum(gap_area),
#     trtmnt_area = median(trtmnt_area),
#     prcnt_area = sum(gap_area)/median(trtmnt_area)) %>%
#   rowwise() %>%
#   mutate(log_bin_to_val(log_bin), .after = 'log_bin')
# 
# ggplot(
#   data = gap_binned,
#   mapping = aes(
#     x = bin_mid,
#     y = prcnt_area,
#     fill = trtmnt,
#     color = trtmnt
#   )
# ) + 
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   labs(y = 'Proportion of study area',
#        x = expression(Gap~size~(m^2)),
#        fill = 'Treatment',
#        color = 'Treatment') +
#   scale_fill_manual(values = fill_pal) +
#   scale_color_manual(values = col_pal) +
#   theme(legend.position = 'bottom') +
#   lemon::facet_rep_wrap(~gap_thresh, ncol = 2) +
#   scale_x_log10()
# 
# ggsave(glue('{figure_output}/gap_area_proportion.png'), 
#        width = 8, height = 6, units = 'in', dpi = 700)
# 
# # ggplot(
# #   data = gap_binned,
# #   mapping = aes(
# #     x = bin_mid,
# #     y = area_sum,
# #     fill = trtmnt,
# #     color = trtmnt
# #   )
# # ) + 
# #   geom_line(linewidth = 1) +
# #   labs(y = 'Cumulative area',
# #        x = expression(Gap~size~(m^2)),
# #        fill = 'Treatment',
# #        color = 'Treatment') +
# #   scale_fill_manual(values = fill_pal) +
# #   scale_color_manual(values = col_pal) +
# #   theme(legend.position = 'bottom') +
# #   lemon::facet_rep_wrap(~gap_thresh, ncol = 2) +
# #   scale_x_log10()
# 
# # ggplot(
# #   data = gap_binned,
# #   mapping = aes(
# #     x = bin_mid,
# #     y = count,
# #     fill = trtmnt,
# #     color = trtmnt
# #   )
# # ) + 
# #   geom_line(linewidth = 1) +
# #   labs(y = 'Count',
# #        x = expression(Gap~size~(m^2)),
# #        fill = 'Treatment',
# #        color = 'Treatment') +
# #   scale_fill_manual(values = fill_pal) +
# #   scale_color_manual(values = col_pal) +
# #   theme(legend.position = 'bottom') +
# #   facet_wrap(~gap_thresh, ncol = 2) +
# #   scale_x_log10()
# 
# ggplot(
#   data = gap_binned,
#   mapping = aes(
#     x = bin_mid,
#     y = scaled_count,
#     fill = trtmnt,
#     color = trtmnt
#   )
# ) + 
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   labs(y = expression(Gaps~per~km^2),
#        x = expression(Gap~size~(m^2)),
#        fill = 'Treatment',
#        color = 'Treatment') +
#   scale_fill_manual(values = fill_pal) +
#   scale_color_manual(values = col_pal) +
#   theme(legend.position = 'bottom') +
#   lemon::facet_rep_wrap(~gap_thresh, ncol = 2) +
#   scale_x_log10()
# 
# ggsave(glue('{figure_output}/gaps_per_sq_km.png'), 
#        width = 8, height = 6, units = 'in', dpi = 700)
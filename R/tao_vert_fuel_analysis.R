# ==============================================================================
#
# TAO Vertical Fuel Continuity
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 16 Dec 2024
# Last commit: 
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

# ==============================================================================

library(tidyverse)
library(sf)
library(quantreg)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(quantreg)
library(Hmisc)
library(permuco)
library(permutes)


# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

metric_file <- 'data/poly_metrics/bldgt_ffs10m_als22_crown_metrics.shp'
study_shp <- 'data/gis/ffs_10m_buffer.shp'

plot_helpers <- 'R/plot_helpers.R'

# ==============================================================================
# ================================== Data prep =================================
# ==============================================================================

source(plot_helpers)

tao_metrics <- metric_file %>%
  st_read() %>%
  mutate(trtmnt = case_match(
    trtmnt,
    'c' ~ 'Control',
    'mf' ~ 'Mech + Fire',
    'f' ~ 'Fire',
    'm' ~ 'Mech'
  ) |>
    fct_relevel(
      'Control',
      'Mech + Fire',
      'Fire', 
      'Mech'
    )) %>%
  mutate(clip_area = st_area(geometry) |>
           units::drop_units(),
         .after = area)  %>%
  mutate(comp = as.integer(comp)) %>%
  filter(z_p95 > 15) %>%                  # Filter to above 15 m TAO height
  filter(lddr_vn > -100 & lddr_vn < 100)  # Filter artifact ladder fuels

comp <- tao_metrics %>%
  st_drop_geometry() %>%
  group_by(comp, trtmnt) %>%
  dplyr::summarize(n = n())

comp_area = st_read(study_shp, quiet = T) %>%
  mutate(comp_area = st_area(geometry)) %>%
  left_join(comp, by = join_by(COMP == comp)) %>%
  rename(comp = COMP) %>%
  st_drop_geometry() %>%
  select(trtmnt, comp, comp_area)

trtmnt_area <- comp_area %>%
  group_by(trtmnt) %>%
  dplyr::summarize(
    trtmnt_area = as.numeric(sum(comp_area))
  )

tao_metrics <- tao_metrics %>%
  left_join(comp_area, by = join_by(trtmnt, comp)) %>%
  left_join(trtmnt_area, by = join_by(trtmnt))%>%
  select(comp, trtmnt, area, clip_area, cc1r_48:trtmnt_area)

# ==============================================================================
# ==============================================================================
# ================ Percent total area vs vertical fuel continuity ============== 
# ==============================================================================
# ==============================================================================

# ------------------------------------------------------------------------------
# --------------------------------- Statistics ---------------------------------
# ------------------------------------------------------------------------------ 

# -------------------------------- LM Mixed FX --------------------------------- 

vertfuel_mixedfx <- lmer(lddr_vn ~ trtmnt + (1|comp),
                 data = tao_metrics,
                 weights = area)

anova(vertfuel_mixedfx, type = 3)

emmeans(vertfuel_mixedfx, pairwise ~ trtmnt)
performance::r2_nakagawa(vertfuel_mixedfx)
performance::icc(vertfuel_mixedfx)

# ----------------------------- Permutation MANOVA ----------------------------- 

bin_width = 0.05

vertfuel_perm_bin <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(bin = cut_width(lddr_vn, width = bin_width)) %>%
  group_by(comp, trtmnt, bin, comp_area) %>%
  dplyr::summarise(
    tao_area = sum(area),
    .groups = "drop") %>%
  group_by(comp, trtmnt) %>%
  mutate(prop_tao = tao_area / sum(tao_area),
         prop_comp = tao_area / median(comp_area) |>
           as.numeric()) %>%
  ungroup()

per_manova <- function(prop_var) {
  
  curve_mat <- tao_perm_bin %>%
    select(comp, trtmnt, bin, !!sym(prop_var)) %>%
    pivot_wider(
      names_from = bin,
      values_from = !!sym(prop_var),
      values_fill = 0
    ) %>%
    arrange(trtmnt)
  
  Y <- curve_mat %>%
    select(-comp, -trtmnt) %>%
    as.matrix
  
  groups  <- curve_mat$trtmnt
  
  # permutation MANOVA (compartment = unit of permutation)
  adon <- adonis2(Y ~ groups, permutations = 9999, method = "euclidean")
  
  return(adon)
  
}

# Proportional coverage of TAO area
'prop_tao' %>%
  per_manova()


# Proportional coverage of full compartment
'prop_comp' %>%
  per_manova()


# -------------------------- Permutation anova of height proportion ---------------------------

mean_vert_fuel = wtd.mean(tao_metrics$lddr_vn, tao_metrics$area)
mean_vert_fuel

# Full canopy

set.seed(942651)
pcomp_perm_full <- lmperm(
  prop_comp ~ trtmnt,
  data = vertfuel_perm_bin,
  np = 9999)

summary(pcomp_perm_full)

# Upper canopy

vertfuel_perm_bin_upper <- tao_metrics %>%
  st_drop_geometry() %>%
  filter(lddr_vn > mean_vert_fuel) %>%
  mutate(bin = cut_width(lddr_vn, width = bin_width)) %>%
  group_by(comp, trtmnt, bin, comp_area) %>%
  dplyr::summarise(
    tao_area = sum(area),
    .groups = "drop") %>%
  group_by(comp, trtmnt) %>%
  mutate(prop_tao = tao_area / sum(tao_area),
         prop_comp = tao_area / median(comp_area) |>
           as.numeric()) %>%
  ungroup()

set.seed(942651)
pcomp_perm_upper <- lmperm(
  prop_comp ~ trtmnt,
  data = vertfuel_perm_bin_upper,
  np = 9999)

summary(pcomp_perm_upper)

# Lower canopy

vertfuel_perm_bin_lower <- tao_metrics %>%
  st_drop_geometry() %>%
  filter(lddr_vn <= mean_vert_fuel) %>%
  mutate(bin = cut_width(lddr_vn, width = bin_width)) %>%
  group_by(comp, trtmnt, bin, comp_area) %>%
  dplyr::summarise(
    tao_area = sum(area),
    .groups = "drop") %>%
  group_by(comp, trtmnt) %>%
  mutate(prop_tao = tao_area / sum(tao_area),
         prop_comp = tao_area / median(comp_area) |>
           as.numeric()) %>%
  ungroup()

set.seed(942651)
pcomp_perm_lower <- lmperm(
  prop_comp ~ trtmnt,
  data = vertfuel_perm_bin_lower,
  np = 9999)

summary(pcomp_perm_lower)

# ---------------------------- Quantile regression ----------------------------- 

tau <- c(0.25, 0.5, 0.75, 0.95)

quant_reg <- rq(
  lddr_vn ~ trtmnt, 
  tau = tau,
  weights = area,
  data = tao_metrics, 
)

set.seed(837294)
bstrap_quant_reg <- summary(
  quant_reg,
  se = 'boot', 
  R = 3000,
  cluster = tao_metrics$comp
)

bstrap_quant_reg

# ------------------------------- IQR Perm anova ------------------------------- 

iqr_vertfuel_comp <- tao_metrics %>%
  group_by(trtmnt, comp) %>%
  dplyr::summarise(
    q25 = wtd.quantile(lddr_vn, area, 0.25),
    q75 = wtd.quantile(lddr_vn, area, 0.75),
    iqr = q75 - q25,
    .groups = "drop"
  )

set.seed(39413)
iqr_vertfuel_perm <- lmperm(
  iqr ~ trtmnt,
  data = iqr_vertfuel_comp,
  np = 9999)

summary(iqr_perm)


# ------------------------------------------------------------------------------
# ------------------------------------ Plots -----------------------------------
# ------------------------------------------------------------------------------


# --------------------- Proportion of total treatment area --------------------- 

bin_vertfuel_trtmnt <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(
    # bin = cut_interval(z_p95, n = 15)) %>%
    bin = cut_width(lddr_vn, width = bin_width)) %>%
  group_by(bin, trtmnt) %>%
  dplyr::summarize(
    area_sum = sum(area),
    trtmnt_area = median(trtmnt_area),
    prcnt_area_trtmnt = sum(area)/median(trtmnt_area)) %>%
  mutate(bin_to_val(bin))

ggplot(
  data = bin_vertfuel_trtmnt,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area_trtmnt,
    fill = trtmnt,
    color = trtmnt
  )
) +
  geom_point(size = 2) +
  geom_line(linewidth = 1,
            alpha = 0.8) +
  labs(y = 'Proportion of study area',
       x = 'Vertical fuel continuity',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent) +
  # lims(x = c(10,60)) +
  theme(legend.position = 'bottom')

ggsave('figures/vertical_fuel_continuity_trtmnt.png', 
       width = 6, height = 4, units = 'in', dpi = 700)


# --------------------- Proportion of per compartment area --------------------- 

comp_vertfuel_bin_height <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(
    across(c(area, comp_area), as.numeric),
    bin = cut_width(lddr_vn, width = bin_width)) %>%
  group_by(bin, trtmnt, comp) %>%
  dplyr::summarize(
    area_sum = sum(area),
    prcnt_area_comp = sum(area)/median(comp_area)) %>%
  mutate(bin_to_val(bin)) %>%
  ungroup()

ggplot(
  data = comp_vertfuel_bin_height,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area_comp,
    fill = trtmnt,
    color = trtmnt,
    group = comp
  )
) +
  geom_line(
    data = comp_vertfuel_bin_height |>
      filter(trtmnt == 'Control') |>
      select(-trtmnt),
    mapping = aes(
      x = bin_mid,
      y = prcnt_area_comp,
      group = comp
    ),
    alpha = 0.8,
    linewidth = 0.7,
    color = 'grey70',
    linetype = 'dashed',
    inherit.aes = FALSE) +
  geom_point(
    size = 2) +
  geom_line(
    linewidth = 1,
    alpha = 0.8) +
  labs(y = 'Proportion of study area',
       x = 'TAO height (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent) +
  # lims(x = c(10,60)) +
  theme(legend.position = 'bottom') + 
  lemon::facet_rep_wrap(~trtmnt)

ggsave('figures/vertfuel_by_treatment_cover_compartment.png', 
       width = 6, height = 4, units = 'in', dpi = 700)




# 
# crown = poly_metrics %>%
#   st_drop_geometry() %>%
#   # filter(lddr_vn > -100 & lddr_vn < 100,
#   #        z_p95 > 15) %>%
#   filter(z_p95 > 15) %>%
#   group_by(trtmnt) %>%
#   summarize(crown_area = sum(clip_area))
# 
# gap = gap_shp %>%
#   # filter(lddr_vn > -100 & lddr_vn < 100,
#   #        z_p95 > 15) %>%
#   filter(gap_thresh == 15) %>%
#   group_by(trtmnt) %>%
#   summarize(gap_area = sum(gap_area))
# 
# x = left_join(crown, gap) %>%
#   mutate(tot = crown_area + gap_area)




ggplot(
  data = ladder_bin_metrics,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area,
    fill = trtmnt,
    color = trtmnt
  )
) +
  geom_point(size = 2) +
  geom_line(linewidth = 1,
            alpha = 0.8) +
  labs(y = 'Proportion of study area',
       x = 'Vertical fuel continuity',
       fill = 'Treatment',
       color = 'Treatment') +
  lims(x = c(0, 1.0005))+
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme(legend.position = 'bottom')

ggsave('figures/vertical_fuels_by_treatment_cover.png', 
       width = 6, height = 4, units = 'in', dpi = 700)






ladder_bin_metrics %>%
  group_by(trtmnt) %>%
  slice_max(prcnt_area) %>%
  select(trtmnt, prcnt_area, bin_mid)





# ggplot(
#   data = ladder_bin_metrics,
#   mapping = aes(
#     x = bin_mid,
#     y = scaled_count,
#     fill = trtmnt,
#     color = trtmnt
#   )
# ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   labs(y = expression(TAO~per~km^2),
#        x = 'Ladder Fuels',
#        fill = 'Treatment',
#        color = 'Treatment') +
#   scale_fill_manual(values = fill_pal) +
#   scale_color_manual(values = col_pal) +
#   theme(legend.position = 'bottom')
# 
# ggsave('figures/ladder_fuels_by_treatment_count.png', 
#        width = 7, height = 5, units = 'in', dpi = 700)

# ladder_bin_metrics %>%
#   group_by(trtmnt) %>%
#   summarize(
#     sum(prcnt_area),
#     sum(area_sum))
# 
# p = poly_metrics %>%
#   st_drop_geometry() %>%
#   # filter(lddr_vn > -100 & lddr_vn < 100,
#   #        z_p95 > 15) %>%
#   group_by(trtmnt) %>%
#   summarize(
#     p_area = sum(clip_area))
# 
# g = gap_binned %>%
#   group_by(trtmnt) %>%
#   filter(gap_thresh == 'd) Height cutoff: 10 m') %>%
#   summarize(
#     g_area = sum(area_sum))
# 
# x = left_join(p, g) %>%
#   rowwise() %>%
#   mutate(tot = p_area + g_area) %>%
#   left_join(study_area) %>%
#   mutate(perc = tot/trtmnt_area)





height_bin_metrics <- poly_metrics %>%
  st_drop_geometry() %>%
  filter(z_p95 > 15) %>%
  left_join(study_area, by = join_by(trtmnt)) %>%
  mutate(
    # bin = cut_interval(z_p95, n = 15)) %>%
    bin = cut_width(z_p95, width = 3)) %>%
  group_by(bin, trtmnt) %>%
  summarize(
    count = n(),
    scaled_count = n()/median(trtmnt_area) * 1000,
    area_sum = sum(area),
    trtmnt_area = median(trtmnt_area),
    prcnt_area = sum(area)/median(trtmnt_area),
    mean_p95 = mean(z_p95)) %>%
  rowwise() %>%
  mutate(bin_to_val(bin))

ggplot(
  data = height_bin_metrics,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area,
    fill = trtmnt,
    color = trtmnt
  )
) +
  geom_point(size = 2) +
  geom_line(linewidth = 1,
            alpha = 0.8) +
  labs(y = 'Proportion of study area',
       x = 'TAO height (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme(legend.position = 'bottom')

ggsave('figures/height_by_treatment_cover.png', 
       width = 6, height = 4, units = 'in', dpi = 700)


height_bin_metrics %>%
  group_by(trtmnt) %>%
  slice_max(prcnt_area) %>%
  select(trtmnt, prcnt_area, bin_mid)



############## Shannon vs cover

filtered_poly_metrics <- poly_metrics %>%
  st_drop_geometry() %>%
  filter(lddr_vn > -100 & lddr_vn < 100,
         z_p95 > 15) %>%
  left_join(study_area, by = join_by(trtmnt))

df = filtered_poly_metrics %>%
  rowwise() %>%
  mutate(
    cc_shannon = {
      x = c_across(cc1r_48:cc1r_32)
      x = x[1:max(which(x > 0))]
      ((sum((x / sum(x)) ^ 2) * length(x)) ^ -1 )
      },
    cc_mean = {
      x = c_across(cc1r_48:cc1r_32)
      x = x[1:max(which(x > 0))]
      mean(x)
    },
    cc_even = cc_shannon * cc_mean
  ) %>%
  ungroup()





col_pal = c(
  'Ladder Evenness' = 'black',
  'Cover 4 - 8 m' = '#D55E00',
  'Cover 8 - 16 m' = '#CC79A7',
  'Cover 16 - 32 m' = '#009E73',
  'Cover > 32 m' = '#0072B2'
)

fill_pal = c(
  'Control' = 'black',
  'Fire' = '#CC6677',
  'Mechanical' = '#537747',
  'Mech + Fire' = '#557EA2'
)


p = ggplot(
  data = df,
  mapping = aes(
    x = cc_mean,
    y = cc_shannon
    # ,
    # color = trtmnt
  )
) +
  stat_density_2d(
    aes(fill = after_stat(level), group = trtmnt), 
    geom = "polygon", 
    color = 'white',
    alpha = 0.5
    ) +
  # stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
  labs(y = 'Cover diversity index',
       x = 'Mean canopy cover',
       fill = 'Kernel Density') +
  lims(x = c(0, 1),
       y = c(0, 1.05)) +
  # scale_fill_manual(values = fill_pal) +
  # scale_color_manual(values = fill_pal) +
  # scale_fill_gradient(low = 'white', high = 'black') +
  theme(legend.position = 'bottom') +
  lemon::facet_rep_wrap(~trtmnt)


# Extract the built data
pbuilt <- ggplot_build(p)
contours <- pbuilt$data[[1]]

# Identify which PANEL is Control
panel_map <- pbuilt$layout$layout
control_panel <- panel_map %>% filter(as.character(trtmnt) == "Control") %>% pull(PANEL)

# Filter to just the Control panel contours
control_contours <- contours %>%
  filter(PANEL == control_panel)

# Convert to numeric for level comparison
control_contours <- control_contours %>%
  mutate(level_num = as.numeric(as.character(level)))

# Find contour levels
min_level <- min(control_contours$level_num, na.rm = TRUE)
max_level <- max(control_contours$level_num, na.rm = TRUE)
control_contours <- control_contours %>%
  filter(level_num %in% seq(min_level, max_level, length.out = 3)) %>%
  mutate(trtmnt = 'Control')


# Final plot
ggplot(df, aes(x = cc_mean, y = cc_shannon)) +
  # 1. Plot Control contours in background
  geom_polygon(
    data = control_contours |>
      select(-trtmnt),
    aes(x = x, y = y, group = group),
    fill = NA,
    color = "grey80",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  # 2. Foreground filled densities by trtmnt
  stat_density_2d(
    aes(fill = after_stat(level), group = trtmnt),
    geom = "polygon",
    color = "white",
    alpha = 0.5
  ) +
  geom_polygon(
    data = control_contours,
    aes(x = x, y = y, group = group),
    fill = NA,
    color = "black",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  facet_rep_wrap(~trtmnt) +
  labs(
    x = "Mean canopy cover",
    y = "Vertical cover diversity index",
    fill = "Kernel Density"
  ) +
  scale_fill_viridis_c(option = 'magma') +
  theme(legend.position = "bottom") 

ggsave('figures/mean_cover_diveristy_comp.png',
       width = 7, height = 7, units = 'in', dpi = 700)


#####################  Cover classes by treatment

filtered_poly_metrics <- poly_metrics %>%
  st_drop_geometry() %>%
  filter(lddr_vn > -100 & lddr_vn < 100,
         z_p95 > 15) %>%
  left_join(study_area, by = join_by(trtmnt))


summarize_binned_metric <- function(df, variable) {
  
  df %>%
    st_drop_geometry() %>%
    mutate(bin = cut(.data[[variable]], 
                     breaks = seq(0, 1, length.out = 21), 
                     include.lowest = TRUE)) %>%
    group_by(bin, trtmnt) %>%
    summarize(
      prcnt_area = sum(clip_area)/median(trtmnt_area), .groups = "drop",
      mean = mean(.data[[variable]])) %>%
    mutate(bin_mid = map_dbl(as.character(bin), bin_to_val)) %>%
    mutate(variable = variable)
  
}

bin_var = c('lddr_vn', 'cc1r_48', 'cc1_816', 'c1_1632', 'cc1r_32')

all_metrics_summary <- bin_var %>%
  map(function(v) {
    summarize_binned_metric(
      df = filtered_poly_metrics,
      variable = v
    )
  }) %>%
  list_rbind() %>%
  mutate(
    variable = variable |>
      as.factor() |>
      case_match(
        'lddr_vn' ~ 'Ladder Evenness',
        'cc1r_48' ~ 'Cover 4 - 8 m',
        'cc1_816' ~ 'Cover 8 - 16 m',
        'c1_1632' ~ 'Cover 16 - 32 m',
        'cc1r_32' ~ 'Cover > 32 m',
      ) |>
      fct_relevel(
        'Ladder Evenness',
        'Cover 4 - 8 m',
        'Cover 8 - 16 m',
        'Cover 16 - 32 m',
        'Cover > 32 m',
      ))



col_pal = c(
  'Ladder Evenness' = 'black',
  'Cover 4 - 8 m' = '#D55E00',
  'Cover 8 - 16 m' = '#CC79A7',
  'Cover 16 - 32 m' = '#009E73',
  'Cover > 32 m' = '#0072B2'
)

ggplot(
  data = all_metrics_summary,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area,
    fill = variable,
    color = variable
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(y = 'Proportion of study area',
       x = 'Ladder Fuels',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = col_pal) +
  scale_color_manual(values = col_pal) +
  theme(legend.position = 'bottom') +
  facet_wrap(~trtmnt)

ggsave('figures/ladder_fuels_by_treatment_cover.png', 
       width = 7, height = 5, units = 'in', dpi = 700)
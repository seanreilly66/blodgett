# ==============================================================================
#
# TAO height analysis
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

bin_width = 3

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
  filter(z_p95 > 15) %>%                       # Filter to above 15 m TAO height
  mutate(z_p95_c = z_p95 - wtd.mean(z_p95, area)) 

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
  left_join(trtmnt_area, by = join_by(trtmnt))


# ==============================================================================
# ==============================================================================
# ============================= TAO area vs height ============================= 
# ==============================================================================
# ==============================================================================


# ------------------------------------------------------------------------------
# --------------------------------- Statistics --------------------------------- 
# ------------------------------------------------------------------------------

mean_h = wtd.mean(tao_metrics$z_p95, tao_metrics$area)
mean_h

tao_avh_mod <- lmer(log(area) ~ trtmnt * z_p95_c + (1 | comp), data = tao_metrics)
summary(tao_avh_mod)
anova(tao_avh_mod, type = 3)

# slope comparison
slopes <- emtrends(tao_avh_mod, ~ trtmnt, var = 'z_p95_c')
slopes
pairs(slopes, adjust = 'tukey')                          # Tukey contrasts
exp(slopes@linfct %*% fixef(tao_avh_mod)) - 1

# intercepts comparison
intercept <- emmeans(tao_avh_mod, ~ trtmnt | z_p95_c,
                     at = list(z_p95_c = c(0, 40 - mean_h)))
pairs(intercept, adjust = "tukey")

# conditional and marginal R2
performance::r2_nakagawa(tao_avh_mod)
performance::icc(tao_avh_mod)

slice_max(tao_metrics,
          order_by = area,
          by = c(comp))

slice_max(tao_metrics,
          order_by = z_p95,
          by = c(comp))

# ------------------------------------------------------------------------------
# ------------------------------------ Plots -----------------------------------
# ------------------------------------------------------------------------------

# Prediction values

h_range <- seq(from = min(tao_metrics$z_p95_c),
               to   = max(tao_metrics$z_p95_c),
               length = 100)


pred_val <- expand.grid(
  z_p95_c = h_range,
  trtmnt = levels(tao_metrics$trtmnt)) %>%
  mutate(
    log_pred_area = predict(
      tao_avh_mod,
      pick(everything()),
      re.form = NA,
      type = 'response'
    )) %>%
  mutate(
    area = exp(log_pred_area),
    z_p95 = z_p95_c + mean_h
  )

obs_range <- tao_metrics %>%
  group_by(trtmnt) %>%
  summarise(obs_min_z = min(z_p95),
            obs_max_z = max(z_p95), 
            .groups = "drop")

pred_val <- pred_val %>%
  left_join(obs_range, by = "trtmnt") %>%
  filter(z_p95 >= obs_min_z, 
         z_p95 <= obs_max_z) %>%
  select(-obs_min_z, -obs_max_z)   


ggplot(mapping = aes(
  x = z_p95,
  y = area,
  colour = trtmnt,
)) +
  geom_point(data = tao_metrics,
             size = 0.5, 
             alpha = 0.3) +
  geom_line(data = pred_val |>
              filter(trtmnt == 'Control') |>
              select(-trtmnt),
            alpha = 0.8,
            linewidth = 0.5,
            color = 'grey40',
            linetype = 'dashed') +
  geom_line(data = pred_val,
            mapping = aes(
              x = z_p95,
              y = area),
            linewidth = 1.8,
            color = 'white',
            alpha = 0.75,
            inherit.aes = FALSE) +
  geom_line(data = pred_val, linewidth = 0.8) + #, color = 'black') +
  labs(y = bquote('TAO area ( m'^2~')'),
       x = 'TAO height ( m )',
       # fill = 'Treatment',
       color = 'Treatment') +
  # scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(
    breaks = seq(0, 400, 100), 
    labels = c(0, '', 200, '', 400)) +
  # lims(x = c(10,60)) +
  theme(legend.position = 'bottom') +
  lemon::facet_rep_wrap(~trtmnt)

ggsave(
  'figures/TAO_height_x_area.png',
  width = 6,
  height = 4,
  units = 'in',
  dpi = 700
)


# ==============================================================================
# ==============================================================================
# ========================= Percent total area vs height ======================= 
# ==============================================================================
# ==============================================================================


# ------------------------------------------------------------------------------
# --------------------------------- Statistics ---------------------------------
# ------------------------------------------------------------------------------ 

# -------------------------------- LM Mixed FX --------------------------------- 

h_mixedfx <- lmer(z_p95_c ~ trtmnt + (1|comp),
                         data = tao_metrics,
                         weights = area)

anova(h_mixedfx, type = 3)

emmeans(h_mixedfx, pairwise ~ trtmnt)
performance::r2_nakagawa(h_mixedfx)
performance::icc(h_mixedfx)


# ----------------------------- Permutation MANOVA ----------------------------- 

tao_perm_bin <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(bin = cut_width(z_p95, width = bin_width)) %>%
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

# Full canopy

set.seed(942651)
pcomp_perm_full <- lmperm(
  prop_comp ~ trtmnt,
  data = tao_perm_bin,
  np = 9999)

summary(pcomp_perm_full)

# Upper canopy

tao_perm_bin_upper <- tao_metrics %>%
  st_drop_geometry() %>%
  filter(z_p95 > mean_h) %>%
  mutate(bin = cut_width(z_p95, width = bin_width)) %>%
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
  data = tao_perm_bin_upper,
  np = 9999)

summary(pcomp_perm_upper)

# Lower canopy

tao_perm_bin_lower <- tao_metrics %>%
  st_drop_geometry() %>%
  filter(z_p95 <= mean_h) %>%
  mutate(bin = cut_width(z_p95, width = bin_width)) %>%
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
  data = tao_perm_bin_lower,
  np = 9999)

summary(pcomp_perm_lower)

# ---------------------------- Quantile regression ----------------------------- 

tau <- c(0.25, 0.5, 0.75, 0.95)

quant_reg <- rq(
  z_p95 ~ trtmnt, 
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

iqr_comp <- tao_metrics %>%
  group_by(trtmnt, comp) %>%
  dplyr::summarise(
    q25 = wtd.quantile(z_p95, area, 0.25),
    q75 = wtd.quantile(z_p95, area, 0.75),
    iqr = q75 - q25,
    .groups = "drop"
  )

set.seed(39413)
iqr_perm <- lmperm(
  iqr ~ trtmnt,
  data = iqr_comp,
  np = 9999)

summary(iqr_perm)


# ------------------------------------------------------------------------------
# ------------------------------------ Plots -----------------------------------
# ------------------------------------------------------------------------------


# --------------------- Proportion of total treatment area --------------------- 

trtmnt_bin_height <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(
    # bin = cut_interval(z_p95, n = 15)) %>%
    bin = cut_width(z_p95, width = bin_width)) %>%
  group_by(bin, trtmnt) %>%
  dplyr::summarize(
    count = n(),
    scaled_count = n()/median(trtmnt_area) * 1000,
    area_sum = sum(area),
    trtmnt_area = median(trtmnt_area),
    prcnt_area_trtmnt = sum(area)/median(trtmnt_area),
    mean_p95 = mean(z_p95)) %>%
  mutate(bin_to_val(bin))

ggplot(
  data = trtmnt_bin_height,
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
       x = 'TAO height (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent) +
  # lims(x = c(10,60)) +
  theme(legend.position = 'bottom')

ggsave('figures/height_by_treatment_cover.png', 
       width = 6, height = 4, units = 'in', dpi = 700)

bin_height %>%
  group_by(trtmnt) %>%
  slice_max(prcnt_area) %>%
  select(trtmnt, prcnt_area, bin_mid)




# --------------------- Proportion of per compartment area --------------------- 

comp_bin_height <- tao_metrics %>%
  st_drop_geometry() %>%
  mutate(
    across(c(area, comp_area), as.numeric),
    bin = cut_width(z_p95, width = bin_width)) %>%
  group_by(bin, trtmnt, comp) %>%
  dplyr::summarize(
    area_sum = sum(area),
    prcnt_area_comp = sum(area)/median(comp_area),
    mean_p95 = mean(z_p95)) %>%
  mutate(bin_to_val(bin)) %>%
  ungroup()

ggplot(
  data = comp_bin_height,
  mapping = aes(
    x = bin_mid,
    y = prcnt_area_comp,
    fill = trtmnt,
    color = trtmnt,
    group = comp
  )
) +
  geom_line(
    data = comp_bin_height |>
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
  labs(y = 'Proportion of compartment area',
       x = 'TAO height (m)',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = trtmnt_fill_pal) +
  scale_color_manual(values = trtmnt_col_pal) +
  scale_y_continuous(labels = scales::percent) +
  # lims(x = c(10,60)) +
  theme(legend.position = 'bottom') + 
  lemon::facet_rep_wrap(~trtmnt)

ggsave('figures/height_by_treatment_cover_compartment.png', 
       width = 6, height = 4, units = 'in', dpi = 700)


# ==============================================================================
# ==============================================================================
# ====================== Cumulative distribution function ====================== 
# ==============================================================================
# ==============================================================================

cdf_height <- tao_metrics %>% 
  arrange(trtmnt, z_p95) %>%
  filter(z_p95 > 15) %>%
  group_by(trtmnt) %>% 
  mutate(cum_area = cumsum(area),
         total_area = sum(area),
         trtmnt_area = median(trtmnt_area),
         cdf_tao_ext = cum_area / total_area,
         cdf_trtmnt_ext = cum_area/trtmnt_area
         ) %>%
  ungroup()

ggplot(data = cdf_height,
       mapping = aes(
         x = z_p95,
         y = cdf_tao_ext,
         colour = trtmnt,
         group = trtmnt
       )) +
  geom_hline(
    yintercept = 0.5,
    linetype = 'dashed',
    linewidth = 1,
    color = 'grey70',
    alpha = 0.7
  ) +
  # geom_hline(
  #   yintercept = 0.95,
  #   linetype = 'dashed',
  #   linewidth = 1,
  #   color = 'grey70',
  #   alpha = 0.7
  # ) +
  geom_step(linewidth = 1,
            direction = "hv",
            alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(15, 50)) +
  labs(x = "TAO height (m)", 
       y = "Cumulative cover\nProportion of TAO extent", 
       colour = "Treatment") +
  scale_color_manual(values = trtmnt_col_pal) +
  theme(legend.position = 'bottom')

ggsave(
  'figures/height_cdf_normalized.png',
  width = 6,
  height = 4,
  units = 'in',
  dpi = 700
)

ggplot(data = cdf_height,
       mapping = aes(
         x = z_p95,
         y = cdf_trtmnt_ext,
         colour = trtmnt,
         group = trtmnt
       )) +
  geom_step(linewidth = 1,
            direction = "hv",
            alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(15, 50)) +
  labs(x = "TAO height (m)", 
       y = "Cumulative cover\nProportion of treatment extent", 
       colour = "Treatment") +
  scale_color_manual(values = trtmnt_col_pal) +
  theme(legend.position = 'bottom')

ggsave(
  'figures/height_cdf_trtmntext.png',
  width = 6,
  height = 4,
  units = 'in',
  dpi = 700
)

# ==============================================================================
# ==============================================================================
# ==============================================================================
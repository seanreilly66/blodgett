# ==============================================================================
#
# Crown attribute analysis
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

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

metric_file <- 'data/poly_metrics/bldgt_ffs10m_als22_crown_metrics.shp'



poly_metrics <- 'data/poly_metrics/bldgt_als22_crown_metrics.shp' %>%
  st_read() %>%
  st_drop_geometry() %>%
  mutate(trtmnt = case_match(
    trtmnt,
    'c' ~ 'Control',
    'mf' ~ 'Mech + Fire',
    'f' ~ 'Fire',
    'm' ~ 'Mechanical'
  ) |>
    fct_relevel(
                'Control',
                'Mech + Fire',
                'Fire', 
                'Mechanical'
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
    title = element_text(size = 12.8)
  )
)

col_pal = c(
  'Control' = 'grey20',
  'Fire' = '#CC6677',
  'Mechanical' = '#537747',
  'Mech + Fire' = '#557EA2'
)

fill_pal = c(
  'Control' = 'white',
  'Fire' = '#CC6677',
  'Mechanical' = '#537747',
  'Mech + Fire' = '#557EA2'
)

# ==============================================================================
# ==================== Density plots - metrics by treatment ==================== 
# ==============================================================================

# Ladder fuels by treatment

ggplot(
  data = poly_metrics |>
    filter(lddr_vn > -100 & lddr_vn < 100),
  mapping = aes(
    x = lddr_vn,
    fill = trtmnt,
    color = trtmnt
  )
) + geom_density(alpha = 0.3, linewidth = 1) +
  labs(y = 'Density',
       x = 'Ladder Fuels',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme(legend.justification.top = "left",
        legend.justification.left = "top",
        legend.justification.bottom = "right",
        legend.justification.inside = c(.9, .9),
        legend.location = 'plot',
        legend.position = 'inside')

ggsave('figures/ladder_fuels_by_treatment.png', 
       width = 5, height = 5, units = 'in', dpi = 700)

lddr_vn = aov(lddr_vn ~ trtmnt, data = poly_metrics|>
                filter(lddr_vn > -100 & lddr_vn < 100))
summary(lddr_vn)

TukeyHSD(lddr_vn)



ggplot(
  data = poly_metrics,
  mapping = aes(
    x = area,
    # y = z_p95,
    fill = trtmnt
  )
) + geom_histogram(position = 'jitter')

area = aov(area ~ trtmnt, data = poly_metrics)
summary(area)

TukeyHSD(area)



ggplot(
  data = poly_metrics,
  mapping = aes(
    x = z_p95,
    fill = trtmnt
  )
) + geom_density(alpha = 0.5) +
  facet_wrap(~trtmnt)

z_p95 = aov(z_p95 ~ trtmnt, data = poly_metrics)
summary(z_p95)

TukeyHSD(z_p95)



# Ladder fuels by treatment

ggplot(
  data = poly_metrics |>
    filter(lddr_vn > -100 & lddr_vn < 100),
  mapping = aes(
    x = lddr_vn,
    fill = trtmnt,
    color = trtmnt
  )
) + geom_density(alpha = 0.3, linewidth = 1) +
  labs(y = 'Density',
       x = 'Ladder Fuels',
       fill = 'Treatment',
       color = 'Treatment') +
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme(legend.justification.top = "left",
        legend.justification.left = "top",
        legend.justification.bottom = "right",
        legend.justification.inside = c(.9, .9),
        legend.location = 'plot',
        legend.position = 'inside')

ggsave('figures/ladder_fuels_by_treatment.png', 
       width = 5, height = 5, units = 'in', dpi = 700)

lddr_vn = aov(lddr_vn ~ trtmnt, data = poly_metrics|>
                filter(lddr_vn > -100 & lddr_vn < 100))
summary(lddr_vn)

TukeyHSD(lddr_vn)
 library(lidR)

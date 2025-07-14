library(tidyverse)
library(ggplot2)


load('data/spatial_pattern/k_quart.RData') # contains k_quart_c, k_quart_f, k_quart_m, and k_quart_mf

tidy_k <- function(k_obj) {
  k_obj$which %>%
    as.data.frame() %>%
    rownames_to_column("mark_i") %>%
    pivot_longer(-mark_i, names_to = "mark_j", values_to = "idx") %>%
    mutate(envelope = map(idx, ~ {as.data.frame(k_obj$fns[[.x]])})) %>%
    select(-idx) %>%
    unnest(envelope)
}

spat_df <- ls(pattern = "^k_quart_") %>%
  map(~ {
    
    k_name <- .x
    
    trtmnt <- str_extract(k_name, "(?<=_)[^_]+$")  
    k_obj  <- get(k_name, envir = .GlobalEnv)
    
    tidy_k(k_obj) %>%
      mutate(trtmnt = trtmnt, .before = 1)
  
    }) %>%
  list_rbind()


spat_df <- spat_df %>%
  mutate(
    across(c(obs, theo, hi, lo),
           .fns = ~ sqrt(.x / pi), 
           .names = '{.col}_L')) %>% 
  mutate(
    across(c(obs_L, theo_L, hi_L, lo_L), 
           .fns = ~ .x - r, 
           .names = "{.col}_norm"))

spat_df <- spat_df %>%
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

spat_df <- spat_df %>%
  mutate(across(c(mark_i, mark_j),
                .fns = ~ case_match(
                  .x,
                  '2' ~ 'Low',
                  '3' ~ 'Med',
                  '4' ~ 'High'
                ) |>
                  fct_relevel(
                    'Low',
                    'Med',
                    'High'
                  )))

pairs_to_plot <- tribble(
  ~mark_i, ~mark_j, ~panel,
  'High',  "High",  "a) High → High",
  "Med",   "Med",   "b) Med → Med",
  "Low",   "Low",   "c) Low → Low",
  "High",  "Med",   "d) High → Med",
  "High",  "Low",   "e) High → Low",
  "Med",   "Low",   "f) Med → Low"
)

plot_df <- spat_df %>% 
  inner_join(pairs_to_plot, by = c("mark_i", "mark_j")) %>% 
  mutate(panel = factor(panel, levels = pairs_to_plot$panel))

seg_data <- tibble(
  panel  = factor("c) Low → Low", levels = pairs_to_plot$panel),
  trtmnt = "Mechanical",
  x     = 65,   y     = 11,
  xend  = 175,  yend  = 11
)

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
    strip.text = element_text(size = 16, 
                              hjust = 0, 
                              margin = margin(l = 0,
                                              b = 2))
  )
)

col_pal = c(
  'Control' = 'grey20',
  'Fire' = '#CC6677',
  'Mechanical' = '#537747',
  'Mech + Fire' = '#557EA2'
)


ggplot(
  data = plot_df, 
  mapping = aes(
    x = r, 
    color = trtmnt,
    fill = trtmnt)) +
  geom_segment(
    data = seg_data,
    mapping = aes(
      x = x,
      y = y,
      xend = xend,
      color = trtmnt
    ),
    linetype = 'dashed',
    linewidth = 0.8,
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray40") +
  geom_ribbon(
    mapping = aes(ymin = lo_L_norm, ymax = hi_L_norm), 
    alpha = 0.2, 
    color = NA) +
  geom_line(
    mapping = aes(y = obs_L_norm), 
    linewidth = 0.8, 
    alpha = 0.8) +
  scale_color_manual(values = col_pal, name = 'Treatment') +
  scale_fill_manual(values = col_pal, name = 'Treatment') +
  lemon::facet_rep_wrap(~ panel, nrow = 2
  ) +
  labs(
    x = 'radius (m)',
    y = expression(L[r] - r)
  ) +
  lims(
    x = c(0, 175),
    y = c(-5, 11)
  ) +
  theme(legend.position = 'bottom')

ggsave('figures/spatial_distribution.png', 
       width = 6, height = 5, units = 'in', dpi = 700)
               
write_csv(spat_df, 'data/spatial_pattern/spatstat_df.csv')      

#
als <- list.files(path = 'data/las/als',
                  pattern = '.las$',
                  full.names = T)

mf180 <- str_subset(als, 'c180') %>%
  readLAS() %>%
  filter_duplicates() %>%
  filter_poi(Classification != 7, Classification != 18) %>%
  normalize_height(tin()) 

plot(mf180, bg = "white", size = 4)

kernel <- matrix(1,3,3)
mf180_chm <- rasterize_canopy(mf180, 0.5, p2r(subcircle = 0.2), pkg = "terra") %>%
  terra::focal(w = kernel, fun = median, na.rm = TRUE)
  

col <- height.colors(50)
plot(mf180_chm, col = col)
ws_f <- function(x) {x * 0.1 + 3}




mf180_ttop <- locate_trees(x, lmf(ws_f))

p <- plot(filter_poi(mf180, Z > 3), bg = "white", size = 4)
add_treetops3d(p, mf180_ttop)

algo <- dalponte2016(mf180_chm, mf180_ttop)
mf180_seg <- segment_trees(mf180, algo)

crowns <- algo()
plot(crowns, col = pastel.colors(200))


plot(filter_poi(mf180_seg, Z > 3), bg = "white", size = 4, color = "treeID")


mf180_mtrc <- crown_metrics(mf180_seg, func = .stdtreemetrics, geom = "convex")

plot(mf180_mtrc["Z"], main = "Height")


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

mf180_mtrc <- mf180_mtrc %>%
  filter(npoints > 1000)

ggplot(
  data = as.data.frame(mf180_mtrc),
  aes(x = Z)) +
  geom_histogram(fill = '#77B68A', color = 'black', binwidth = 5, boundary = 0) +
  labs(y = 'Count',
       x = 'Height (m)') +
  scale_x_continuous(breaks = seq(0,50,5))

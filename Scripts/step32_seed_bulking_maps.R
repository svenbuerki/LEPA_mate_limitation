# step32_seed_bulking_maps.R — reproductive status of every population on the
# talk background (same theme as step30_talk_maps.R).
#
#   seed_bulking_status_map  every population at its centroid, sized by
#                            flowering plants (two-year average), coloured by
#                            reproductive type from the seed-by-size test:
#                            mate-limited / reproductively intact /
#                            over-producing / not yet clear; insets: type
#                            icons and a zoom on EO27 (one EO, five populations
#                            of different types, one polygon per population)
#   seed_bulking_model_map   the reproductively intact populations (the model
#                            populations for seed bulking) highlighted and
#                            labelled; all others in light grey
#
# Types and colours come from tables/step32_seed_bulking_populations.tsv
# (written by step32_seed_bulking.py), so maps and plots share one palette.
# Outputs go to figures_seed_bulking/ (kept local: gitignored).
#
# Run from the repository root:  Rscript Scripts/step32_seed_bulking_maps.R

OUT <- "figures_seed_bulking"
pdf(NULL)   # null device: building grobs (zoom inset) must not write Rplots.pdf
source("Scripts/srk_map_background.R")

# order of Figures 2-3 and of the icon inset
TYPE_ORDER <- c("Reproductively intact", "Mate-limited", "Over-producing", "Not yet clear",
                "One season only")
st <- read.delim("tables/step32_seed_bulking_populations.tsv")
st$candidate <- st$candidate %in% c("True", TRUE)
TYPE_COLORS <- setNames(st$color, st$type)[!duplicated(st$type)][TYPE_ORDER]

cent <- do.call(rbind, lapply(sort(unique(ev_sf$populationID)), function(p) {
  xy <- st_coordinates(st_centroid(st_union(ev_sf[ev_sf$populationID == p, ])))
  data.frame(populationID = p, X = xy[1], Y = xy[2])
}))
cent <- merge(merge(cent, pops[, c("populationID", "avg_fertile")], by = "populationID"),
              st[, c("populationID", "type", "candidate")], by = "populationID")
cent$type <- factor(cent$type, levels = TYPE_ORDER)
cent <- cent[order(as.integer(cent$type), decreasing = TRUE), ]   # greys drawn first
n_of <- table(cent$type)
lab_of <- setNames(sprintf("%s (%d)", TYPE_ORDER, as.integer(n_of[TYPE_ORDER])), TYPE_ORDER)

# BL footprints as dashed outlines in the BL colours (as on population_EO_map_BL)
bls <- sort(unique(ev_sf$BL))
bl_hulls <- do.call(rbind, lapply(bls, function(bl) {
  st_sf(BL = bl, geometry = st_buffer(st_convex_hull(st_union(ev_sf[ev_sf$BL == bl, ])), 1500))
}))
bl_layers <- function() list(
  geom_sf(data = bl_hulls, aes(color = BL), fill = NA, linetype = "longdash", linewidth = 0.8,
          key_glyph = "path"),
  scale_color_manual(values = BL_COLORS, name = "Bottleneck lineage",
                     guide = guide_legend(order = 3, override.aes = list(linewidth = 1.2))),
  new_scale_color())

size_scale <- scale_size_area(name = "Flowering plants\n(two-year average)", max_size = 14,
                              breaks = c(10, 100, 500, 1000), labels = scales::comma)

# ---------------------------------------------------------------- 1. all types
typed <- cent[cent$type != "One season only", ]
typed <- typed[order(typed$avg_fertile, decreasing = TRUE), ]   # small circles on top (P13 under P16)
p_status <- ggplot() + base_layers(dem_reg, bb_reg) + bl_layers() +
  geom_point(data = cent[cent$type == "One season only", ], aes(X, Y, size = avg_fertile,
             fill = type), shape = 21, color = "grey60", stroke = 0.4, alpha = 0.85) +
  geom_point(data = typed[!typed$candidate, ], aes(X, Y, size = avg_fertile, fill = type),
             shape = 21, color = "grey10", stroke = 0.6, alpha = 0.95) +
  geom_point(data = typed[typed$candidate, ], aes(X, Y, size = avg_fertile),
             shape = 21, fill = "white", color = TYPE_COLORS[["Mate-limited"]], stroke = 2) +
  scale_fill_manual(values = TYPE_COLORS, labels = lab_of, drop = FALSE,
                    name = "Seed set for plant size\n(populations; hollow =\ncandidate mate-limited)") +
  size_scale +
  city_layers(bb_reg) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 6, color = "grey10",
                                                            stroke = 0.6)),
         size = guide_legend(order = 2, override.aes = list(fill = "grey75", color = "grey10",
                                                            stroke = 0.6))) +
  label_layer(data.frame(X = typed$X, Y = typed$Y, label = paste0("P", typed$populationID),
                         kind = "pop"), city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
# zoom on EO27: one EO, five populations of different types (Seed_bulking § 6.4) —
# populations as polygons (convex hull of their slick spots + 60 m, to stay visible), the EO as a
# dashed hull (as on population_EO_map_BL), a box marking the zoom on the main map
ZOOM_EO <- "EO27"
zev <- ev_sf[ev_sf$EO == ZOOM_EO, ]
zpop <- do.call(rbind, lapply(sort(unique(zev$populationID)), function(p) {
  g <- ev_sf[ev_sf$populationID == p, ]
  st_sf(populationID = p, geometry = st_buffer(st_convex_hull(st_union(g)), 60))
}))
zpop <- merge(zpop, st[, c("populationID", "type", "candidate")], by = "populationID")
zpop$type <- factor(zpop$type, levels = TYPE_ORDER)
zeo <- st_buffer(st_convex_hull(st_union(zev)), 150)
zbb <- st_bbox(st_buffer(zeo, 250))
zlab <- cbind(st_drop_geometry(zpop), st_coordinates(st_point_on_surface(zpop)))
zlab$label <- paste0("P", zlab$populationID)
sb_x <- zbb[["xmin"]] + 0.05 * (zbb[["xmax"]] - zbb[["xmin"]])     # 1 km scale bar
sb_y <- zbb[["ymin"]] + 0.07 * (zbb[["ymax"]] - zbb[["ymin"]])
p_zoom <- ggplot() +
  geom_sf(data = zeo, fill = NA, color = "grey20", linewidth = 0.5, linetype = "22") +
  geom_sf(data = zpop, aes(fill = type), color = "grey10", linewidth = 0.5, alpha = 0.9) +
  scale_fill_manual(values = TYPE_COLORS, guide = "none") +
  scale_color_manual(values = TYPE_COLORS, guide = "none") +
  geom_label_repel(data = zlab, aes(X, Y, label = label), size = 3.2, fontface = "bold",
                   label.size = 0.2, min.segment.length = 0, box.padding = 0.9,
                   point.padding = 0.4, force = 3, max.overlaps = Inf, seed = 2) +
  annotate("segment", x = sb_x, xend = sb_x + 1000, y = sb_y, yend = sb_y, linewidth = 1) +
  annotate("text", x = sb_x + 500, y = sb_y, label = "1 km", vjust = -0.6, size = 3) +
  coord_sf(crs = UTM, xlim = zbb[c("xmin", "xmax")], ylim = zbb[c("ymin", "ymax")],
           expand = FALSE, datum = NA) +
  labs(title = sprintf("%s (dashed): one EO, %d populations", ZOOM_EO, nrow(zpop))) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", color = "grey30", linewidth = 0.6),
        panel.background = element_rect(fill = "#FAF7F0", color = NA),
        plot.margin = margin(4, 4, 4, 4))
# zoom drawn in map coordinates under the icon inset; the box marks its area
ZOOM_POS <- c(left = 0.465, bottom = 0.49, right = 0.955, top = 0.785)   # map fractions
pp <- ggplot_build(p_status)$layout$panel_params[[1]]      # the map's real extent
at_x <- function(f) pp$x_range[1] + f * diff(pp$x_range)
at_y <- function(f) pp$y_range[1] + f * diff(pp$y_range)
p_status <- p_status + geom_sf(data = st_as_sfc(zbb), fill = NA, color = "grey10",
                               linewidth = 0.6, inherit.aes = FALSE) +
  annotation_custom(ggplotGrob(p_zoom), xmin = at_x(ZOOM_POS[["left"]]),
                    xmax = at_x(ZOOM_POS[["right"]]), ymin = at_y(ZOOM_POS[["bottom"]]),
                    ymax = at_y(ZOOM_POS[["top"]])) +
  coord_sf(crs = UTM, expand = FALSE)
# reproductive-type icons where population_EO_map_BL has its inset (type_icons.py)
icons <- grid::rasterGrob(png::readPNG(file.path(OUT, "type_icons_strip.png")))
icons <- wrap_elements(full = icons) +
  theme(plot.background = element_rect(fill = NA, colour = NA))   # no white box behind the panel
p_status <- p_status + inset_element(icons, left = 0.53, bottom = 0.755, right = 0.97,
                                     top = 0.98, align_to = "panel")   # ~12 % larger than Fig 7's
save_map(p_status, "seed_bulking_status_map")

# ---------------------------------------------------------------- 2. model populations
model <- cent[cent$type == "Reproductively intact", ]
other <- cent[cent$type != "Reproductively intact", ]
p_model <- ggplot() + base_layers(dem_reg, bb_reg) + bl_layers() +
  geom_point(data = other, aes(X, Y, size = avg_fertile), shape = 21, fill = "grey85",
             color = "grey55", stroke = 0.4, alpha = 0.8) +
  geom_point(data = model, aes(X, Y, size = avg_fertile), shape = 21,
             fill = TYPE_COLORS[["Reproductively intact"]], color = "grey10", stroke = 0.7) +
  size_scale +
  city_layers(bb_reg) +
  guides(size = guide_legend(order = 1,
                             override.aes = list(fill = TYPE_COLORS[["Reproductively intact"]]))) +
  label_layer(data.frame(X = model$X, Y = model$Y, label = paste0("P", model$populationID),
                         kind = "pop"), city_anchors(bb_reg)) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("Reproductively intact over two seasons (%d) — candidate sources for seed bulking; gray = other populations",
                          nrow(model))) +
  coord_sf(crs = UTM, expand = FALSE) + map_theme +
  theme(plot.subtitle = element_text(size = 12, face = "bold",
                                     color = TYPE_COLORS[["Reproductively intact"]]))
save_map(p_model, "seed_bulking_model_map")

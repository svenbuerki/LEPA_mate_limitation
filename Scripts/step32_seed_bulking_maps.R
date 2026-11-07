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
#   seed_bulking_status_map_talk  the same without the EO27 zoom (talk; the
#                            close-ups are BL1_BL3_seed_bulking_status_map)
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
# populations seen in one season with no seed data yet are not in the type table;
# they join "One season only" (gray), so every population of an EO appears
TYPE_ORDER <- c("Reproductively intact", "Mate-limited", "Over-producing", "Not yet clear",
                "One season only")
st <- read.delim("tables/step32_seed_bulking_populations.tsv")
st$candidate <- st$candidate %in% c("True", TRUE)
TYPE_COLORS <- setNames(st$color, st$type)[!duplicated(st$type)][TYPE_ORDER]
with_type <- function(d) {                     # attach type; no seed data yet → one season only
  d <- merge(d, st[, c("populationID", "type", "candidate")], by = "populationID", all.x = TRUE)
  d$type[is.na(d$type)] <- "One season only"
  d$candidate[is.na(d$candidate)] <- FALSE
  d$type <- factor(d$type, levels = TYPE_ORDER)
  d
}

cent <- do.call(rbind, lapply(sort(unique(ev_sf$populationID)), function(p) {
  xy <- st_coordinates(st_centroid(st_union(ev_sf[ev_sf$populationID == p, ])))
  data.frame(populationID = p, X = xy[1], Y = xy[2])
}))
cent <- with_type(merge(cent, pops[, c("populationID", "avg_fertile")], by = "populationID"))
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
untyped <- cent$type == "One season only"
# candidate sources: one season so far (2025), significantly within half – twice
# (equivalence test), and 2026 seed collected but not yet cleaned (Seed_bulking § 7.4.4)
cand_src <- st$populationID[st$type == "One season only" & st$n_plants_2026 == 0 &
                              !is.na(st$qeq_2025) & st$qeq_2025 < 0.05 &
                              !is.na(st$mothers_collected_2026) & st$mothers_collected_2026 > 0]
cand <- cent[cent$populationID %in% cand_src, ]
typed <- cent[!untyped, ]
typed <- typed[order(typed$avg_fertile, decreasing = TRUE), ]   # small circles on top (P13 under P16)
p_status <- ggplot() + base_layers(dem_reg, bb_reg) + bl_layers() +
  geom_point(data = cent[untyped, ], aes(X, Y, size = avg_fertile,
             fill = type), shape = 21, color = "grey60", stroke = 0.4, alpha = 0.85) +
  geom_point(data = cand, aes(X, Y, size = avg_fertile), shape = 21, fill = "grey85",
             color = TYPE_COLORS[["Reproductively intact"]], stroke = 2.2) +
  geom_point(data = typed[!typed$candidate, ], aes(X, Y, size = avg_fertile, fill = type),
             shape = 21, color = "grey10", stroke = 0.6, alpha = 0.95) +
  geom_point(data = typed[typed$candidate, ], aes(X, Y, size = avg_fertile),
             shape = 21, fill = "white", color = TYPE_COLORS[["Mate-limited"]], stroke = 2) +
  scale_fill_manual(values = TYPE_COLORS, labels = lab_of, drop = FALSE,
                    name = "Seed set for plant size\n(populations; hollow =\ncandidate mate-limited;\nteal ring = candidate source,\n2026 seed being cleaned)") +
  size_scale +
  city_layers(bb_reg) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 6, color = "grey10",
                                                            stroke = 0.6)),
         size = guide_legend(order = 2, override.aes = list(fill = "grey75", color = "grey10",
                                                            stroke = 0.6))) +
  label_layer(data.frame(X = c(typed$X, cand$X), Y = c(typed$Y, cand$Y),
                         label = paste0("P", c(typed$populationID, cand$populationID)),
                         kind = "pop"), city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
# reproductive-type icons where population_EO_map_BL has its inset (type_icons.py)
icons <- grid::rasterGrob(png::readPNG(file.path(OUT, "type_icons_strip.png")))
icons <- wrap_elements(full = icons) +
  theme(plot.background = element_rect(fill = NA, colour = NA))   # no white box behind the panel
add_icons <- function(p) p + inset_element(icons, left = 0.53, bottom = 0.755, right = 0.97,
                                           top = 0.98, align_to = "panel")   # ~12 % larger than Fig 7's
# talk version: no EO27 zoom (the close-ups are their own slide, step32_BL_pair_map.R)
save_map(add_icons(p_status), "seed_bulking_status_map_talk")

# zoom on EO27: one EO, five populations of different types (Seed_bulking § 6.4) —
# populations as polygons (convex hull of their slick spots + 60 m, to stay visible), the EO as a
# dashed hull (as on population_EO_map_BL), a box marking the zoom on the main map
ZOOM_EO <- "EO27"
zev <- ev_sf[ev_sf$EO == ZOOM_EO, ]
zpop <- do.call(rbind, lapply(sort(unique(zev$populationID)), function(p) {
  g <- ev_sf[ev_sf$populationID == p, ]
  st_sf(populationID = p, geometry = st_buffer(st_convex_hull(st_union(g)), 60))
}))
zpop <- with_type(zpop)
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
  labs(title = sprintf("%s (dashed): %d populations", ZOOM_EO, nrow(zpop))) +
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
p_status <- add_icons(p_status)
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

# ---------------------------------------------------------------- 3. BL1 close-up (talk)
# BL1 at slide format (16:9): every population as a polygon (slick spots,
# convex hull + 80 m so the smallest stay visible) filled by reproductive type,
# EOs as dashed hulls; shows EO27 holding the source P13 next to mate-limited
# populations. Fine-scale per-BL map: kept local (gitignored), like
# BL*_populations_EO_map.
ZOOM_BL <- "BL1"
sub <- ev_sf[ev_sf$BL == ZOOM_BL, ]
bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(sub)), 1500))
w <- bb[["xmax"]] - bb[["xmin"]]; h <- bb[["ymax"]] - bb[["ymin"]]
cx <- mean(bb[c("xmin", "xmax")]); cy <- mean(bb[c("ymin", "ymax")])
h <- max(h, w / 1.45); w <- h * 1.45                    # panel ~1.45:1, legend to the right
bb[c("xmin", "xmax")] <- cx + c(-1, 1) * w / 2
bb[c("ymin", "ymax")] <- cy + c(-1, 1) * h / 2
bb_wgs <- st_bbox(st_transform(st_as_sfc(bb), 4326))
dem_bl <- crop(get_dem(bb_wgs, 11, ZOOM_BL), ext(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]]))
bl_hull <- function(g, buf) st_buffer(st_convex_hull(st_union(g)), buf)
eo_h <- do.call(rbind, lapply(unique(sub$EO), function(eo) {
  st_sf(EO = eo, geometry = bl_hull(ev_sf[ev_sf$EO == eo, ], 150))
}))
pop_h <- do.call(rbind, lapply(sort(unique(sub$populationID)), function(p) {
  st_sf(populationID = p, geometry = bl_hull(sub[sub$populationID == p, ], 80))
}))
pop_h <- with_type(pop_h)
pop_lab <- cbind(st_drop_geometry(pop_h), st_coordinates(st_point_on_surface(pop_h)))
pop_lab$label <- sprintf("P%d", pop_lab$populationID)
eo_lab <- cbind(st_drop_geometry(eo_h),
                as.data.frame(t(sapply(st_geometry(eo_h), function(g) {
                  b <- st_bbox(g); c(X = unname(b["xmax"]), Y = unname(b["ymax"]))
                }))))
n_bl <- table(pop_h$type)
lab_bl <- setNames(sprintf("%s (%d)", TYPE_ORDER, as.integer(n_bl[TYPE_ORDER])), TYPE_ORDER)
p_bl <- ggplot() + base_layers(dem_bl, bb, show_dem_legend = FALSE) +
  geom_sf(data = eo_h, aes(linetype = "Element occurrence (EO)"), fill = NA, color = "grey15",
          linewidth = 0.6) +
  scale_linetype_manual(values = c("Element occurrence (EO)" = "22"), name = NULL,
                        guide = guide_legend(order = 2)) +
  geom_sf(data = pop_h[!pop_h$candidate, ], aes(fill = type), color = "grey10",
          linewidth = 0.5, alpha = 0.95) +
  geom_sf(data = pop_h[pop_h$candidate, ], fill = "white",
          color = TYPE_COLORS[["Mate-limited"]], linewidth = 1.3) +
  scale_fill_manual(values = TYPE_COLORS, labels = lab_bl, drop = TRUE,
                    name = "Population\n(hollow = candidate\nmate-limited)",
                    guide = guide_legend(order = 1)) +
  city_layers(bb, nudge = h / 40, size = 4.2) +
  label_layer(rbind(data.frame(X = eo_lab$X, Y = eo_lab$Y, label = eo_lab$EO, kind = "eo"),
                    data.frame(X = pop_lab$X, Y = pop_lab$Y, label = pop_lab$label, kind = "pop")),
              city_anchors(bb, nudge = h / 40)) +
  coord_sf(crs = UTM, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], expand = FALSE) +
  labs(x = NULL, y = NULL,
       title = sprintf("%s — reproductive type of its %d populations in %d EOs", ZOOM_BL,
                       nrow(pop_h), sum(unique(sub$EO) != "New 2026 population"))) +
  map_theme +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0, margin = margin(b = 8)))
save_map(p_bl, paste0(ZOOM_BL, "_seed_bulking_status_map"), w = 16, h = 9)

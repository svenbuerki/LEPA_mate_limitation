# step32_BL_pair_map.R — populations and EOs of the two bottleneck lineages
# holding EO24 (BL3) and EO27 (BL1), side by side for the seed-bulking talk
# (population vs EO, Seed_bulking § 6). Each panel shows one whole BL:
# populations as polygons (convex hull of their slick spots + 80 m, so the
# smallest stay visible) in the BL colour, EOs as dashed convex hulls (+ 150 m),
# both labelled (EO labels give their number of populations). Fine-scale per-BL map: kept local (gitignored), like
# BL*_populations_EO_map.
#
# Run from the repository root:  Rscript Scripts/step32_BL_pair_map.R

OUT <- "figures_seed_bulking"
pdf(NULL)
source("Scripts/srk_map_background.R")

PAIR <- c("BL1", "BL3")                       # BLs holding EO27 and EO24
PAD <- 1500                                   # m around the slick spots
hull <- function(g, buf) st_buffer(st_convex_hull(st_union(g)), buf)

bl_panel <- function(bl) {
  sub <- ev_sf[ev_sf$BL == bl, ]
  bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(sub)), PAD))
  bb_wgs <- st_bbox(st_transform(st_as_sfc(bb), 4326))
  dem <- crop(get_dem(bb_wgs, 12, paste0(bl, "_pair")),
              ext(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]]))
  eo_h <- do.call(rbind, lapply(unique(sub$EO), function(eo) {
    st_sf(EO = eo, geometry = hull(ev_sf[ev_sf$EO == eo, ], 150))
  }))
  pop_h <- do.call(rbind, lapply(sort(unique(sub$populationID)), function(p) {
    st_sf(populationID = p, BL = bl, geometry = hull(sub[sub$populationID == p, ], 80))
  }))
  pop_lab <- cbind(st_drop_geometry(pop_h), st_coordinates(st_point_on_surface(pop_h)))
  # EO label anchored on the hull outline (point nearest its top-right corner),
  # so the leader line reaches the dashed polygon
  eo_lab <- cbind(st_drop_geometry(eo_h),
                  as.data.frame(t(sapply(st_geometry(eo_h), function(g) {
                    b <- st_bbox(g)
                    corner <- st_sfc(st_point(c(b[["xmax"]], b[["ymax"]])), crs = UTM)
                    xy <- st_coordinates(st_cast(st_nearest_points(
                      st_boundary(st_sfc(g, crs = UTM)), corner), "POINT"))[1, ]
                    c(X = unname(xy[1]), Y = unname(xy[2]))
                  }))))
  n_in <- sapply(eo_lab$EO, function(eo) length(unique(ev_sf$populationID[ev_sf$EO == eo])))
  eo_lab$label <- sprintf("%s\n%d population%s", eo_lab$EO, n_in, ifelse(n_in == 1, "", "s"))
  n_eo <- sum(unique(sub$EO) != "New 2026 population")
  n_pop <- nrow(pop_h)
  p <- ggplot() + base_layers(dem, bb, show_dem_legend = FALSE) +
    geom_sf(data = eo_h, aes(linetype = "Element occurrence (EO)"), fill = NA,
            color = "grey15", linewidth = 0.6) +
    scale_linetype_manual(values = c("Element occurrence (EO)" = "22"), name = NULL) +
    geom_sf(data = pop_h, aes(fill = BL), color = "grey15", linewidth = 0.4,
            alpha = 0.8) +
    scale_fill_manual(values = BL_COLORS[PAIR], guide = "none") +
    city_layers(bb, nudge = 600, size = 4.2) +
    label_layer(rbind(data.frame(X = eo_lab$X, Y = eo_lab$Y, label = eo_lab$label, kind = "eo"),
                      data.frame(X = pop_lab$X, Y = pop_lab$Y,
                                 label = sprintf("P%d", pop_lab$populationID), kind = "pop")),
                city_anchors(bb, nudge = 600)) +
    coord_sf(crs = UTM, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
             expand = FALSE) +
    labs(x = NULL, y = NULL,
         title = sprintf("%s — %d EO%s · %d populations", bl, n_eo,
                         ifelse(n_eo == 1, "", "s"), n_pop)) +
    map_theme +
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0,
                                    color = BL_COLORS[[bl]], margin = margin(b = 6)))
  list(plot = p, aspect = (bb[["xmax"]] - bb[["xmin"]]) / (bb[["ymax"]] - bb[["ymin"]]))
}

panels <- lapply(PAIR, bl_panel)
p <- wrap_plots(lapply(panels, `[[`, "plot"), nrow = 1,
                widths = sapply(panels, `[[`, "aspect")) +
  plot_annotation(caption = "Filled polygons = populations (P), in the color of their bottleneck lineage") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13),
        plot.caption = element_text(size = 13, hjust = 0.5))
save_map(p, "BL1_BL3_populations_EO_map", w = 16, h = 6.4)

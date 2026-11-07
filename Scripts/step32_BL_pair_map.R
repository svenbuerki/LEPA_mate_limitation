# step32_BL_pair_map.R — populations and EOs of the two bottleneck lineages
# holding EO24 (BL3) and EO27 (BL1), side by side for the seed-bulking talk
# (population vs EO, Seed_bulking § 6). Each panel shows one whole BL:
# populations as polygons (convex hull of their slick spots + 80 m, so the
# smallest stay visible) in the BL colour, EOs as dashed convex hulls (+ 150 m),
# both labelled (EO labels give their number of populations). Two versions:
#   BL1_BL3_populations_EO_map       populations in the BL colour
#   BL1_BL3_seed_bulking_status_map  populations by reproductive type (talk
#                                    close-ups of Figure 11; hollow gold =
#                                    candidate mate-limited, teal outline =
#                                    candidate source)
# Fine-scale per-BL maps: kept local (gitignored), like BL*_populations_EO_map.
#
# Run from the repository root:  Rscript Scripts/step32_BL_pair_map.R

OUT <- "figures_seed_bulking"
pdf(NULL)
source("Scripts/srk_map_background.R")

PAIR <- c("BL3", "BL1")                       # BLs holding EO24 and EO27, west to east
PAD <- 1500                                   # m around the slick spots
hull <- function(g, buf) st_buffer(st_convex_hull(st_union(g)), buf)

# reproductive types and colors, as on Figure 11 (step32_seed_bulking_maps.R)
TYPE_ORDER <- c("Reproductively intact", "Mate-limited", "Over-producing", "Not yet clear",
                "One season only")
st <- read.delim("tables/step32_seed_bulking_populations.tsv")
st$candidate <- st$candidate %in% c("True", TRUE)
TYPE_COLORS <- setNames(st$color, st$type)[!duplicated(st$type)][TYPE_ORDER]
cand_src <- st$populationID[st$type == "One season only" & st$n_plants_2026 == 0 &
                              !is.na(st$qeq_2025) & st$qeq_2025 < 0.05 &
                              !is.na(st$mothers_collected_2026) & st$mothers_collected_2026 > 0]
with_type <- function(d) {
  d <- merge(d, st[, c("populationID", "type", "candidate")], by = "populationID", all.x = TRUE)
  d$type[is.na(d$type)] <- "One season only"
  d$candidate[is.na(d$candidate)] <- FALSE
  d$source <- d$populationID %in% cand_src
  d$type <- factor(d$type, levels = TYPE_ORDER)
  d
}

# types present in the two BLs (same legend in both panels)
PAIR_TYPES <- intersect(TYPE_ORDER, as.character(with_type(
  data.frame(populationID = unique(ev$populationID[ev$BL %in% PAIR])))$type))

pop_layers <- function(pop_h, by_type, legend = TRUE) {
  if (!by_type) return(list(
    geom_sf(data = pop_h, aes(fill = BL), color = "grey15", linewidth = 0.4, alpha = 0.8),
    scale_fill_manual(values = BL_COLORS[PAIR], guide = "none")))
  pop_h <- with_type(pop_h)
  plain <- pop_h[!pop_h$candidate & !pop_h$source, ]
  list(
    geom_sf(data = plain, aes(fill = type), color = "grey15", linewidth = 0.4, alpha = 0.95,
            key_glyph = "polygon"),
    geom_sf(data = pop_h[pop_h$source, ], fill = "grey85",
            color = TYPE_COLORS[["Reproductively intact"]], linewidth = 1.4),
    geom_sf(data = pop_h[pop_h$candidate, ], fill = "white",
            color = TYPE_COLORS[["Mate-limited"]], linewidth = 1.3),
    scale_fill_manual(values = TYPE_COLORS, limits = PAIR_TYPES, drop = FALSE,
                      name = "Seed set for plant size",
                      guide = if (legend) guide_legend(order = 1) else "none"))
}

bl_panel <- function(bl, by_type = FALSE) {
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
    scale_linetype_manual(values = c("Element occurrence (EO)" = "22"), name = NULL,
                          guide = guide_legend(order = 2)) +
    pop_layers(pop_h, by_type, legend = bl == "BL1") +   # legend once: BL1 has every type
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

pair_map <- function(by_type, caption, name) {
  panels <- lapply(PAIR, bl_panel, by_type = by_type)
  p <- wrap_plots(lapply(panels, `[[`, "plot"), nrow = 1,
                  widths = sapply(panels, `[[`, "aspect")) +
    plot_annotation(caption = caption) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 13),
          plot.caption = element_text(size = 13, hjust = 0.5))
  save_map(p, name, w = 16, h = 6.4)
}
pair_map(FALSE, "Filled polygons = populations (P), in the color of their bottleneck lineage",
         "BL1_BL3_populations_EO_map")
pair_map(TRUE, paste("Populations (P) by seed set for plant size; hollow gold = candidate",
                     "mate-limited; teal outline = candidate source (2026 seed being cleaned)"),
         "BL1_BL3_seed_bulking_status_map")

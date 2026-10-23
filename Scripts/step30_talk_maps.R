# step30_talk_maps.R — presentation maps on the topographic background of
# LEPA_EO_spatial_clustering (EO_BL_fragmentation_network_map.png):
# elevation from AWS Terrain Tiles (elevatr) in a warm-sepia ramp, the Snake
# River (Natural Earth, scale 10) and four reference cities.
#
#   1. BL fragmentation network — the five population-based BLs (Step 05):
#      nodes at BL centroids sized by flowering plants (two-year average),
#      edges = shortest event-to-event distance between BLs (km), inset
#      lower-triangle distance matrix.
#   2. EO sampling map — one circle per element occurrence, coloured by BL
#      (EO_sampling_map_BL), and the same map with black circles
#      (EO_sampling_map_black) for before clustering,
#      sized by its average number of flowering plants over 2025 + 2026
#      (mean over the years with plants above ground).
#   2b. Population map — every population at its real extent, coloured by
#      BL, EO convex hulls, BL outlines; inset violin of population size.
#   3. One map per BL — populations as polygons (convex hull of their slick
#      spots, buffered 25 m) and each EO as a dashed convex hull (buffered
#      50 m), on the same background zoomed to the BL.
#
# Outputs go to talk_figures/ (kept local: gitignored). Background and data
# setup are shared with the project maps (srk_map_background.R).
#
# Run from the repository root:  Rscript Scripts/step30_talk_maps.R

OUT <- "talk_figures"
source("Scripts/srk_map_background.R")

# ---------------------------------------------------------------- 1. BL network
bl_pts <- ev_sf %>% group_by(BL) %>% summarise(geometry = st_union(geometry), .groups = "drop")
bl_nodes <- cbind(st_drop_geometry(bl_pts),
                  st_coordinates(st_centroid(bl_pts)))
bl_size <- pops %>% group_by(BL) %>%
  summarise(n_pop = n(), avg_fertile = sum(avg_fertile), .groups = "drop")
bl_nodes <- merge(bl_nodes, bl_size, by = "BL")
# Node labels count element occurrences (EOs), not populations: the map
# opens the talk, before populations are introduced. The 2026 site without
# an EO number is counted as a new site.
bl_eo <- ev %>% distinct(BL, EO) %>% group_by(BL) %>%
  summarise(n_eo = sum(EO != "New 2026 population"),
            n_new = sum(EO == "New 2026 population"), .groups = "drop")
bl_nodes <- merge(bl_nodes, bl_eo, by = "BL")
bl_nodes$label <- sprintf("%s\n%d EO%s%s\n~%s plants / yr", bl_nodes$BL,
                          bl_nodes$n_eo, ifelse(bl_nodes$n_eo == 1, "", "s"),
                          ifelse(bl_nodes$n_new > 0, " + 1 new site", ""),
                          trimws(format(round(bl_nodes$avg_fertile, -1), big.mark = ",")))
pairs <- t(combn(sort(bl_nodes$BL), 2))
bl_edges <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  d <- st_distance(ev_sf[ev_sf$BL == a, ], ev_sf[ev_sf$BL == b, ])
  na <- bl_nodes[bl_nodes$BL == a, ]; nb <- bl_nodes[bl_nodes$BL == b, ]
  data.frame(BL_A = a, BL_B = b, dist_km = min(as.numeric(d)) / 1000,
             x_A = na$X, y_A = na$Y, x_B = nb$X, y_B = nb$Y,
             mid_x = (na$X + nb$X) / 2, mid_y = (na$Y + nb$Y) / 2)
}))
bl_edges$closeness <- 1 / bl_edges$dist_km
bl_edges$closeness_n <- (bl_edges$closeness - min(bl_edges$closeness)) /
  diff(range(bl_edges$closeness))

lv <- paste0("BL", 1:5)
hm <- transform(bl_edges, BL_row = factor(BL_B, levels = lv[2:5]),
                BL_col = factor(BL_A, levels = lv[1:4]))
p_heat <- ggplot(hm, aes(BL_col, BL_row, fill = dist_km)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", dist_km)), color = "grey15", size = 2.8,
            fontface = "bold") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, guide = "none") +
  scale_x_discrete(limits = lv[1:4], drop = FALSE) +
  scale_y_discrete(limits = rev(lv[2:5]), drop = FALSE) +
  coord_fixed() + labs(title = "Inter-BL distance (km)") +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(color = BL_COLORS[lv[1:4]], face = "bold", size = 8),
        axis.text.y = element_text(color = BL_COLORS[rev(lv[2:5])], face = "bold", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.5),
        plot.margin = margin(4, 4, 4, 4))

p_net <- ggplot() + base_layers(dem_reg, bb_reg) +
  geom_segment(data = bl_edges, aes(x = x_A, y = y_A, xend = x_B, yend = y_B,
                                    alpha = closeness_n, linewidth = closeness_n),
               color = "grey20") +
  scale_alpha_continuous(range = c(0.25, 0.95), guide = "none") +
  scale_linewidth_continuous(range = c(0.4, 2.8), guide = "none") +
  geom_point(data = bl_nodes, aes(X, Y, size = avg_fertile, fill = BL), shape = 21,
             color = "grey15", stroke = 0.9, alpha = 0.95) +
  scale_size_area(name = "Flowering plants\n(two-year average)", max_size = 22,
                  breaks = c(500, 1000, 2500, 5000), labels = scales::comma) +
  scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage") +
  city_layers(bb_reg) +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 6,
                                                            color = "grey15")),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey75",
                                                            color = "grey15"))) +
  label_layer(rbind(data.frame(X = bl_nodes$X, Y = bl_nodes$Y, label = bl_nodes$label,
                               kind = "node"),
                    data.frame(X = bl_edges$mid_x, Y = bl_edges$mid_y,
                               label = sprintf("%.1f km", bl_edges$dist_km), kind = "edge")),
              city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
p_net <- p_net + inset_element(p_heat, left = 0.53, bottom = 0.66, right = 0.89,
                               top = 0.99, align_to = "panel")
save_map(p_net, "BL_fragmentation_network_map")

# ---------------------------------------------------------------- 2. EO map
eo_pts <- ev_sf %>% group_by(EO) %>% summarise(geometry = st_union(geometry), .groups = "drop")
eo_nodes <- cbind(st_drop_geometry(eo_pts), st_coordinates(st_centroid(eo_pts)))
eo_nodes <- merge(eo_nodes, eo_tab, by = "EO")
p_eo <- ggplot() + base_layers(dem_reg, bb_reg) +
  geom_point(data = eo_nodes, aes(X, Y, size = avg_fertile, fill = BL), shape = 21,
             color = "grey15", stroke = 0.7, alpha = 0.9) +
  scale_size_area(name = "Flowering plants per EO\n(two-year average)", max_size = 16,
                  breaks = c(10, 100, 500, 1000, 2000), labels = scales::comma) +
  scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage") +
  city_layers(bb_reg) +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 6)),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey75"))) +
  label_layer(data.frame(X = eo_nodes$X, Y = eo_nodes$Y, label = eo_nodes$EO, kind = "eo"),
              city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
save_map(p_eo, "EO_sampling_map_BL")

# Same map before clustering: black circles, no BL colours (sampling comes
# first in the talk); EO_sampling_map_BL keeps the BL colours for after clustering
p_eo_plain <- ggplot() + base_layers(dem_reg, bb_reg) +
  geom_point(data = eo_nodes, aes(X, Y, size = avg_fertile), shape = 21,
             fill = "black", color = "white", stroke = 0.8, alpha = 1) +
  scale_size_area(name = "Flowering plants per EO\n(two-year average)", max_size = 16,
                  breaks = c(10, 100, 500, 1000, 2000), labels = scales::comma) +
  city_layers(bb_reg) +
  guides(size = guide_legend(order = 1, override.aes = list(shape = 21, fill = "black",
                                                            color = "white"))) +
  label_layer(data.frame(X = eo_nodes$X, Y = eo_nodes$Y, label = eo_nodes$EO, kind = "eo"),
              city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
save_map(p_eo_plain, "EO_sampling_map_black")

# ---------------------------------------------------------------- 2b. population map
# All populations on the regional background at their real extent (convex
# hull of their slick spots, buffered 25 m), coloured by BL; a thin outline
# in the same colour keeps the smallest visible. EOs as dashed convex hulls,
# each BL's footprint as a dashed outline. Population size (flowering plants,
# two-year average) is shown in an inset violin plot per BL.
pop_poly <- do.call(rbind, lapply(sort(unique(ev_sf$populationID)), function(p) {
  g <- ev_sf[ev_sf$populationID == p, ]
  st_sf(populationID = p, BL = g$BL[1],
        geometry = st_buffer(st_convex_hull(st_union(g)), 25))
}))
eo_hulls <- do.call(rbind, lapply(sort(unique(ev_sf$EO)), function(eo) {
  st_sf(EO = eo, geometry = st_buffer(st_convex_hull(st_union(ev_sf[ev_sf$EO == eo, ])), 400))
}))
bl_hulls <- do.call(rbind, lapply(lv, function(bl) {
  st_sf(BL = bl, geometry = st_buffer(st_convex_hull(st_union(ev_sf[ev_sf$BL == bl, ])), 1500))
}))
eo_hl <- cbind(st_drop_geometry(eo_hulls),
               as.data.frame(t(sapply(st_geometry(eo_hulls), function(g) {
                 b <- st_bbox(g); c(X = unname(b["xmax"]), Y = unname(b["ymax"]))
               }))))
p_violin <- ggplot(pops, aes(BL, avg_fertile, fill = BL)) +
  geom_violin(scale = "width", color = "grey25", linewidth = 0.3, alpha = 0.85, trim = TRUE) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, color = "grey15", alpha = 0.8) +
  scale_fill_manual(values = BL_COLORS, guide = "none") +
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma) +
  labs(title = "Population size", x = NULL,
       y = "Flowering plants\n(two-year average)") +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text.x = element_text(color = BL_COLORS[lv], face = "bold", size = 8),
        axis.title.y = element_text(size = 7.5),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.5),
        plot.margin = margin(4, 6, 4, 4))
p_pop <- ggplot() + base_layers(dem_reg, bb_reg) +
  geom_sf(data = bl_hulls, aes(color = BL), fill = NA, linetype = "longdash",
          linewidth = 0.7, show.legend = FALSE) +
  scale_color_manual(values = BL_COLORS, guide = "none") +
  new_scale_color() +
  geom_sf(data = eo_hulls, aes(linetype = "Element occurrence (convex hull)"),
          fill = NA, color = "grey20", linewidth = 0.45) +
  scale_linetype_manual(values = c("Element occurrence (convex hull)" = "22"),
                        name = NULL, guide = guide_legend(order = 3)) +
  geom_sf(data = pop_poly, aes(fill = BL, color = BL), linewidth = 0.6) +
  scale_fill_manual(values = BL_COLORS, name = "Population, by\nBottleneck lineage") +
  scale_color_manual(values = BL_COLORS, guide = "none") +
  city_layers(bb_reg) +
  guides(fill = guide_legend(order = 1)) +
  label_layer(data.frame(X = eo_hl$X, Y = eo_hl$Y, label = eo_hl$EO, kind = "eo"),
              city_anchors(bb_reg)) +
  coord_sf(crs = UTM, expand = FALSE) + labs(x = NULL, y = NULL) + map_theme
p_pop <- p_pop + inset_element(p_violin, left = 0.58, bottom = 0.60, right = 0.97,
                               top = 0.98, align_to = "panel")
save_map(p_pop, "population_EO_map_BL")

# ---------------------------------------------------------------- 3. per-BL maps
hull <- function(g, buf) st_buffer(st_convex_hull(st_union(g)), buf)
for (bl in lv) {
  sub <- ev_sf[ev_sf$BL == bl, ]
  bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(sub)), 1500))
  side <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"], 4000)
  cx <- mean(bb[c("xmin", "xmax")]); cy <- mean(bb[c("ymin", "ymax")])
  bb[c("xmin", "xmax")] <- cx + c(-1, 1) * side / 2 * 1.35   # wider for the legend
  bb[c("ymin", "ymax")] <- cy + c(-1, 1) * side / 2
  bb_wgs <- st_bbox(st_transform(st_as_sfc(bb), 4326))
  z <- if (side > 40000) 10 else if (side > 15000) 11 else 12
  dem <- crop(get_dem(bb_wgs, z, bl), ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
  eo_h <- do.call(rbind, lapply(unique(sub$EO), function(eo) {
    st_sf(EO = eo, geometry = hull(ev_sf[ev_sf$EO == eo, ], 50))
  }))
  pop_h <- do.call(rbind, lapply(sort(unique(sub$populationID)), function(p) {
    st_sf(populationID = p, geometry = hull(sub[sub$populationID == p, ], 25))
  }))
  pop_h <- merge(pop_h, pops[, c("populationID", "avg_fertile")], by = "populationID")
  pop_lab <- cbind(st_drop_geometry(pop_h), st_coordinates(st_point_on_surface(pop_h)))
  pop_lab$label <- sprintf("P%d", pop_lab$populationID)
  eo_lab <- cbind(st_drop_geometry(eo_h),
                  as.data.frame(t(sapply(st_geometry(eo_h), function(g) {
                    b <- st_bbox(g); c(X = unname(b["xmax"]), Y = unname(b["ymax"]))
                  }))))
  p <- ggplot() + base_layers(dem, bb) +
    geom_sf(data = eo_h, aes(linetype = "Element occurrence (convex hull)"),
            fill = NA, color = "grey20", linewidth = 0.6) +
    scale_linetype_manual(values = c("Element occurrence (convex hull)" = "22"),
                          name = NULL, guide = guide_legend(order = 2)) +
    geom_sf(data = pop_h, aes(fill = "Population"), color = BL_COLORS[[bl]],
            linewidth = 0.5, alpha = 0.75) +
    scale_fill_manual(values = setNames(BL_COLORS[[bl]], "Population"), name = bl,
                      guide = guide_legend(order = 1)) +
    city_layers(bb, nudge = side / 40, size = 4.2) +
    label_layer(rbind(data.frame(X = eo_lab$X, Y = eo_lab$Y, label = eo_lab$EO, kind = "eo"),
                      data.frame(X = pop_lab$X, Y = pop_lab$Y, label = pop_lab$label,
                                 kind = "pop")),
                city_anchors(bb, nudge = side / 40)) +
    coord_sf(crs = UTM, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
             expand = FALSE) +
    labs(x = NULL, y = NULL) + map_theme
  n_eo  <- sum(unique(sub$EO) != "New 2026 population")
  n_new <- sum(unique(sub$EO) == "New 2026 population")
  n_pop <- length(unique(sub$populationID))
  header <- sprintf("%s — %d EO%s%s · %d population%s", bl, n_eo,
                    ifelse(n_eo == 1, "", "s"), ifelse(n_new > 0, " + 1 new site", ""),
                    n_pop, ifelse(n_pop == 1, "", "s"))
  p <- p + labs(title = header) +
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0,
                                    color = BL_COLORS[[bl]], margin = margin(b = 8)))
  save_map(p, paste0(bl, "_populations_EO_map"), w = 13, h = 9)
}
message("[step30] talk maps written to ", OUT, "/")

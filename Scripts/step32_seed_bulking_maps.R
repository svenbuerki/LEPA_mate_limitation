# step32_seed_bulking_maps.R — reproductive status of every population on the
# talk background (same theme as step30_talk_maps.R).
#
#   seed_bulking_status_map  every population at its centroid, sized by
#                            flowering plants (two-year average), coloured by
#                            reproductive type from the seed-by-size test:
#                            mate-limited / reproductively intact /
#                            over-producing / not yet clear
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
source("Scripts/srk_map_background.R")

TYPE_ORDER <- c("Mate-limited", "Reproductively intact", "Over-producing", "Not yet clear",
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
       subtitle = sprintf("Reproductively intact over two seasons (%d) — candidate sources for seed bulking; grey = other populations",
                          nrow(model))) +
  coord_sf(crs = UTM, expand = FALSE) + map_theme +
  theme(plot.subtitle = element_text(size = 12, face = "bold",
                                     color = TYPE_COLORS[["Reproductively intact"]]))
save_map(p_model, "seed_bulking_model_map")

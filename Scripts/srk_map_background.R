# srk_map_background.R — shared map setup for the R map scripts
# (step28_model_population_design_map.R, step30_talk_maps.R).
#
# Background of LEPA_EO_spatial_clustering's EO_BL_fragmentation_network_map:
# elevation from AWS Terrain Tiles (elevatr) in a warm-sepia ramp on one
# regional colour scale, the Snake River (Natural Earth, scale 10) and four
# reference cities. Also loads the slick spots (event coordinates, in memory
# only) with their population, deme, BL and EO, and the label/theme helpers.
# Callers set OUT before sourcing.

suppressPackageStartupMessages({
  library(sf); library(terra); library(elevatr); library(tidyterra)
  library(rnaturalearth); library(ggplot2); library(ggrepel)
  library(ggnewscale); library(patchwork); library(dplyr)
})
sf_use_s2(FALSE)
if (Sys.getenv("PROJ_LIB") == "") {
  sf_proj <- system.file("proj", package = "sf")
  if (dir.exists(sf_proj)) Sys.setenv(PROJ_LIB = sf_proj)
}

# Downloaded elevation tiles and river data (cached, gitignored)
CACHE <- "map_cache"
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
UTM <- 32611

BL_COLORS <- c(BL1 = "#984EA3", BL2 = "#377EB8", BL3 = "#E41A1C",
               BL4 = "#FF7F00", BL5 = "#4DAF4A")       # = srk_bl_constants.py
DEM_COLOURS <- c("#fff8e7", "#f0d9a8", "#d4a574", "#a87545", "#6b4423", "#3d2817")

# ---------------------------------------------------------------- data
cw <- read.delim("tables/step29a_population_crosswalk.tsv")          # local only
units <- read.delim("tables/step29a_event_to_deme.tsv")
pops <- read.delim("tables/step30g_populations_classified.tsv")
ev <- merge(units[, c("event_year", "eventID", "populationID", "BL", "demeID",
                      "locationCode", "n_fertile")],
            cw[, c("event_year", "eventID", "lat", "lon")],
            by = c("event_year", "eventID"))
ev$EO <- ifelse(grepl("^EO[0-9]", ev$locationCode),
                sub("^(EO[0-9]+).*$", "\\1", ev$locationCode), ev$locationCode)
ev$EO[ev$EO == "EO_NEW2026"] <- "New 2026 population"
ev_sf <- st_transform(st_as_sf(ev, coords = c("lon", "lat"), crs = 4326), UTM)

# two-year average size: mean over the years with plants above ground
avg_size <- function(n25, n26) {
  v <- cbind(n25, n26); v[v == 0] <- NA
  rowMeans(v, na.rm = TRUE)
}
pops$avg_fertile <- avg_size(pops$n_fertile_2025, pops$n_fertile_2026)

eo_year <- ev %>% group_by(EO, BL, event_year) %>%
  summarise(n = sum(n_fertile), .groups = "drop")
eo_tab <- eo_year %>% group_by(EO, BL) %>%
  summarise(avg_fertile = mean(n[n > 0]), years = n(), .groups = "drop")

cities <- data.frame(
  name = c("Boise", "Mountain Home", "Glenns Ferry", "New Plymouth"),
  lon  = c(-116.2023, -115.6912, -115.3022, -116.8195),
  lat  = c(  43.6150,   43.1326,   42.9551,   43.9685))
cities_utm <- st_transform(st_as_sf(cities, coords = c("lon", "lat"), crs = 4326), UTM)

# ---------------------------------------------------------------- background
get_dem <- function(bbox_wgs, z, tag) {
  f <- file.path(CACHE, paste0("dem_", tag, "_z", z, ".tif"))
  if (file.exists(f)) return(rast(f))
  bb <- st_as_sf(data.frame(id = 1L),
                 geometry = st_sfc(st_as_sfc(bbox_wgs)[[1]], crs = 4326))
  d <- get_elev_raster(locations = bb, z = z, clip = "bbox", prj = "EPSG:4326")
  d <- project(rast(d), paste0("EPSG:", UTM), method = "bilinear")
  writeRaster(d, f, overwrite = TRUE)
  d
}
snake_file <- file.path(CACHE, "snake_river.gpkg")
if (!file.exists(snake_file)) {
  rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines",
                        category = "physical", returnclass = "sf")
  snake <- rivers[grepl("Snake", rivers$name, ignore.case = TRUE), ]
  st_write(st_transform(snake, UTM), snake_file, quiet = TRUE)
}
snake_utm <- st_read(snake_file, quiet = TRUE)

DEM_LIMITS <- NULL   # set from the regional DEM so every map shares one colour scale
base_layers <- function(dem, bbox_utm, show_dem_legend = TRUE) {
  list(
    geom_spatraster(data = dem, maxcell = 5e5),
    scale_fill_gradientn(colours = DEM_COLOURS, name = "Elevation (m)", na.value = NA,
                         limits = DEM_LIMITS, oob = scales::squish,
                         guide = if (show_dem_legend)
                           guide_colorbar(barheight = unit(3.2, "cm"),
                                          barwidth = unit(0.4, "cm"), order = 4)
                         else "none"),
    geom_sf(data = st_crop(snake_utm, bbox_utm), aes(color = "Snake River"),
            linewidth = 0.8, alpha = 0.85, show.legend = "line"),
    scale_color_manual(values = c("Snake River" = "#1f77b4"), name = "Map features",
                       guide = guide_legend(order = 5,
                                            override.aes = list(linewidth = 1.2))),
    new_scale_fill(), new_scale_color()
  )
}
city_layers <- function(bbox_utm, nudge = 3500, size = 4.8) {
  inside <- st_crop(cities_utm, bbox_utm)
  if (nrow(inside) == 0) return(list())
  cc <- cbind(as.data.frame(st_coordinates(inside)), name = inside$name)
  # Glenns Ferry sits at the edge of the EO8 / EO26 cluster: label it below
  cc$dy <- ifelse(cc$name == "Glenns Ferry", -nudge, nudge)
  list(geom_sf(data = inside, shape = 22, fill = "black", color = "white",
               size = 3.2, stroke = 0.5),
       geom_text_repel(data = cc, aes(x = X, y = Y, label = name), fontface = "bold",
                       size = size, color = "grey15", bg.color = "white", bg.r = 0.15,
                       nudge_y = cc$dy, point.padding = unit(0.5, "lines"),
                       box.padding = unit(0.45, "lines"), min.segment.length = 0.2,
                       segment.color = "grey35", segment.size = 0.3,
                       max.overlaps = Inf, seed = 42))
}
city_anchors <- function(bbox_utm, nudge = 3500) {
  inside <- st_crop(cities_utm, bbox_utm)
  if (nrow(inside) == 0) return(data.frame(X = numeric(), Y = numeric()))
  xy <- as.data.frame(st_coordinates(inside))
  dy <- ifelse(inside$name == "Glenns Ferry", -nudge, nudge)
  rbind(xy, transform(xy, Y = Y + dy))
}
# One repel layer for every label of a map, so labels avoid each other
label_layer <- function(lab, anchors) {
  pts <- rbind(lab[, c("X", "Y", "label", "kind")],
               if (nrow(anchors)) data.frame(anchors, label = "", kind = "anchor"))
  sizes <- c(node = 3.4, edge = 3.0, eo = 3.6, pop = 3.0, anchor = 1)
  faces <- c(node = "bold", edge = "plain", eo = "bold", pop = "bold", anchor = "plain")
  fills <- c(node = "white", edge = "white", eo = "white", pop = NA, anchor = NA)
  bords <- c(node = 0, edge = 0.25, eo = 0.2, pop = 0, anchor = 0)
  list(new_scale_fill(),
       geom_label_repel(data = pts,
                        aes(X, Y, label = label, fill = I(fills[kind]),
                            size = I(sizes[kind]), fontface = faces[kind],
                            label.size = I(bords[kind])),
                        color = "grey15", alpha = 0.92, lineheight = 0.95,
                        box.padding = 0.45, point.padding = 0.3,
                        min.segment.length = 0.1, segment.color = "grey35",
                        segment.size = 0.3, max.overlaps = Inf, seed = 42,
                        show.legend = FALSE))
}
map_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 7, color = "grey50"),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(14, 14, 14, 14))
save_map <- function(p, name, w = 13, h = 9, out_dir = OUT) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  message("  wrote ", file.path(out_dir, paste0(name, ".png")))
}

# Region extent: all slick spots + the four cities, padded ~0.1 degree
reg_wgs <- st_bbox(rbind(st_as_sf(ev, coords = c("lon", "lat"), crs = 4326)["geometry"],
                         st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)["geometry"]))
reg_wgs[c("xmin", "ymin")] <- reg_wgs[c("xmin", "ymin")] - 0.1
reg_wgs[c("xmax", "ymax")] <- reg_wgs[c("xmax", "ymax")] + 0.1
dem_reg <- get_dem(reg_wgs, 9, "region")
bb_reg <- st_bbox(dem_reg)
DEM_LIMITS <- as.numeric(global(dem_reg, range, na.rm = TRUE))


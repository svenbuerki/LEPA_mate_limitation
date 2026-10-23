# step28_model_population_design_map.R — maps of the demes selected for the
# model-population sampling plan (Step 28, tables/step30_D_model_design_demes.tsv):
# one figure per population, 2025 and 2026 side by side, on the shared
# topographic background. Selected demes in the BL colour (label: deme ·
# mothers), fit but not selected in grey, unfit demes as dashed outlines.
#
# The maps show slick spots at full resolution, so they stay local
# (figures/step30_D_model_design_demes_P*.png|pdf are gitignored).
#
# Run from the repository root:  Rscript Scripts/step28_model_population_design_map.R

OUT <- "figures"
source("Scripts/srk_map_background.R")

# Selected demes of the model populations (Step 28 design): one figure per
# population, 2025 and 2026 side by side. Selected = BL colour (label: deme
# and mothers), fit but not selected = grey, unfit = dashed outline.
design_f <- "tables/step30_D_model_design_demes.tsv"
if (file.exists(design_f)) {
  des <- read.delim(design_f)
  for (pop in unique(des$populationID)) {
    dp <- des[des$populationID == pop, ]
    sub_all <- ev_sf[ev_sf$populationID == pop, ]
    bl <- sub_all$BL[1]
    bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(sub_all)), 300))
    side <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"], 600)
    cx <- mean(bb[c("xmin", "xmax")]); cy <- mean(bb[c("ymin", "ymax")])
    bb[c("xmin", "xmax")] <- cx + c(-1, 1) * side / 2
    bb[c("ymin", "ymax")] <- cy + c(-1, 1) * side / 2
    bb_wgs <- st_bbox(st_transform(st_as_sfc(bb), 4326))
    dem <- crop(get_dem(bb_wgs, if (side > 5000) 12 else 13, paste0("P", pop)),
                ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
    panels <- lapply(sort(unique(dp$year)), function(yr) {
      ey <- sub_all[sub_all$event_year == yr, ]
      dy <- dp[dp$year == yr, ]
      dh <- do.call(rbind, lapply(sort(unique(ey$demeID)), function(dm) {
        st_sf(demeID = dm, geometry = st_buffer(st_convex_hull(st_union(ey[ey$demeID == dm, ])), 25))
      }))
      dh <- merge(dh, dy, by = "demeID")
      dh$status <- ifelse(dh$selected == "True" | dh$selected == TRUE, "Selected",
                   ifelse(dh$fit == "True" | dh$fit == TRUE, "Fit, not selected", "Unfit"))
      lab <- cbind(st_drop_geometry(dh), st_coordinates(st_point_on_surface(dh)))
      lab <- lab[lab$status == "Selected", ]
      lab$label <- sprintf("d%d · %d", lab$demeID, lab$mother_target)
      lab$kind <- "pop"
      n_sel <- sum(dh$status == "Selected"); m_sel <- sum(dh$mother_target[dh$status == "Selected"])
      ggplot() + base_layers(dem, bb, show_dem_legend = FALSE) +
        geom_sf(data = dh, aes(fill = status, linetype = status), color = "grey15",
                linewidth = 0.4, alpha = 0.85) +
        scale_fill_manual(values = c("Selected" = BL_COLORS[[bl]], "Fit, not selected" = "grey70",
                                     "Unfit" = NA), name = "Deme", drop = FALSE) +
        scale_linetype_manual(values = c("Selected" = "solid", "Fit, not selected" = "solid",
                                         "Unfit" = "22"), name = "Deme", drop = FALSE) +
        label_layer(lab[, c("X", "Y", "label", "kind")], data.frame(X = numeric(), Y = numeric())) +
        coord_sf(crs = UTM, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
                 expand = FALSE) +
        labs(title = sprintf("%d — %d demes selected, %d mothers", yr, n_sel, m_sel),
             x = NULL, y = NULL) +
        map_theme + theme(plot.title = element_text(face = "bold", size = 13))
    })
    p <- wrap_plots(panels, nrow = 1, guides = "collect") +
      plot_annotation(title = sprintf("P%d (%s) — demes sampled for the model-population test (label: deme · mothers)",
                                      pop, bl),
                      theme = theme(plot.title = element_text(face = "bold", size = 15,
                                                              color = BL_COLORS[[bl]])))
    save_map(p, sprintf("step30_D_model_design_demes_P%d", pop), w = 16, h = 8)
  }
}
message("[step28] model-design deme maps written to figures/")

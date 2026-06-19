# analysis/01_population_map.R ------------------------------------------
#
# Figure 1: Map of collection localities, overlaid on the proposed ranges
# of Cx. pipiens, the hybrid zone, and Cx. quinquefasciatus (Barr 1957).
#
# Now includes Baltimore and SERC as mid-latitude populations, used in
# the hybrid-zone extension (analysis/08).
# ------------------------------------------------------------------------

library(usmap)

# --- Build population points ------------------------------------------
# Union of main and hybrid-zone populations — all 7 get mapped.
pops_all <- union(POPS_MAIN, POPS_HZ)
pop_info_map <- pop_info[pop_info$pop_id %in% pops_all, ]

pop_pts <- sf::st_as_sf(
  pop_info_map,
  coords = c("lon", "lat"),
  crs    = 4326
)
pop_pts_tx <- usmap::usmap_transform(pop_pts)


# --- Barr (1957) range lines ------------------------------------------
lon_grid <- seq(-108, -72, by = 0.5)

north_line <- usmap::usmap_transform(
  data.frame(lon = lon_grid, lat = 39)
)
south_line <- usmap::usmap_transform(
  data.frame(lon = lon_grid, lat = 36)
)

labels_df <- usmap::usmap_transform(
  data.frame(
    lon      = rep(-83, 3),
    lat      = c(40, 37.5, 35),
    label    = c("Culex pipiens", "Hybrid Zone", "Culex quinquefasciatus"),
    fontface = c("italic", "plain", "italic"),
    stringsAsFactors = FALSE
  )
)


# --- Plot --------------------------------------------------------------
p_map <- usmap::plot_usmap(
  exclude = c(.mountain, .pacific),
  fill    = "lightyellow",
  color   = "#888888"
) +
  ggplot2::geom_sf(data = pop_pts_tx,
                   ggplot2::aes(fill = pop_id),
                   shape = 23, size = 8, stroke = 0.8) +
  ggplot2::geom_sf(data = north_line, size = 1) +
  ggplot2::geom_sf(data = south_line, size = 1) +
  ggplot2::geom_sf_text(
    data = labels_df,
    ggplot2::aes(label = label, fontface = fontface),
    size = 6, angle = 12
  ) +
  ggplot2::scale_fill_manual(values = pop_colors, name = "Population",
                             breaks = POP_LEVELS[POP_LEVELS %in% pops_all]) +
  ggplot2::theme_void(base_size = 20) +
  ggplot2::theme(
    legend.position   = c(0.18, 0.7),
    legend.background = ggplot2::element_rect(
      fill = "white", color = "black", linewidth = 0.5
    ),
    legend.margin     = ggplot2::margin(0.3, 0.3, 0.3, 0.3, "cm")
  )

save_fig(p_map, "Fig1_population_map",
         width = 10, height = 8)

#' Function to create plot of DTM and density structure
#'
#' @param sw_rad_2D Dataframe with shortwave radiative values
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @import ggplot2
#' @export
#'
plots_dtm_struct <- function(dtm, grid){

  #########
  # PLOTS #
  #########

  # Plot DTM

  dtm_plot = ggplot2::ggplot(dtm, ggplot2::aes(x = X - min(X), y = Y - min(Y), fill = Z - min(Z))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "magma",
                                  guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) + # Mooie kleurenschaal
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, margin = ggplot2::margin(r = 10))) + # Labels horizontaal maken
    ggplot2::labs(title = "Digital terrain model (DTM), top view",
                  x = "west     ——     east (m)",
                  y = "north\n \n |\n \n south (m)",
                  fill = "Ground\nelevation (m)")
  print(dtm_plot)

  ggplot2::ggsave(paste0('Output/TLS_dtm.png'), plot = dtm_plot, width = 9, height = 3, dpi = 300)


  # Plot structural grid

  # Ensure output directory exists
  if (!dir.exists('Output')) dir.create('Output', recursive = TRUE)

  # Calculate the average of the density over all y-slices for each combination of x and z
  average_density <- grid |>
    dplyr::group_by(X, Z) |>
    dplyr::summarise(mean_density = mean(density, na.rm = TRUE), .groups = 'keep')

  # Voeg een lichtblauwe kleur toe voor mean_density == 0 en verwijder uit de legende
  average_density$color <- ifelse(average_density$mean_density == 0, "lightblue", NA)


  density_plot = ggplot2::ggplot(average_density, ggplot2::aes(x = X, y = Z, fill = mean_density)) +
    ggplot2::geom_tile(ggplot2::aes(fill = ifelse(mean_density == 0, NA, mean_density)), na.rm = TRUE) +
    ggplot2::scale_fill_gradientn(
      colors = c("lightblue", "darkgoldenrod1", "darkorange", "darkolivegreen2", "darkgreen"),
      values = scales::rescale(c(0, 0.1, 0.4, 0.6, 1)),
      na.value = "lightblue",
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")
    ) +
    ggplot2::labs(title = paste("Normalised density, averaged across Y-slices"),
                  x = "X (m)", y = "Z (m)", fill = "Density \n(unitless)")+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw()

  print(density_plot)

  ggplot2::ggsave(paste0('Output/TLS_grid_density.png'), plot = density_plot, width = 9, height = 3, dpi = 300)

}

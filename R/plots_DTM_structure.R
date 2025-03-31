#' Function to create plot of DTM and density structure
#'
#' @param sw_rad_2D Dataframe with shortwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_dtm_struct <- function(dtm, grid, output_path){

  #########
  # PLOTS #
  #########

  # Plot DTM

  dtm_plot = ggplot(dtm, aes(x = X - min(X), y = Y - min(Y), fill = Z - min(Z))) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma",
                                  guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) + # Mooie kleurenschaal
    coord_fixed() +
    theme_bw() +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 10))) + # Labels horizontaal maken
    labs(title = "Digital terrain model (DTM), top view",
                  x = "west     ——     east (m)",
                  y = "north\n \n |\n \n south (m)",
                  fill = "Ground\nelevation (m)")
  print(dtm_plot)

  ggsave(paste0(output_path, '/TLS_dtm.png'), plot = dtm_plot, width = 9, height = 3, dpi = 300)


  # Plot structural grid

  # Ensure output directory exists
  if (!dir.exists('Output')) dir.create('Output', recursive = TRUE)

  # Calculate the average of the density over all y-slices for each combination of x and z
  average_density <- grid |>
    group_by(X, Z) |>
    summarise(mean_density = mean(density, na.rm = TRUE), .groups = 'keep')

  # Voeg een lichtblauwe kleur toe voor mean_density == 0 en verwijder uit de legende
  average_density$color <- ifelse(average_density$mean_density == 0, "lightblue", NA)


  density_plot = ggplot(average_density, aes(x = X, y = Z, fill = mean_density)) +
    geom_tile(aes(fill = ifelse(mean_density == 0, NA, mean_density)), na.rm = TRUE) +
    scale_fill_gradientn(
      colors = c("lightblue", "darkgoldenrod1", "darkorange", "darkolivegreen2", "darkgreen"),
      values = rescale(c(0, 0.1, 0.4, 0.6, 1)),
      na.value = "lightblue",
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")
    ) +
    labs(title = paste("Normalised density, averaged across Y-slices"),
                  x = "X (m)", y = "Z (m)", fill = "Density \n(unitless)")+
    coord_fixed(ratio = 1) +
    theme_bw()

  print(density_plot)

  ggsave(paste0(output_path, '/TLS_grid_density.png'), plot = density_plot, width = 9, height = 3, dpi = 300)

}

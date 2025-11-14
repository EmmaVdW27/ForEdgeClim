#' Function to create several longwave radiation plots
#'
#' @param lw_rad_2D Dataframe with longwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_lw <- function(lw_rad_2D, output_path){

  # 2D RTM plots
  ##############

  # Calculate averages across Y-slices
  final_avg_results_2D <- lw_rad_2D |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE), .groups = 'keep')

  # Caption to be plotted below the plots
  caption = paste0(
    "Diffuse LW radiation = ", round(F_sky_lw,2), ' W/m²'
  )


  # F_d_down plot:
  F_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse LW downward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "LW down (W/m²)",
         caption = caption) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(F_d_down)

  ggsave(paste0(output_path, '/LW_2D_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 300)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse LW upward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "LW up (W/m²)",
         caption = caption) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(F_d_up)

  ggsave(paste0(output_path, '/LW_2D_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 300)




}

#' Function to create several longwave radiation plots
#'
#' @param sw_rad_2D Dataframe with shortwave radiative values
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @import ggplot2
#' @export
#'
plots_lw <- function(lw_rad_2D){

  # 2D RTM plots
  ##############

  # Calculate averages across Y-slices
  final_avg_results_2D <- lw_rad_2D |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE), .groups = 'keep')


  # F_d_down plot:
  F_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse LW downward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_d_down (W/m2)",
         caption = paste0(
           "Diffuse LW radiation = ", F_sky_lw, ' W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_d_down)

  ggsave(paste0('Output/LW_2D_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 300)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse LW upward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_d_up (W/m2)",
         caption = paste0(
           "Diffuse LW radiation = ", F_sky_lw, ' W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_d_up)

  ggsave(paste0('Output/LW_2D_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 300)




}

#' Function to create several shortwave radiation plots
#'
#' @param sw_rad_2D Dataframe with shortwave radiative values
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @import ggplot2
#' @export
#'
plots_sw <- function(sw_rad_2D){

  # 2D RTM plots
  ##############

  # Calculate averages across Y-slices
  final_avg_results_2D <- sw_rad_2D |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE),
              avg_F_b_down = mean(F_b_down, na.rm = TRUE), .groups = 'keep')


  # F_d_down plot:
  F_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse downward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_d_down (W/m2)",
         caption = paste(
           "Direct radiation vertical =", F_sky_dir_v, 'W/m2',
           "| Direct radiation horizontal =", F_sky_dir_h, 'W/m2',
           "| Diffuse radiation =", F_sky_diff_init, 'W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_d_down)

  ggsave(paste0('Output/SW_2D_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 300)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral diffuse upward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_d_up (W/m2)",
         caption = paste(
           "Direct radiation vertical =", F_sky_dir_v, 'W/m2',
           "| Direct radiation horizontal =", F_sky_dir_h, 'W/m2',
           "| Diffuse radiation =", F_sky_diff_init, 'W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_d_up)

  ggsave(paste0('Output/SW_2D_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 300)

  # F_b_down plot:
  F_b_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_b_down)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral direct beam downward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_b_down (W/m2)",
         caption = paste(
           "Direct radiation vertical =", F_sky_dir_v, 'W/m2',
           "| Direct radiation horizontal =", F_sky_dir_h, 'W/m2',
           "| Diffuse radiation =", F_sky_diff_init, 'W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_b_down)

  ggsave(paste0('Output/SW_2D_F_b_down.png'), plot = F_b_down, width = 9, height = 3, dpi = 300)

  # F_d_down + F_b_down plot:
  F_b_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down + avg_F_b_down)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("darkgreen", "seagreen", "yellowgreen", "orange"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Vertical & lateral, direct & diffuse downward radiation, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "F_b_d_down (W/m2)",
         caption = paste(
           "Direct radiation vertical =", F_sky_dir_v, 'W/m2',
           "| Direct radiation horizontal =", F_sky_dir_h, 'W/m2',
           "| Diffuse radiation =", F_sky_diff_init, 'W/m2'
         )) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "seagreen")
    )

  print(F_b_d_down)

  ggsave(paste0('Output/SW_2D_F_b_d_down.png'), plot = F_b_d_down, width = 9, height = 3, dpi = 300)


}

#' Function to create several temperature plots
#'
#' @param micro_grid Grid with micro variables
#' @param T_air_vec Vector with air temperatures
#' @param net_radiaiton Vector with net radiation
#' @param sensible_flux Vector with sensible heat fluxes
#' @param latent_flux Vector with latent heat fluxes
#' @param G Vector with ground heat fluxes
#' @return Several temperature plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @import ggplot2
#' @export
#'
plots_temp_flux <- function(micro_grid, T_air_vec, net_radiation, sensible_flux, latent_flux, G){
  # Define grid with surface temperatures
  temp_surf_grid = micro_grid
  # Convert surface temperature values from Kelvin to degrees Celsius
  temp_surf_grid$temperature = micro_grid$temperature - 273.15
  temp_surf_grid$T_soil = ifelse(micro_grid$z == 1, micro_grid$T_soil - 273.15, NA)
  temp_surf_grid$density = den

  # Define grid with air temperatures
  temp_air_grid = micro_grid
  # Convert air temperature values from Kelvin to degrees Celsius
  temp_air_grid$temperature = T_air_vec - 273.15

  # Define grid with fluxes
  fluxes_grid = micro_grid
  # Convert surface temperature values from Kelvin to degrees Celsius
  fluxes_grid$net_radiation = net_radiation
  fluxes_grid$sensible_flux = sensible_flux
  fluxes_grid$latent_flux = latent_flux
  fluxes_grid$ground_flux = G

  # Calculate the average of the net radiation over all y-slices for each combination of x and z
  average_radiation <- fluxes_grid |>
    dplyr::filter(net_radiation != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_radiation = mean(net_radiation, na.rm = TRUE), .groups = 'keep')

  # Check whether 'Output' directory already exists, if not create it.
  if (!dir.exists('Output')) {
    dir.create('Output', recursive = TRUE)
    message("'Output/' directory was created.")
  }

  # 2D tile plot:
  rad_plot = ggplot2::ggplot(average_radiation, ggplot2::aes(x = x, y = z, fill = mean_radiation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Netto radiation, averaged across Y-slices"),
         subtitle = "+ value ⮕ flux entering surface",
         x = "X (m)", y = "Z (m)", fill = "Rn (W/m2)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(rad_plot)

  ggplot2::ggsave(paste0('Output/flux_net_radiation.png'), plot = rad_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the sensible heat flux over all y-slices for each combination of x and z
  average_H <- fluxes_grid |>
    dplyr::filter(sensible_flux != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_H = mean(sensible_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  H_plot = ggplot2::ggplot(average_H, ggplot2::aes(x = x, y = z, fill = mean_H)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Sensible heat surface flux, averaged across Y-slices"),
         subtitle = "+ value ⮕ flux lost to air",
         x = "X (m)", y = "Z (m)", fill = "H (W/m2)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(H_plot)

  ggplot2::ggsave(paste0('Output/flux_sensible.png'), plot = H_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the latent heat flux over all y-slices for each combination of x and z
  average_LE <- fluxes_grid |>
    dplyr::filter(latent_flux != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_LE = mean(latent_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  LE_plot = ggplot2::ggplot(average_LE, ggplot2::aes(x = x, y = z, fill = mean_LE)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Latent heat surface flux, averaged across Y-slices"),
         subtitle = "+ value ⮕ flux lost to air",
         x = "X (m)", y = "Z (m)", fill = "LE (W/m2)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(LE_plot)

  ggplot2::ggsave(paste0('Output/flux_latent.png'), plot = LE_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the ground heat flux over all y-slices for each combination of x and z
  average_G <- fluxes_grid |>
    dplyr::filter(ground_flux != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_G = mean(ground_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  ground_plot = ggplot2::ggplot(average_G, ggplot2::aes(x = x, y = z, fill = mean_G)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Ground heat flux, averaged across Y-slices"),
         subtitle = "+ value ⮕ flux entering soil",
         x = "X (m)", y = "Z (m)", fill = "G (W/m2)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue"),
      axis.text.y = ggplot2::element_blank()
    )

  print(ground_plot)

  ggplot2::ggsave(paste0('Output/flux_ground.png'), plot = ground_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the surface temperature over all y-slices for each combination of x and z
  average_surface_temperature <- temp_surf_grid |>
    dplyr::filter(den != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_surf_plot = ggplot2::ggplot(average_surface_temperature, ggplot2::aes(x = x, y = z, fill = mean_temperature)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Surface temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(temp_surf_plot)

  ggplot2::ggsave(paste0('Output/temp_surface.png'), plot = temp_surf_plot, width = 9, height = 3, dpi = 300)


  # Calculate the average of the soil temperature over all y-slices for each combination of x and z
  average_soil_temperature <- temp_surf_grid |>
    dplyr::filter(T_soil != 0) |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_temperature = mean(T_soil, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_soil_plot = ggplot2::ggplot(average_soil_temperature, ggplot2::aes(x = x, y = z, fill = mean_temperature)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Soil temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue"),
      axis.text.y = ggplot2::element_blank()
    )

  print(temp_soil_plot)

  ggplot2::ggsave(paste0('Output/temp_soil.png'), plot = temp_soil_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the air temperature over all y-slices for each combination of x and z
  average_air_temperature <- temp_air_grid |>
    dplyr::group_by(x, z) |>
    dplyr::summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_air_plot = ggplot2::ggplot(average_air_temperature, ggplot2::aes(x = x, y = z, fill = mean_temperature)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "lightblue", "orange", "red"),
      values = scales::rescale(c(0.25, 0.5, 0.75, 1)),
      guide = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    ggplot2::labs(title = paste("Air temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         ))+
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(temp_air_plot)

  ggplot2::ggsave(paste0('Output/temp_air.png'), plot = temp_air_plot, width = 9, height = 3, dpi = 300)


  # Average air temperature along the x-axis over the xy-plane at Z = ...
  reqhgt <- temp_air_grid |>
    dplyr::filter(z == 5) |>
    dplyr::group_by(x) |>
    dplyr::summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # 1D graph:
  Temp_height <- ggplot2::ggplot(reqhgt, ggplot2::aes(x = x, y = mean_temperature)) +
    ggplot2::geom_line(color = 'cornflowerblue') +
    ggplot2::geom_point(color = 'cornflowerblue') +
    ggplot2::geom_smooth(color = "blue", method = "loess", formula = y ~ x, se = FALSE) +
    ggplot2::labs(title = "Air temperature at 5m height from forest core to edge, averaged across horizontal plane", x = "Distance (m)", y = " Mean temperature (°C)",
         caption = paste0(
           "Macro T = ", macro_temp - 273.15, ' °C',
           " | Date-time = ", datetime, '\n',
           "Direct radiation vertical = ", F_sky_dir_v, ' W/m2',
           " | Direct radiation horizontal = ", F_sky_dir_h, ' W/m2',
           " | Diffuse radiation = ", F_sky_diff_init, ' W/m2'
         )) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, color = "cornflowerblue")
    )

  print(Temp_height)

  ggplot2::ggsave(paste0('Output/temp_air_reqhgt.png'), plot = Temp_height, width = 10, height = 6, dpi = 300)


}

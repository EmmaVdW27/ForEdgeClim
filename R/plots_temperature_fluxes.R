#' Function to create several temperature plots
#'
#' @param micro_grid Grid with micro variables
#' @param T_air_vec Vector with air temperatures
#' @param output_path The output_path to store the plots
#' @return Several temperature plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_temp <- function(micro_grid, T_air_vec, output_path){

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

  # Calculate the average of the surface temperature over all y-slices for each combination of x and z
  average_surface_temperature <- temp_surf_grid |>
    filter(den != 0) |>
    group_by(x, z) |>
    summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # Caption to be plotted below the plots
  caption = substitute(
    atop(
      "Macrotemperature = "*T*"°C | Date-time = "*DT*" h UTC",
      "Direct radiation vertical = "*V*" W/m"^2*
        " | Direct radiation horizontal = "*H*" W/m"^2*
        " | Diffuse radiation = "*D*" W/m"^2
    ),
    list(
      T = round(macro_temp, 2) - 273.15,
      DT = format(datetime, "%Y-%m-%d %H"),
      V = round(F_sky_dir_v, 2),
      H = round(F_sky_dir_h, 2),
      D = round(F_sky_diff_init, 2)
    )
  )



  # Colors to be used in the plots
  colors = c("blue", "lightblue", "orange", "red")

  # 2D tile plot:
  temp_surf_plot = ggplot(average_surface_temperature, aes(x = x, y = z, fill = mean_temperature)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Surface temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(temp_surf_plot)

  ggsave(paste0(output_path, '/temp_surface.png'), plot = temp_surf_plot, width = 9, height = 3, dpi = 300)


  # Calculate the average of the soil temperature over all y-slices for each combination of x and z
  average_soil_temperature <- temp_surf_grid |>
    filter(T_soil != 0) |>
    group_by(x, z) |>
    summarise(mean_temperature = mean(T_soil, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_soil_plot = ggplot(average_soil_temperature, aes(x = x, y = z, fill = mean_temperature)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Soil temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.text.y = element_blank()
    )

  print(temp_soil_plot)

  ggsave(paste0(output_path, '/temp_soil.png'), plot = temp_soil_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the air temperature over all y-slices for each combination of x and z
  average_air_temperature <- temp_air_grid |>
    group_by(x, z) |>
    summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_air_plot = ggplot(average_air_temperature, aes(x = x, y = z, fill = mean_temperature)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Air temperature, averaged across Y-slices"),
         x = "X (m)", y = "Z (m)", fill = "Temperature (°C)",
         caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(temp_air_plot)

  ggsave(paste0(output_path, '/temp_air.png'), plot = temp_air_plot, width = 9, height = 3, dpi = 300)


  # Modelled air temperature at DTS position
  reqhgt <- temp_air_grid |>
    filter(z == 1, y == 15, x <= length_transect)

  # DTS observations
  DTS <- read.csv("Data/DTS_filtered_distance_temp.csv")

  # Add dataset lable
  reqhgt$model <- "Modelled (every 1m)"
  DTS$model <- "DTS observations (every 0.25m)"

  # 1D graph:
  Temp_height <- ggplot() +
    geom_line(data = reqhgt, aes(x = x, y = temperature, color = "Modelled (every 1m)")) +
    geom_line(data = DTS, aes(x = distance, y = temp, color = "DTS observations (every 0.25m)")) +
    labs(
      title = "Air temperature at 1m height from forest core to edge",
      x = "Distance (m)",
      y = "Temperature (°C)",
      caption = caption) +
    scale_color_manual(
      values = c("Modelled (every 1m)" = "cornflowerblue", "DTS observations (every 0.25m)" = "red")
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(Temp_height)

  ggsave(paste0(output_path, '/temp_air_reqhgt.png'), plot = Temp_height, width = 10, height = 6, dpi = 300)


}


#' Function to create several flux plots
#'
#' @param micro_grid Grid with micro variables
#' @param net_radiaiton Vector with net radiation
#' @param sensible_flux Vector with sensible heat fluxes
#' @param latent_flux Vector with latent heat fluxes
#' @param G Vector with ground heat fluxes
#' @param output_path The output_path to store the plots
#' @return Several flux plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_flux <- function(micro_grid, net_radiation, sensible_flux, latent_flux, G, output_path){

  # Define grid with fluxes
  fluxes_grid = micro_grid
  # Convert surface temperature values from Kelvin to degrees Celsius
  fluxes_grid$net_radiation = net_radiation
  fluxes_grid$sensible_flux = sensible_flux
  fluxes_grid$latent_flux = latent_flux
  fluxes_grid$ground_flux = G

  # Calculate the average of the net radiation over all y-slices for each combination of x and z
  average_radiation <- fluxes_grid |>
    filter(net_radiation != 0) |>
    group_by(x, z) |>
    summarise(mean_radiation = mean(net_radiation, na.rm = TRUE), .groups = 'keep')

  # Caption to be plotted below the plots
  caption = substitute(
    atop(
      "Macrotemperature = "*T*"°C | Date-time = "*DT*" h UTC",
      "Direct radiation vertical = "*V*" W/m"^2*
        " | Direct radiation horizontal = "*H*" W/m"^2*
        " | Diffuse radiation = "*D*" W/m"^2
    ),
    list(
      T = round(macro_temp, 2) - 273.15,
      DT = format(datetime, "%Y-%m-%d %H"),
      V = round(F_sky_dir_v, 2),
      H = round(F_sky_dir_h, 2),
      D = round(F_sky_diff_init, 2)
    )
  )


  # Colors to be used in the plots
  colors = c("blue", "lightblue", "orange", "red")

  # 2D tile plot:
  rad_plot = ggplot(average_radiation, aes(x = x, y = z, fill = mean_radiation)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Netto radiation, averaged across Y-slices"),
                  subtitle = "+ value ⮕ flux entering surface",
                  x = "X (m)", y = "Z (m)", fill = "Rn (W/m2)",
                  caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(rad_plot)

  ggsave(paste0(output_path, '/flux_net_radiation.png'), plot = rad_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the sensible heat flux over all y-slices for each combination of x and z
  average_H <- fluxes_grid |>
    filter(sensible_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_H = mean(sensible_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  H_plot = ggplot(average_H, aes(x = x, y = z, fill = mean_H)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Sensible heat surface flux, averaged across Y-slices"),
                  subtitle = "+ value ⮕ flux lost to air",
                  x = "X (m)", y = "Z (m)", fill = "H (W/m2)",
                  caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(H_plot)

  ggsave(paste0(output_path, '/flux_sensible.png'), plot = H_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the latent heat flux over all y-slices for each combination of x and z
  average_LE <- fluxes_grid |>
    filter(latent_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_LE = mean(latent_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  LE_plot = ggplot(average_LE, aes(x = x, y = z, fill = mean_LE)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Latent heat surface flux, averaged across Y-slices"),
                  subtitle = "+ value ⮕ flux lost to air",
                  x = "X (m)", y = "Z (m)", fill = "LE (W/m2)",
                  caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  print(LE_plot)

  ggsave(paste0(output_path, '/flux_latent.png'), plot = LE_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the ground heat flux over all y-slices for each combination of x and z
  average_G <- fluxes_grid |>
    filter(ground_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_G = mean(ground_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  ground_plot = ggplot(average_G, aes(x = x, y = z, fill = mean_G)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Ground heat flux, averaged across Y-slices"),
                  subtitle = "+ value ⮕ flux entering soil",
                  x = "X (m)", y = "Z (m)", fill = "G (W/m2)",
                  caption = caption)+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.text.y = element_blank()
    )

  print(ground_plot)

  ggsave(paste0(output_path, '/flux_ground.png'), plot = ground_plot, width = 9, height = 3, dpi = 300)


}

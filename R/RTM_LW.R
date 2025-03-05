#' Two-stream RTM for longwave radiation and for a single column existing of layers/voxels with density values (%).
#' This model is inspired by the longwave RTM of the Ecosystem Demography model.
#' @param soil_reflect Soil albedo (%)
#' @param density Density of each layer (%)
#' @param temperature Temperature of each layer (K)
#' @param F_sky_diff Diffuse radiation (W/m2)
#' @param Kd Diffuse extinction coefficient, per unit density
#' @param omega Scattering coefficient
#' @param beta Fraction of scattering in the backward direction for diffuse radiation (%)
#' @return Dataframe containing the combined results of the 2D shortwave RTM
#' @export

lw_two_stream <- function(soil_reflect, density, F_sky_lw, temperature, T_soil,
                          T_macro,
                          Kd,
                          omega, beta) {
  omega_g = soil_reflect
  F_sky_diff = F_sky_lw

  # Calculate black body emission at the temperature of the layer
  # Black body emission of atmosphere is defined to be 0 because F_sky_lw has already been defined as above canopy input
  # Black body emission is per surface unit, so a density correction must be applied.
  # Like this, zero-density values with no structure or surface don't produce BB emission (as eg the atmosphere).
  black = c(density*e_forest*sigma_SB * temperature^4,0)

  # Number of (cohort) layers
  n <- length(density)
  # Add atmosphere layer
  density <- c(density, 0)
  # Inverse optical depths
  # mu atmosphere (density = 0, no forest structure) is defined to be 1
  mu <- ifelse(density == 0, 1, 1 / (Kd * density))  # Inverse optical depth diffuse radiation

  # Initialisation matrices
  S <- matrix(0, nrow = 2 * n + 2, ncol = 2 * n + 2)
  y <- numeric(2 * n + 2)

  # black body emission of the soil
  # emissivity of soil = absorption of soil for lw radiation => emissivity = 1 - omega_g
  black_g = (1-omega_g)*sigma_SB*T_soil^4

  # Auxiliary variables
  # Gamma
  gamma_plus <- rep(0.5 * (1 + sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_plus[n+1] = 1
  gamma_minus <- rep(0.5 * (1 - sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_minus[n+1] = 0
  # Lambda
  lambda <- sqrt((1 - (1 - 2 * beta) * omega) * (1 - omega) / (mu)^2)
  lambda[n+1] = 0
  # no Kappa and Delta as in shortwave RTM

  # Boundary conditions
  y[2*n+2] <- F_sky_diff - black[n+1] # Incoming diffuse radiation at top canopy
  y[1] <- (1 - omega_g) * black_g - (1 - omega_g)*black[1]

  S[2 * n + 2, 2 * n + 1] <- gamma_plus[n+1]
  S[2 * n + 2, 2 * n + 2] <- gamma_minus[n+1]
  S[1, 1] <- (gamma_minus[1] - omega_g*gamma_plus[1])*exp(-lambda[1]*density[1])
  S[1, 2] <- (gamma_plus[1] - omega_g*gamma_minus[1])*exp(lambda[1]*density[1])

  # Filling of matrix S and vector y
  for (i in 1:n) {
    S[2 * i, 2 * i - 1] <- gamma_plus[i]
    S[2 * i, 2 * i] <- gamma_minus[i]
    S[2 * i, 2 * i + 1] <- -gamma_plus[i+1] * exp(-lambda[i+1] * (density[i+1]))
    S[2 * i, 2 * i + 2] <- -gamma_minus[i+1] * exp(lambda[i+1] * (density[i+1]))

    S[2 * i + 1, 2 * i - 1] <- gamma_minus[i]
    S[2 * i + 1, 2 * i] <- gamma_plus[i]
    S[2 * i + 1, 2 * i + 1] <- -gamma_minus[i+1] * exp(-lambda[i+1] * (density[i+1]))
    S[2 * i + 1, 2 * i + 2] <- -gamma_plus[i+1] * exp(lambda[i+1] * (density[i+1]))

    y[2 * i] <- black[i+1] - black[i]
    y[2 * i + 1] <- black[i+1] - black[i]
  }

  # Solving matrix equation
  x <- solve(S, y)

  # Calculate diffuse radiative fluxes
  # Radiative fluxes refer to the rate of energy transfer by electromagnetic radiation through a surface or area,
  # typically measured in W/m², encompassing all or specific wavelengths (e.g., shortwave or longwave radiation).
  # These values are mainly of interest when modelling e.g. energy balance, heat fluxes...
  # Every down[k] and up[k] respectively refer to the downward and upward flux below layer k.
  down <- numeric(n + 1)
  up <- numeric(n + 1)
  for (k in seq_len(n + 1)) {
    k2 <- 2 * k
    k2m1 <- k2 - 1
    down[k] <- x[k2m1] * gamma_plus[k] * exp(-lambda[k]*density[k]) +
      x[k2] * gamma_minus[k] * exp(lambda[k]*density[k]) +
      black[k]
    up[k] <- x[k2m1] * gamma_minus[k] * exp(-lambda[k]*density[k]) +
      x[k2] * gamma_plus[k] * exp(lambda[k]*density[k]) +
      black[k]
  }

  # Calculate light levels and netto absorbed or emitted radiation
  # Light levels refer to the intensity of photosynthetically active radiation (PAR, 400-700 nm) available for plants.
  # These values are mainly of interest when modelling e.g. photosyntheses, plant growth, schadow structure...
  light_diff_level <- numeric(n)
  net_lw <- numeric(n)
  for (k in seq_len(n)) {
    kp1 <- k + 1
    light_diff_level[k] <- 0.5 * (down[k] + down[kp1])
    net_lw[k] <- (down[kp1] - down[k]) + (up[k] - up[kp1])
  }
  # There are no light levels for atmosphere
  light_diff_level[n+1] = NA


  # Voeg de geabsorbeerde straling toe aan de output
  result <- data.frame(
    layer = c(paste("layer", 1:n), "atmosphere"),
    density = density,
    temperature = c(temperature, NA),
    F_d_down = round(down, 0),
    F_d_up = round(up, 0),
    BB_emission = black,
    net_lw = c(round(net_lw, 0), 0)  # No netto radiation in the atmosphere (no structure there)
  )

  return(result)

}

#' 2D Longwave Two-Stream Radiative Transfer Model (RTM)
#'
#' This function performs a two-directional (vertical & horizontal) longwave
#' two-stream RTM to estimate radiation fluxes in a forest voxel grid.
#'
#' @param voxel_grid Dataframe containing voxel density data (columns: X, Y, Z, density).
#' @param micro_grid Dataframe containing microclimate temperature data (columns: x, y, z, temperature, T_soil).
#' @param lw_two_stream Function that calculates longwave two-stream radiation for a single column.
#' @param F_sky_lw Longwave sky flux.
#' @param omega_g_lw_v Ground reflectance for vertical LW RTM.
#' @param omega_g_lw_h Ground reflectance for horizontal LW RTM.
#' @param macro_temp Macroclimate temperature outside the forest.
#' @param Kd_lw_v Diffuse extinction coefficient for vertical LW RTM.
#' @param Kd_lw_h Diffuse extinction coefficient for horizontal LW RTM.
#' @param omega_lw Single-scattering albedo for longwave RTM.
#' @param beta_lw Phase function parameter for longwave RTM.
#' @return A dataframe with computed longwave fluxes and net longwave radiation.
#' A lot of these parameters are global parameters and actually do not need to be stated as arguments.
#' @export
longwave_two_stream_RTM <- function(voxel_grid, micro_grid, lw_two_stream,
                             F_sky_lw, omega_g_lw_v, omega_g_lw_h, macro_temp,
                             Kd_lw_v, Kd_lw_h, omega_lw, beta_lw) {

  ##########################
  # Vertical shortwave RTM #
  ##########################

  # Unique XY-combinations
  xy_combinations <- unique(voxel_grid[, c("X", "Y")])

  # Save results
  results_list_lw_v <- list()

  # Iteration over unique XY-combinations
  for (row in seq_len(nrow(xy_combinations))) {
    xy <- xy_combinations[row, ]
    x_pos <- xy$X
    y_pos <- xy$Y

    # Density values and temperature values of current column
    current_density <- voxel_grid$density[voxel_grid$X == x_pos & voxel_grid$Y == y_pos]
    current_temperature <- micro_grid$temperature[micro_grid$x == x_pos & micro_grid$y == y_pos]
    soil_temp <- micro_grid$T_soil[micro_grid$x == x_pos & micro_grid$y == y_pos & micro_grid$z == 1]

    # Number of layers
    n = length(current_density)

    fluxes = lw_two_stream(soil_reflect = omega_g_lw_v, density = current_density, F_sky_lw = F_sky_lw,
                           temperature = current_temperature, T_soil = soil_temp,
                           T_macro = macro_temp,
                           Kd = Kd_lw_v,
                           omega = omega_lw, beta = beta_lw)

    # Save results without atmosphere layer
    result <- data.frame(
      X = rep(x_pos, n),  # Fill X for every layer
      Y = rep(y_pos, n),  # Fill Y for every layer
      Z = 1:n,
      density = current_density,
      temperature = current_temperature,
      F_d_down = fluxes$F_d_down[-(n+1)],
      F_d_up = fluxes$F_d_up[-(n+1)],
      net_lw = fluxes$net_lw[-(n+1)]
    )

    results_list_lw_v[[row]] <- result
  }

  # Combine all results
  final_results_lw_v <- do.call(rbind, results_list_lw_v)
  # Sort by X and Y
  final_results_lw_vertical <- final_results_lw_v[order(final_results_lw_v$X, final_results_lw_v$Y), ]

  ############################
  # Horizontal shortwave RTM #
  ############################

  # Unique YZ-combinations
  yz_combinations <- unique(voxel_grid[, c("Y", "Z")])

  # Save results
  results_list_lw_h <- list()

  # Iteration over unique YZ-combinations
  for (row in seq_len(nrow(yz_combinations))) {
    yz <- yz_combinations[row, ]
    y_pos <- yz$Y
    z_pos <- yz$Z

    # Density values and temperature values of current YZ column
    current_density <- voxel_grid$density[voxel_grid$Y == y_pos & voxel_grid$Z == z_pos]
    current_temperature <- micro_grid$temperature[micro_grid$y == y_pos & micro_grid$z == z_pos]

    # Number of layers
    n = length(current_density)

    fluxes = lw_two_stream(soil_reflect = omega_g_lw_h, density = current_density, F_sky_lw = F_sky_lw/4,
                           temperature = current_temperature, T_soil = 0, # no BB emission from the forest interior
                           T_macro = macro_temp,
                           Kd = Kd_lw_h,
                           omega = omega_lw, beta = beta_lw)

    # Save results without atmosphere layer
    result <- data.frame(
      Y = rep(y_pos, n),  # Fill Y for every layer
      Z = rep(z_pos, n),  # Fill Z for every layer
      X = 1:n,
      density = current_density,
      temperature = current_temperature,
      F_d_down = fluxes$F_d_down[-(n+1)],
      F_d_up = fluxes$F_d_up[-(n+1)],
      net_lw = fluxes$net_lw[-(n+1)]
    )

    results_list_lw_h[[row]] <- result

  }

  # Combine all results
  final_results_lw_h <- do.call(rbind, results_list_lw_h)
  # Sort by Y en Z
  final_results_lw_horizontal <- final_results_lw_h[order(final_results_lw_h$Y, final_results_lw_h$Z), ]

  #########################################
  # Combining both directions into 2D RTM #
  #########################################

  # Make sure both dataframes have same order of rows
  final_results_lw_horizontal <- final_results_lw_horizontal[order(final_results_lw_horizontal$X, final_results_lw_horizontal$Y, final_results_lw_horizontal$Z), ]
  final_results_lw_vertical <- final_results_lw_vertical[order(final_results_lw_vertical$X, final_results_lw_vertical$Y, final_results_lw_vertical$Z), ]

  # Calculate sum of corresponding columns and create new dataframe
  final_results_lw_2D <- data.frame(
    X = final_results_lw_horizontal$X,
    Y = final_results_lw_horizontal$Y,
    Z = final_results_lw_horizontal$Z,
    density = final_results_lw_horizontal$density,
    temperature = final_results_lw_horizontal$temperature,
    F_d_down = final_results_lw_horizontal$F_d_down + final_results_lw_vertical$F_d_down,
    F_d_up = final_results_lw_horizontal$F_d_up + final_results_lw_vertical$F_d_up,
    net_lw = final_results_lw_horizontal$net_lw + final_results_lw_vertical$net_lw
  )

  return(final_results_lw_2D)

}

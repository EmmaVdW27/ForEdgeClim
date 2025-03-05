#' Two-stream RTM for shortwave radiation and for a single column existing of layers/voxels with density values (%).
#' This model is inspired by the shortwave RTM of the Ecosystem Demography model.
#' @param soil_reflect Soil albedo (%)
#' @param density Density of each layer (%)
#' @param F_sky_diff Diffuse radiation (W/m2)
#' @param F_sky_dir Direct beam radiation (W/m2)
#' @param Kd Diffuse extinction coefficient, per unit density
#' @param Kb Direct beam extinction coefficient, per unit density
#' @param omega Scattering coefficient
#' @param beta Fraction of scattering in the backward direction for diffuse radiation (%)
#' @param beta0 Fraction of scattering in the backward direction for direct beam radiation (%)
#' @return Dataframe containing the combined results of the 2D shortwave RTM
#' @export

sw_two_stream <- function(soil_reflect, density,
                          F_sky_diff, F_sky_dir,
                          Kd, Kb,
                          omega, beta, beta0) {

  omega_g = soil_reflect
  # Number of (cohort) layers
  n <- length(density)
  # Add atmosphere layer
  density <- c(density, 0)
  # Inverse optical depths
  # mu atmosphere (density = 0, no forest structure) is defined to be 1
  mu0 <- ifelse(density == 0, 1, 1 / (Kb * density)) # Inverse optical depth direct beam radiation
  mu <- ifelse(density == 0, 1, 1 / (Kd * density))  # Inverse optical depth diffuse radiation

  # Initialisation matrices
  S <- matrix(0, nrow = 2 * n + 2, ncol = 2 * n + 2)
  y <- numeric(2 * n + 2)

  # Calculate direct radiation profile (exp decay)
  direct_radiation = numeric(n+1)
  direct_radiation[n+1] = F_sky_dir
  for (i in seq(n, 1)) {
    direct_radiation[i] <- direct_radiation[i+1] * exp(-density[i]/mu0[i])
  }

  # Auxiliary variables
  # Gamma
  gamma_plus <- rep(0.5 * (1 + sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_plus[n+1] = 1
  gamma_minus <- rep(0.5 * (1 - sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_minus[n+1] = 0
  # Lambda
  lambda <- sqrt((1 - (1 - 2 * beta) * omega) * (1 - omega) / (mu)^2)
  lambda[n+1] = 0
  # Kappa
  kappa_plus <- -((1 - (1 - 2 * beta) * omega)/ mu[1:n] + (1 - 2 * beta0) / mu0[1:n]) * omega * direct_radiation[2:(n+1)] / mu0[1:n]
  kappa_plus[n+1] =  3*F_sky_dir
  kappa_minus = - ( (1-2*beta0)*(1-omega)/mu[1:n] + 1/mu0[1:n]) * omega*direct_radiation[2:(n+1)]/mu0[1:n]
  kappa_minus[n+1] = -F_sky_dir
  # Delta
  delta_plus <- (kappa_plus + kappa_minus) * (mu0)^2 / (2 * (1 - lambda^2 * (mu0)^2))
  delta_plus[n+1] = 0.5*(kappa_plus[n+1]+kappa_minus[n+1])
  delta_minus <- (kappa_plus - kappa_minus) * (mu0)^2 / (2 * (1 - lambda^2 * (mu0)^2))
  delta_minus[n+1] = 0.5*(kappa_plus[n+1]-kappa_minus[n+1])

  # Boundary conditions
  y[2*n+2] <- F_sky_diff - delta_plus[n+1] # Incoming diffuse radiation at top canopy
  y[1] <- omega_g*direct_radiation[1] - (delta_minus[1] - omega_g*delta_plus[1])*exp(-density[1]/mu0[1])

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

    y[2 * i] <- delta_plus[i+1] * exp(-density[i + 1] / mu0[i + 1]) - delta_plus[i]
    y[2 * i + 1] <- delta_minus[i+1] * exp(-density[i + 1] / mu0[i + 1]) - delta_minus[i]
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
      delta_plus[k] * exp(-density[k]/mu0[k])
    up[k] <- x[k2m1] * gamma_minus[k] * exp(-lambda[k]*density[k]) +
      x[k2] * gamma_plus[k] * exp(lambda[k]*density[k]) +
      delta_minus[k] * exp(-density[k]/mu0[k])
  }

  # Calculate light levels and netto absorbed or emitted radiation
  # Light levels refer to the intensity of photosynthetically active radiation (PAR, 400-700 nm) available for plants.
  # These values are mainly of interest when modelling e.g. photosyntheses, plant growth, schadow structure...
  light_level <- numeric(n)
  light_beam_level <- numeric(n)
  light_diff_level <- numeric(n)
  net_sw <- numeric(n)
  for (k in seq_len(n)) {
    kp1 <- k + 1
    light_level[k] <- 0.5 * (down[k] + down[kp1] + direct_radiation[k] + direct_radiation[kp1])
    light_beam_level[k] <- 0.5 * (direct_radiation[k] + direct_radiation[kp1])
    light_diff_level[k] <- 0.5 * (down[k] + down[kp1])
    net_sw[k] <- (down[kp1] - down[k]) + (up[k] - up[kp1]) + (direct_radiation[kp1] - direct_radiation[k] )
  }
  # There are no light levels for atmosphere
  light_level[n+1] = NA
  light_beam_level[n+1] = NA
  light_diff_level[n+1] = NA

  # Voeg de geabsorbeerde straling toe aan de output
  result <- data.frame(
    layer = c(paste("layer", 1:n), "atmosphere"),
    density = density,
    F_d_down = round(down, 0),
    F_d_up = round(up, 0),
    F_b_down = round(direct_radiation, 0),
    net_sw = c(round(net_sw, 0), 0),  # No netto radiation in atmosphere
    albedo = c(rep(' ', n), round( up[length(up)]/(F_sky_diff + F_sky_dir), 3 ) )
  )

  return(result)

}


#' 2D (vertical & horizontal) shortwave two-stream RTM
#'
#' @param datetime Datetime object representing the current simulation time
#' @param lat Latitude of the location
#' @param lon Longitude of the location
#' @param voxel_grid Dataframe containing the voxel grid with density values
#' @param F_sky_diff_init Initial diffuse shortwave radiation (W/m2)
#' @param F_sky_dir_init Initial direct shortwave radiation (W/m2)
#' @param omega_g_v Ground reflectance for vertical RTM
#' @param omega_g_h Ground reflectance for horizontal RTM
#' @param Kd_v Diffuse extinction coefficient for vertical RTM
#' @param Kb_v Direct extinction coefficient for vertical RTM
#' @param Kd_h Diffuse extinction coefficient for horizontal RTM
#' @param Kb_h Direct extinction coefficient for horizontal RTM
#' @param omega Scattering coefficient
#' @param betad Fraction of scattered diffuse radiation in the backward direction
#' @param beta0 Fraction of scattered direct radiation in the backward direction
#' @return Dataframe containing the combined results of the 2D shortwave RTM
#' A lot of these parameters are global parameters and actually do not need to be stated as arguments.
#' @importFrom suncalc getSunlightPosition
#' @export

shortwave_two_stream_RTM <- function(datetime, lat, lon, voxel_grid,
                                     F_sky_diff_init, F_sky_dir_init,
                                     omega_g_v, omega_g_h,
                                     Kd_v, Kb_v, Kd_h, Kb_h,
                                     omega, betad, beta0) {

  # Function to determine solar angles
  calculate_sun_angle <- function(datetime, lat, lon) {
    sun_angle <- suncalc::getSunlightPosition(date = datetime, lat = lat, lon = lon)
    sun_altitude = sun_angle$altitude
    sun_azimuth <- sun_angle$azimuth
    # Is sun in eastern half?
    # This is 1 if azimuth is between 0 and -pi radians (direction along the horizon, measured from south to west).
    in_eastern_half <- ifelse(sun_azimuth < 0, 1, 0)
    # Is sun above the horizon?
    # This is 1 if altitude is between 0 and pi/2 radians
    above_horizon <- ifelse(sun_altitude >= 0, 1, 0)
    return(list(in_eastern_half = in_eastern_half, above_horizon = above_horizon,
                sun_altitude = sun_altitude, sun_azimuth = sun_azimuth))
  }

  sun_angles <- calculate_sun_angle(datetime, lat, lon)

  ##########################
  # Vertical shortwave RTM #
  ##########################

  # Radiative values
  F_sky_diff_v = F_sky_diff_init
  F_sky_dir_v <<- F_sky_dir_init * sun_angles$above_horizon

  # Unique XY-combinations
  xy_combinations <- unique(voxel_grid[, c("X", "Y")])

  # Save results
  results_list_v <- list()

  # Iteration over unique XY-combinations
  for (row in seq_len(nrow(xy_combinations))) {
    xy <- xy_combinations[row, ]
    x_pos <- xy$X
    y_pos <- xy$Y

    # Density values of current column
    current_density <- voxel_grid$density[voxel_grid$X == x_pos & voxel_grid$Y == y_pos]
    # Number of layers
    n = length(current_density)

    fluxes = sw_two_stream(soil_reflect = omega_g_v, density = current_density,
                           F_sky_diff = F_sky_diff_v, F_sky_dir = F_sky_dir_v,
                           Kd = Kd_v, Kb = Kb_v,
                           omega = omega, beta = betad, beta0 = beta0)

    # Save results without atmosphere layer
    result <- data.frame(
      X = rep(x_pos, n),  # Fill X for every layer
      Y = rep(y_pos, n),  # Fill Y for every layer
      Z = 1:n,
      density = current_density,
      F_d_down = fluxes$F_d_down[-(n+1)],
      F_d_up = fluxes$F_d_up[-(n+1)],
      F_b_down = fluxes$F_b_down[-(n+1)],
      net_sw = fluxes$net_sw[-(n+1)]
    )

    results_list_v[[row]] <- result
  }

  # Combine all results
  final_results_v <- do.call(rbind, results_list_v)
  # Sort by X and Y
  final_results_vertical <- final_results_v[order(final_results_v$X, final_results_v$Y), ]

  ############################
  # Horizontal shortwave RTM #
  ############################

  # Radiative values
  F_sky_dir_h <<- round( F_sky_dir_init / sin(sun_angles$sun_altitude) *
                         cos(sun_angles$sun_altitude) * abs(cos(sun_angles$sun_azimuth)) *
                         sun_angles$above_horizon * sun_angles$in_eastern_half, 0)     # Direct radiation from eastern edge
  F_sky_diff_h <- F_sky_diff_init/4                                 # Diffuse radiation from eastern edge

  # Unique YZ-combinations
  yz_combinations <- unique(voxel_grid[, c("Y", "Z")])

  # Save results
  results_list_h <- list()

  # Iteration over unique YZ-combinations
  for (row in seq_len(nrow(yz_combinations))) {
    yz <- yz_combinations[row, ]
    y_pos <- yz$Y
    z_pos <- yz$Z

    # Density values of current YZ column
    current_density <- voxel_grid$density[voxel_grid$Y == y_pos & voxel_grid$Z == z_pos]

    # Number of layers
    n = length(current_density)

    fluxes = sw_two_stream(soil_reflect = omega_g_h, density = current_density,
                           F_sky_diff = F_sky_diff_h, F_sky_dir = F_sky_dir_h,
                           Kd = Kd_h, Kb = Kb_h,
                           omega = omega, beta = betad, beta0 = beta0)

    # Save results without atmosphere layer
    result <- data.frame(
      Y = rep(y_pos, n),  # Fill Y for every layer
      Z = rep(z_pos, n),  # Fill Z for every layer
      X = 1:n,
      density = current_density,
      F_d_down = fluxes$F_d_down[-(n+1)],
      F_d_up = fluxes$F_d_up[-(n+1)],
      F_b_down = fluxes$F_b_down[-(n+1)],
      net_sw = fluxes$net_sw[-(n+1)]
    )

    results_list_h[[row]] <- result

  }

  # Combine all results
  final_results_h <- do.call(rbind, results_list_h)
  # Sort by Y en Z
  final_results_horizontal <- final_results_h[order(final_results_h$Y, final_results_h$Z), ]

  #########################################
  # Combining both directions into 2D RTM #
  #########################################

  # Make sure both dataframes have same order of rows
  final_results_horizontal <- final_results_horizontal[order(final_results_horizontal$X, final_results_horizontal$Y, final_results_horizontal$Z), ]
  final_results_vertical <- final_results_vertical[order(final_results_vertical$X, final_results_vertical$Y, final_results_vertical$Z), ]

  # Calculate sum of corresponding columns and create new dataframe
  final_results_2D <- data.frame(
    X = final_results_horizontal$X,
    Y = final_results_horizontal$Y,
    Z = final_results_horizontal$Z,
    density = final_results_horizontal$density,
    F_d_down = final_results_horizontal$F_d_down + final_results_vertical$F_d_down,
    F_d_up = final_results_horizontal$F_d_up + final_results_vertical$F_d_up,
    F_b_down = final_results_horizontal$F_b_down + final_results_vertical$F_b_down,
    net_sw = final_results_horizontal$net_sw + final_results_vertical$net_sw
  )

  return(final_results_2D)

}


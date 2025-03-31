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

  # Calculate direct radiation profile (exp decay from top to bottom)
  direct_radiation <- F_sky_dir * rev(cumprod(rev(exp(-density / mu0))))

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

  # **Vectorized matrix filling**
  indices <- 1:n
  idx_even <- 2 * indices
  idx_odd <- idx_even - 1
  idx_next <- indices + 1

  # Calculate exponential terms
  exp_lambda_density_neg <- exp(-lambda[idx_next] * density[idx_next])
  exp_lambda_density_pos <- exp(lambda[idx_next] * density[idx_next])

  exp_density_mu0_neg <- exp(-density[idx_next] / mu0[idx_next])

  # Populate the matrix S
  S[cbind(idx_even, idx_odd)] <- gamma_plus[indices]
  S[cbind(idx_even, idx_even)] <- gamma_minus[indices]
  S[cbind(idx_even, idx_even + 1)] <- -gamma_plus[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even, idx_even + 2)] <- -gamma_minus[idx_next] * exp_lambda_density_pos

  S[cbind(idx_even + 1, idx_odd)] <- gamma_minus[indices]
  S[cbind(idx_even + 1, idx_even)] <- gamma_plus[indices]
  S[cbind(idx_even + 1, idx_even + 1)] <- -gamma_minus[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even + 1, idx_even + 2)] <- -gamma_plus[idx_next] * exp_lambda_density_pos

  # Populate the vector y
  y[idx_even] <- delta_plus[idx_next] * exp_density_mu0_neg - delta_plus[indices]
  y[idx_even + 1] <- delta_minus[idx_next] * exp_density_mu0_neg - delta_minus[indices]

  # Solving matrix equation using arma::solve from C++
  x <- solve_rtm(S, y)

  # Calculate diffuse radiative fluxes
  # Radiative fluxes refer to the rate of energy transfer by electromagnetic radiation through a surface or area,
  # typically measured in W/m², encompassing all or specific wavelengths (e.g., shortwave or longwave radiation).
  # These values are mainly of interest when modelling e.g. energy balance, heat fluxes...
  # Every down[k] and up[k] respectively refer to the downward and upward flux below layer k.
  indices <- seq_len(n + 1)
  k2 <- 2 * indices
  k2m1 <- k2 - 1

  down <- x[k2m1] * gamma_plus[indices] * exp(-lambda[indices]*density[indices]) +
      x[k2] * gamma_minus[indices] * exp(lambda[indices]*density[indices]) +
      delta_plus[indices] * exp(-density[indices]/mu0[indices])

  up <- x[k2m1] * gamma_minus[indices] * exp(-lambda[indices]*density[indices]) +
      x[k2] * gamma_plus[indices] * exp(lambda[indices]*density[indices]) +
      delta_minus[indices] * exp(-density[indices]/mu0[indices])

  # Calculate light levels and netto absorbed or emitted radiation
  # Light levels refer to the intensity of photosynthetically active radiation (PAR, 400-700 nm) available for plants.
  # These values are mainly of interest when modelling e.g. photosyntheses, plant growth, schadow structure...
  indices2 <- seq_len(n)
  kp1 <- indices2 + 1

  light_level <- 0.5 * (down[indices2] + down[kp1] + direct_radiation[indices2] + direct_radiation[kp1])
  light_beam_level <- 0.5 * (direct_radiation[indices2] + direct_radiation[kp1])
  light_diff_level <- 0.5 * (down[indices2] + down[kp1])
  net_sw <- (down[kp1] - down[indices2]) + (up[indices2] - up[kp1]) + (direct_radiation[kp1] - direct_radiation[indices2] )

  # There are no light levels for atmosphere
  light_level = c(light_level, NA)
  light_beam_level = c(light_beam_level, NA)
  light_diff_level = c(light_diff_level, NA)

  result <- data.frame(
    layer = c(paste("layer", 1:n), "atmosphere"),
    density = density,
    F_d_down = round(down, 0),
    F_d_up = round(up, 0),
    F_b_down = round(direct_radiation, 0),
    net_sw = c(round(net_sw, 0), 0),  # No net radiation in atmosphere (no structure here)
    albedo = c(rep(' ', n), round( up[length(up)]/(F_sky_diff + F_sky_dir), 3 ) )
  )

  return(result)

}


#' 2D (vertical & horizontal) shortwave two-stream RTM
#'
#' This function performs a two-directional (vertical & horizontal) shortwave
#' two-stream RTM to estimate radiation fluxes in a forest voxel grid.
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
#' @importFrom data.table as.data.table setorder
#' @export


shortwave_two_stream_RTM <- function(datetime, lat, lon, voxel_grid,
                                     F_sky_diff_init, F_sky_dir_init,
                                     omega_g_v, omega_g_h,
                                     Kd_v, Kb_v, Kd_h, Kb_h,
                                     omega, betad, beta0) {

  voxel_grid <- as.data.table(voxel_grid)

  # Calculate solar angles
  calculate_sun_angle <- function(datetime, lat, lon) {
    sun_angle <- getSunlightPosition(date = datetime, lat = lat, lon = lon)
    sun_altitude <- sun_angle$altitude
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

  # Diffuse radiation from canopy top
  F_sky_diff_v <- F_sky_diff_init
  # Direct radiation from canopy top
  F_sky_dir_v <<- F_sky_dir_init * sun_angles$above_horizon

  # Calculate fluxes per unique XY combination
  final_results_vertical <- voxel_grid[, {
    # Number of layers
    n_layers <- .N

    fluxes <- sw_two_stream(
      soil_reflect = omega_g_v, density = density,
      F_sky_diff = F_sky_diff_v, F_sky_dir = F_sky_dir_v,
      Kd = Kd_v, Kb = Kb_v,
      omega = omega, beta = betad, beta0 = beta0
    )

    .(X, Y, Z = seq_len(n_layers), density,
      F_d_down = fluxes$F_d_down[1:n_layers],
      F_d_up = fluxes$F_d_up[1:n_layers],
      F_b_down = fluxes$F_b_down[1:n_layers],
      net_sw = fluxes$net_sw[1:n_layers])
  }, by = .(X, Y)]


  setorder(final_results_vertical, X, Y, Z)

  ############################
  # Horizontal shortwave RTM #
  ############################

  # Direct radiation from eastern edge
  F_sky_dir_h <<- round(
    F_sky_dir_init / sin(sun_angles$sun_altitude) *
      cos(sun_angles$sun_altitude) * abs(cos(sun_angles$sun_azimuth)) *
      sun_angles$above_horizon * sun_angles$in_eastern_half, 0
  )
  # Diffuse radiation from eastern edge
  F_sky_diff_h <- F_sky_diff_init

  # Calculate fluxes per unique YZ combination
  final_results_horizontal <- voxel_grid[, {
    # Number of layers
    n_layers <- .N

    fluxes <- sw_two_stream(
      soil_reflect = omega_g_h, density = density,
      F_sky_diff = F_sky_diff_h, F_sky_dir = F_sky_dir_h,
      Kd = Kd_h, Kb = Kb_h,
      omega = omega, beta = betad, beta0 = beta0
    )

    .(Y, Z, X = seq_len(n_layers), density,
      F_d_down = fluxes$F_d_down[1:n_layers],
      F_d_up = fluxes$F_d_up[1:n_layers],
      F_b_down = fluxes$F_b_down[1:n_layers],
      net_sw = fluxes$net_sw[1:n_layers])
  }, by = .(Y, Z)]

  setorder(final_results_horizontal, Y, Z, X)

  #########################################
  # Combining both directions into 2D RTM #
  #########################################

  final_results_2D <- merge(
    final_results_horizontal[, .(X, Y, Z, density,
                                 F_d_down_h = F_d_down,
                                 F_d_up_h = F_d_up,
                                 F_b_down_h = F_b_down,
                                 net_sw_h = net_sw)],
    final_results_vertical[, .(X, Y, Z, density,
                               F_d_down_v = F_d_down,
                               F_d_up_v = F_d_up,
                               F_b_down_v = F_b_down,
                               net_sw_v = net_sw)],
    by = c("X", "Y", "Z", "density"),
    all = TRUE
  )

  # Sum fluxes
  final_results_2D[, `:=`(
    F_d_down = F_d_down_h + F_d_down_v,
    F_d_up = F_d_up_h + F_d_up_v,
    F_b_down = F_b_down_h + F_b_down_v,
    net_sw = net_sw_h + net_sw_v
  )]

  # Select relevant columns
  final_results_2D <- final_results_2D[, .(X, Y, Z, density, F_d_down, F_d_up, F_b_down, net_sw)]

  return(final_results_2D)
}



###############################################################################
#
#                               ForEdgeClim
#
# A microclimate model that solves the energy balance equation in fragmented
# forest areas: Rn = LE + H + G, where,
# Rn = Net radiation (shortwave & longwave)
# LE = Latent heat flux (~ evapotranspiration)
# H = sensible heat flux
# G = Ground heat flux
###############################################################################


#' Execute the microclimate model ForEdgeClim
#'
#' @param structure_grid Dataframe with structural density values
#' @param datetime Datetime object representing the current simulation time
#' @return Dataframe with simulated microclimate temperatures and fluxes
#' @importFrom dplyr group_by ungroup mutate left_join
#' @export
run_foredgeclim <- function(structure_grid, datetime) {


  #print('🚩 𝙁𝙤𝙧𝙀𝙙𝙜𝙚𝘾𝙡𝙞𝙢')
  #start = Sys.time()
  #print('════════════════════════')
  #############################################
  # 2D SHORTWAVE RADIATIVE TRANSFER MODELLING #
  #############################################
  #print('SW RTM ☀️')
  voxel_grid = structure_grid
  final_results_2D <<- shortwave_two_stream_RTM(datetime, lat, lon, voxel_grid,
                           F_sky_diff_init, F_sky_dir_init,
                           omega_g_v, omega_g_h,
                           Kd_v, Kb_v, Kd_h, Kb_h,
                           omega, betad, beta0) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.
  # saveRDS(final_results_2D, 'Output/SW_DATA.rds')

  den <<- final_results_2D$density
  net_SW <- final_results_2D$net_sw
  SW_down <- final_results_2D$F_d_down + final_results_2D$F_b_down
  SW_up <- final_results_2D$F_d_up
  x_coords <- voxel_grid$X
  y_coords <- voxel_grid$Y
  z_coords <<- voxel_grid$Z

  # dimensions grid
  x_dim <<- max(x_coords)  # Length in x-direction
  y_dim <<- max(y_coords)  # Length in y-direction
  z_dim <<- max(z_coords)  # Length in z-direction

  #############################
  # INITIALIZING TEMPERATURES #
  #############################

  #print('Init temp 🌡️')

  # First longwave RTM for initializing temperatures
  micro_grid <- data.frame(
    x = x_coords,
    y = y_coords,
    z = z_coords,
    temperature = macro_temp,
    T_soil = macro_temp
  )
  final_results_lw_2D = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                                F_sky_lw, omega_g_lw_v, omega_g_lw_h, macro_temp,
                                                Kd_lw_v, Kd_lw_h, omega_lw, beta_lw)
  net_LW <- final_results_lw_2D$net_lw
  LW_down <- final_results_lw_2D$F_d_down
  LW_up <- final_results_lw_2D$F_d_up
  net_rad_ground <- ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

  # Net radiation forest surface
  net_radiation_init <- net_SW + net_LW

  # Define ground flux as a % of net radiation at ground
  G = calculate_G(net_rad_ground)

  # Define T_soil at ground surface from ground flux with ground conductance k_soil
  micro_grid$T_soil = G * stable_soil_depth / k_soil + T_soil_deep

  # Define micro_grid$temperature based on macro temp and initial net radiation + random noise
  micro_grid$temperature =  macro_temp + 0.01*net_radiation_init + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)

  # First approximation of air temperature, based on forest surface temperature and random noise
  T_air_vec = micro_grid$temperature + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)

  # Define T_ground for every voxel in a vertical xy-column as T_soil
  micro_grid <- micro_grid |>
    group_by(x, y) |>
    mutate(T_ground = T_soil[z == 1][1])  |>
    ungroup()

  # Determine weights for each component in the linearisation via exponential weighting
  x_threshold = length_transect # x-value from where macro influence is being felt, ie, forest edge value
  z_threshold = height_canopy  # z-value from where macro influence is being felt, ie, canopy top value
  dis_macro_x = ifelse(x_coords > x_threshold, 0, x_threshold - x_coords + 1)
  dis_macro_z = ifelse(z_coords > z_threshold, 0, z_threshold - z_coords + 1)
  dis_soil = z_coords
  dis_forest = 0.5 # size of a voxel edge = 1m, so distance from structure is on average 0.5m
  alpha_macro = log(0.5) / infl_macro
  alpha_soil = log(0.5)/infl_soil
  alpha_forest = log(0.5)/infl_forest
  w_macro_x = exp(alpha_macro*dis_macro_x)
  w_macro_z = exp(alpha_macro*dis_macro_z)
  w_soil = exp(alpha_soil*dis_soil)
  w_forest = exp(alpha_forest*dis_forest)

  # Air temperature update based on a linearization as is done in microclimc (Maclean and Klinges, 2021) ad its successor microclimf.
  # The updated air temperature is weighted by the convection coefficient and the distance to each boundary.
  # This linearization simplifies the complex nonlinear interactions between leaves and air, making the calculations more efficient.
  # This linearization is mostly valid when there are small temperature differences between T air en T leaf, when there is
  # sufficient air streaming, when the net radiation or the humidity is not extreme, when the forest structure is quite homogeneous.
  # The linearization calculates the air temperature as influenced by the macro and soil boundary temperature and the surface temperature.
  T_air_vec = ( (w_macro_x + w_macro_z)*g_macro * macro_temp  + w_soil*g_soil * micro_grid$T_ground + w_forest*g_forest * micro_grid$temperature ) / ( (w_macro_x + w_macro_z)*g_macro + w_soil*g_soil + w_forest*g_forest)

  #########################################
  # ITERATION TO CLOSE THE ENERGY BALANCE #
  #########################################

  #print('E bal ⚖️')

  # Start iteration
  max_iter = 1000
  W = 1  # Weighting step, start value
  W_min = 0.01  # Weighting step, min value
  error_prev = Inf

  for (iter in 1:max_iter) {

    #print(paste0('-> Iteration ', iter))


    #################
    # NET RADIATION #
    #################

    #print('   LW RTM ⛅ ')
    final_results_lw_2D = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                                  F_sky_lw, omega_g_lw_v, omega_g_lw_h, macro_temp,
                                                  Kd_lw_v, Kd_lw_h, omega_lw, beta_lw) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.


    net_LW = final_results_lw_2D$net_lw
    LW_down = final_results_lw_2D$F_d_down
    LW_up = final_results_lw_2D$F_d_up
    net_rad_ground = ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

    # Net radiation for structure
    net_radiation = net_SW + net_LW


    ###############
    # Ground flux #
    ###############

    #print('   G 🟫')

    # G is positive if E is entering the soil
    G = calculate_G(net_rad_ground)
    micro_grid$T_soil = G*stable_soil_depth/k_soil + T_soil_deep # micro_grid$T_soil_deep
    micro_grid <- micro_grid |>
      group_by(x, y) |>
      mutate(T_ground = T_soil[z == 1][1])  |>
      ungroup()


    #############################
    # AIR TO AIR HEAT DIFFUSION #
    #############################

    #print('   D 💨')
    dt = 1 # diffusion step is 1, ie, the iteration step taken for temp convergence
    T_air_vec = T_air_vec - calculate_D(T_air_vec,micro_grid$T_ground)*dt/(Cp*V*rho)


    ######################
    # SENSIBLE HEAT FLUX #
    ######################

    #print('   H ♨️')
    sensible_flux = calculate_H(micro_grid$temperature,T_air_vec)


    ####################
    # LATENT HEAT FLUX #
    ####################

    #print('   LE 💧')
    latent_flux <- calculate_LE(micro_grid$temperature-273.15,net_radiation) # temp in calculate_LE must be in °C


    ##################
    # ENERGY BALANCE #
    ##################

    energy_balance_surf <- net_radiation - sensible_flux - latent_flux
    #print(paste0('   Max E_bal error = ', round(max(abs(energy_balance_surf)), 2)))
    error_current = max(abs(energy_balance_surf))
    # Check oscillations, if so, reduce W
    if (iter > 1 && error_current > error_prev) {
      W = max(W * 0.8, W_min)  # Reduce W by 20%, but never smaller than W_min
    }
    error_prev = error_current  # Update previous error

    # print(summary(net_radiation))
    # print(summary(sensible_flux))
    # print(summary(latent_flux))
    # print(summary(G))

    # Check convergence
    if (max(abs(energy_balance_surf)) < energy_balance_tolerance){
      #print(paste0('Convergence is reached after ', iter, ' iterations.'))
      break
    }

    # Use Newton's method to update surface temperature as is done in the SCOPE 2.0 model (Yang et al. 2021)
    micro_grid$temperature <- micro_grid$temperature - ifelse(energy_balance_surf == 0, 0, W*energy_balance_surf/de_dT(micro_grid$temperature,net_radiation) )

    # Determine planes-averaged surface temperature (to be used in T_air_Vec update for zero-density voxels).
    micro_grid <- micro_grid |>
      group_by(z) |>
      mutate(mean_temp_z = mean(temperature[den>0], na.rm = TRUE)) |>
      ungroup()
    micro_grid <- micro_grid |>
      group_by(y) |>
      mutate(mean_temp_y = mean(temperature[den>0], na.rm = TRUE)) |>
      ungroup()
    micro_grid <- micro_grid |>
      group_by(x) |>
      mutate(mean_temp_x = mean(temperature[den>0], na.rm = TRUE)) |>
      ungroup()
    micro_grid$mean_temp = (micro_grid$mean_temp_x + micro_grid$mean_temp_y + micro_grid$mean_temp_z)/3

    # Air temperature update based on a linearization as is done in microclimc (Maclean & Klinges, 2021).
    T_air_vec <- ( (w_macro_x + w_macro_z)*g_macro * macro_temp + w_soil*g_soil * micro_grid$T_ground + ifelse(den == 0, w_forest*g_forest * micro_grid$mean_temp, w_forest*g_forest *micro_grid$temperature) )/ ( (w_macro_x + w_macro_z)*g_macro + w_soil*g_soil + w_forest*g_forest)

    }

  res <- list(
    sw_rad_2D = final_results_2D, # SW RTM output
    lw_rad_2D = final_results_lw_2D, # LW RTM output
    micro_grid = micro_grid,
    air_temperature = T_air_vec,
    net_radiation = net_radiation,
    sensible_flux = sensible_flux,
    latent_flux = latent_flux,
    ground_flux = G
  )

  # Results
  # cat("Final surface temperature distribution (°C):\n")
  # print(summary(micro_grid$temperature - 273.15))
  # cat("Final air temperature distribution (°C):\n")
  # print(summary(T_air_vec- 273.15))
  # cat("Final soil temperature distribution (°C):\n")
  # print(summary(micro_grid$T_soil- 273.15))


  #finish = Sys.time()
  #print(paste0('Run time = ', round(as.numeric(finish - start, units = "secs"), 2), ' s'))
  #print('════════════════════════')
  #print('🏁 𝙁𝙤𝙧𝙀𝙙𝙜𝙚𝘾𝙡𝙞𝙢')

  return(res)

}




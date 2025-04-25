library(ForEdgeClim)

####################
# INPUT TIMESERIES #
####################

start_timeseries = Sys.time()
start_time <- as.POSIXct("2023-07-08 20:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-07-08 20:00:00", tz = "UTC")
datetime_series <- seq(start_time, end_time, by = "hour")

for (current_datetime in datetime_series) {

  # Make sure the current datetime is a POSIXct object
  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  print(paste0('Running model for ', current_datetime, ' ...'))

  ##########################
  # INPUT MODEL PARAMETERS #
  ##########################

  create_input_drivers(datetime = current_datetime)
  create_physical_constants()
  create_model_parameters()
  # parameters for shortwave RTM
  betad <<- 0.5          # Fraction of scattered diffuse radiation in backward direction
  beta0 <<- 0.5         # Fraction of scattered direct beam radiation in backward direction
  omega <<- 0.5        # Scattering coefficient
  # -> Vertical RTM
  Kd_v <<- 0.2          # Diffuse extinction coefficient for vertical radiation, per unit density
  Kb_v <<- 0.774         # Direct beam extinciton coefficient for vertical radiation, per unit density
  omega_g_v <<- 0.276     # Ground scattering
  # -> Horizontal RTM
  Kd_h <<- 0.15          # Diffuse extinciton coefficient for lateral radiation, per unit density
  Kb_h <<- 0.2          # Direct beam extinction coefficient for lateral radiation, per unit density
  omega_g_h <<- 0.15       # This is not ground scattering, but scattering by the inner forest. This is a balanced
  # value between reflection by leaves and transmission by open spaces.

  # parameters for longwave RTM
  e_atm <<- 0.842        # Emissivity atmosphere (0.9 on cloudy days)
  e_forest <<- 0.937     # Emissivity forest (leaves and wood)
  beta_lw <<- 0      # Fraction of scattered longwave radiation in backward direction
  omega_lw <<- 0.1      # Scattering coefficient for longwave radiation
  # -> Vertical RTM
  Kd_lw_v <<- 0.1       # Longwave extinction coefficient for vertical radiation, per unit density
  omega_g_lw_v <<- 0.05 # Ground scattering for longwave radiation
  # -> Horizontal RTM
  Kd_lw_h <<- 0.2       # Longwave extinciton coefficient for lateral radiation, per unit density
  omega_g_lw_h <<- 0.05    # This is not ground scattering, but scattering by the inner forest. This is a balanced
  # value between reflection by leaves and transmission by open spaces.

  # parameters for air to air heat convection (C)
  h <<- 10.86              # Convection coefficient of air (W/m2/K), higher value means more wind/more turbulence

  # parameters for update linearisation air temperature
  g_macro <<- 11.43          # Convective heat transfer coefficient between (voxel) air and (macro) air (W/m2/K)
  infl_macro <<- 50.03       # Distance over which the influence of macro temp on air temp is reduced by 50% (m)
  infl_soil <<- 6.028         # Distance over which the influence of soil temp on air temp is reduced by 50% (m)
  infl_forest <<- 6.169       # Distance over which the influence of forest temp on air temp is reduced by 50% (m)

  # parameters for sensible heat flux (H)
  g_forest <<- 10.188   # Combined conductive & convective heat transfer coefficient between (voxel) air and structure (leaf) (W/m2/K)

  # parameters for ground heat flux (G) and soil temperature
  p_ground <<- 0.303          # Fraction of net ground radiation to define ground flux
  g_soil <<- 9.233              # Convective heat transfer coefficient between (voxel) air of air layer just above the ground and ground surface (W/m2/K)
  k_soil <<- 0.629             # Thermal conductance soil (W/m/K)

  # Create output path
  output_path <- paste0('Output/', format(current_datetime, "%Y-%m-%d_%H"))
  # Make sure the output path exists
  if (!dir.exists(output_path)) dir.create(output_path)

  ######################
  # INPUT OBSERVATIONS #
  ######################

  # Import observations as input variables and as variables to compare the model with
  import_DTS_observations()
  import_RMI_observations()
  import_pyr_observations()
  import_soil_temperature()
  import_PAR_observations()


  #######################
  # Creation of 3D grid #
  #######################

  # print('3D grid 🌲🌳🌲🌳🌲')
  # voxel_TLS = generate_DTM_grid_TLS(las_file = TLS_input_file, voxel_size = 1)
  # saveRDS(voxel_TLS, TLS_filtered_file)
  voxel_TLS = readRDS(TLS_filtered_file)

  #############
  # Run model #
  #############

  res = run_foredgeclim(voxel_TLS$grid)
  saveRDS(res, paste0(output_path, '/model_results.rds'))

  ############
  # Plotting #
  ############

  # Digital Terrain Model & structural grid plots
  # plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid, output_path)
  #
  # Shortwave radiation plots
  # plots_sw(sw_rad_2D = res$sw_rad_2D, output_path)
  #
  # # Longwave radiation plots
  # plots_lw(lw_rad_2D = res$lw_rad_2D, output_path)
  #
  # Flux plots
  # plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path)
  #
  # # Rad at reqhgt plot
  # plot_rad(res$micro_grid, res$sw_rad_2D$F_d_down + res$sw_rad_2D$F_b_down, output_path)

  # Temperature plots
  # plots_temp(res$micro_grid, res$air_temperature, output_path)

}
end_timeseries = Sys.time()
print(paste0('Total running time timeseries = ', round(as.numeric(end_timeseries - start_timeseries, units = "secs"), 2), ' s'))



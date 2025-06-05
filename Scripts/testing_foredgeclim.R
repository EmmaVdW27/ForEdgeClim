library(ForEdgeClim)

####################
# INPUT TIMESERIES #
####################

start_timeseries = Sys.time()
start_time <- as.POSIXct("2023-07-08 12:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-07-08 12:00:00", tz = "UTC")
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

  #print('3D grid 🌲🌳🌲🌳🌲')
  # Structure from TLS las file
  #voxel_TLS = generate_DTM_grid_TLS(las_file = TLS_input_file, voxel_size = voxel_length)
  #saveRDS(voxel_TLS, TLS_filtered_file)
  voxel_TLS = readRDS(TLS_filtered_file)

  # Structure from TLS vox file
  # voxel_TLS = generate_DTM_grid_Vox(vox_file = vox_input_file)
  # saveRDS(voxel_TLS, vox_filtered_file)
  #voxel_TLS = readRDS(vox_filtered_file)

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
  # #
  # # # Longwave radiation plots
  # plots_lw(lw_rad_2D = res$lw_rad_2D, output_path)
  # #
  # # Flux plots
  # plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path)
  # #
  # # # Rad at reqhgt plot
  # plot_rad(res$micro_grid, res$sw_rad_2D$F_d_down + res$sw_rad_2D$F_b_down, output_path)
  #
  # # Temperature plots
  # plots_temp(res$micro_grid, res$air_temperature, output_path)

}
end_timeseries = Sys.time()
print(paste0('Total running time timeseries = ', round(as.numeric(end_timeseries - start_timeseries, units = "secs"), 2), ' s'))



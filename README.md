# ForEdgeClim

ForEdgeClim is an R package for modelling microclimates in forests, including radiation processes and heat transfer.
The model is initially written to simulate microclimate gradients along transect lines from a forest’s core towards its edge.
A more detailed overview with the applied physical equations can be found in the attached PDF file 'Formulae ForEdgeClim.pdf'.

## Installation
Install the package directly from GitHub:
```r
# Using devtools
install.packages("devtools")
devtools::install_github("EmmaVdW27/ForEdgeClim")
```

## Example functions
```r
library(ForEdgeClim)

####################
# INPUT TIMESERIES #
####################

start_timeseries = Sys.time()
start_time <- as.POSIXct("2023-07-08 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-07-08 23:00:00", tz = "UTC")
datetime_series <- seq(start_time, end_time, by = "hour")

for (current_datetime in datetime_series) {

  # Make sure the current datetime is a POSIXct object
  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  print(paste0('Running model for ', current_datetime, ' ...'))

  ##########################
  # INPUT MODEL PARAMETERS #
  ##########################

  create_input_data(datetime = current_datetime)
  energy_balance_tolerance <<- 2

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

  #######################
  # Creation of 3D grid #
  #######################

  print('3D grid 🌲🌳🌲🌳🌲')
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
  plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid, output_path)

  # Shortwave radiation plots
  plots_sw(sw_rad_2D = res$sw_rad_2D, output_path)

  # Longwave radiation plots
  plots_lw(lw_rad_2D = res$lw_rad_2D, output_path)

  # Flux plots
  plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path)

  # Temperature plots
  plots_temp(res$micro_grid, res$air_temperature,output_path)

}
end_timeseries = Sys.time()
print(paste0('Total running time timeseries = ', round(as.numeric(end_timeseries - start_timeseries, units = "secs"), 2), ' s'))



```

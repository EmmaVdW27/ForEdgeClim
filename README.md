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
###############################################################################
#
# This is an example script to show you how the microclimate model ForEdgeClim
# can be run. To run the model, you need to define several variables, input and
# output files. The code below is similar to the script 'testing_foredgeclim.R'.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################


library(ForEdgeClim)

####################
# INPUT TIMESERIES #
####################

start_timeseries = Sys.time()
start_time <- as.POSIXct("2023-07-08 12:00:00", tz = "UTC") # first hour you want to run the model for
end_time <- as.POSIXct("2023-07-08 12:00:00", tz = "UTC")   # last hour you want to run the model for
datetime_series <- seq(start_time, end_time, by = "hour")   # ForEdgeClim is run for all the hours between the first and last hour
TLS_input_file <- 'Data/2025-04-08_ForSe_Gontrode_5cm_transect_emma.las' # original TLS las file
TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_July2023.rds' # Once the structural voxel grid has been made using the
# function generate_DTM_grid_TLS() (see further below), you can also immediately use the TLS_filtered_file.

for (current_datetime in datetime_series) {

  # Make sure the current datetime is a POSIXct object
  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  print(paste0('Running model for ', current_datetime, ' ...'))

  #################################################
  # INPUT DRIVERS, CONSTANTS AND MODEL PARAMETERS #
  #################################################

  create_input_drivers()
  create_physical_constants()
  create_model_parameters()

  # Create output path
  output_path <- paste0('Output/', format(current_datetime, "%Y-%m-%d_%H"))
  # Make sure the output path exists
  if (!dir.exists(output_path)) dir.create(output_path)

  ########################################################
  # EXTRACT THE OBSERVATIONS FROM THE INPUT DRIVER FILES #
  ########################################################

  # Import observations as input variables and as variables to compare the model with
  import_RMI_observations(current_datetime)
  # In case you don't have longwave downward radiation values from the sky, you can use the
  # macrotemperature as a proxy:
  # if(is_empty(F_sky_lw)){
  #   assign("F_sky_lw", sigma_SB * 0.75 * macro_temp^4, envir = .GlobalEnv)
  # }
  import_pyr_observations(current_datetime)
  import_soil_temperature(current_datetime)

  #######################
  # CREATION OF 3D GRID #
  #######################

  #print('3D grid 🌲🌳🌲🌳🌲')
  # Structure from TLS las file
  voxel_TLS = generate_DTM_grid_TLS(las_file = TLS_input_file, voxel_size = voxel_length)
  saveRDS(voxel_TLS, TLS_filtered_file)
  #voxel_TLS = readRDS(TLS_filtered_file)

  #############
  # RUN MODEL #
  #############

  res = run_foredgeclim(voxel_TLS$grid, current_datetime)
  saveRDS(res, paste0(output_path, '/model_results.rds'))

  ############
  # PLOTTING #
  ############

  # Digital Terrain Model & structural grid plots
  plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid, output_path)

  # Shortwave radiation plots
   plots_sw(sw_rad_2D = res$sw_rad_2D, output_path)

  # Longwave radiation plots
  plots_lw(lw_rad_2D = res$lw_rad_2D, output_path)

  # Flux plots
  plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path, current_datetime)

  # Temperature plots
   plots_temp(res$micro_grid, res$air_temperature, output_path, current_datetime)

}
end_timeseries = Sys.time()
print(paste0('Total running time timeseries = ', round(as.numeric(end_timeseries - start_timeseries, units = "secs"), 2), ' s'))






```

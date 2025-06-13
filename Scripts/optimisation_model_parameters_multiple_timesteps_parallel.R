############################################################################################
# This is a script to optimise the model parameters of the microclimate model ForEdgeClim.
# For this, we will use an optimisation algorithm.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


library(ForEdgeClim)
library(ggplot2)
library(dplyr)
library(cmaes)
library(future.apply)

plan(multisession, workers = 10) # of multicore op Linux

start_analysis = Sys.time()

#########
# INPUT #
#########

# List of time steps for which we want to optimise simultaneously
all_datetimes <- as.POSIXct(c("2023-07-08 08:00:00",
                   "2023-07-08 12:00:00",
                   "2023-07-08 16:00:00",
                   "2023-07-08 20:00:00",
                   "2023-07-08 00:00:00",
                   "2023-07-08 04:00:00"), tz = "UTC")


logfile <- "Output/parameter_optimisation/optimisation_log_CMA-ES_20230708-6timesteps_ngen5_iteration3.csv"
results_output_file = "Output/parameter_optimisation/calibration_results_iteration3.rds"
max_it = 5 # max amount of generations in the CMA-ES algorithm; set min to 5 for minimal optimisation
#lambda = 2 # amount of offspring for every generation; around 13 for 26 parameters
stop_fitness = 1 # 1°C as convergence criterium for RMSE

# During CMA-ES run, current fitness output = the RMSE of the parent vector (params)
# that is passed to the next generation. The parent vector is the mean vector of the
# mu (here around 6 if lambda is around 13) best fits of the lambda offsprings.
# Run time CMA-ES is around no_timesteps x max_it x lambda x model_runtime.
# Total run time for 2 timesteps = 1u 51min.


#################
# START PROGRAM #
#################

iteration_counter <- 0

# Make log file with column names
param_names = c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h", "e_forest",
                "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h", "h", "g_macro", "infl_macro", "infl_soil", "infl_forest",
                "g_forest", "p_ground", "g_soil", "k_soil")

log_header <- data.frame(
  iteration = integer(),
  rmse = numeric()
)

for (p in param_names) {
  log_header[[p]] <- numeric()
}

write.table(log_header, file = logfile, sep = ",", row.names = FALSE, col.names = TRUE)


##################################
# INITIAL VALUE MODEL PARAMETERS #
##################################

initial_model_parameters = c(0.45, 0.25, 0.2218234, 0.45, 0.7, 0.08, 0.6, 1, 0.001,
                             0.99, 0.001, 0.001, 0.951, 0.001, 0.95, 0.001, 9.938226, 11.78616, 50.99152, 4.205386,
                             6.992731, 12.24181, 0.35, 7.587941, 0.3140566)



######################
# OBJECTIVE FUNCITON #
######################

# This is the function to be minimized or maximized


cost_function <- function(par) {


  # Parallellisation over time steps
  rmse_vector <- future_lapply(all_datetimes, function(dt, par_values) {
    tryCatch({
      datetime <- as.POSIXct(unlist(dt), origin = "1970-01-01", tz = "UTC")

      create_input_drivers()
      req_height <- 1
      length_transect <- 135

      create_physical_constants()
      create_model_parameters()
      for(i in seq_along(param_names)) {
        assign(param_names[i], par_values[i], envir = .GlobalEnv)
      }


      import_DTS_observations(datetime)
      import_RMI_observations(datetime)
      import_pyr_observations(datetime)
      import_soil_temperature(datetime)

      #############
      # RUN MODEL #
      #############

      voxel_TLS = readRDS(TLS_filtered_file)
      res = run_foredgeclim(voxel_TLS$grid, datetime)
      micro_grid = res$micro_grid
      air_temp = res$air_temperature


      if (any(is.na(air_temp)) || any(is.infinite(air_temp))) {
        cat("air_temp is NA or inf")
        return(1e6)
      }

      ##################################################
      # CALCULATE ERROR BETWEEN MODEL AND OBSERVATIONS #
      ##################################################

      # TOMST horizontal observations
      TOMST_hor <- read.csv(paste0("Data/TOMST_filtered_distance_temp_", format(datetime, "%Y%m%d_%H%M"), ".csv")) |>
        mutate(D_edge = 135 - D_edge) |>
        mutate(D_edge = ifelse(D_edge == 0, 1, D_edge)) |> # Like this, 10 TOMST measurements can be used
        rename(position_X_or_Z = D_edge) |>
        arrange(position_X_or_Z)
      if (any(is.na(TOMST_hor)) || any(is.infinite(TOMST_hor))) {
        cat("TOMST_hor is NA or inf \n")
      }

      # TOMST vertical observations
      TOMST_ver <- read.csv(paste0("Data/TOMST_filtered_height_temp_", format(datetime, "%Y%m%d_%H%M"), ".csv")) |>
        rename(position_X_or_Z = height) |>
        filter(position_X_or_Z != 0) |> # remove sensor at ground level since this one is already included in TOMST_hor
        arrange(position_X_or_Z)
      if (any(is.na(TOMST_ver)) || any(is.infinite(TOMST_ver))) {
        cat("TOMST_ver is NA or inf \n")
      }

      # Combine TOMST datasets
      TOMST_all <- bind_rows(TOMST_hor, TOMST_ver)
      if (any(is.na(TOMST_all)) || any(is.infinite(TOMST_all))) {
        cat("TOMST_all is NA or inf \n")
      }


      # Define grid with air temperatures
      temp_air_grid = micro_grid
      # Convert air temperature values from Kelvin to degrees Celsius
      temp_air_grid$temperature = air_temp - 273.15

      # Modelled air temperature at reqhgt, horizontally
      reqhgt <- temp_air_grid |>
        filter(z == req_height, y == 15, x <= length_transect) |>
        filter(x %in% TOMST_hor$position_X_or_Z) |>
        arrange(x)

      # Modelled air temperature along the vertical line
      vertical <- temp_air_grid |>
        filter(x == 135 - 75, y == 15) |> # tower is positioned at 75m from eastern side transect
        filter(z %in% TOMST_ver$position_X_or_Z) |>
        arrange(z)

      # Combine model datasets
      model_all <- bind_rows(reqhgt, vertical)
      if (any(is.na(model_all)) || any(is.infinite(model_all))) {
        cat("model_all is NA or inf \n")
      }


      # Calculate RMSE
      rmse <- sqrt(mean((model_all$temperature - TOMST_all$Tair)^2, na.rm = TRUE))
      if (is.na(rmse) || is.infinite(rmse)) rmse <- 1e6
      return(rmse)

    }, error = function(e) {
      cat("Error at datetime:", datetime, "-->", e$message, "\n")
      return(1e6)
    })
  }, par_values = par, future.seed = TRUE)

  final_rmse <- mean(rmse_vector)

  # write to log file
  iteration_counter <<- iteration_counter + 1
  log_entry <- data.frame(
    iteration = iteration_counter,
    rmse = final_rmse,
    t(par)
  )
  write.table(log_entry, file = logfile, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)

  return(final_rmse)



}


##########################
# OPTIMISATION ALGORITHM #
##########################

# Initial values for model parameters
start_vals <- initial_model_parameters

# Optimisation algorithm = CMA-ES = Covariance matrix adaptation evolution strategy
# Evolutionary strategy: tries to find solutions like evolution does (mutations, selections)
# Covariance matrix adaptation: learns which directions in the parameter space bring the greatest gains,
# and cleverly adjusts its search accordingly.
opt_result <<- cma_es(
  par = start_vals,
  fn = cost_function,
  lower = c(0.45, 0.25, 0.2, 0.25, 0.4, 0.08, 0.25, 0.4, 0,
            0.94, 0, 0, 0.95, 0, 0.95, 0, 2, 5, 20,
            2, 1, 5, 0.1, 2, 0.1),
  upper = c(0.55, 0.4, 0.5, 0.45, 0.7, 0.18, 0.6, 1, 0.001,
            0.99, 0.001, 0.001, 0.951, 0.001, 0.951, 0.001, 15, 30, 60,
            10, 10, 25, 0.35, 10, 1.5),
  control = list(maxit = max_it, trace = TRUE, stopfitness = stop_fitness)
)

saveRDS(opt_result, file = results_output_file)

# Get optimisation results
print(opt_result$par)     # optimised parameters
print(opt_result$value)   # minimal RMSE
log_df <- read.csv(logfile)
best <- log_df[which.min(log_df$rmse), ]
print(paste("Best RMSE: ", best$rmse))
print(best)


end_analysis = Sys.time()
print(paste0('Total running time optimisation algorithm = ', round(as.numeric(end_analysis - start_analysis, units = "secs"), 2), ' s'))

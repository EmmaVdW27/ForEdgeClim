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


start_analysis = Sys.time()

#########
# INPUT #
#########

current_datetime = "2023-07-08 08:00:00" # runnen voor 8h, 12h and 23h
logfile <- "Output/parameter_optimisation/optimisation_log_CMA-ES_20230708-08_maxit5.csv"
max_it = 5


#################
# START PROGRAM #
#################

iteration_counter <- 0
param_names = c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h", "e_atm", "e_forest",
                "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h", "h", "g_macro", "infl_macro", "infl_soil", "infl_forest",
                "g_forest", "p_ground", "g_soil", "k_soil")
# Make log file with column names
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

initial_model_parameters = c(0.3, 0.2, 0.85, 0.3, 0.75, 0.2, 0.2, 0.3, 0, 0.65,
                             0.95, 0.03, 0.1, 0.2, 0.02, 0.1, 0, 10, 10, 50, 5,
                             6, 10, 0.10, 10, 0.8)

######################
# OBJECTIVE FUNCITON #
######################

# This is the function to be minimized or maximized

cost_function <- function(par) {

  ####################
  # MODEL PARAMETERS #
  ####################
  betad <<- par[1]; beta0 <<- par[2]; omega <<- par[3]; Kd_v <<- par[4];
  Kb_v <<- par[5]; omega_g_v <<- par[6]; Kd_h <<- par[7];  Kb_h <<- par[8];
  omega_g_h <<- par[9];  e_atm <<- par[10];  e_forest <<- par[11];
  beta_lw <<- par[12];  omega_lw <<- par[13]  Kd_lw_v <<- par[14];
  omega_g_lw_v <<- par[15];  Kd_lw_h <<- par[16];  omega_g_lw_h <<- par[17];
  h <<- par[18];  g_macro <<- par[19];  infl_macro <<- par[20];
  infl_soil <<- par[21];  infl_forest <<- par[22];  g_forest <<- par[23];
  p_ground <<- par[24];  g_soil <<- par[25];  k_soil <<- par[26]

  tryCatch({
    #################
    # INPUT DRIVERS #
    #################
    create_input_drivers(datetime = current_datetime)

    ######################
    # PHYSICAL CONSTANTS #
    ######################
    create_physical_constants()

    ######################
    # INPUT OBSERVATIONS #
    ######################

    # Import observations as input variables and as variables to compare the model with
    #import_DTS_observations()
    import_RMI_observations()
    import_pyr_observations()
    import_soil_temperature()
    #import_PAR_observations()


    #######################
    # Creation of 3D grid #
    #######################

    print('3D grid 🌲🌳🌲🌳🌲')
    # voxel_TLS = generate_DTM_grid_TLS(las_file = TLS_input_file, voxel_size = 1)
    # saveRDS(voxel_TLS, TLS_filtered_file)
    voxel_TLS = readRDS(TLS_filtered_file)

    #############
    # RUN MODEL #
    #############

    res = run_foredgeclim(voxel_TLS$grid)
    micro_grid = res$micro_grid
    air_temp = res$air_temperature



    # Check for NaN/Inf in the output
    if (any(is.na(air_temp)) || any(is.infinite(air_temp))) {
      warning("NaN or Inf detected")
      return(1e6)  # RMSE gets a high value
    }

    ##################################################
    # CALCULATE ERROR BETWEEN MODEL AND OBSERVATIONS #
    ##################################################

    # TOMST observations
    TOMST <- read.csv("Data/TOMST_filtered_distance_temp.csv") |>
      mutate(D_edge = 135 - D_edge) |>
      mutate(D_edge = ifelse(D_edge == 0, 1, D_edge)) |> # Like this, 10 TOMST measurements can be used
      arrange(D_edge)

    # Define grid with air temperatures
    temp_air_grid = micro_grid
    # Convert air temperature values from Kelvin to degrees Celsius
    temp_air_grid$temperature = air_temp - 273.15

    # Modelled air temperature at reqhgt
    reqhgt <- temp_air_grid |>
      filter(z == req_height, y == 15, x <= length_transect) |>
      filter(x %in% TOMST$D_edge) |>
      arrange(x)

    # Calculate RMSE
    rmse <- sqrt(mean((reqhgt$temperature - TOMST$Tair)^2, na.rm = TRUE))
    if (is.na(rmse) || is.infinite(rmse)) return(1e6)

    iteration_counter <<- iteration_counter + 1
    log_entry <- data.frame(
      iteration = iteration_counter,
      rmse = rmse,
      t(par)
    )
    write.table(log_entry, file = logfile, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)


    return(rmse)

  }, error = function(e) {
    warning(paste("Crash during model run:", e$message))
    return(1e6)  # high RMSE if model crashes
  })
}


############################
# OPTIMISATION ALGORITHM #
############################

# Initial values for model parameters
start_vals <- initial_model_parameters

# Optimisation algorithm = CMA-ES = Covariance matrix adaptation evolution strategy
# Evolutionary strategy: tries to find solutions like evolution does (mutations, selections)
# Covariance matrix adaptation: learns which directions in the parameter space bring the greatest gains,
# and cleverly adjusts its search accordingly.
# 15 min with maxit = 2
opt_result <<- cma_es(
  par = start_vals,
  fn = cost_function,
  lower = c(0.15, 0.25, 0.5, 0.2, 0.6, 0.1, 0.15, 0.2, 0,
            0.6, 0.92, 0, 0, 0.1, 0, 0.05, 0, 5, 5, 10,
            1, 1, 5, 0.1, 5, 0.2),
  upper = c(0.5, 0.5, 0.95, 0.5, 0.9, 0.3, 0.3, 0.4, 0.15,
            0.9, 0.98, 0.2, 0.1, 0.3, 0.05, 0.2, 0.05, 15, 15, 60,
            10, 15, 15, 0.35, 15, 2),
  control = list(maxit = max_it, trace = TRUE, stopfitness = 0.3)
  # for minimal optimisation maxit should be 5
  # 5 generations (each with around 13 evaluations for 26 parameters)
  # running time around 5 x 13 x 40s = around 43 min
  # for timestep 2023/07/08-12h this took 35 min
)

# Get optimisation results
print(opt_result$par)     # optimised parameters
print(opt_result$value)   # minimal RMSE
log_df <- read.csv(logfile)
best <- log_df[which.min(log_df$rmse), ]
print(paste("Best RMSE: ", best$rmse))
print(best)


end_analysis = Sys.time()
print(paste0('Total running time optimisation algorithm = ', round(as.numeric(end_analysis - start_analysis, units = "secs"), 2), ' s'))

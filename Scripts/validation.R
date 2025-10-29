############################################################################################
# This script validates ForEdgeClim's model output with observational TOMST TMS sensors.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClim)
  library(dplyr)
  library(Metrics) # voor r2, nse etc.
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(future)
  library(rlang) # voor is_empty
  library(future.apply)
})

# Parallelisation on HPC (supercomputer)
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

params = list(
  # Calibrated parameters (HEAT)
  g_macro = 25.28, infl_macro = 32.5, infl_soil = 5,
  # uncalibrated, fixed values: g_macro = 25, infl_macro = 32.5, infl_soil = 5,

  # Fixed parameters (SW RTM)
  betad = 0.325, beta0 = 0.325, omega = 0.52, Kd_v = 0.775, Kb_v = 1.25, omega_g_v = 0.13,
  Kd_h = 0.725, Kb_h = 1.15, omega_g_h = 0.15,
  # Fixed parameters (LW RTM)
  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035, Kd_lw_v = 0.3, omega_g_lw_v = 0.055, Kd_lw_h = 0.3, omega_g_lw_h = 0.035,
  # Fixed parameters (HEAT)
  h = 10, infl_forest = 5, g_forest = 12.5, p_ground = 0.225, g_soil = 10, k_soil = 1.225
)
# Dates for which you want to validate the model

season_string = "summer_uncalibrated" # summer spring autumn winter |_uncalibrated

all_datetimes <- seq(
    from = as.POSIXct("2023-07-01 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2023-07-31 23:00:00", tz = "UTC"),
    by   = "hour"
  )

# Seasonal structure
structure <- "Data/TLS_scaled_DTM_and_grid_July2023.rds"
structure_PAI_scaling_factor <- 5.91 / 6.18



# -----------------
# VALIDATION: STATISTICAL COMPARISON MODEL VS OBSERVATIONS
# -----------------

# Parallelize
parts <- future_lapply(
  all_datetimes,
  FUN = function(dt) {

    # make params accessible in de worker-GlobalEnv
    list2env(params, envir = .GlobalEnv)

    # model inputs
    create_input_drivers()
    voxel_TLS <- readRDS(structure)
    voxel_TLS$grid$density = voxel_TLS$grid$density * structure_PAI_scaling_factor
    create_physical_constants()

    datetime <- as.POSIXct(dt, tz = "UTC")

    out <- tryCatch({

      # Import observations
      import_RMI_observations(datetime)
      if(is_empty(F_sky_lw)){
        assign("F_sky_lw", sigma_SB * 0.75 * macro_temp^4, envir = .GlobalEnv)
      }
      import_pyr_observations(datetime)
      import_soil_temperature(datetime)

      # Run model
      res <- run_foredgeclim(voxel_TLS$grid, datetime)
      micro_grid <- res$micro_grid
      air_temp   <- res$air_temperature
      if (any(!is.finite(air_temp))) stop("NaN/Inf in air_temp")

      # Calculate error between observations and model
      # Observations
      key <- format(datetime, "%Y%m%d_%H%M")

      TOMST_hor <- read.csv(paste0("Data/TOMST_filtered_distance_temp_", key, ".csv"))
      TOMST_hor <- TOMST_hor %>%
        mutate(D_edge = 135 - D_edge, D_edge = ifelse(D_edge == 0, 1, D_edge)) %>%
        rename(position_X_or_Z = D_edge) %>%
        arrange(position_X_or_Z)

      TOMST_ver <- read.csv(paste0("Data/TOMST_filtered_height_temp_", key, ".csv")) %>%
        rename(position_X_or_Z = height) %>%
        filter(position_X_or_Z != 0) %>%
        arrange(position_X_or_Z)

      TOMST_all <- bind_rows(TOMST_hor, TOMST_ver)

      # Model
      temp_air_grid <- micro_grid
      temp_air_grid$temperature <- air_temp - 273.15

      reqhgt <- temp_air_grid %>%
        filter(z == req_height, y == 15, x <= length_transect) %>%
        filter(x %in% TOMST_hor$position_X_or_Z) %>%
        arrange(x)

      vertical <- temp_air_grid %>%
        filter(x == 135 - 75, y == 15) %>%
        filter(z %in% TOMST_ver$position_X_or_Z) %>%
        arrange(z)

      model_all <- bind_rows(reqhgt, vertical)

      # error
      if (nrow(model_all) != nrow(TOMST_all)) stop("Model-observation lengths mismatch")
      diffs <- model_all$temperature - TOMST_all$Tair
      valid <- is.finite(diffs)
      if (!any(valid)) stop("No valid differences")

      list(
        sse = sum(diffs[valid]^2),
        n   = sum(valid),
        data = tibble::tibble(
          datetime = datetime,
          obs = TOMST_all$Tair[valid],
          mod = model_all$temperature[valid],
          resid = model_all$temperature[valid] - TOMST_all$Tair[valid],
          TMS_position = TOMST_all$position_X_or_Z[valid]
        )
      )
    }, error = function(e) {
      warning(paste("Crash during model run:", format(datetime), "-", e$message))
      list(sse = NA_real_, n = 0L)
    })

    out
  },
  future.seed = TRUE  # reproducible results over workers (cores)
)

# combine data
val_data <- dplyr::bind_rows(lapply(parts, `[[`, "data"))

# save validated data
out_dir  <- file.path("Output", "validation")
out_file <- file.path(out_dir, paste0("val_data_", season_string, ".rds"))
saveRDS(val_data, out_file)

# calculate statistics
mean_error <- mean(val_data$resid, na.rm = TRUE)
rmse       <- sqrt(mean(val_data$resid^2, na.rm = TRUE))
r2         <- cor(val_data$obs, val_data$mod)^2
nse        <- 1 - sum(val_data$resid^2) / sum((val_data$obs - mean(val_data$obs))^2)


plan(sequential)

# Console summaries
cat("\nValidation statistics:\n")
cat(sprintf("Mean Error (Bias): %.2f °C\n", mean_error))
cat(sprintf("RMSE: %.2f °C\n", rmse))
cat(sprintf("R²: %.2f\n", r2))
cat(sprintf("NSE: %.2f\n", nse))



end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n', as.numeric(end_analysis - start_analysis, units = 'mins')))

############################################################################################
# This script optimizes ForEdgeClim parameters with the CMA-ES algorithm (via cmaesr).
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClim)
  library(cmaesr)
  library(smoof)
  library(ParamHelpers)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(future)
  library(future.apply)
  library(lubridate)
})

# Parallelization on HPC (supercomputer)
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

# We calibrate on 3 days per season: the most sunny day, the most cloudy day
# and the day with most solar fluctuations.

date_string <- 'spring'

all_datetimes <- c(
  seq(
    from = as.POSIXct("2025-04-30 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2025-04-30 23:00:00", tz = "UTC"),
    by   = "hour"
  ),
  seq(
    from = as.POSIXct("2025-04-23 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2025-04-23 23:00:00", tz = "UTC"),
    by   = "hour"
  ),
  seq(
    from = as.POSIXct("2024-04-02 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2024-04-02 23:00:00", tz = "UTC"),
    by   = "hour"
  )
)

# For spring and autumn, there are 2 years to be taken into account,
# for summer and winter, there is only 1 year.
structures <- list(
  `2024` = readRDS("Data/TLS_scaled_DTM_and_grid_April2024.rds"),
  `2025` = readRDS("Data/TLS_scaled_DTM_and_grid_April2025.rds")
)

scales <- c(
  `2024` = 2.98 / 6.18,
  `2025` = 2.24 / 6.18
)

structures_scaled <- lapply(names(structures), function(yr){
  vox <- structures[[yr]]
  vox$grid$density <- vox$grid$density * scales[[yr]]
  vox
})
names(structures_scaled) <- names(structures)

param_set <- "top_3"   # "all" | "focused" | "top_3"

max_it        <- 39     # generations (there will be max_it + 1 generations)
stop_fitness  <- 1      # RMSE target (°C)
lambda        <- NULL   # if NULL we set a heuristic later

# Output paths
log_dir     <- file.path("Output", "calibration")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
logfile_offspring <- file.path(log_dir, paste0("logging_offspring_", date_string, ".csv"))
logfile_gener     <- file.path(log_dir, paste0("logging_generation_", date_string, ".csv"))
results_output_file <- file.path(log_dir, paste0("calibration_results_", date_string, ".rds"))

# -----------------
# PARAMETER NAMES & DEFAULTS
# -----------------
if (param_set == "all") {
  param_names <- c(
    "betad","beta0","omega","Kd_v","Kb_v","omega_g_v","Kd_h","Kb_h","omega_g_h",
    "e_forest","beta_lw","omega_lw","Kd_lw_v","omega_g_lw_v","Kd_lw_h","omega_g_lw_h",
    "h","g_macro","infl_macro","infl_soil","infl_forest","g_forest","p_ground","g_soil","k_soil"
  )
} else if (param_set == "focused") {
  param_names <- c("h","g_macro","infl_macro","infl_soil","infl_forest","g_forest","p_ground","g_soil","k_soil")
  # Fixed parameters (SW RTM)
  betad <<- 0.325; beta0 <<- 0.325; omega <<- 0.52; Kd_v <<- 0.775; Kb_v <<- 1.25; omega_g_v <<- 0.13;
  Kd_h <<- 0.725; Kb_h <<- 1.15; omega_g_h <<- 0.15
  # Fixed parameters (LW RTM)
  e_forest <<- 0.965; beta_lw <<- 0.325; omega_lw <<- 0.035; Kd_lw_v <<- 0.3; omega_g_lw_v <<- 0.055; Kd_lw_h <<- 0.3; omega_g_lw_h <<- 0.035
} else { # top_3
  param_names <- c("g_macro", "g_forest", "infl_macro")
  # Fixed parameters (SW RTM)
  betad <<- 0.325; beta0 <<- 0.325; omega <<- 0.52; Kd_v <<- 0.775; Kb_v <<- 1.25; omega_g_v <<- 0.13;
  Kd_h <<- 0.725; Kb_h <<- 1.15; omega_g_h <<- 0.15
  # Fixed parameters (LW RTM)
  e_forest <<- 0.965; beta_lw <<- 0.325; omega_lw <<- 0.035; Kd_lw_v <<- 0.3; omega_g_lw_v <<- 0.055; Kd_lw_h <<- 0.3; omega_g_lw_h <<- 0.035
  # Fixed parameters (HEAT)
  h <<- 10; infl_soil <<- 5; infl_forest <<- 5; g_forest <<- 12.5; p_ground <<- 0.225; g_soil <<- 10; k_soil <<- 1.225
}

# Initial values (mean value from uniform parameter distributions)
if (param_set == "all") {
  start_vals <- c(0.325,0.325,0.52,0.775,1.25,0.13,0.725,1.15,0.15,
                  0.965,0.325,0.035,0.3,0.055,0.3,0.035,
                  10,25,32.5,5,5,12.5,0.225,10,1.225)
} else if (param_set == "focused") {
  start_vals <- c(10,25,32.5,5,5,12.5,0.225,10,1.225)
} else { # top_3
  start_vals <- c(25,32.5,5)
}

# Bounds
if (param_set == "all") {
  lower_bounds <- c(0.3,0.2,0.43,0.6,0.5,0.08,0.5,0.3,0.1, 0.94,0.3,0.01,0.2,0.04,0.2,0.01, 0,10,5,0,0,5,0.1,5,0.25)
  upper_bounds <- c(0.35,0.45,0.61,0.95,2,0.18,0.95,2,0.2, 0.99,0.35,0.06,0.4,0.07,0.4,0.06, 20,40,60,10,10,20,0.35,15,2.2)
} else if (param_set == "focused") {
  lower_bounds <- c(0,10,5,0,0,5,0.1,5,0.25)
  upper_bounds <- c(20,40,60,10,10,20,0.35,15,2.2)
} else { # top_3
  lower_bounds <- c(10,5,0)
  upper_bounds <- c(40,60,10)
}

# -----------------
# MODEL INPUTS
# -----------------
create_input_drivers()

create_physical_constants()

# -----------------
# COST = OBJECTIVE FUNCTION (RMSE)
# -----------------
compute_rmse <- function(par) {
  # Parallelize over all 24h (each core for 1h)
  parts <- future_lapply(
    all_datetimes,
    FUN = function(dt) {

      # Assign parameters to globals used in ForEdgeClim
      if (param_set == "all") {
        betad <<- par[1];  beta0 <<- par[2];  omega <<- par[3];  Kd_v <<- par[4];  Kb_v <<- par[5];  omega_g_v <<- par[6]
        Kd_h <<- par[7];   Kb_h <<- par[8];   omega_g_h <<- par[9]; e_forest <<- par[10]
        beta_lw <<- par[11]; omega_lw <<- par[12]; Kd_lw_v <<- par[13]; omega_g_lw_v <<- par[14]; Kd_lw_h <<- par[15]; omega_g_lw_h <<- par[16]
        h <<- par[17]; g_macro <<- par[18]; infl_macro <<- par[19]; infl_soil <<- par[20]; infl_forest <<- par[21]; g_forest <<- par[22]
        p_ground <<- par[23]; g_soil <<- par[24]; k_soil <<- par[25]
      } else if (param_set == "focused") {
        h <<- par[1]; g_macro <<- par[2]; infl_macro <<- par[3]; infl_soil <<- par[4]; infl_forest <<- par[5]; g_forest <<- par[6]
        p_ground <<- par[7]; g_soil <<- par[8]; k_soil <<- par[9]
      } else { # top_3
        g_macro <<- par[1]; infl_macro <<- par[2]; infl_soil <<- par[3]
      }

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
        yr <- as.character(lubridate::year(dt))
        voxel_TLS = structures_scaled[[yr]]
        res <- run_foredgeclim(voxel_TLS$grid, datetime)
        micro_grid <- res$micro_grid
        air_temp   <- res$air_temperature
        if (any(!is.finite(air_temp))) stop("NaN/Inf in air_temp")

        # Calculate error between observations and model
        # Observations
        key <- format(datetime, "%Y%m%d_%H%M")

        TOMST_hor <- read.csv(paste0("Data/TOMST_filtered_distance_temp_", key, ".csv"))
        TOMST_hor <- TOMST_hor %>%
          dplyr::mutate(D_edge = 135 - D_edge, D_edge = ifelse(D_edge == 0, 1, D_edge)) %>%
          dplyr::rename(position_X_or_Z = D_edge) %>%
          dplyr::arrange(position_X_or_Z)

        TOMST_ver <- read.csv(paste0("Data/TOMST_filtered_height_temp_", key, ".csv")) %>%
          dplyr::rename(position_X_or_Z = height) %>%
          dplyr::filter(position_X_or_Z != 0) %>%
          dplyr::arrange(position_X_or_Z)

        # Double the vertical points since there are twice as many horizontal points
        TOMST_ver <- rbind(TOMST_ver, TOMST_ver)

        TOMST_all <- dplyr::bind_rows(TOMST_hor, TOMST_ver)


        # Model
        temp_air_grid <- micro_grid
        temp_air_grid$temperature <- air_temp - 273.15

        reqhgt <- temp_air_grid %>%
          dplyr::filter(z == req_height, y == 15, x <= length_transect) %>%
          dplyr::filter(x %in% TOMST_hor$position_X_or_Z) %>%
          dplyr::arrange(x)

        vertical <- temp_air_grid %>%
          dplyr::filter(x == 135 - 75, y == 15) %>%
          dplyr::filter(z %in% TOMST_ver$position_X_or_Z) %>%
          dplyr::arrange(z)

        # Double the vertical points since there are twice as many horizontal points
        vertical <- rbind(vertical, vertical)

        model_all <- dplyr::bind_rows(reqhgt, vertical)

        # error
        if (nrow(model_all) != nrow(TOMST_all)) stop("Model/observation lengths mismatch")
        diffs <- model_all$temperature - TOMST_all$Tair
        valid <- is.finite(diffs)
        if (!any(valid)) stop("No valid differences")

        list(sse = sum(diffs[valid]^2), n = sum(valid))
      }, error = function(e) {
        warning(paste("Crash during model run:", format(datetime), "-", e$message))
        list(sse = NA_real_, n = 0L)
      })

      out
    },
    future.seed = TRUE  # reproducible results over workers (cores)
  )

  # aggregation over 24h
  sse_total <- sum(vapply(parts, `[[`, numeric(1), "sse"), na.rm = TRUE)
  n_total   <- sum(vapply(parts, `[[`, integer(1), "n"))
  if (n_total == 0) return(1e6)
  sqrt(sse_total / n_total)  # RMSE
}



# State for logging without monitor callbacks
.eval_counter <- 0L # variable accessible in multiple functions '.' is a convention: 'do not touch'


# -------------------
# Objective wrapper to log each offspring evaluation
# -------------------
# Wrapper around objective function so logging is performed
wrapped_obj_fn <- function(x) {
  rmse <- compute_rmse(x)
  .eval_counter <<- .eval_counter + 1L
  # infer generation and offspring id from eval count and lambda
  gen <- 1L + floor((.eval_counter - 1L) / lambda)
  off <- 1L + ((.eval_counter - 1L) %% lambda)
  # build row
  row <- tibble::tibble(
    eval_id = .eval_counter,
    iter = gen,
    offspring_id = off,
    rmse = rmse,
    time = Sys.time()
  ) %>% dplyr::bind_cols(tibble::as_tibble_row(stats::setNames(as.list(as.numeric(x)), param_names)))
  offspring_log <<- dplyr::bind_rows(offspring_log, row)
  rmse
}

# -----------------
# SMOOF WRAPPER with bounds
# -----------------
# Specific Smoof Wrapper around wrapped ojective function
# Necessary to set the objective function (compute_rmse) into a readable format for cmaesr, ie, cmaesr needs a smoof object
dim_par <- length(start_vals)
obj_logged <- smoof::makeSingleObjectiveFunction(
  name = "ForEdgeClim_RMSE_logged",
  fn   = function(x) wrapped_obj_fn(x), # function to be optimized
  has.simple.signature = TRUE, # fn accepts a simple vector as input, not a list
  par.set = ParamHelpers::makeNumericParamSet(len = dim_par, lower = lower_bounds, upper = upper_bounds), # number of params and search space (lower/upper bounds)
  minimize = TRUE # minimize RMSE, direction of optimalization
)


# -----------------
# DEFINING LOGGING STRUCTURES
# -----------------
# offspring log: one row per objective evaluation (i.e. per offspring)
offspring_log <- tibble::tibble(
  eval_id = integer(),
  iter    = integer(),           # generation index (1-based)
  offspring_id = integer(),      # lambda within generation
  rmse    = numeric(),
  time    = as.POSIXct(NA, tz = "UTC")
) %>%
  dplyr::bind_cols(tibble::as_tibble_row(stats::setNames(as.list(rep(NA_real_, dim_par)), param_names))) %>%
  dplyr::slice(0) # no data yet

# generation log (to be derived post-run)
generation_log <- tibble::tibble(
  iter = integer(),
  best_value = numeric(),
  sigma_proxy = numeric(),
  spread_max = numeric(),
  spread_mean = numeric(),
  mean_step = numeric()
)

# -----------------
# RUN CMA-ES
# -----------------
# Choose lambda if not provided
if (is.null(lambda)) lambda <- 4 + floor(3 * log(dim_par))

# initial sigma heuristic (initial global step size)
init_sigma <- 0.3 * mean(upper_bounds - lower_bounds)

# Build a robust stop list
# Stop after max_it iterations
stop_list <- list(cmaesr::stopOnMaxIters(max_it))
# Add default stopping conditions if available in your cmaesr version
# These include e.g. too small step sizes, no variance in population any more...
if ("getDefaultStoppingConditions" %in% getNamespaceExports("cmaesr")) {
  stop_list <- c(stop_list, cmaesr::getDefaultStoppingConditions())
}
# Target fitness as additional criterion if supported
# Stop if best RMSE =< stop_fitness
if ("stopOnOptFitness" %in% getNamespaceExports("cmaesr")) {
  stop_list <- c(stop_list, list(cmaesr::stopOnOptFitness(stop_fitness)))
}

# make results reproducible
set.seed(1)
res <- cmaesr::cmaes(
  obj_logged, # Smoof object with logging
  monitor = cmaesr::makeSimpleMonitor(), # simple progressing output in console
  control = list(
    lambda   = lambda,
    sigma    = init_sigma,
    stop.ons = stop_list
  )
)

saveRDS(res, file = results_output_file)

# -----------------
# DERIVE GENERATION-LEVEL LOG FROM OFFSPRING LOG
# -----------------

# Function to get spread values to investigate convergence
# M is a matrix of dimension lambda x d (number of params) with the offspring params for several iterations.
# Convergence would lead to small sigma_proxy, spread_max and spread_mean values.
get_spreads <- function(M) {
  sds <- apply(M, 2, stats::sd) # get sd for each param (= each column)
  # Covariance matrix representing the multivariate spread (including correlations between params)
  # The eigenvalues of this matrix are the variances along the main axes (principal components, PC).
  # The diagonal of a covariance matrix contains the variance of 1 param, the other elements the correlations between 2 params.
  # CMA-ES learns correlations between parameters: the algorithm eg detects improvements when changing combinations of parameters,
  # then the matrix is adapted, meaning the search directions (= PCs) are updated. Main search is in the direction of the main PC,
  # since it has the highest eigenvalue.
  # There are as many searching directions as parameters.
  covM <- stats::cov(M)
  ev <- tryCatch(eigen(covM, symmetric = TRUE)$values, error = function(e) rep(NA_real_, ncol(M)))
  list(
    sigma_proxy = mean(sds, na.rm = TRUE), # mean of param sds is a proxy for the global step size (= distribution, spread)
    spread_max  = sqrt(max(ev, na.rm = TRUE)), # largest deviation (sd) along main PC = spread in main searching direction
    spread_mean = sqrt(mean(ev, na.rm = TRUE)) # mean deviation in the PC-space = mean spread over all PCs/searching directions = general value for spread
  )
}

param_mat_cols <- param_names

generation_log <- offspring_log %>%
  dplyr::group_by(iter) %>%
  dplyr::summarise(
    best_value = min(rmse, na.rm = TRUE), # get lowest RMSE within each generation
    params = list(as.matrix(dplyr::pick(dplyr::all_of(param_mat_cols)))), # make a list with for every generation a lambda x d matrix with paramvalues
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    spreads = purrr::map(params, get_spreads), # get spread values per generation
    sigma_proxy = purrr::map_dbl(spreads, "sigma_proxy"), # get sigma_proxy from spreads
    spread_max  = purrr::map_dbl(spreads, "spread_max"),
    spread_mean = purrr::map_dbl(spreads, "spread_mean"),
    mean_param  = purrr::map(params, ~ colMeans(.x, na.rm = TRUE)) # mean param value per generation
  ) %>%
  # Get change in population mean between generations
  # This measures how strong the population center shifts from one generation to the next.
  # This change is defined as the Euclidean distance between consecutive generation means.
  # A low value means convergence.
  dplyr::mutate(mean_step = c(NA_real_,
                              purrr::map2_dbl(mean_param[-length(mean_param)], mean_param[-1],
                                              ~ sqrt(sum((.y - .x)^2))))) %>%
  # Keep columns of interest
  dplyr::select(iter, best_value, sigma_proxy, spread_max, spread_mean, mean_step)

# -----------------
# SAVE LOGS
# -----------------
readr::write_csv(offspring_log, logfile_offspring)
readr::write_csv(generation_log, logfile_gener)

plan(sequential)

# Console summaries
cat("\nOptimised parameters (final):\n")
print(res$best.param)
cat("\nMinimal RMSE observed:\n")
print(res$best.fitness)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f h\n', as.numeric(end_analysis - start_analysis, units = 'hours')))

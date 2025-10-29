###############################################################################
# This is a script to perform a Sobol sensitivity analysis on the ForEdgeClim
# model. Sobol indices per metric of interest are calculated and saved as well.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(ForEdgeClim)
library(dplyr)
library(ggplot2)
library(scales)
library(lhs)
library(future.apply)
library(future)
library(tidyr)
library(forcats)
library(purrr)
library(ggridges)
library(sensitivity)
library(ggtext)
library(progressr)

# Parallelisation on HPC (supercomputer)
handlers("progress")
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

#########
# INPUT #
#########

datetime <- as.POSIXct("2025-04-30 12:00:00", tz = "UTC")
structure <- "Data/TLS_scaled_DTM_and_grid_April2025.rds"
structure_PAI_scaling_factor <- 2.24 / 6.18
n <- 400
n_boot <- 1000

output_path <- "Output/sensitivity_analysis/"

# Parameter ranges
param_ranges <- data.frame(
  parameter = c(
    "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
    "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h",
    "h", "g_macro", "infl_macro", "infl_soil", "infl_forest", "g_forest", "p_ground", "g_soil", "k_soil"
  ),
  min = c(
    0.3, 0.2, 0.43, 0.6, 0.5, 0.08, 0.5, 0.3, 0.1,
    0.94, 0.3, 0.01, 0.2, 0.04, 0.2, 0.01,
    0, 10, 5, 0, 0, 5, 0.1, 5, 0.25
  ),
  max = c(
    0.35, 0.45, 0.61, 0.95, 2, 0.18, 0.95, 2, 0.2,
    0.99, 0.35, 0.06, 0.4, 0.07, 0.4, 0.06,
    20, 40, 60, 10, 10, 20, 0.35, 15, 2.2
  ),
  stringsAsFactors = FALSE
)

param_names <- param_ranges$parameter
k <- nrow(param_ranges)

formatted_datetime <- format(datetime, "%Hh_%d%m%Y")

output_name <- paste0(n, "samples_", k, "parameters_", formatted_datetime)

##################
# SOBOL SAMPLING #
##################

# Sobol design
set.seed(42)
X1 <- randomLHS(n, k)
X2 <- randomLHS(n, k)
colnames(X1) <- param_names
colnames(X2) <- param_names
sobol_design <- soboljansen(NULL, X1 = X1, X2 = X2, nboot = n_boot, conf = 0.95)

# Scaling
scale_sobol <- function(X) {
  for (i in seq_len(k)) {
    X[, i] <- param_ranges$min[i] + X[, i] * (param_ranges$max[i] - param_ranges$min[i])
  }
  colnames(X) <- param_ranges$parameter
  as.data.frame(X)
}
scaled_input <- scale_sobol(sobol_design$X)

#############################################
# MODEL DEFINITION WITH METRIC CALCULATIONS #
#############################################

model_function <- function(param_values) {
  create_input_drivers()
  create_physical_constants()
  create_model_parameters()
  for (i in seq_along(param_names)) {
    assign(param_names[i], param_values[i], envir = .GlobalEnv)
  }

  import_RMI_observations(datetime)
  if(is_empty(F_sky_lw)){
    assign("F_sky_lw", sigma_SB * 0.75 * macro_temp^4, envir = .GlobalEnv)
  }
  import_pyr_observations(datetime)
  import_soil_temperature(datetime)

  voxel_TLS <- readRDS(structure)
  voxel_TLS$grid$density = voxel_TLS$grid$density * structure_PAI_scaling_factor
  res <- run_foredgeclim(voxel_TLS$grid, datetime)

  # horizontal gradient
  sel_h <- res$micro_grid$z == 1 & res$micro_grid$y == 15
  xs_h <- res$micro_grid$x[sel_h]
  ts_h <- res$air_temperature[sel_h] - 273.15

  avT_h <- mean(ts_h, na.rm = TRUE)
  gradT_h <- (ts_h[xs_h == max(xs_h)][1] - ts_h[xs_h == min(xs_h)][1]) / (max(xs_h) - min(xs_h))
  SDT_h <- sd(ts_h, na.rm = TRUE)

  # vertical gradient
  sel_v <- res$micro_grid$x == 68 & res$micro_grid$y == 15
  zs_v <- res$micro_grid$z[sel_v]
  ts_v <- res$air_temperature[sel_v] - 273.15

  avT_v <- mean(ts_v, na.rm = TRUE)
  gradT_v <- (ts_v[zs_v == max(zs_v)][1] - ts_v[zs_v == min(zs_v)][1]) / (max(zs_v) - min(zs_v))
  SDT_v <- sd(ts_v, na.rm = TRUE)

  return(c(avT_h = avT_h, gradT_h = gradT_h, SDT_h = SDT_h, avT_v = avT_v, gradT_v = gradT_v, SDT_v = SDT_v))
}


##############
# SOBOL RUNS #
##############

sobol_outputs_all <- with_progress({
  p <- progressor(steps = nrow(scaled_input))
  future_lapply(seq_len(nrow(scaled_input)), function(i) {
    metrics <- model_function(unlist(scaled_input[i, ]))
    p()
    return(metrics)
  }, future.seed = TRUE)
})

plan(sequential)

# Metric extraction
# horizontal gradient
sobol_outputs_metric1 <- map_dbl(sobol_outputs_all, "avT_h")
sobol_outputs_metric2 <- map_dbl(sobol_outputs_all, "gradT_h")
sobol_outputs_metric3 <- map_dbl(sobol_outputs_all, "SDT_h")

sobol_metric1 <- tell(sobol_design, sobol_outputs_metric1)
sobol_metric2 <- tell(sobol_design, sobol_outputs_metric2)
sobol_metric3 <- tell(sobol_design, sobol_outputs_metric3)

# vertical gradient
sobol_outputs_metric4 <- map_dbl(sobol_outputs_all, "avT_v")
sobol_outputs_metric5 <- map_dbl(sobol_outputs_all, "gradT_v")
sobol_outputs_metric6 <- map_dbl(sobol_outputs_all, "SDT_v")

sobol_metric4 <- tell(sobol_design, sobol_outputs_metric4)
sobol_metric5 <- tell(sobol_design, sobol_outputs_metric5)
sobol_metric6 <- tell(sobol_design, sobol_outputs_metric6)

##########
# SAVING #
##########

# horizontal gradient
saveRDS(sobol_metric1, file = file.path(output_path, paste0(output_name, "_avT_h.rds")))
saveRDS(sobol_metric2, file = file.path(output_path, paste0(output_name, "_gradT_h.rds")))
saveRDS(sobol_metric3, file = file.path(output_path, paste0(output_name, "_SDT_h.rds")))

# vertical gradient
saveRDS(sobol_metric4, file = file.path(output_path, paste0(output_name, "_avT_v.rds")))
saveRDS(sobol_metric5, file = file.path(output_path, paste0(output_name, "_gradT_v.rds")))
saveRDS(sobol_metric6, file = file.path(output_path, paste0(output_name, "_SDT_v.rds")))


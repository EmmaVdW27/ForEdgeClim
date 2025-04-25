###############################################################################
# This is a script to run a OAT (one at a time) global (or local) sensitivity analysis on the ForEdgeClim microclimate model.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(ForEdgeClim)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(future.apply)

#########
# INPUT #
#########

# PLAN PARALLELISATION
cores = 10 # 20 available (best max 10)
plan(multisession, workers = cores)

start_analysis = Sys.time()

datetime = "2023-07-08 12:00:00"

params = "non-RTM" # "RTM" or "non-RTM" or "input_drivers"

####################################################
# MODEL PARAMETERS TO RUN SENSITIVITY ANALYSIS FOR #
####################################################

# Define ranges for parameters
if(params == "non-RTM"){
  # non-RTM parameters (4.2 min)
  param_ranges <- data.frame(
    parameter = c("g_macro", "g_forest", "g_soil", "infl_macro", "infl_soil", "infl_forest", "p_ground", "k_soil", "h"),
    min = c(5, 5, 5, 10, 1, 1, 0.1, 0.2, 5),
    max = c(15, 15, 15, 60, 10, 15, 0.35, 2, 15),
    steps = rep(5,9)
  )
} else if(params == "RTM"){
  # RTM parameters (6.7 min)
  param_ranges <- data.frame(
    parameter = c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
                  "e_atm", "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h"),
    min = c(0.15, 0.25, 0.5, 0.2, 0.6, 0.1, 0.15, 0.2, 0,
            0.6, 0.92, 0, 0, 0.1, 0, 0.05, 0),
    max = c(0.5, 0.5, 0.95, 0.5, 0.9, 0.3, 0.3, 0.4, 0.15,
            0.9, 0.98, 0.2, 0.1, 0.3, 0.05, 0.2, 0.05),
    steps = rep(5,17)
  )
} else {
  # input driver variables (2.1 min)
  param_ranges <- data.frame(
    parameter = c("T_soil_deep", "F_sky_dir_init", "F_sky_diff_init", "macro_temp"),
    min = c(12+273.15, 0, 0, 12+273.15),
    max = c(20+273.15, 1000, 400, 32+273.15),
    steps = rep(5,4)
  )
}

# Define title and output path for plot
if(params == "non-RTM"){
  title = "Non-RTM parameters"
  output_plot = 'Output/sensitivity_analysis/OAT_sensitivity_analysis_non-RTM.png'
} else if(params == "RTM"){
  title = "RTM parameters"
  output_plot = 'Output/sensitivity_analysis/OAT_sensitivity_analysis_RTM.png'
} else {
  title = "Input driver variables"
  output_plot = 'Output/sensitivity_analysis/OAT_sensitivity_analysis_drivers.png'
}

# Generate values for the parameters based on the ranges
params_to_test <- lapply(1:nrow(param_ranges), function(i) {
  seq(param_ranges$min[i], param_ranges$max[i], length.out = param_ranges$steps[i])
})
names(params_to_test) <- param_ranges$parameter


####################################
# RUN THE SIMULATION PER PARAMETER #
####################################

# list of all (parameter, value)-combinations
param_grid <- param_ranges %>%
  rowwise() %>%
  mutate(value = list(seq(min, max, length.out = steps))) %>%
  unnest(value) %>%
  select(parameter, value)

param_list <- split(param_grid, seq_len(nrow(param_grid)))

# parallel iteration with future_lapply
results_list <- future_lapply(param_list, function(df_row) {
  param <- df_row$parameter
  val   <- df_row$value

  # import & init
  current_datetime <- as.POSIXct(datetime, tz = "UTC")

  create_input_drivers(current_datetime)
  create_model_parameters()
  create_physical_constants()
  if(params != "input_drivers") {assign(param, val, envir = .GlobalEnv)}

  import_DTS_observations()
  import_RMI_observations()
  import_pyr_observations()
  import_soil_temperature()
  if(params == "input_drivers") {assign(param, val, envir = .GlobalEnv)}

  # run model
  voxel_TLS <- readRDS(TLS_filtered_file)
  res <- run_foredgeclim(voxel_TLS$grid)

  metric <- mean(res$air_temperature[res$micro_grid$z == 1], na.rm = TRUE) - 273.15 # K to °C

  # result data.frame
  data.frame(parameter = param, value = val, metric = metric)
}, future.seed = TRUE)

# bind all results into one data.frame
results <- do.call(rbind, results_list)

#########
# PLOTS #
#########

results_scaled <- results %>%
  group_by(parameter) %>%
  mutate(
    scaled_value = (value - min(value)) / (max(value) - min(value))
  ) %>%
  ungroup()

# add temperature in °C
results_labeled <- results_scaled %>%
  mutate(temp_C = round(metric, 2))

# get last point per parameter (for labeling)
label_points <- results_labeled %>%
  group_by(parameter) %>%
  filter(scaled_value == max(scaled_value)) %>%
  ungroup()

if(params == "non-RTM"){
  # parameters curve specifics
  my_colors <- c(
    "g_macro" = "steelblue",
    "infl_macro" = "steelblue",
    "g_forest" = "forestgreen",
    "infl_forest" = "forestgreen",
    "h" = "forestgreen",
    "g_soil" = "saddlebrown",
    "infl_soil" = "saddlebrown",
    "k_soil" = "saddlebrown",
    "p_ground" = "saddlebrown"
  )

  # annotate specifics
  y_base   <- max(results_labeled$temp_C)
  y_offset <- c(1.2, 0.6, 0)

  labels <- c("Macro environment",
              "Forest environment",
              "Ground environment")
  colors <- c("steelblue",
              "forestgreen",
              "saddlebrown")
} else if(params == "RTM"){
  my_colors <- c("betad" = "goldenrod", "beta0"  = "goldenrod", "omega"  = "goldenrod", "Kd_v"  = "goldenrod",
                 "Kb_v"  = "goldenrod", "omega_g_v"  = "goldenrod", "Kd_h"  = "goldenrod", "Kb_h"  = "goldenrod",
                 "omega_g_h"  = "goldenrod",
                 "e_atm" = "firebrick", "e_forest" = "firebrick", "beta_lw" = "firebrick", "omega_lw" = "firebrick",
                 "Kd_lw_v" = "firebrick", "omega_g_lw_v" = "firebrick", "Kd_lw_h" = "firebrick", "omega_g_lw_h" = "firebrick")

  y_base   <- max(results_labeled$temp_C)
  y_offset <- c(0.05, 0)

  labels <- c("Shortwave RTM",
              "Longwave RTM")
  colors <- c("goldenrod",
              "firebrick")
} else {
  my_colors <- c("T_soil_deep" = "saddlebrown", "F_sky_dir_init" = "goldenrod",
                 "F_sky_diff_init" = "goldenrod", "macro_temp" = "steelblue")

  y_base   <- max(results_labeled$temp_C)
  y_offset <- c(1.2, 0.6, 0)

  labels <- c("Ground environment",
              "Sun environment",
              "Macro environment")
  colors <- c("saddlebrown",
              "goldenrod",
              "steelblue")
}

# Plot
sensitivity_plot <- ggplot(results_labeled, aes(x = scaled_value, y = temp_C, colour = parameter)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_text_repel(
    data = label_points,
    aes(label = parameter),
    nudge_x = 0.05,
    direction = "y",
    hjust = 0,
    show.legend = FALSE,
    size = 4,
    max.overlaps = Inf
  ) +
  labs(title = title,
       x = "Scaled parameter value", y = "Air temperature (°C)",
       colour = "Parameter") +
  theme_bw() +
  theme(
    legend.position      = "none",
    plot.title = element_text(face = "bold", size = 20), # title
    axis.title.x         = element_text(size = 16),  # axis name x
    axis.title.y         = element_text(size = 16),  # axis name y
    axis.text.x          = element_text(size = 14),  # tick labels x
    axis.text.y          = element_text(size = 14)   # tick labels y
  ) +
  scale_colour_manual(values = my_colors) +
  annotate("text",
           x       = rep(0.5, length(labels)),
           y       = y_base + y_offset,
           label   = labels,
           colour  = colors,
           size    = 10,
           fontface= "bold",
           hjust   = 0.5)

print(sensitivity_plot)



ggsave(paste0(output_plot), plot = sensitivity_plot, width = 10, height = 6, dpi = 300)


# reset PLAN
plan(sequential)

end_analysis = Sys.time()
print(paste0('Total running time sensitivity analysis = ', round(as.numeric(end_analysis - start_analysis, units = "mins"), 2), ' min'))

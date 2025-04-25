###############################################################################
# This is a script to run an AAT (all at a time) global sensitivity analysis on the ForEdgeClim microclimate model.
# By this, we mean a variance decomposition to see how input parameter variance weights on the
# variance of the output metric.
# It uses the SOBOL technique to analyse the contributions of each parameter.
# The input parameter distributions are currently uniform distributions defined with a min and max value.
# We use Latin Hypercube Sampling to perform a stratified sampling from these distributions.
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
library(ggtext)

#########
# INPUT #
#########

# progress bar
handlers("progress")

# PLAN PARALLELISATION
cores = 10 # 20 available (best max 10)
plan(multisession, workers = cores)


start_analysis = Sys.time()

datetime = as.POSIXct("2023-07-08 12:00:00", tz = "UTC")

n = 400 # no of samples. Ideal around 250/500

n_boot = 1000 # no of bootstrap-resamples to build confidence intervals (necessary for CI)
           # (ideal around 500/1000, 200 for decent testing); > 0 for output plot

# SOBOL runtime is around (k + 2) * n * ~30s / cores with k = no of parameters
# for n = 5 and cores = 10 and nboot = 0: 4.5 min
# for n = 5 and cores = 10 and nboot = 1000: 8.15 min
# for n = 50 and cores = 13: 46.1 min
# for n = 100 and cores = 10: 102.3 min
# for n = 100 and cores = 10 and n_boot = 200: 91.7 min
# for n = 250 and cores = 10 and n_boot = 200: 230.6 min
# for n = 250 and cores = 10 and n_boot = 1000: 232.2 min
# for n = 400 and cores = 10 and n_boot = 1000: 400.8 min

###############################################################
# MODEL PARAMETERS TO RUN AAT GLOBAL SENSITIVITY ANALYSIS FOR #
###############################################################

# Define parameters & ranges
param_ranges <- data.frame(
  parameter = c("g_macro", "g_forest", "g_soil", "infl_macro", "infl_forest", "infl_soil", "p_ground", "k_soil", "h"),
  min = c(5, 5, 5, 10, 1, 1, 0.1, 0.2, 5), # real possible values in temperate forests
  max = c(15, 15, 15,  60, 10, 15, 0.35, 2, 15),
  # min = c(-100, -100, -100, -100, -100, -100, 0, -100, -100),
  # max = c(100, 100, 100,  100, 100, 100, 1, 100, 100),
  stringsAsFactors = FALSE
)

param_names <- param_ranges$parameter
k <- nrow(param_ranges)   # no of parameters


#############
# MODEL RUN #
#############

# Function to calculate output metric on which you want to run the analysis for.
# Here we define a sum of the averaged temperature at the position of the transect line
# and the grandient along this transect.
compute_metric <- function(res) {
  # 1) subset z=1, y=15 (transect line position)
  sel <- res$micro_grid$z == 1 & res$micro_grid$y == 15
  xs  <- res$micro_grid$x[sel]
  ts  <- res$air_temperature[sel] - 273.15 # K to °C

  # 2) mean air temperature along transect line
  mean_t <- mean(ts, na.rm = TRUE)

  # 3) core to edge gradient along transect line
  x_min <- min(xs, na.rm = TRUE)
  x_max <- max(xs, na.rm = TRUE)
  t_min <- ts[xs == x_min][1]   # value at x_min
  t_max <- ts[xs == x_max][1]   # value at x_max
  grad <- (t_max - t_min) / (x_max - x_min)

  # 4) combine
  metric <- mean_t + 5*grad

  return(metric)

}

# Model function that returns the metric
model_function <- function(param_values) {

    current_datetime <- datetime

    create_input_drivers(datetime)
    create_physical_constants()
    create_model_parameters()
    for(i in seq_along(param_names)) {
      assign(param_names[i], param_values[i], envir = .GlobalEnv) # overwrite parameters in global environment
    }

    import_DTS_observations()
    import_RMI_observations()
    import_pyr_observations()
    import_soil_temperature()


    voxel_TLS <- readRDS(TLS_filtered_file)
    res <- run_foredgeclim(voxel_TLS$grid)

    compute_metric(res)

}


##################
# SOBOL SAMPLING #
##################

X1 <- randomLHS(n, k)  # n × k matrix
X2 <- randomLHS(n, k)

# set column names so sobol() can use them
colnames(X1) <- param_names
colnames(X2) <- param_names

# make sobol design object
sobol_design <- soboljansen(NULL, X1 = X1, X2 = X2, nboot = n_boot, conf = 0.95)

# function to scale parameters
scale_sobol <- function(X) {
  rownames(X) <- NULL
  for (i in seq_len(k)) {
    X[, i] <- param_ranges$min[i] + X[, i] * (param_ranges$max[i] - param_ranges$min[i])
  }
  colnames(X) <- param_ranges$parameter
  as.data.frame(X, stringsAsFactors = FALSE)

}

# scale combined samples
sobol_input <- sobol_design$X # (k+2)*n * k columns
scaled_input <- scale_sobol(sobol_input)

# run model on all parameter combinations (like with MC)
# sobol_outputs <- future_apply(as.matrix(scaled_input), 1, model_function, future.seed = TRUE)
sobol_outputs <- with_progress({          # ➊ start progressbar
  p <- progressor(steps = nrow(scaled_input))  # ➋ total = no of rows

  future_apply(
    as.matrix(scaled_input),
    1,
    function(param_values) {
      result <- model_function(param_values)
      p()                           # ➌ +1 step on bar
      result
    },
    future.seed = TRUE
  )
})                                        # ➍ done


# no names for output-vector
sobol_outputs <- unname(sobol_outputs)

# calculate sobol indices/ tell results to sobol_design
sobol_design <- tell(sobol_design, sobol_outputs)


#####################
# RESULTS AND PLOTS #
#####################

# First order SOBOL index S: S_i = var(E[Y|X_i])/var(Y)
# Direct contribution of parameter X_i on the variance of Y, without taking interaction with other parameters into account.
# eg: S_i = 0.3 -> X_i explains 30% of the variance in Y.
# S < 0.005 can be a fixed parameter with not much influence on the output variance.
# Total order SOBOL index T: T_i = 1 - var(E[Y|X~i])/var(Y) with X~i all parameters except X_i
# This is the total contribution of X_i in the variance of Y, including all interaction effects with other parameters.
# eg: T_i = 0.5 -> X_i explains 50% from itself and in combination with others of the variance in Y.
# T_i - S_i gives the contribution of the interactions in which X_i participates.

# wide-table with consistant column names *_original, *_lower, *_upper
df_sobol <- tibble(
  parameter       = param_names,
  first_original  = sobol_design$S$original,
  first_lower     = sobol_design$S$`min. c.i.`,
  first_upper     = sobol_design$S$`max. c.i.`,
  total_original  = sobol_design$T$original,
  total_lower     = sobol_design$T$`min. c.i.`,
  total_upper     = sobol_design$T$`max. c.i.`
) %>%
  # set all negative values to 0, SOBOL indices cannot be negative, if that is the case it is due to the sampling strategy (LHS)
  mutate(
    across(starts_with("first_"), ~ pmax(0, .)),
    across(starts_with("total_"), ~ pmax(0, .))
  ) %>%
  # calculate interactions
  mutate(
    interactions_original = total_original - first_original,
    # no negative lower boundaries
    interactions_lower    = pmax(0, total_lower  - first_upper),
    interactions_upper    = total_upper     - first_lower,
    across(starts_with("interactions_"), ~ pmax(0, .)),
    # zorg dat total_upper minimaal first_upper is
    total_original = pmax(total_original, first_original)
  )

# pivot longer per (parameter, order) so you get original value + lower/upper
df_sobol_long <- df_sobol %>%
  pivot_longer(
    cols      = -parameter,
    names_to  = c("order", "stat"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from  = stat,
    values_from = value
  ) %>%
  mutate(
    order = factor(order,
                   levels = c("first", "total", "interactions"),
                   labels = c("First order", "Total order", "Interactions")),
    parameter = factor(parameter, levels = param_ranges$parameter)
  )

parameter_order <- rev(c("g_macro", "infl_macro", "g_forest", "infl_forest","h",
                         "g_soil", "infl_soil", "p_ground", "k_soil"))
param_colors <- c(
  "g_macro"     = "steelblue",
  "g_forest"    = "forestgreen",
  "g_soil"      = "saddlebrown",
  "infl_macro"  = "steelblue",
  "infl_forest" = "forestgreen",
  "infl_soil"   = "saddlebrown",
  "p_ground"    = "saddlebrown",
  "k_soil"      = "saddlebrown",
  "h"           = "forestgreen"
)

# make HTML‑labels with color
colored_labels <- parameter_order %>%
  setNames(., .) %>%
  imap_chr(~ paste0("<span style='color:", param_colors[.x], ";'>", .x, "</span>"))

df_sobol_long <- df_sobol_long %>%
  mutate(parameter = factor(parameter, levels = parameter_order))

# bar plot with error bars
sobol_indices = ggplot(df_sobol_long, aes(x = parameter, y = original, fill = order)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  coord_flip() +
  scale_x_discrete(labels = colored_labels) +
  scale_fill_manual(
    values = c(
      "First order"  = '#17becf',
      "Total order"  = "#bcbd22",
      "Interactions" = "orange"
    )
  ) +
  labs(
    title = "SOBOL indices with 95% CI",
    x     = "Parameter",
    y     = "Fraction of explained variance",
    fill  = "Order"
  ) +
  theme_bw() +
  theme(
    # render HTML in y‑labels
    plot.title  = element_text(face = "bold", size = 18),
    # as‑titels
    axis.title.x   = element_text(size = 18, face = "bold"),
    axis.title.y   = element_text(size = 18, face = "bold"),
    # tick‑labels
    axis.text.x    = element_text(size = 12),
    axis.text.y = element_markdown(size = 12),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )
print(sobol_indices)


ggsave('Output/sensitivity_analysis/global_sensitivity_analysis_SOBOL_indices.png', plot = sobol_indices, width = 6, height = 9, dpi = 300)


# reset PLAN
plan(sequential)
end_analysis = Sys.time()
print(paste0('Total running time AAT global sensitivity analysis = ', round(as.numeric(end_analysis - start_analysis, units = "mins"), 2), ' min'))

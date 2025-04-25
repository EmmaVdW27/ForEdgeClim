###############################################################################
# This is a script to run an uncertainty analysis on the ForEdgeClim microclimate model.
# It uses the Latin Hypercube sampling technique to extract stratified samples from input parameter distributions.
# These distributions are currently uniform distributions defined with a min and max value.
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
library(ggh4x)
library(ggforce)

#########
# INPUT #
#########

# PLAN PARALLELISATION
cores = 15 # 20 available (best max 15)
plan(multisession, workers = cores)


start_analysis = Sys.time()

datetime = as.POSIXct("2023-07-08 12:00:00", tz = "UTC")

n_samples = 1000 # number of parameter samples for MC (ideal around 1000)
                # MC simulation will take n_samples x model_run_time (~30s) / cores
                # 9min if n_samples = 100 and cores = 6
                # 7min if n_samples = 100 and cores = 15
                # 33min if n_samples = 500 and cores = 15
                # 81.6min if n_samples = 1000 and cores = 15


####################################################
# MODEL PARAMETERS TO RUN UNCERTAINTY ANALYSIS FOR #
####################################################

# Define parameters & ranges
param_ranges <- data.frame(
  parameter = c("g_macro", "g_forest", "g_soil", "infl_macro", "infl_forest", "infl_soil", "p_ground", "k_soil", "h"),
  min = c(5, 5, 5, 10, 1, 1, 0.1, 0.2, 5),
  max = c(15, 15, 15,  60, 10, 15, 0.35, 2, 15)
)

param_names <- param_ranges$parameter
k = nrow(param_ranges)

# generate parameter samples
set.seed(123)
n <- n_samples
# Latin Hypercube Sample (better spread over parameter domain wrt random sampling -> statified sampling)
# n = number of trials (model runs)
# k = number of parameters
param_samples <- randomLHS(n, k)

# scale matrix to real parameter values
scale_params <- function(X) {
  X_scaled <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  for (i in 1:ncol(X)) {
    X_scaled[, i] <- param_ranges$min[i] + X[, i] * (param_ranges$max[i] - param_ranges$min[i])
  }
  colnames(X_scaled) <- param_names
  return(X_scaled)
}

param_matrix <- scale_params(param_samples)


#############
# MODEL RUN #
#############

# Model function that returns mean temperature at 1m height
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

    mean(res$air_temperature[res$micro_grid$z == 1], na.rm = TRUE) - 273.15 # K to °C
}


mc_outputs <- future_apply(param_matrix, 1, model_function, future.seed = TRUE)


############
# PLOTTING #
############

# input + output in 1 dataframe
df <- data.frame(g_macro = param_matrix[, "g_macro"],
                 g_forest = param_matrix[, "g_forest"],
                 g_soil = param_matrix[, "g_soil"],
                 infl_macro = param_matrix[, "infl_macro"],
                 infl_forest = param_matrix[, "infl_forest"],
                 infl_soil = param_matrix[, "infl_soil"],
                 p_ground = param_matrix[, "p_ground"],
                 k_soil = param_matrix[, "k_soil"],
                 h = param_matrix[, "h"],
                 output = mc_outputs)

# # 1) distribution model output
# print(summary(mc_outputs))
# print(quantile(mc_outputs, probs = c(0.05, 0.5, 0.95))) # uncertainty intervals on model output
#
# # density plot
# print(ggplot(df, aes(x = output)) +
#   geom_density(fill = "lightblue") +
#   labs(title = "Probability distribution of mean air temperature at 1m height", x = "Temperature (°C)", y = "Probability") +
#   theme_bw()
#   )


# calculations
mu    <- mean(mc_outputs)
sdv   <- sd(mc_outputs)
dens  <- density(mc_outputs)
dens_df <- data.frame(x = dens$x, y = dens$y)

# Y-values on curve for line segments
y_mu        <- approx(dens_df$x, dens_df$y, xout = mu)$y
y_1sd_low   <- approx(dens_df$x, dens_df$y, xout = mu - sdv)$y
y_1sd_high  <- approx(dens_df$x, dens_df$y, xout = mu + sdv)$y
y_2sd_low   <- approx(dens_df$x, dens_df$y, xout = mu - 2*sdv)$y
y_2sd_high  <- approx(dens_df$x, dens_df$y, xout = mu + 2*sdv)$y

# set small x‐offset and y-positions for labels
offset  <- sdv * 0.1
y_label <- max(dens_df$y) * 0.02  # 2% above x‐axis

output_dist = ggplot(dens_df, aes(x = x, y = y)) +
  # density
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = y), fill = "lightblue") +

  # segments for μ and σ‐lines
  annotate("segment",
           x = mu, xend = mu, y = 0, yend = y_mu,
           colour = "blue", size = 1) +
  annotate("segment",
           x = mu - sdv, xend = mu - sdv, y = 0, yend = y_1sd_low,
           linetype = "dashed", colour = "blue", size = 0.8) +
  annotate("segment",
           x = mu + sdv, xend = mu + sdv, y = 0, yend = y_1sd_high,
           linetype = "dashed", colour = "blue", size = 0.8) +
  annotate("segment",
           x = mu - 2*sdv, xend = mu - 2*sdv, y = 0, yend = y_2sd_low,
           linetype = "dotdash", colour = "blue", size = 0.8) +
  annotate("segment",
           x = mu + 2*sdv, xend = mu + 2*sdv, y = 0, yend = y_2sd_high,
           linetype = "dotdash", colour = "blue", size = 0.8) +

  # labels right of vertical lines, horizontally above y=0
  annotate("text",
           x = mu + offset, y = y_label,
           label = "mu", parse = TRUE,
           hjust = 0, vjust = 0, size = 3, colour = "blue") +
  annotate("text",
           x = mu - sdv + offset, y = y_label,
           label = "mu - sigma", parse = TRUE,
           hjust = 0, vjust = 0, size = 3) +
  annotate("text",
           x = mu + sdv + offset, y = y_label,
           label = "mu + sigma", parse = TRUE,
           hjust = 0, vjust = 0, size = 3) +
  annotate("text",
           x = mu - 2*sdv + offset, y = y_label,
           label = "mu - 2*sigma", parse = TRUE,
           hjust = 0, vjust = 0, size = 3) +
  annotate("text",
           x = mu + 2*sdv + offset, y = y_label,
           label = "mu + 2*sigma", parse = TRUE,
           hjust = 0, vjust = 0, size = 3) +

  labs(
    title = "Probability distribution of mean air temperature at 1m height",
    x     = "Temperature (°C)",
    y     = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(face = "bold", size = 18),
    # as‑titels
    axis.title.x   = element_text(size = 18, face = "bold"),
    axis.title.y   = element_text(size = 18, face = "bold"),
    # tick‑labels
    axis.text.x    = element_text(size = 12),
    axis.text.y = element_markdown(size = 12)
    )
print(output_dist)
ggsave('Output/uncertainty_analysis/uncertainty_analysis_output_dist.png', plot = output_dist, width = 10, height = 6, dpi = 300)





# 2) correlations between input and output (~informal global sensitivity analysis)
# This gives indication of sensitivity without formal method like SOBOL
# Eg linear relationship or other ...

df_scaled <- df %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

# parameter order for facet plot
parameter_order <- c("g_macro", "infl_macro", "g_forest", "infl_forest", "h",
                     "g_soil", "infl_soil", "p_ground", "k_soil")

# long format (excl. 'output' as parameter)
df_long <- df_scaled %>%
  pivot_longer(cols = -output, names_to = "parameter", values_to = "value")%>%
  mutate(parameter = factor(parameter, levels = parameter_order))


# calculate correlation per parameter
cor_info <- df_long %>%
  group_by(parameter) %>%
  summarise(
    cor_val = cor(value, output, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label_text = paste0(
      "cor = ", round(cor_val, 2)
    )
  )

# add labels to df_long (join via parameter)
df_long <- df_long %>%
  left_join(cor_info, by = "parameter")

# colors for plot
param_colors <- c(
  "g_macro"     = "steelblue",
  "infl_macro"  = "steelblue",
  "g_forest"    = "forestgreen",
  "infl_forest" = "forestgreen",
  "h"           = "forestgreen",
  "g_soil"      = "saddlebrown",
  "infl_soil"   = "saddlebrown",
  "p_ground"    = "saddlebrown",
  "k_soil"      = "saddlebrown"

)

correlations <- ggplot(df_long, aes(x = value, y = output)) +
  geom_point(alpha = 0.6, colour = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "orange") +
  facet_wrap2(
    ~ parameter,
    scales = "free",
    axes = "margins",
    remove_labels = "all",
    strip = strip_themed(
      # color per level
      background_x = elem_list_rect(
        fill = param_colors
      ),
      # text color per level, eg white for contrast
      text_x = elem_list_text(
        colour = ifelse(
          grDevices::col2rgb(param_colors)[1, ] +
            grDevices::col2rgb(param_colors)[2, ] +
            grDevices::col2rgb(param_colors)[3, ] < 380,
          "white", "white" # first color is text color on dark background, second one on light background
        )
      )
    )
  ) +
  geom_label(
    aes(x = -Inf, y = -Inf, label = label_text),
    hjust = -0.05, vjust = -0.5,
    fill   = "white", colour = "darkblue", family = "mono",
    size   = 8
  ) +facetted_pos_scales(
    y = list(
      COL != 1 ~ scale_y_continuous(guide = "none")  # only on 1st column y-scale
    ),
    x = list(
      ROW != 3 ~ scale_x_continuous(guide = "none")# only on last row x-scale
    )
  ) +
  labs(
    title = element_blank(),
    x     = "Parameter value (z-score)",
    y     = "Air temperature (z-score)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    # strip‑titels (facet headers)
    strip.text.x   = element_text(size = 16, face = "bold"),
    # as‑titels
    axis.title.x   = element_text(size = 18, face = "bold"),
    axis.title.y   = element_text(size = 18, face = "bold"),
    # tick‑labels
    axis.text.x    = element_text(size = 12),
    axis.text.y    = element_text(size = 12)
  )

#print(correlations)
ggsave('Output/uncertainty_analysis/uncertainty_analysis_correlations.png', plot = correlations, width = 10, height = 10, dpi = 300)



# reset PLAN
plan(sequential)
end_analysis = Sys.time()
print(paste0('Total running time LHS uncertainty analysis = ', round(as.numeric(end_analysis - start_analysis, units = "mins"), 2), ' min'))

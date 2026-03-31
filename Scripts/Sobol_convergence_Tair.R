###############################################################################
# Sobol convergence analysis (200 vs 400 LHS samples)
# Uses ALL conditions (seasons × moments × metrics)
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(forcats)
library(glue)
library(patchwork)

#########
# INPUT #
#########

input_path <- "Output/sensitivity_analysis/Sobol_QoI/data/"
output_path_plots <- "Output/sensitivity_analysis/Sobol_QoI/plots_output/"

direction <- 'h'  # 'h' or 'v'

##########
# FILES
##########

files_info <- tribble(
  ~file, ~label,

  # winter
  glue("400samples_25parameters_01h_13012025_avTa_{direction}.rds"), "Wi_Ni_≈T",
  glue("400samples_25parameters_01h_13012025_SDTa_{direction}.rds"), "Wi_Ni_σT",
  glue("400samples_25parameters_01h_13012025_gradTa_{direction}.rds"), "Wi_Ni_∇T",

  glue("400samples_25parameters_08h_13012025_avTa_{direction}.rds"), "Wi_Mo_≈T",
  glue("400samples_25parameters_08h_13012025_SDTa_{direction}.rds"), "Wi_Mo_σT",
  glue("400samples_25parameters_08h_13012025_gradTa_{direction}.rds"), "Wi_Mo_∇T",

  glue("400samples_25parameters_12h_13012025_avTa_{direction}.rds"), "Wi_No_≈T",
  glue("400samples_25parameters_12h_13012025_SDTa_{direction}.rds"), "Wi_No_σT",
  glue("400samples_25parameters_12h_13012025_gradTa_{direction}.rds"), "Wi_No_∇T",

  # spring
  glue("400samples_25parameters_01h_30042025_avTa_{direction}.rds"), "Sp_Ni_≈T",
  glue("400samples_25parameters_01h_30042025_SDTa_{direction}.rds"), "Sp_Ni_σT",
  glue("400samples_25parameters_01h_30042025_gradTa_{direction}.rds"), "Sp_Ni_∇T",

  glue("400samples_25parameters_05h_30042025_avTa_{direction}.rds"), "Sp_Mo_≈T",
  glue("400samples_25parameters_05h_30042025_SDTa_{direction}.rds"), "Sp_Mo_σT",
  glue("400samples_25parameters_05h_30042025_gradTa_{direction}.rds"), "Sp_Mo_∇T",

  glue("400samples_25parameters_12h_30042025_avTa_{direction}.rds"), "Sp_No_≈T",
  glue("400samples_25parameters_12h_30042025_SDTa_{direction}.rds"), "Sp_No_σT",
  glue("400samples_25parameters_12h_30042025_gradTa_{direction}.rds"), "Sp_No_∇T",

  # summer
  glue("400samples_25parameters_01h_07072023_avTa_{direction}.rds"), "Su_Ni_≈T",
  glue("400samples_25parameters_01h_07072023_SDTa_{direction}.rds"), "Su_Ni_σT",
  glue("400samples_25parameters_01h_07072023_gradTa_{direction}.rds"), "Su_Ni_∇T",

  glue("400samples_25parameters_05h_07072023_avTa_{direction}.rds"), "Su_Mo_≈T",
  glue("400samples_25parameters_05h_07072023_SDTa_{direction}.rds"), "Su_Mo_σT",
  glue("400samples_25parameters_05h_07072023_gradTa_{direction}.rds"), "Su_Mo_∇T",

  glue("400samples_25parameters_12h_07072023_avTa_{direction}.rds"), "Su_No_≈T",
  glue("400samples_25parameters_12h_07072023_SDTa_{direction}.rds"), "Su_No_σT",
  glue("400samples_25parameters_12h_07072023_gradTa_{direction}.rds"), "Su_No_∇T",

  # autumn
  glue("400samples_25parameters_01h_01102023_avTa_{direction}.rds"), "Au_Ni_≈T",
  glue("400samples_25parameters_01h_01102023_SDTa_{direction}.rds"), "Au_Ni_σT",
  glue("400samples_25parameters_01h_01102023_gradTa_{direction}.rds"), "Au_Ni_∇T",

  glue("400samples_25parameters_06h_01102023_avTa_{direction}.rds"), "Au_Mo_≈T",
  glue("400samples_25parameters_06h_01102023_SDTa_{direction}.rds"), "Au_Mo_σT",
  glue("400samples_25parameters_06h_01102023_gradTa_{direction}.rds"), "Au_Mo_∇T",

  glue("400samples_25parameters_12h_01102023_avTa_{direction}.rds"), "Au_No_≈T",
  glue("400samples_25parameters_12h_01102023_SDTa_{direction}.rds"), "Au_No_σT",
  glue("400samples_25parameters_12h_01102023_gradTa_{direction}.rds"), "Au_No_∇T"
)

##########
# CREATE 200 + 400 DATASETS
##########

files_info_all <- bind_rows(
  files_info %>%
    mutate(file = gsub("400samples", "200samples", file),
           n_samples = 200),

  files_info %>%
    mutate(n_samples = 400)
)

##########
# PARAMETERS
##########

param_order <- c(
  "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
  "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h",
  "h", "g_macro", "infl_macro", "infl_soil", "infl_forest", "g_forest", "p_ground", "g_soil", "k_soil"
)

top_params <- c("infl_macro", "infl_soil", "k_soil")

##########
# LOAD + NORMALISE
##########

sobol_df_all <- files_info_all %>%
  mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
  mutate(S = map(data, ~ .x$T$original)) %>%
  select(label, S, n_samples) %>%
  unnest(S) %>%
  mutate(parameter = rep(param_order, times = nrow(.) / length(param_order))) %>%
  group_by(label, n_samples) %>%
  mutate(norm_value = pmax(0, S) / sum(pmax(0, S))) %>%
  ungroup()

###############################################################################
# CONVERGENCE ANALYSIS
###############################################################################

# Individual parameters
conv_params <- sobol_df_all %>%
  filter(parameter %in% top_params) %>%
  group_by(n_samples, parameter) %>%
  summarise(
    mean_index = mean(norm_value),
    sd_index   = sd(norm_value),
    se         = sd_index / sqrt(n()),
    .groups = "drop"
  )

# Combined contribution (KEY RESULT)
conv_sum <- sobol_df_all %>%
  filter(parameter %in% top_params) %>%
  group_by(n_samples, label) %>%
  summarise(sum_index = sum(norm_value), .groups = "drop") %>%
  group_by(n_samples) %>%
  summarise(
    mean_sum = mean(sum_index),
    sd_sum   = sd(sum_index),
    se       = sd_sum / sqrt(n()),
    .groups = "drop"
  )

###############################################################################
# PLOT
###############################################################################

p1 <- ggplot(conv_params, aes(x = n_samples, y = mean_index, colour = parameter)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_index - se, ymax = mean_index + se), width = 20) +
  scale_x_continuous(breaks = c(200, 400)) +
  labs(
    title = "(a) Individual parameters",
    x = "Number of LHS samples",
    y = "Mean normalised Sobol index",
    colour = "Parameter"
  ) +
  theme_bw(base_size = 14)

p2 <- ggplot(conv_sum, aes(x = n_samples, y = mean_sum)) +
  geom_line(linewidth = 1.2, colour = "black") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_sum - se, ymax = mean_sum + se), width = 20) +
  scale_x_continuous(breaks = c(200, 400)) +
  labs(
    title = "(b) Combined contribution",
    x = "Number of LHS samples",
    y = "Sum of Sobol indices"
  ) +
  theme_bw(base_size = 14)

p_convergence <- p1 / p2

ggsave(
  file.path(output_path_plots, paste0("Sobol_convergence_airT_", direction, ".png")),
  p_convergence,
  width = 8,
  height = 10,
  dpi = 400
)

###############################################################################
# OPTIONAL: CHECK DIFFERENCES
###############################################################################

diff_check <- conv_params %>%
  pivot_wider(names_from = n_samples, values_from = mean_index) %>%
  mutate(diff = abs(`400` - `200`))

print(diff_check)

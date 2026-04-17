###############################################################################
# Sobol convergence analysis (200 vs 400 LHS samples)
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(glue)
library(latex2exp)

#########
# INPUT #
#########

input_path <- "Output/sensitivity_analysis/Sobol_QoI/data/"
output_path_plots <- "Output/sensitivity_analysis/Sobol_QoI/plots_output/"

##########
# PARAMETERS
##########

param_order <- c(
  "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
  "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h",
  "h", "g_macro", "infl_macro", "infl_soil", "infl_forest", "g_forest", "p_ground", "g_soil", "k_soil"
)

top_params <- c("infl_macro", "infl_soil", "k_soil")

###############################################################################
# FUNCTION: run one case
###############################################################################

run_case <- function(temp, direction) {

  files_info <- tribble(
    ~file, ~label,

    # winter
    glue("400samples_25parameters_01h_13012025_avT{temp}_{direction}.rds"), "Wi_Ni_≈T",
    glue("400samples_25parameters_01h_13012025_SDT{temp}_{direction}.rds"), "Wi_Ni_σT",
    glue("400samples_25parameters_01h_13012025_gradT{temp}_{direction}.rds"), "Wi_Ni_∇T",

    glue("400samples_25parameters_08h_13012025_avT{temp}_{direction}.rds"), "Wi_Mo_≈T",
    glue("400samples_25parameters_08h_13012025_SDT{temp}_{direction}.rds"), "Wi_Mo_σT",
    glue("400samples_25parameters_08h_13012025_gradT{temp}_{direction}.rds"), "Wi_Mo_∇T",

    glue("400samples_25parameters_12h_13012025_avT{temp}_{direction}.rds"), "Wi_No_≈T",
    glue("400samples_25parameters_12h_13012025_SDT{temp}_{direction}.rds"), "Wi_No_σT",
    glue("400samples_25parameters_12h_13012025_gradT{temp}_{direction}.rds"), "Wi_No_∇T",

    # spring
    glue("400samples_25parameters_01h_30042025_avT{temp}_{direction}.rds"), "Sp_Ni_≈T",
    glue("400samples_25parameters_01h_30042025_SDT{temp}_{direction}.rds"), "Sp_Ni_σT",
    glue("400samples_25parameters_01h_30042025_gradT{temp}_{direction}.rds"), "Sp_Ni_∇T",

    glue("400samples_25parameters_05h_30042025_avT{temp}_{direction}.rds"), "Sp_Mo_≈T",
    glue("400samples_25parameters_05h_30042025_SDT{temp}_{direction}.rds"), "Sp_Mo_σT",
    glue("400samples_25parameters_05h_30042025_gradT{temp}_{direction}.rds"), "Sp_Mo_∇T",

    glue("400samples_25parameters_12h_30042025_avT{temp}_{direction}.rds"), "Sp_No_≈T",
    glue("400samples_25parameters_12h_30042025_SDT{temp}_{direction}.rds"), "Sp_No_σT",
    glue("400samples_25parameters_12h_30042025_gradT{temp}_{direction}.rds"), "Sp_No_∇T",

    # summer
    glue("400samples_25parameters_01h_07072023_avT{temp}_{direction}.rds"), "Su_Ni_≈T",
    glue("400samples_25parameters_01h_07072023_SDT{temp}_{direction}.rds"), "Su_Ni_σT",
    glue("400samples_25parameters_01h_07072023_gradT{temp}_{direction}.rds"), "Su_Ni_∇T",

    glue("400samples_25parameters_05h_07072023_avT{temp}_{direction}.rds"), "Su_Mo_≈T",
    glue("400samples_25parameters_05h_07072023_SDT{temp}_{direction}.rds"), "Su_Mo_σT",
    glue("400samples_25parameters_05h_07072023_gradT{temp}_{direction}.rds"), "Su_Mo_∇T",

    glue("400samples_25parameters_12h_07072023_avT{temp}_{direction}.rds"), "Su_No_≈T",
    glue("400samples_25parameters_12h_07072023_SDT{temp}_{direction}.rds"), "Su_No_σT",
    glue("400samples_25parameters_12h_07072023_gradT{temp}_{direction}.rds"), "Su_No_∇T",

    # autumn
    glue("400samples_25parameters_01h_01102023_avT{temp}_{direction}.rds"), "Au_Ni_≈T",
    glue("400samples_25parameters_01h_01102023_SDT{temp}_{direction}.rds"), "Au_Ni_σT",
    glue("400samples_25parameters_01h_01102023_gradT{temp}_{direction}.rds"), "Au_Ni_∇T",

    glue("400samples_25parameters_06h_01102023_avT{temp}_{direction}.rds"), "Au_Mo_≈T",
    glue("400samples_25parameters_06h_01102023_SDT{temp}_{direction}.rds"), "Au_Mo_σT",
    glue("400samples_25parameters_06h_01102023_gradT{temp}_{direction}.rds"), "Au_Mo_∇T",

    glue("400samples_25parameters_12h_01102023_avT{temp}_{direction}.rds"), "Au_No_≈T",
    glue("400samples_25parameters_12h_01102023_SDT{temp}_{direction}.rds"), "Au_No_σT",
    glue("400samples_25parameters_12h_01102023_gradT{temp}_{direction}.rds"), "Au_No_∇T"
  )

  ##########
  # ONLY TRUE 200–400 PAIRS
  ##########

  valid_cases <- files_info %>%
    mutate(file_200 = gsub("400samples", "200samples", file)) %>%
    filter(file.exists(file.path(input_path, file_200)))

  files_info_all <- bind_rows(
    valid_cases %>% transmute(label, file = file_200, n_samples = 200),
    valid_cases %>% transmute(label, file = file,     n_samples = 400)
  )

  ##########
  # LOAD + NORMALISE
  ##########

  sobol_df <- files_info_all %>%
    mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
    mutate(S = map(data, ~ .x$T$original)) %>%
    select(label, S, n_samples) %>%
    unnest(S) %>%
    mutate(parameter = rep(param_order, times = nrow(.) / length(param_order))) %>%
    group_by(label, n_samples) %>%
    mutate(norm_value = pmax(0, S) / sum(pmax(0, S))) %>%
    ungroup()

  ##########
  # CONVERGENCE (combined)
  ##########

  if (temp == "a") {
    top_params <- c("infl_macro", "infl_soil", "k_soil")
    labels_latex <- c(
      "infl_macro" = TeX("$i_m$"),
      "infl_soil"  = TeX("$i_s$"),
      "k_soil"     = TeX("$k_s$")
    )
  } else if (temp == "f") {
    top_params <- c("infl_macro", "infl_soil", "g_forest")
    labels_latex <- c(
      "infl_macro" = TeX("$i_m$"),
      "infl_soil"  = TeX("$i_s$"),
      "g_forest"   = TeX("$g_f$")
    )
  }

  # Individual parameters
  conv_params <- sobol_df %>%
    filter(parameter %in% top_params) %>%
    group_by(n_samples, parameter) %>%
    summarise(
      mean_value = mean(norm_value),
      se = sd(norm_value) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(type = parameter)

  # Combined contribution
  conv_sum <- sobol_df %>%
    filter(parameter %in% top_params) %>%
    group_by(n_samples, label) %>%
    summarise(sum_index = sum(norm_value), .groups = "drop") %>%
    group_by(n_samples) %>%
    summarise(
      mean_value = mean(sum_index),
      se = sd(sum_index) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(type = "sum")

  # Combine + RETURN
  result <- bind_rows(conv_params, conv_sum) %>%
    mutate(
      temp = temp,
      direction = direction,
      case = paste0(
        ifelse(temp == "a", "Air temperature", "Forest surface temperature"),
        " \n ",
        ifelse(direction == "v", "vertical", "horizontal")
      )
    )


  return(result)
}

###############################################################################
# RUN ALL 4 CASES
###############################################################################

all_results <- expand_grid(
  temp = c("a", "f"),
  direction = c("v", "h")
) %>%
  pmap_dfr(~ run_case(..1, ..2))

all_results$type <- factor(
  all_results$type,
  levels = c("sum", "infl_macro", "infl_soil", "k_soil", "g_forest")
)

###############################################################################
# FINAL PLOT (FACETTED)
###############################################################################

p_final <- ggplot(all_results,
                  aes(x = n_samples, y = mean_value, colour = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_value - se, ymax = mean_value + se),
    width = 20
  ) +
  facet_wrap(~ case, ncol = 2) +
  scale_x_continuous(breaks = c(200, 400)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_colour_manual(
    values = c(
      "sum" = "black",
      "infl_macro" = "blue",
      "infl_soil"  = "red",
      "k_soil"     = "orange",
      "g_forest"   = "darkgreen"
    ),
    labels = c(
      "sum"          = "Sum",
      "infl_macro"   = expression(italic(i)[m]),
      "infl_soil"    = expression(italic(i)[s]),
      "k_soil"       = expression(italic(k)[s]),
      "g_forest"     = expression(italic(g)[f])
    ),
    breaks = c("sum", "infl_macro", "infl_soil", "k_soil", "g_forest"),
    name = NULL
  ) +
  labs(
    x = "Number of LHS samples",
    y = "Normalised Sobol index"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),

    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(fill = "grey90"),

    legend.text = element_text(size = 18),
    legend.title = element_text(size = 14),

    legend.position = "bottom"
  ) +
  geom_text(
    data = all_results %>% filter(type == "sum"),
    aes(label = sprintf("%.2f", mean_value)),
    nudge_y = 0.1,
    size = 5,
    colour = "black"
  )


print(p_final)

ggsave(
  file.path(output_path_plots, "Sobol_convergence_all_cases.png"),
  p_final,
  width = 8,
  height = 8,
  dpi = 300
)




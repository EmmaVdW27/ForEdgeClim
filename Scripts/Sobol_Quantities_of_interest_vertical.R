###############################################################################
# In this script all 'Sobol_all_at_once.R' rds outputfiles for the several
# seasons and time points are merged together to get an overview of the
# Quantities of Interest (4 seasons x 3 time points x 3 metrics).
# This script runs for the vertical line, along the central Y- and X-line.
###############################################################################


library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(forcats)

#########
# INPUT #
#########

input_path <- "Output/sensitivity_analysis/Sobol_QoI/data/"
output_path_plots <- "Output/sensitivity_analysis/Sobol_QoI/plots_output/"
output_path_numbers <- "Output/sensitivity_analysis/Sobol_QoI/numbers_output/"

output_plot_normalized <- "Sobol_QoI_vertical_normalized.png"
output_plot_unnormalized <- "Sobol_QoI_vertical_unnormalized.png"

output_plot_focused_normalized <- "Sobol_QoI_focused_vertical_normalized.png"
output_plot_focused_unnormalized <- "Sobol_QoI_focused_vertical_unnormalized.png"

output_plot_parameters_by_condition_normalized <- "Sobol_parameter_by_condition_vertical_normalized.png"
output_plot_parameters_by_condition_unnormalized <- "Sobol_parameter_by_condition_vertical_unnormalized.png"

output_plot_focused_parameters_by_condition_normalized <- "Sobol_focused_parameter_by_condition_vertical_normalized.png"
output_plot_focused_parameters_by_condition_unnormalized <- "Sobol_focused_parameter_by_condition_vertical_unnormalized.png"

# File list with labels
files_info <- tribble(
  ~file, ~label,
  # winter files
    # night
  "400samples_25parameters_01h_13012025_avT_v.rds",     "Wi_Ni_≈T",
  "400samples_25parameters_01h_13012025_SDT_v.rds",     "Wi_Ni_σT",
  "400samples_25parameters_01h_13012025_gradT_v.rds",   "Wi_Ni_∇T",
    # morning
  "400samples_25parameters_08h_13012025_avT_v.rds",     "Wi_Mo_≈T",
  "400samples_25parameters_08h_13012025_SDT_v.rds",     "Wi_Mo_σT",
  "400samples_25parameters_08h_13012025_gradT_v.rds",   "Wi_Mo_∇T",
    # noon
  "400samples_25parameters_12h_13012025_avT_v.rds",     "Wi_No_≈T",
  "400samples_25parameters_12h_13012025_SDT_v.rds",     "Wi_No_σT",
  "400samples_25parameters_12h_13012025_gradT_v.rds",   "Wi_No_∇T",
  # spring
    # night
  "400samples_25parameters_01h_30042025_avT_v.rds",     "Sp_Ni_≈T",
  "400samples_25parameters_01h_30042025_SDT_v.rds",     "Sp_Ni_σT",
  "400samples_25parameters_01h_30042025_gradT_v.rds",   "Sp_Ni_∇T",
    # morning
  "400samples_25parameters_05h_30042025_avT_v.rds",     "Sp_Mo_≈T",
  "400samples_25parameters_05h_30042025_SDT_v.rds",     "Sp_Mo_σT",
  "400samples_25parameters_05h_30042025_gradT_v.rds",   "Sp_Mo_∇T",
    # noon
  "400samples_25parameters_12h_30042025_avT_v.rds",     "Sp_No_≈T",
  "400samples_25parameters_12h_30042025_SDT_v.rds",     "Sp_No_σT",
  "400samples_25parameters_12h_30042025_gradT_v.rds",   "Sp_No_∇T",
  # summer
    # night
  "400samples_25parameters_01h_07072023_avT_v.rds",     "Su_Ni_≈T",
  "400samples_25parameters_01h_07072023_SDT_v.rds",     "Su_Ni_σT",
  "400samples_25parameters_01h_07072023_gradT_v.rds",   "Su_Ni_∇T",
    # morning
  "400samples_25parameters_05h_07072023_avT_v.rds",     "Su_Mo_≈T",
  "400samples_25parameters_05h_07072023_SDT_v.rds",     "Su_Mo_σT",
  "400samples_25parameters_05h_07072023_gradT_v.rds",   "Su_Mo_∇T",
    # noon
  "400samples_25parameters_12h_07072023_avT_v.rds",     "Su_No_≈T",
  "400samples_25parameters_12h_07072023_SDT_v.rds",     "Su_No_σT",
  "400samples_25parameters_12h_07072023_gradT_v.rds",   "Su_No_∇T",
  # autumn files
    # night
  "400samples_25parameters_01h_01102023_avT_v.rds",     "Au_Ni_≈T",
  "400samples_25parameters_01h_01102023_SDT_v.rds",     "Au_Ni_σT",
  "400samples_25parameters_01h_01102023_gradT_v.rds",   "Au_Ni_∇T",
    # morning
  "400samples_25parameters_06h_01102023_avT_v.rds",     "Au_Mo_≈T",
  "400samples_25parameters_06h_01102023_SDT_v.rds",     "Au_Mo_σT",
  "400samples_25parameters_06h_01102023_gradT_v.rds",   "Au_Mo_∇T",
    # noon
  "400samples_25parameters_12h_01102023_avT_v.rds",     "Au_No_≈T",
  "400samples_25parameters_12h_01102023_SDT_v.rds",     "Au_No_σT",
  "400samples_25parameters_12h_01102023_gradT_v.rds",   "Au_No_∇T"

)

##########
# OUTPUT #
##########

# Parameter sequence
param_order <- c(
  "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
  "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h",
  "h", "g_macro", "infl_macro", "infl_soil", "infl_forest", "g_forest", "p_ground", "g_soil", "k_soil"
)

# Color palette
param_colors <- c(
  colorRampPalette(c("gold", "darkorange"))(9),
  colorRampPalette(c("mediumpurple", "darkviolet"))(7),
  colorRampPalette(c("darkseagreen", "darkgreen"))(9)
)
names(param_colors) <- param_order

# Focused parameter sequence
focus_params <- c("h", "g_macro", "infl_macro", "infl_soil", "infl_forest",
                  "g_forest", "p_ground", "g_soil", "k_soil")

# Focused color palette: contrast rich
focus_colors <- setNames(RColorBrewer::brewer.pal(n = 9, name = "Paired"), focus_params)

# Load, combine and normalise
sobol_df <- files_info %>%
  mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
  mutate(S = map(data, ~ .x$T$original)) %>%
  select(label, S) %>%
  unnest(S) %>%
  mutate(parameter = rep(param_order, times = nrow(.) / length(param_order))) %>%
  group_by(label) %>%
  mutate(norm_value = pmax(0, S) / sum(pmax(0, S))) %>%
  ungroup() %>%
  mutate(parameter = factor(parameter, levels = param_order),
         label = fct_rev(fct_inorder(label))) %>%
  # define conditions
  mutate(
    season = substr(label, 1, 2),
    moment = substr(label, 4, 5),
    metric = case_when(
      grepl("≈T", label) ~ "avT",
      grepl("σT", label) ~ "sigmaT",
      grepl("∇T", label) ~ "gradT",
      TRUE ~ NA_character_
    )
  )

# Define to each parameter a process group
param_sets <- tribble(
  ~parameter,        ~process,
  "betad","SW","beta0","SW","omega","SW","Kd_v","SW","Kb_v","SW","omega_g_v","SW","Kd_h","SW","Kb_h","SW","omega_g_h","SW",
  "e_forest","LW","beta_lw","LW","omega_lw","LW","Kd_lw_v","LW","omega_g_lw_v","LW","Kd_lw_h","LW","omega_g_lw_h","LW",
  "h","HEAT","g_macro","HEAT","infl_macro","HEAT","infl_soil","HEAT","infl_forest","HEAT","g_forest","HEAT","p_ground","HEAT","g_soil","HEAT","k_soil","HEAT"
) %>% mutate(process = factor(process, levels = c("SW","LW","HEAT")))

sobol_df2 <- sobol_df %>%
  left_join(param_sets, by = "parameter")

###############################################################################
# PLOTTING RESULTS
###############################################################################


# Plot QoI all parameters

# normalized
p_n = ggplot(sobol_df, aes(x = norm_value, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Parameter contribution to QoI variance along the vertical line",
    x = "Normalised total-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = param_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
        )
ggsave(paste0(output_path_plots, output_plot_normalized), plot = p_n, width = 10, height = 10, dpi = 500)

# unnormalized
p_u = ggplot(sobol_df, aes(x = S, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Parameter contribution to QoI variance along the vertical line",
    x = "Absolute total-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = param_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
  )
ggsave(paste0(output_path_plots, output_plot_unnormalized), plot = p_u, width = 10, height = 10, dpi = 500)

#------------------------------------------------------------------------------

# Plot QoI focused parameters

sobol_focus_df <- sobol_df %>%
  filter(parameter %in% focus_params)

# normalized
p_focus_n <- ggplot(sobol_focus_df, aes(x = norm_value, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Key parameter contribution to QoI variance along the vertical line",
    x = "Normalised total-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = focus_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_focused_normalized), plot = p_focus_n, width = 10, height = 10, dpi = 500)

# unnormalized
p_focus_u <- ggplot(sobol_focus_df, aes(x = S, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Key parameter contribution to QoI variance along the vertical line",
    x = "Absolute total-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = focus_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_focused_unnormalized), plot = p_focus_u, width = 10, height = 10, dpi = 500)

#------------------------------------------------------------------------------

# Plot mean parameter contribution per condition (metric, moment, season)

# normalized Sobol index

# Dataframe per condition
# Mean per metric
mean_metric_n <- sobol_df %>%
  filter(!is.na(metric)) %>%
  group_by(metric, parameter) %>%
  summarise(mean_index = mean(norm_value), .groups = "drop") %>%
  rename(condition = metric)

# Mean per moment
mean_moment_n <- sobol_df %>%
  filter(!is.na(moment)) %>%
  group_by(moment, parameter) %>%
  summarise(mean_index = mean(norm_value), .groups = "drop") %>%
  rename(condition = moment)

# Mean per season
mean_season_n <- sobol_df %>%
  filter(!is.na(season)) %>%
  group_by(season, parameter) %>%
  summarise(mean_index = mean(norm_value), .groups = "drop") %>%
  rename(condition = season)

# Combine dataframes
mean_combined_n <- bind_rows(
  mutate(mean_metric_n, group = "Metric"),
  mutate(mean_moment_n, group = "Moment"),
  mutate(mean_season_n, group = "Season")
) %>%
  mutate(
    condition = as.character(condition),
    condition = recode(condition,
                       "avT" = "≈T",
                       "sigmaT" = "σT",
                       "gradT" = "∇T"
    ),
    condition = factor(condition, levels = rev(c("≈T", "σT", "∇T", "Ni", "Mo", "No", "Wi", "Sp", "Su", "Au"))),
    group = factor(group, levels = c("Metric", "Moment", "Season")),
    parameter = factor(parameter, levels = param_order)
  )

# Plot
p_condition_n <- ggplot(mean_combined_n, aes(x = mean_index, y = condition, fill = parameter)) +
  geom_col(width = 0.5) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = param_colors) +
  labs(
    title = "Parameter contribution to condition variance along the vertical line",
    x = "Average normalised total-order Sobol-index",
    y = "Condition",
    fill = "Parameter",
    caption = "Sp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_parameters_by_condition_normalized), plot = p_condition_n, width = 12, height = 10, dpi = 500)

#------------------------------------------------------------------------------

# Plot mean focused parameter contribution per condition (metric, moment, season)

mean_combined_focused_n <- mean_combined_n %>%
  filter(parameter %in% focus_params) %>%
  mutate(parameter = factor(parameter, levels = focus_params))

# Plot
p_condition_focused_n <- ggplot(mean_combined_focused_n, aes(x = mean_index, y = condition, fill = parameter)) +
  geom_col(width = 0.5) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = focus_colors) +
  labs(
    title = "Key parameter contribution to condition variance along the vertical line",
    x = "Average normalised total-order Sobol-index",
    y = "Condition",
    fill = "Parameter",
    caption = "Sp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_focused_parameters_by_condition_normalized), plot = p_condition_focused_n, width = 12, height = 10, dpi = 500)

#------------------------------------------------------------------------------

# Unnormalized Sobol index

# Dataframe per condition
# Mean per metric
mean_metric_u <- sobol_df %>%
  filter(!is.na(metric)) %>%
  group_by(metric, parameter) %>%
  summarise(mean_index = mean(S), .groups = "drop") %>%
  rename(condition = metric)

# Mean per moment
mean_moment_u <- sobol_df %>%
  filter(!is.na(moment)) %>%
  group_by(moment, parameter) %>%
  summarise(mean_index = mean(S), .groups = "drop") %>%
  rename(condition = moment)

# Mean per season
mean_season_u <- sobol_df %>%
  filter(!is.na(season)) %>%
  group_by(season, parameter) %>%
  summarise(mean_index = mean(S), .groups = "drop") %>%
  rename(condition = season)

# Combine dataframes
mean_combined_u <- bind_rows(
  mutate(mean_metric_u, group = "Metric"),
  mutate(mean_moment_u, group = "Moment"),
  mutate(mean_season_u, group = "Season")
) %>%
  mutate(
    condition = as.character(condition),
    condition = recode(condition,
                       "avT" = "≈T",
                       "sigmaT" = "σT",
                       "gradT" = "∇T"
    ),
    condition = factor(condition, levels = rev(c("≈T", "σT", "∇T", "Ni", "Mo", "No", "Wi", "Sp", "Su", "Au"))),
    group = factor(group, levels = c("Metric", "Moment", "Season")),
    parameter = factor(parameter, levels = param_order)
  )

# Plot
p_condition_u <- ggplot(mean_combined_u, aes(x = mean_index, y = condition, fill = parameter)) +
  geom_col(width = 0.5) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = param_colors) +
  labs(
    title = "Parameter contribution to condition variance along the vertical line",
    x = "Average absolute total-order Sobol-index",
    y = "Condition",
    fill = "Parameter",
    caption = "Sp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_parameters_by_condition_unnormalized), plot = p_condition_u, width = 12, height = 10, dpi = 500)

#------------------------------------------------------------------------------

# Plot mean focused parameter contribution per condition (metric, moment, season)

mean_combined_focused_u <- mean_combined_u %>%
  filter(parameter %in% focus_params) %>%
  mutate(parameter = factor(parameter, levels = focus_params))

# Plot
p_condition_focused_u <- ggplot(mean_combined_focused_u, aes(x = mean_index, y = condition, fill = parameter)) +
  geom_col(width = 0.5) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = focus_colors) +
  labs(
    title = "Key parameter contribution to condition variance along the vertical line",
    x = "Average absolute total-order Sobol-index",
    y = "Condition",
    fill = "Parameter",
    caption = "Sp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(paste0(output_path_plots, output_plot_focused_parameters_by_condition_unnormalized), plot = p_condition_focused_u, width = 12, height = 10, dpi = 500)



###############################################################################
# QUANTITATIVE RESULTS (only for normalized Sobol analysis (yet))
###############################################################################

# Top 3 parameter contributors per condition

top3_per_label <- sobol_df %>%
  group_by(label) %>%
  arrange(desc(norm_value), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  summarise(
    top1_param = parameter[1], top1_val = norm_value[1],
    top2_param = parameter[2], top2_val = norm_value[2],
    top3_param = parameter[3], top3_val = norm_value[3],
    top3_sum   = sum(norm_value[1:3], na.rm = TRUE),
    .groups = "drop"
  )

# Herfindahl–Hirschman (concentration) and entropy as measure for spread of influence

concentration_per_label <- sobol_df %>%
  group_by(label) %>%
  summarise(
    # HHI: square each normalized value and sum them
    # If all influence is concentrated in one category (monoculture), norm_value = 1 for one, 0 for rest => HHI = 1
    # If influence is evenly distributed across n categories: HHI = 1/n (so the lower the value, the greater the spread)
    HHI = sum(norm_value^2),
    # Shannon entropy measures uncertainty/spread (-sum(norm_value * log(norm_value)))
    # Maximum entropy occurs when all categories are equal in size (maximum spread)
    # Entropy = 0 when one category holds 100% of the influence (no spread)
    entropy = -sum(ifelse(norm_value>0, norm_value*log(norm_value), 0)),
    .groups = "drop"
  )


# Share per process per condition

process_share_per_label <- sobol_df2 %>%
  group_by(label, process) %>%
  summarise(process_share = sum(norm_value), .groups = "drop")


#------------------------------------------------------------------------------

# Mean and sd parameter contribution per metric, moment and season
avg_param_by_metric  <- sobol_df %>% filter(!is.na(metric)) %>%
  group_by(metric, parameter) %>%
  summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")

avg_param_by_moment  <- sobol_df %>% filter(!is.na(moment)) %>%
  group_by(moment, parameter) %>%
  summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")

avg_param_by_season  <- sobol_df %>% filter(!is.na(season)) %>%
  group_by(season, parameter) %>%
  summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")

# get ranges on the contributions for top 3 parameters across conditions
target_params <- c("infl_macro", "infl_soil", "g_macro")

range_metric <- avg_param_by_metric %>%
  filter(parameter %in% target_params) %>%
  group_by(metric) %>%
  summarise(sum_share = sum(mean_index), .groups = "drop") %>%
  summarise(min_sum = min(sum_share),
            max_sum = max(sum_share))

range_moment <- avg_param_by_moment %>%
  filter(parameter %in% target_params) %>%
  group_by(moment) %>%
  summarise(sum_share = sum(mean_index), .groups = "drop") %>%
  summarise(min_sum = min(sum_share),
            max_sum = max(sum_share))

range_season <- avg_param_by_season %>%
  filter(parameter %in% target_params) %>%
  group_by(season) %>%
  summarise(sum_share = sum(mean_index), .groups = "drop") %>%
  summarise(min_sum = min(sum_share),
            max_sum = max(sum_share))

print(range_metric)
print(range_moment)
print(range_season)


# Mean and sd process contribution per metric, moment and season
process_by_label <- sobol_df2 %>%
  group_by(label, metric, moment, season, process) %>%
  summarise(process_share = sum(norm_value), .groups = "drop")

avg_process_by_metric <- process_by_label %>% filter(!is.na(metric)) %>%
  group_by(metric, process) %>%
  summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")

avg_process_by_moment <- process_by_label %>% filter(!is.na(moment)) %>%
  group_by(moment, process) %>%
  summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")

avg_process_by_season <- process_by_label %>% filter(!is.na(season)) %>%
  group_by(season, process) %>%
  summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")

#------------------------------------------------------------------------------

# Stability : rank correlation and variation
# Rank per condition (1 = most important)
ranks_long <- sobol_df %>%
  group_by(label) %>%
  mutate(rank = dense_rank(desc(norm_value))) %>%
  ungroup()

# Consistency: mean rang + sd per parameter for every condition
rank_consistency_moment <- ranks_long %>%
  group_by(moment, parameter) %>%
  summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")

rank_consistency_season <- ranks_long %>%
  group_by(season, parameter) %>%
  summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")

rank_consistency_metric <- ranks_long %>%
  group_by(metric, parameter) %>%
  summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")

# get ranges on the rank and sd on rank for top 3 parameters across conditions
# Moment
range_moment_rank <- rank_consistency_moment %>%
  filter(parameter %in% target_params) %>%
  summarise(
    min_mean_rank = min(mean_rank),
    max_mean_rank = max(mean_rank),
    min_sd_rank   = min(sd_rank),
    max_sd_rank   = max(sd_rank)
  )

# Season
range_season_rank <- rank_consistency_season %>%
  filter(parameter %in% target_params) %>%
  summarise(
    min_mean_rank = min(mean_rank),
    max_mean_rank = max(mean_rank),
    min_sd_rank   = min(sd_rank),
    max_sd_rank   = max(sd_rank)
  )

# Metric
range_metric_rank <- rank_consistency_metric %>%
  filter(parameter %in% target_params) %>%
  summarise(
    min_mean_rank = min(mean_rank),
    max_mean_rank = max(mean_rank),
    min_sd_rank   = min(sd_rank),
    max_sd_rank   = max(sd_rank)
  )

print(range_metric_rank)
print(range_moment_rank)
print(range_season_rank)

#------------------------------------------------------------------------------

# export tables
write_csv(top3_per_label,          file.path(output_path_numbers, "top3_per_condition_v.csv"))
write_csv(concentration_per_label, file.path(output_path_numbers, "concentration_per_condition_v.csv"))
write_csv(process_share_per_label, file.path(output_path_numbers, "process_share_per_condition_v.csv"))

write_csv(avg_param_by_metric,     file.path(output_path_numbers, "avg_param_by_metric_v.csv"))
write_csv(avg_param_by_moment,     file.path(output_path_numbers, "avg_param_by_moment_v.csv"))
write_csv(avg_param_by_season,     file.path(output_path_numbers, "avg_param_by_season_v.csv"))

write_csv(avg_process_by_metric,   file.path(output_path_numbers, "avg_process_by_metric_v.csv"))
write_csv(avg_process_by_moment,   file.path(output_path_numbers, "avg_process_by_moment_v.csv"))
write_csv(avg_process_by_season,   file.path(output_path_numbers, "avg_process_by_season_v.csv"))

write_csv(rank_consistency_metric, file.path(output_path_numbers, "rank_consistency_metric_v.csv"))
write_csv(rank_consistency_moment, file.path(output_path_numbers, "rank_consistency_moment_v.csv"))
write_csv(rank_consistency_season, file.path(output_path_numbers, "rank_consistency_season_v.csv"))



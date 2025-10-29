############################################################################################
# This script validates ForEdgeClim's model output (calibrated & uncalibrated)
# against observational TOMST TMS sensors.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClim)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
  library(purrr)
  library(stringr)
})

# -----------------
# INPUT
# -----------------

# Set the base season. The script will also look for "<base>_uncalibrated".
base_season <- "spring"   # summer | spring | autumn | winter

# Build file map: names = model label, values = season-string used in filenames
season_map <- c(
  calibrated   = base_season,
  uncalibrated = paste0(base_season, "_uncalibrated")
)

# -----------------
# OUTPUT
# -----------------

output_path <- "Output/validation/plots/"

# -----------------
# LOAD DATA (both calibrated & uncalibrated)
# -----------------

read_val <- function(season_string, model_label){
  fp <- file.path("Output", "validation/data", paste0("val_data_", season_string, ".rds"))
  if(!file.exists(fp)){
    stop(sprintf("File not found: %s", fp))
  }
  readRDS(fp) %>%
    mutate(model = model_label) %>%
    # ensure expected columns exist
    relocate(model, .before = 1)
}

val_list <- imap(season_map, ~ read_val(.x, .y))  # .y = model label
val_data_all <- bind_rows(val_list)

# For convenience keep the calibrated subset as before
val_data <- val_data_all %>% filter(model == "calibrated")

# -----------------
# STATISTICS (per model)
# -----------------

stats_df <- val_data_all %>%
  group_by(model) %>%
  summarise(
    ME   = mean(resid, na.rm = TRUE), # mean error (bias)
    RMSE = sqrt(mean(resid^2, na.rm = TRUE)), # root mean square error (spread)
    R2   = cor(obs, mod, use = "complete.obs")^2, # R2 (model fit)
    NSE  = 1 - sum(resid^2, na.rm = TRUE) / # Nash–Sutcliffe model efficiency coefficient (overall performance)
      sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE),
    .groups = "drop"
  )

# Labels for annotations per model (for the scatter/facets)
lab_df <- stats_df %>%
  transmute(
    model,
    lab_scatter = sprintf("R\u00B2 = %.2f\nNSE = %.2f", R2, NSE),
    lab_resid   = sprintf("RMSE = %.2f \u00B0C\nME = %.2f \u00B0C", RMSE, ME)
  )

# -----------------
# MODEL vs OBSERVATION (all timesteps x positions) — facet per model
# -----------------

# Position of the 1:1 label per facet
one2one_pos <- val_data_all %>%
  group_by(model) %>%
  summarise(
    x_min = min(obs, na.rm = TRUE),
    y_min = min(obs, na.rm = TRUE),
    .groups = "drop"
  )

m_o_time_and_position <- ggplot(val_data_all, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.6, colour = "cornflowerblue") +
  geom_abline(slope = 1, intercept = 0, colour = "#f38ba8", linewidth = 1) +
  labs(x = "Observed Tair (°C)", y = "Modelled Tair (°C)",
       title = "Model vs observation (timesteps × loggers datapoints)") +
  theme_bw() +
  facet_wrap(~ model, ncol = 2) +
  # per-facet stats label
  geom_label(
    data = lab_df,
    aes(x = -Inf, y = Inf, label = lab_scatter),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.2,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white", colour = "black"
  ) +
  # 1:1 annotation per facet
  geom_text(
    data = one2one_pos,
    aes(x = x_min, y = y_min, label = "1:1"),
    inherit.aes = FALSE,
    colour = "#f38ba8", hjust = -2, vjust = -1, fontface = "bold"
  )

print(m_o_time_and_position)
ggsave(file.path(output_path, paste0("model_vs_obs_over_time_and_position_", base_season, "_both.png")),
       plot = m_o_time_and_position, width = 10, height = 6, dpi = 300)

# -----------------
# MODEL vs OBS per position for specific time points (1, 5, 12h)
# -----------------

# If there are only 24 unique datetimes, duplicate rows (as in your original code)
val_data_all_2 <- val_data_all %>%
  group_by(model) %>%
  group_modify(~ {
    x <- .x
    if (length(unique(x$datetime)) == 24) rbind(x, x) else x
  }) %>%
  ungroup()

# Mean + 95% CI per TMS_position & moment & model
val_summary <- val_data_all_2 %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
         datetime = ymd_hms(datetime),
         hour = hour(datetime)) %>%
  filter(hour %in% c(1, 5, 12),
         !(TMS_position %in% c(7, 14, 21, 28, 35))) %>% # exclude vertical ones here
  mutate(
    moment = case_when(
      hour == 1  ~ "Night (01:00)",
      hour == 5  ~ "Morning (05:00)",
      hour == 12 ~ "Day (12:00)"
    ),
    moment = factor(moment, levels = c("Night (01:00)", "Morning (05:00)", "Day (12:00)"))
  ) %>%
  group_by(model, TMS_position, moment) %>%
  summarise(
    n        = dplyr::n(),
    mean_obs = mean(obs, na.rm = TRUE),
    sd_obs   = sd(obs, na.rm = TRUE),
    mean_mod = mean(mod, na.rm = TRUE),
    sd_mod   = sd(mod, na.rm = TRUE),
    .groups  = "drop_last"
  ) %>%
  mutate(
    se_obs   = sd_obs / sqrt(pmax(n, 1)),
    se_mod   = sd_mod / sqrt(pmax(n, 1)),
    tcrit    = qt(0.975, df = pmax(n - 1, 1)), # 95% t-CI
    ci_obs_low    = mean_obs - tcrit * se_obs,
    ci_obs_high   = mean_obs + tcrit * se_obs,
    ci_mod_low    = mean_mod - tcrit * se_mod,
    ci_mod_high   = mean_mod + tcrit * se_mod
  ) %>%
  ungroup()

# Ranges per moment (inclusive CIs) across both models
ranges <- val_summary %>%
  group_by(moment) %>%
  summarise(
    rmin = min(ci_obs_low,  ci_mod_low,  na.rm = TRUE),
    rmax = max(ci_obs_high, ci_mod_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pad  = 0.2 * (rmax - rmin + 1e-9),
    xmin = rmin - pad,
    xmax = rmax + pad,
    ymin = rmin - pad,
    ymax = rmax + pad
  )

# Helper to extract a single row of limits as a named list
get_limits <- function(rng_tbl, moment_label){
  rng_tbl %>% filter(moment == moment_label) %>% as.list()
}

# Plot function with model colour + dodged error bars
make_plot <- function(df, limits_row) {
  ggplot(df, aes(x = mean_mod, y = mean_obs, label = TMS_position, colour = model)) +
    geom_errorbar(aes(ymin = ci_obs_low, ymax = ci_obs_high),
                  width = 0, position = position_dodge(width = 0.2)) +
    geom_errorbarh(aes(xmin = ci_mod_low, xmax = ci_mod_high),
                   height = 0, position = position_dodge(width = 0.2)) +
    geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#f38ba8") +
    geom_text(nudge_x = 0.05, nudge_y = 0.05, size = 3, show.legend = FALSE) +
    labs(x = "Observed Tair (°C)",
         y = "Modelled Tair (°C)",
         title = unique(df$moment)) +
    theme_bw() +
    coord_fixed(
      xlim = c(limits_row$xmin, limits_row$xmax),
      ylim = c(limits_row$ymin, limits_row$ymax),
      expand = FALSE
    )
}

df_morning <- val_summary %>% filter(moment == "Morning (05:00)")
df_day     <- val_summary %>% filter(moment == "Day (12:00)")
df_night   <- val_summary %>% filter(moment == "Night (01:00)")

lim_morning <- get_limits(ranges, "Morning (05:00)")
lim_day     <- get_limits(ranges, "Day (12:00)")
lim_night   <- get_limits(ranges, "Night (01:00)")

p_morning <- make_plot(df_morning, lim_morning)
p_day     <- make_plot(df_day,     lim_day)
p_night   <- make_plot(df_night,   lim_night)

m_o_position <- p_night + p_morning + p_day + plot_annotation(
  title = "Modelled vs observed temperature per TMS position (calibrated vs uncalibrated)",
  theme = theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12))
)
print(m_o_position)
ggsave(file.path(output_path, paste0("model_vs_obs_over_position_", base_season, "_both.png")),
       plot = m_o_position, width = 11, height = 6, dpi = 300)

# -----------------
# TIME SERIES: temperature fluctuations (obs + both model variants)
# -----------------

# Keep one copy of the observations (avoid plotting them twice)
obs_df <- val_data_all %>%
  distinct(datetime, TMS_position, obs) %>%
  transmute(datetime, TMS_position, series = "obs", Tair = obs)

mods_df <- val_data_all %>%
  transmute(datetime, TMS_position,
            series = paste0("model_", model),  # "model_calibrated", "model_uncalibrated"
            Tair = mod)

val_data_long <- bind_rows(obs_df, mods_df)

time_fluctuations <- ggplot(val_data_long, aes(x = datetime, y = Tair, colour = series)) +
  geom_point(alpha = 0.6) +
  labs(y = "Air temperature (°C)",
       title = "Time series of temperature fluctuations (obs vs calibrated & uncalibrated)") +
  theme_bw()
print(time_fluctuations)
ggsave(file.path(output_path, paste0("time_series_fluctuations_", base_season, "_both.png")),
       plot = time_fluctuations, width = 10, height = 6, dpi = 300)

# -----------------
# TIME SERIES: residuals (model - obs) per model
# -----------------

# Join stats label to use per-model in a facet
lab_resid_df <- lab_df %>% select(model, lab_resid)

x_left_df <- val_data_all %>%
  group_by(model) %>%
  summarise(x_left = min(datetime, na.rm = TRUE), .groups = "drop")

time_residual <- ggplot(val_data_all, aes(x = datetime, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#f38ba8", linewidth = 1) +
  geom_point(colour = "cornflowerblue", alpha = 0.6) +
  labs(y = "Residual (°C)",
       title = "Time series of temperature residuals (model - observation)") +
  theme_bw() +
  facet_wrap(~ model, ncol = 2) +
  # per-facet stats label
  geom_label(
    data = x_left_df %>% left_join(lab_resid_df, by = "model"),
    aes(x = x_left, y = Inf, label = lab_resid),
    inherit.aes = FALSE,
    hjust = 0.15, vjust = 1.2,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white",
    colour = "black"
  )

print(time_residual)
ggsave(file.path(output_path, paste0("time_series_residual_", base_season, "_both.png")),
       plot = time_residual, width = 10, height = 6, dpi = 300)

# -----------------
# CONSOLE SUMMARIES (per model)
# -----------------

cat("\nValidation statistics per model:\n")
stats_df %>%
  mutate(
    ME   = sprintf("%.2f °C", ME),
    RMSE = sprintf("%.2f °C", RMSE),
    R2   = sprintf("%.2f",   R2),
    NSE  = sprintf("%.2f",   NSE)
  ) %>%
  rename(
    `Model` = model,
    `ME (bias)` = ME,
    `RMSE (spread)` = RMSE,
    `R2 (model fit)` = R2,
    `NSE (overall performance)` = NSE
  ) %>%
  print(row.names = FALSE)

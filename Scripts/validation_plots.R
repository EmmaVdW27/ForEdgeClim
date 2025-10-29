############################################################################################
# This script validates ForEdgeClim's model output with observational TOMST TMS sensors.
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

})


# -----------------
# INPUT
# -----------------

season_string = "summer" # summer spring autumn winter |_uncalibrated

# -----------------
# OUTPUT
# -----------------

# load data
val_data <- readRDS(file.path("Output", "validation/data", paste0("val_data_", season_string, ".rds")))
output_path <- "Output/validation/plots/"

# calculate statistics
mean_error <- mean(val_data$resid, na.rm = TRUE)
rmse       <- sqrt(mean(val_data$resid^2, na.rm = TRUE))
r2         <- cor(val_data$obs, val_data$mod)^2
# Nash–Sutcliffe Efficiency (NSE, total model performance in 1 number):
# NSE = 1 → perfect fit,
# NSE = 0 → model not better than average observation,
# NSE < 0 → model worse than average observation.
nse        <- 1 - sum(val_data$resid^2) / sum((val_data$obs - mean(val_data$obs))^2)

# labels for annotations
lab_scatter <- sprintf("R\u00B2 = %.2f\nNSE = %.2f", r2, nse)
lab_resid   <- sprintf("RMSE = %.2f \u00B0C\nME = %.2f \u00B0C", rmse, mean_error)

# Model vs observation over time and position
# -------------------------------------------

m_o_time_and_position = ggplot(val_data, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.6, colour = "cornflowerblue") +
  geom_abline(slope = 1, intercept = 0, colour = "#f38ba8", linewidth = 1) +
  labs(x = "Observed Tair (°C)", y = "Modelled Tair (°C)",
       title = "Model vs observation (timesteps x loggers datapoints)") +
  theme_bw() +
  annotate("label", x = -Inf, y = Inf, label = lab_scatter,
           hjust = -0.1, vjust = 1.2,
           label.size = 0.3,      # box thickness
           label.r = unit(0.15, "lines"), # curved corner
           fill = "white",
           colour = "black")  +
  # 1:1 annotation
  annotate("text", x = min(val_data$obs, na.rm=TRUE),
           y = min(val_data$obs, na.rm=TRUE),
           label = "1:1", colour = "#f38ba8",
           hjust = -2, vjust = -1, fontface = "bold")

print(m_o_time_and_position)
ggsave(paste0(output_path, "model_vs_obs_over_time_and_position_", season_string, ".png"), plot = m_o_time_and_position, width = 10, height = 6, dpi = 300)



# Model vs observations over position for specific time points
# ------------------------------------------------------------

if(length(unique(val_data$datetime)) == 24){
  val_data_2 = rbind(val_data, val_data)
} else{
  val_data_2 = val_data
}

# Mean + 95% CI per TMS_position and moment
val_summary <- val_data_2 %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  mutate(datetime = ymd_hms(datetime), hour = hour(datetime)) %>%
  filter(hour %in% c(1, 5, 12),
         !(TMS_position %in% c(7,14,21,28,35))) %>% # don't take the vertical ones into account here
  mutate(
    moment = case_when(
      hour == 1  ~ "Night (01:00)",
      hour == 5  ~ "Morning (05:00)",
      hour == 12 ~ "Day (12:00)"
    ),
    moment = factor(moment, levels = c("Night (01:00)", "Morning (05:00)", "Day (12:00)"))
  ) %>%
  group_by(TMS_position, moment) %>%
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
    tcrit    = qt(0.975, df = pmax(n - 1, 1)),     # 95% t-CI (better for small n)
    ci_obs_low    = mean_obs - tcrit * se_obs,
    ci_obs_high   = mean_obs + tcrit * se_obs,
    ci_mod_low    = mean_mod - tcrit * se_mod,
    ci_mod_high   = mean_mod + tcrit * se_mod
  ) %>%
  ungroup()

# Range per moment (inclusive CI’s) + small marge
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

# plot function
make_plot <- function(df, limits_row) {
  ggplot(df, aes(x = mean_mod, y = mean_obs, label = TMS_position)) +
    # CI axes
    geom_errorbar(aes(ymin = ci_obs_low, ymax = ci_obs_high), width = 0) +
    geom_errorbarh(aes(xmin = ci_mod_low, xmax = ci_mod_high), height = 0) +
    geom_point(size = 2, colour = "cornflowerblue", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#f38ba8") +
    geom_text(nudge_x = 0.05, nudge_y = 0.05, size = 3) +
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

# Data per moment + plots
df_morning <- val_summary %>% filter(moment == "Morning (05:00)")
df_day     <- val_summary %>% filter(moment == "Day (12:00)")
df_night   <- val_summary %>% filter(moment == "Night (01:00)")

lim_morning <- ranges %>% filter(moment == "Morning (05:00)") %>% as.list()
lim_day     <- ranges %>% filter(moment == "Day (12:00)")     %>% as.list()
lim_night   <- ranges %>% filter(moment == "Night (01:00)")   %>% as.list()

p_morning <- make_plot(df_morning, lim_morning)
p_day     <- make_plot(df_day,     lim_day)
p_night   <- make_plot(df_night,   lim_night)

# Final, combined plot
m_o_position <- p_night + p_morning + p_day + plot_annotation(
    title = "Modelled vs observed temperature per TMS position",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12))
)
print(m_o_position)
ggsave(paste0(output_path, "model_vs_obs_over_position_", season_string, ".png"), plot = m_o_position, width = 10, height = 6, dpi = 300)




# Time series plot: fluctuations model and obs
# -------------------------------------------

val_data_long <- pivot_longer(
  val_data, cols = c("obs","mod"), names_to = "series", values_to = "Tair"
)

time_fluctuations = ggplot(val_data_long, aes(x = datetime, y = Tair, colour = series)) +
  geom_point(alpha = 0.6) +
  labs(y = "Air temperature (°C)", title = "Time series of temperature fluctuations") +
  theme_bw()
print(time_fluctuations)
ggsave(paste0(output_path, "time_series_fluctuations_", season_string, ".png"), plot = time_fluctuations, width = 10, height = 6, dpi = 300)


# Time series plot : residual model - obs
# ---------------------------------------

# horizontal position on datetime axis
x_left <- min(val_data$datetime, na.rm = TRUE)

time_residual = ggplot(val_data, aes(x = datetime, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#f38ba8", linewidth = 1) +
  geom_point(colour = "cornflowerblue", alpha = 0.6) +
  labs(y = "Residual (°C)", title = "Time series of temperature residuals (model - observation)") +
  theme_bw() +
  annotate("label", x = x_left, y = Inf, label = lab_resid,
           hjust = 0.15, vjust = 1.2,
           label.size = 0.3,      # box thickness
           label.r = unit(0.15, "lines"), # curved corner
           fill = "white",
           colour = "black")
print(time_residual)
ggsave(paste0(output_path, "time_series_residual_", season_string, ".png"), plot = time_residual, width = 10, height = 6, dpi = 300)


# Console summaries
cat("\nValidation statistics:\n")
cat(sprintf("ME (bias): %.2f °C\n", mean_error))
cat(sprintf("RMSE (spread): %.2f °C\n", rmse))
cat(sprintf("R² (model fit): %.2f\n", r2))
cat(sprintf("NSE (overal model performance): %.2f\n", nse))




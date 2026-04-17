###############################################################################
# This scripts creates a plot that shows wind speed vs mean error over all time
# steps and differentiates between night and day and between distances to forest
# edge.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################


library(lubridate)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(tidyr)


RMI_data <- read.csv("Data/RMI_Melle.csv")
model_output <- readRDS("Output/timeseries/data/timeseries_year_yearly_calibrated.rds")


# RMI_data timestamp (UTC)
RMI_data$timestamp <- ymd_hms(RMI_data$timestamp, tz = "UTC")


merged_data <- inner_join(
  RMI_data,
  model_output,
  by = c("timestamp" = "datetime")
)


merged_data$D_edge <- factor(
  merged_data$D_edge,
  levels = sort(unique(merged_data$D_edge))
)


# DEFINE DAY AND NIGHT DATA POINTS
merged_data <- merged_data %>%
  mutate(
    hour_utc = hour(timestamp),
    period = case_when(
      hour_utc >= 5 & hour_utc < 17 ~ "day",
      TRUE ~ "night"
    )
  )

merged_data <- merged_data %>%
  mutate(
    box_x = ifelse(period == "day", -1, -2)
  )%>%
  mutate(
    period = factor(period, levels = c("night", "day"))
  )

# SLOPE PER D_edge
slopes <- merged_data %>%
  group_by(D_edge, period) %>%
  summarise(
    model = list(lm(resid ~ wind_speed)),
    slope = coef(model[[1]])[2],
    p_value = summary(model[[1]])$coefficients[2,4],
    .groups = "drop"
  ) %>%
  mutate(
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    label = paste0(signif, formatC(slope, width = 6, format = "f", digits = 3))
  )

y_min <- min(merged_data$resid, na.rm = TRUE)
y_range <- diff(range(merged_data$resid, na.rm = TRUE))

slopes <- slopes %>%
  mutate(
    y_pos = case_when(
      period == "day"   ~ y_min + 0.05 * y_range,
      period == "night" ~ y_min + 0.15 * y_range
    )
  )


slope_title <- slopes %>%
  group_by(D_edge) %>%
  summarise(
    y_pos = max(y_pos),
    .groups = "drop"
  ) %>%
  mutate(
    label = "slope:"
  )


# PLOT: residual vs wind_speed per D_edge
p <- ggplot(merged_data, aes(x = wind_speed, y = resid)) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linetype = "solid",
    linewidth = 0.6
  ) +
  geom_point(aes(colour = period), alpha = 0.3, size = 0.6) +
  geom_boxplot(
    aes(x = box_x, y = resid, fill = period, colour = period),
    width = 0.4,
    alpha = 0.6,
    outlier.size = 0.6,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_smooth(aes(colour = period), method = "lm", se = FALSE) +
  scale_colour_manual(
    values = c(
      "day" = "orange",
      "night" = "blue"
    )
  ) +
  scale_fill_manual(
    values = c(
      "day" = "orange",
      "night" = "blue"
    )
  ) +
  xlim(-3, max(merged_data$wind_speed, na.rm = TRUE)) +
  geom_text(
    data = slope_title,
    aes(
      x = Inf,
      y = y_pos + 0.1 * y_range,
      label = label
    ),
    colour = "black",
    hjust = 1.1,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = slopes,
    aes(
      x = Inf,
      y = y_pos,
      label = label,
      colour = period
    ),
    hjust = 1.1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ D_edge,
    ncol = 5,
    labeller = as_labeller(function(x) {
      sapply(x, function(val) {
        paste0("D[e]==", val, "*' m'")
      })
    }, label_parsed)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  labs(
    x = "Wind speed (m/s)", # * 3.6 to convert to km/h
    y = "Residual (°C)",
    title = expression("Model–observation residual (mean error) in function of wind speed and distance to edge (" * D[e] * ")" )
  ) +
  guides(
    fill = "none",
    colour = guide_legend(
      override.aes = list(
        shape = NA,
        linetype = 1,
        linewidth = 1.2
      )
    )
  )

print(p)


ggsave(
  filename = "Output/residue_vs_windspeed_per_Dedge.png",
  plot = p,
  width = 9,
  height = 5
)





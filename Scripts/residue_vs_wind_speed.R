###############################################################################
# LOAD LIBRARIES
###############################################################################
library(lubridate)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(tidyr)

###############################################################################
# LOAD DATA
###############################################################################
RMI_data <- read.csv("Data/RMI_Melle.csv")
model_output <- readRDS("Output/timeseries/data/timeseries_year.rds")

###############################################################################
# FORMAT TIME COLUMNS
###############################################################################

# RMI_data timestamp (UTC)
RMI_data$timestamp <- ymd_hms(RMI_data$timestamp, tz = "UTC")

# # model_output datetime (ook UTC zetten!)
# model_output$datetime <- ymd_hms(model_output$datetime, tz = "UTC")

###############################################################################
# MERGE DATA (alleen overlappende timestamps)
###############################################################################

merged_data <- inner_join(
  RMI_data,
  model_output,
  by = c("timestamp" = "datetime")
)

###############################################################################
# OPTIONAL: D_edge netjes sorteren
###############################################################################

merged_data$D_edge <- factor(
  merged_data$D_edge,
  levels = sort(unique(merged_data$D_edge))
)

###############################################################################
# DEFINE DAY AND NIGHT DATA POINTS
###############################################################################

merged_data <- merged_data %>%
  mutate(
    hour_utc = hour(timestamp),
    # period = case_when(
    #   hour_utc >= 5 & hour_utc < 9  ~ "morning",
    #   hour_utc >= 9 & hour_utc < 16 ~ "day",
    #   hour_utc >= 16 & hour_utc < 20 ~ "evening",
    #   TRUE ~ "night"
    # )
    period = case_when(
      hour_utc >= 5 & hour_utc < 17 ~ "day",
      TRUE ~ "night"
    )
  )

merged_data <- merged_data %>%
  mutate(
    box_x = ifelse(period == "day", -1, -2)  # links van je data
  )

###############################################################################
# SLOPE BEREKENEN PER D_edge
###############################################################################
slopes <- merged_data %>%
  group_by(D_edge) %>%
  summarise(
    slope = coef(lm(abs(resid) ~ wind_speed))[2],
    .groups = "drop"
  )

# labels maken voor in de plot
slopes <- slopes %>%
  mutate(label = paste0("slope: ", round(slope, 3)))

###############################################################################
# PLOT: residue vs wind_speed per D_edge (in één figuur)
###############################################################################

p <- ggplot(merged_data, aes(x = wind_speed, y = resid, colour = period)) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linetype = "solid",
    linewidth = 0.6
  ) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_boxplot(
    aes(x = box_x, y = resid, fill = period, colour = period),
    width = 0.4,
    alpha = 0.6,
    outlier.size = 0.6,
    inherit.aes = FALSE
  ) +
  geom_smooth(method = "lm", colour = "red") +
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
    data = slopes,
    aes(x = Inf, y = -Inf, label = label),
    hjust = 1.1,
    vjust = -0.5,
    colour = "red",
    inherit.aes = FALSE
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
  labs(
    x = "Wind speed (m/s)", # * 3.6 to convert to km/h
    y = "Residue (°C)",
    title = expression("Model–observation residue (mean error) in function of wind speed and distance to edge (" * D[e] * ")" )
  )

print(p)


###############################################################################
# OPTIONAL: SAVE FIGURE
###############################################################################

ggsave(
  filename = "Output/residue_vs_windspeed_per_Dedge.png",
  plot = p,
  width = 9,
  height = 5
)





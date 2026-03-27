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
# PLOT: residu vs wind_speed per D_edge (in één figuur)
###############################################################################

p <- ggplot(merged_data, aes(x = wind_speed, y = abs(resid))) +
  geom_point(alpha = 0.1, size = 0.6, colour = "blue") +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(
    data = slopes,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.5,
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
    y = "Residu (°C)",
    title = expression("Model-observation residu (mean absolute error) in function of wind speed and distance to edge (" * D[e] * ")" )
  )

print(p)


###############################################################################
# OPTIONAL: SAVE FIGURE
###############################################################################

ggsave(
  filename = "Output/residu_vs_windspeed_per_Dedge.png",
  plot = p,
  width = 9,
  height = 5
)





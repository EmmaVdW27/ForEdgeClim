############################################################################################
# This script creates two model performance plots (modelled vs observed air temperature gradient and
# observed gradient vs % of cases where error < observed gradient).
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

library(ggplot2)
library(dplyr)
library(ggnewscale)
library(grid)
library(lubridate)
library(ggchicklet)

df <- readRDS("Output/timeseries/data/timeseries_year_yearly_calibrated.rds")

# DATA SET MANIPULATION:

df_summary <- df %>%
  group_by(datetime) %>%
  summarise(
    # observed (2 points per zone)
    T_edge_obs = mean(Tair[D_edge %in% c(0, 15)], na.rm = TRUE), # c(0, 15, 30)
    T_core_obs = mean(Tair[D_edge %in% c(120, 135)], na.rm = TRUE), # c(105, 120, 135)
    gradient_obs = T_edge_obs - T_core_obs,

    # modelled (2 points per zone)
    T_edge_mod = mean(temperature[D_edge %in% c(0, 15)], na.rm = TRUE),
    T_core_mod = mean(temperature[D_edge %in% c(120, 135)], na.rm = TRUE),
    gradient_mod = T_edge_mod - T_core_mod,

    # error
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    gradient_AME = abs(gradient_mod - gradient_obs)
  ) %>%
  mutate(hour = as.numeric(format(datetime, "%H")),
         month = as.numeric(format(datetime, "%m"))) %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )


# R2 VALUES:

R2_year <- df_summary %>%
  summarise(
    R2 = cor(gradient_obs, gradient_mod, use = "complete.obs")^2
  )

print(R2_year)

R2_season <- df_summary %>%
  group_by(season) %>%
  summarise(
    R2 = cor(gradient_obs, gradient_mod, use = "complete.obs")^2
  )

print(R2_season)


# PLOT: MODDELED VS OBSERVED GRADIENT
# -----------------------------------

season_cols <- c(
  winter = "brown",
  spring = "yellowgreen",
  summer = "green4",
  autumn = "orange"
)

# R² labels with position
R2_labels <- R2_season %>%
  mutate(season = factor(season,
                         levels = c("summer", "autumn", "spring", "winter"))) %>%
  arrange(season) %>%
  mutate(
    label = paste0(season, ": R² = ", round(R2, 2)),
    x = max(df_summary$gradient_obs, na.rm = TRUE),
    y = max(df_summary$gradient_mod, na.rm = TRUE) - 0.2 -(row_number() - 1) * 0.55
  )

# plot
gradient_plot <- ggplot(df_summary,
                        aes(x = gradient_obs,
                            y = gradient_mod)) +

  geom_point(aes(colour = month),
             alpha = 0.6,
             size = 0.5) +

  scale_colour_gradientn(
    colours = c("brown", "limegreen", "limegreen",
                "green4", "green4",
                "orange", "orange", "brown"),
    limits = c(1, 12),
    breaks = seq(1, 12, by = 1),
    guide = guide_colourbar(
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black"
    )
  ) +

  ggnewscale::new_scale_colour() +

  # regression line per season
  geom_smooth(aes(colour = season),
              method = "lm",
              se = FALSE,
              linewidth = 0.8) +

  scale_colour_manual(values = season_cols, guide = "none") +

  # 1:1 line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +

  annotate("text",
           x = -1.2, y = -1.2,
           label = "1:1",
           angle = 45,
           hjust = -0.5,
           vjust = -0.5,
           size = 5) +

  # R² labels
  geom_label(data = R2_labels,
             aes(x = x, y = y, label = label, fill = season),
             colour = "white",
             size = 5,
             label.size = 0.3,
             label.r = unit(0.2, "lines"),
             hjust = 1,
             nudge_x = -0.2,
             nudge_y = 0.4,
             inherit.aes = FALSE,
             show.legend = FALSE) +

  scale_fill_manual(values = season_cols) +

  scale_x_continuous(breaks = seq(-2, 6, by = 1)) +
  scale_y_continuous(breaks = seq(-1, 3, by = 1)) +

  labs(
    x = "Observed gradient (°C)",
    y = "Modelled gradient (°C)",
    colour = "month",
    title = "a) Modelled vs observed air temperature gradient"
  ) +

  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

print(gradient_plot)

ggsave(
  "Output/timeseries/plots/Observed_vs_modelled_gradient.png",
  plot = gradient_plot,
  width = 6,
  height = 3,
  dpi = 300
)




# PLOT: BINNED PERFORMANCE
# ------------------------

df_binned <- df_summary %>%
  mutate(
    abs_grad = abs(gradient_obs),
    rel_error = gradient_AME / abs_grad,
    below = rel_error < 1,

    warm = factor(T_edge_obs > 20, levels = c(FALSE, TRUE)),

    # bins
    grad_bin = cut(abs_grad,
                   breaks = seq(0, 8, by = 1),
                   include.lowest = TRUE,
                   right = TRUE)
  ) %>%

  filter(!is.na(grad_bin)) %>%

  mutate(
    grad_bin = gsub("\\(", "[", grad_bin),
    grad_bin = gsub("\\)", "]", grad_bin)
  ) %>%

  group_by(grad_bin, warm) %>%
  summarise(
    perc_below = 100 * mean(below, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(grad_bin, warm, fill = list(perc_below = NA)) %>%

  ggplot(aes(x = grad_bin, y = perc_below, fill = warm)) +
  geom_chicklet(position = position_dodge2(preserve = "single"),
                radius = grid::unit(3, "pt")) +

  scale_fill_manual(
    values = c("blue", "orange"),
    labels = c("< 20°C", "> 20°C"),
    name = "Edge temperature"
  ) +

  labs(
    x = "Observed gradient (°C)",
    y = "% cases where error < gradient",
    title = "b) Model ability to resolve spatial temperature gradients"
  ) +

  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

print(df_binned)

ggsave(
  "Output/timeseries/plots/binned_performance.png",
  plot = df_binned,
  width = 6,
  height = 3,
  dpi = 300
)

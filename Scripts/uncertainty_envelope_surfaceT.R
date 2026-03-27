###############################################################################
# UNCERTAINTY ENVELOPE FOR FOREST SURFACE TEMPERATURE
#
# This script extracts model outputs from Sobol runs (soboljansen objects)
# and computes uncertainty envelopes (Bell curves) for forest surface
# temperature across parameter ensembles.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(glue)
library(forcats)
library(ggh4x)

#########
# INPUT #
#########

input_path  <- "Output/sensitivity_analysis/Sobol_QoI/data/"
output_path <- "Output/sensitivity_analysis/Sobol_QoI/plots_output/"

direction <- "h"   # 'h' or 'v'
N <- 400           # number of Sobol base samples

# File list (forest surface temperature version!)
files_info <- tribble(
  ~file, ~label,

  #winter files
  # night
  glue("400samples_25parameters_01h_13012025_avTf_{direction}.rds"),     "Wi_Ni_≈T",
  glue("400samples_25parameters_01h_13012025_SDTf_{direction}.rds"),     "Wi_Ni_σT",
  glue("400samples_25parameters_01h_13012025_gradTf_{direction}.rds"),   "Wi_Ni_∇T",
  # morning
  glue("400samples_25parameters_08h_13012025_avTf_{direction}.rds"),     "Wi_Mo_≈T",
  glue("400samples_25parameters_08h_13012025_SDTf_{direction}.rds"),     "Wi_Mo_σT",
  glue("400samples_25parameters_08h_13012025_gradTf_{direction}.rds"),   "Wi_Mo_∇T",
  # noon
  glue("400samples_25parameters_12h_13012025_avTf_{direction}.rds"),     "Wi_No_≈T",
  glue("400samples_25parameters_12h_13012025_SDTf_{direction}.rds"),     "Wi_No_σT",
  glue("400samples_25parameters_12h_13012025_gradTf_{direction}.rds"),   "Wi_No_∇T",
  # spring files
  # night
  glue("400samples_25parameters_01h_30042025_avTf_{direction}.rds"),     "Sp_Ni_≈T",
  glue("400samples_25parameters_01h_30042025_SDTf_{direction}.rds"),     "Sp_Ni_σT",
  glue("400samples_25parameters_01h_30042025_gradTf_{direction}.rds"),   "Sp_Ni_∇T",
  # morning
  glue("400samples_25parameters_05h_30042025_avTf_{direction}.rds"),     "Sp_Mo_≈T",
  glue("400samples_25parameters_05h_30042025_SDTf_{direction}.rds"),     "Sp_Mo_σT",
  glue("400samples_25parameters_05h_30042025_gradTf_{direction}.rds"),   "Sp_Mo_∇T",
  # noon
  glue("400samples_25parameters_12h_30042025_avTf_{direction}.rds"),     "Sp_No_≈T",
  glue("400samples_25parameters_12h_30042025_SDTf_{direction}.rds"),     "Sp_No_σT",
  glue("400samples_25parameters_12h_30042025_gradTf_{direction}.rds"),   "Sp_No_∇T",
  # summer files
  # night
  glue("400samples_25parameters_01h_07072023_avTf_{direction}.rds"),     "Su_Ni_≈T",
  glue("400samples_25parameters_01h_07072023_SDTf_{direction}.rds"),     "Su_Ni_σT",
  glue("400samples_25parameters_01h_07072023_gradTf_{direction}.rds"),   "Su_Ni_∇T",
  # morning
  glue("400samples_25parameters_05h_07072023_avTf_{direction}.rds"),     "Su_Mo_≈T",
  glue("400samples_25parameters_05h_07072023_SDTf_{direction}.rds"),     "Su_Mo_σT",
  glue("400samples_25parameters_05h_07072023_gradTf_{direction}.rds"),   "Su_Mo_∇T",
  # noon
  glue("400samples_25parameters_12h_07072023_avTf_{direction}.rds"),     "Su_No_≈T",
  glue("400samples_25parameters_12h_07072023_SDTf_{direction}.rds"),     "Su_No_σT",
  glue("400samples_25parameters_12h_07072023_gradTf_{direction}.rds"),   "Su_No_∇T",
  # autumn files
  # night
  glue("400samples_25parameters_01h_01102023_avTf_{direction}.rds"),     "Au_Ni_≈T",
  glue("400samples_25parameters_01h_01102023_SDTf_{direction}.rds"),     "Au_Ni_σT",
  glue("400samples_25parameters_01h_01102023_gradTf_{direction}.rds"),   "Au_Ni_∇T",
  # morning
  glue("400samples_25parameters_06h_01102023_avTf_{direction}.rds"),     "Au_Mo_≈T",
  glue("400samples_25parameters_06h_01102023_SDTf_{direction}.rds"),     "Au_Mo_σT",
  glue("400samples_25parameters_06h_01102023_gradTf_{direction}.rds"),   "Au_Mo_∇T",
  # noon
  glue("400samples_25parameters_12h_01102023_avTf_{direction}.rds"),     "Au_No_≈T",
  glue("400samples_25parameters_12h_01102023_SDTf_{direction}.rds"),     "Au_No_σT",
  glue("400samples_25parameters_12h_01102023_gradTf_{direction}.rds"),   "Au_No_∇T"
)

##########
# LOAD OUTPUTS
##########

output_df <- files_info %>%
  mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
  mutate(Y = map(data, ~ .x$y[1:(2*N)])) %>%
  select(label, Y) %>%
  unnest(Y)

##########
# ADD METADATA
##########

output_df <- output_df %>%
  mutate(
    season = substr(label, 1, 2),
    moment = substr(label, 4, 5),
    metric = case_when(
      grepl("≈T", label) ~ "Mean temperature",
      grepl("σT", label) ~ "Temperature SD",
      grepl("∇T", label) ~ "Temperature gradient"
    )
  )

output_df <- output_df %>%
  mutate(
    season = recode(season,
                    "Wi" = "Winter",
                    "Sp" = "Spring",
                    "Su" = "Summer",
                    "Au" = "Autumn"
    ),
    moment = recode(moment,
                    "Ni" = "Night",
                    "Mo" = "Morning",
                    "No" = "Noon"
    )
  )

output_df <- output_df %>%
  mutate(
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    moment = factor(moment, levels = c("Night", "Morning", "Noon"))
  )

##########
# MEDIAN PER PANEL (for vertical line)
##########

median_df <- output_df %>%
  group_by(season, moment, metric) %>%
  summarise(median = median(Y), .groups = "drop")

##########
# FACET DENSITY PLOT
##########

p <- ggplot(output_df, aes(x = Y)) +

  geom_density(
    aes(y = after_stat(scaled)),
    fill = "steelblue",
    alpha = 0.7
  ) +

  geom_vline(
    data = median_df,
    aes(xintercept = median),
    colour = "black",
    size = 0.6
  ) +

  ggh4x::facet_nested(
    rows = vars(season, moment),
    cols = vars(metric),
    scales = "free_x"
  ) +

  scale_y_continuous(
    breaks = c(0, 1),
    limits = c(0, 1)
  ) +

  theme_bw(base_size = 14) +

  labs(
    title = "Uncertainty of forest surface temperature",
    x = "Forest surface temperature metric (°C)",
    y = "Scaled density"
  ) +

  theme(
    strip.text = element_text(size = 13),
    ggh4x.facet.nestline = element_line(size = 1),
    strip.text.y = element_text(angle = 0),
    axis.title = element_text(size = 16)       # x- en y-as titels
  )

print(p)


##########
# SAVE
##########

ggsave(
  paste0(output_path, "Density_surfaceT_facet_", direction, ".png"),
  plot = p,
  width = 8,
  height = 8,
  dpi = 300
)

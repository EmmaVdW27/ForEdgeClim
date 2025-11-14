############################################################################################
# This script visualizes the parameter optimization process, using the CMA-ES algorithm (cmaesr).
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(ggh4x)
  library(grid)
})

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------
season = "spring" # winter spring summer or autumn

number_of_params = 3

# -----------------
# PLOTS — Convergence visuals
# -----------------

offspring_log <- read.csv(paste0("Output/calibration/data/logging_offspring_", season, ".csv"))
generation_log <- read.csv(paste0("Output/calibration/data/logging_generation_", season, ".csv"))
res <- readRDS(paste0("Output/calibration/data/calibration_results_", season, ".rds"))
output_path <- "Output/calibration/plots/"
lambda <- 4 + floor(3 * log(number_of_params))


# PLOT 1: Offspring RMSE between observed and modelled temperature
# ------

# Best offspring per generation + jitter-x
best_points <- offspring_log %>%
  group_by(iter) %>%
  slice_min(rmse, n = 1, with_ties = FALSE) %>%  # pick 1 value
  ungroup() %>%
  mutate(jitter_x = iter + (offspring_id - 1) / lambda)

p_RMSE <- ggplot() +
  # All offspring
  geom_point(data = offspring_log,
             aes(x = iter + (offspring_id - 1) / lambda, y = rmse, colour = "All offsprings")) +
  geom_smooth(data = offspring_log,
              aes(x = iter + (offspring_id - 1) / lambda, y = rmse, colour = "All offsprings"),
              se = FALSE, method = "loess", formula = y ~ x) +

  # Best offspring per generation is highlighted
  geom_point(data = best_points,
             aes(x = jitter_x, y = rmse, colour = "Best offspring per generation"),
             size = 2) +
  geom_line(data = best_points,
            aes(x = jitter_x, y = rmse, colour = "Best offspring per generation"),
            linewidth = 1, linetype = "dashed") +

  geom_vline(xintercept = unique(offspring_log$iter),
             linetype = "dashed", colour = "orange") +

  scale_colour_manual(values = c("All offsprings" = "blue",
                                 "Best offspring per generation" = alpha("red", 0.5))) +
  # Captions
  labs(x = "Generation",
       y = "RMSE (°C)",
       title = "a) RMSE evolution over generations",
       colour = NULL) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40)
  )

print(p_RMSE)
ggsave(paste0(output_path, "RMSE_", season, ".png"), plot = p_RMSE, width = 16, height = 10, dpi = 500)




# PLOT 2: 2D plot in PCA space
# ------

param_cols <- setdiff(names(offspring_log), c("eval_id","iter","offspring_id","rmse","time"))

# Data matrix (all generations together)
X <- offspring_log %>% select(all_of(param_cols)) %>% as.matrix()
keep <- complete.cases(X) # only keep complete parameter setes
Xc  <- scale(X[keep, ], center = TRUE, scale = TRUE)  # standardize/rescale parameters: mu = 0, variance = 1

iters_all <- offspring_log$iter[keep]

# PCA fit on all data (so axes stay consistent) (quite similar to PCs from CMA-ES algorithm: PCA on final population would be similar to final PCA in CMA-ES)
pca <- prcomp(Xc, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:2]) # scores of the 2 main PCs
names(scores) <- c("PC1", "PC2")
scores$iter <- iters_all

# Facet plot per generation
p_PCA_space <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "blue") +
  stat_ellipse(type = "norm", level = 0.95, colour = "red", linewidth = 1) +
  coord_equal() +
  facet_wrap(~iter, ncol = 8) +
  labs(
    title = "b) PCA on parameter space",
    #subtitle = "Points = offspring,\nEllipse = 95% spread",
    x = "PC1", y = "\nPC2"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    strip.text = element_text(size = 40) # facet labels
  ) +
  scale_x_continuous(labels = function(x) rep(" ", length(x))) +
  scale_y_continuous(labels = function(x) rep(" ", length(x)))


print(p_PCA_space)
ggsave(paste0(output_path, "PCA_space_", season, ".png"), plot = p_PCA_space, width = 16, height = 14, dpi = 500)

#
# # PLOT 3: Global step size over generations
# # ------
#
# p_step_size <- generation_log %>%
#   ggplot(aes(x = iter, y = sigma_proxy)) +
#   geom_line(colour = "cornflowerblue") + geom_point(colour = "cornflowerblue") +
#   scale_x_continuous(breaks = seq(min(generation_log$iter), max(generation_log$iter), by = 1)) +
#   labs(x = "Generation", y = "Step size (unitless)", title = "Global step size over generations",
#        subtitle = "(Mean parameter SD is used as a proxy for global step size.)") +
#   theme_bw(base_size = 12)
# print(p_step_size)
# ggsave(paste0(output_path, "step_size_", season, ".png"), plot = p_step_size, width = 10, height = 6, dpi = 300)
#
#
# # PLOT 4: Population spread over generations
# # ------
#
# p_population_spread <- generation_log %>%
#   ggplot(aes(x = iter, y = spread_mean)) +
#   geom_line(colour = "cornflowerblue") + geom_point(colour = "cornflowerblue") +
#   scale_x_continuous(breaks = seq(min(generation_log$iter), max(generation_log$iter), by = 1)) +
#   labs(x = "Generation", y = "Population spread (unitless)",
#        title = "Population spread over generations",
#        subtitle = "(Square root over mean eigenvalue of covariance matrix is used as a proxy for population spread.\n An eigenvalue is the weight of a search direction within the parameter space.)") +
#   theme_bw(base_size = 12)
# print(p_population_spread)
# ggsave(paste0(output_path, "population_spread_", season, ".png"), plot = p_population_spread, width = 10, height = 6, dpi = 300)
#
#
# # PLOT 5: Change in population center between generations
# # ------
# generation_log_plot <- generation_log %>%
#   mutate(iter_mid = iter - 0.5)  # zet het punt tussen huidige en vorige generatie
#
# p_center_change <- generation_log_plot %>%
#   ggplot(aes(x = iter_mid, y = mean_step)) +
#   geom_line(colour = "cornflowerblue") + geom_point(colour = "cornflowerblue") +
#   scale_x_continuous(breaks = seq(min(generation_log$iter), max(generation_log$iter), by = 1)) +
#   labs(x = "Generation", y = "Change", title = "Change in population center from one generation to the next",
#        subtitle = "(Euclidean distance between consecutive generation means is used as a proxy for this change.)") +
#   theme_bw(base_size = 12)
# print(p_center_change)
# ggsave(paste0(output_path, "center_change_", season, ".png"), plot = p_center_change, width = 10, height = 6, dpi = 300)




end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f s\n', as.numeric(end_analysis - start_analysis, units = 'secs')))


# Console summaries
cat("\nOptimised parameters (final):\n")
print(res$best.param)
cat("\nMinimal RMSE observed:\n")
print(res$best.fitness)

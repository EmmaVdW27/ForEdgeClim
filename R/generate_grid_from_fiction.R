#' Generate a virtual sctructural grid with trees
#'
#' @param season Season as string: "summer", "winter", "spring", "autumn"
#' @param output_folder Map to save 3D-plot (Default "Output/")
#' @param phi Viewing angle for 3D-plot (Default -80)
#' @param trunk_size Stem size in plot (Default 1)
#' @param canopy_size Canopy/leave size in plot (Default 1)
#' @return A dataframe with the generated voxel structure
#' @importFrom dplyr mutate bind_rows
#' @export
generate_fictive_grid <- function(season) {

  # Typical canopy height range for a tree (global variables)
  canopy_min_height <- 10
  canopy_max_height <- 18

  set.seed(123)

  # Grid dimensies
  x_dim <- 135  # Length of transect in x-direction (m)
  y_dim <- 30   # Width of transect in y-direction (m)
  z_dim <- 30   # Maximum height of trees (m)

  # Create an empty voxel grid with density set to 0
  voxel_grid <- expand.grid(X = 1:x_dim, Y = 1:y_dim, Z = 1:z_dim) |>
    dplyr::mutate(density = 0)

  # Function to add a single tree or shrub
  add_single_tree <- function(grid, x_pos, y_pos, trunk_height, canopy_max, season, is_shrub = FALSE) {
    if (!is_shrub) {
      # Add trunk
      canopy_min <- trunk_height
      grid <- grid |>
        dplyr::mutate(
          distance_from_trunk = sqrt((X - x_pos)^2 + (Y - y_pos)^2),
          density = ifelse(distance_from_trunk < 1 & Z < trunk_height, runif(nrow(grid), min = 0.8, max = 1), density)
        )
    } else {
      canopy_min <- 0  # Shrubs start from ground level
      canopy_max <- 0
    }

    # Seasonal adjustments for canopy density
    density_range <- switch(season,
                            "summer" = c(0.6, 0.9),
                            "winter" = c(0.3, 0.6),  # Sparser leaves in winter
                            "spring" = c(0.4, 0.7),
                            "autumn" = c(0.4, 0.7))

    # # Sparse canopy in winter, connected to trunk
    # if (season == "winter") {
    #   canopy_min <- trunk_height
    #   canopy_max <- trunk_height + 6  # Keep small canopy near the top of the trunk
    # }

    # Add canopy leaves
    horizontal_radius <- ifelse(is_shrub, 4, 5)
    grid <- grid |>
      dplyr::mutate(
        horizontal_distance = sqrt((X - x_pos)^2 + (Y - y_pos)^2),
        random_fringe_factor = runif(nrow(grid)),
        density = ifelse(horizontal_distance < horizontal_radius & random_fringe_factor > 0.8 & Z >= canopy_min & Z <= canopy_max,
                         runif(nrow(grid), min = density_range[1], max = density_range[2]), density)
      )

    return(grid)
  }

  # Add multiple trees and shrubs at random positions
  num_trees <- 60  # Default number of trees
  tree_positions <- data.frame(
    x_pos = sample(1:x_dim, num_trees, replace = TRUE),
    y_pos = sample(1:y_dim, num_trees, replace = TRUE),
    trunk_height = sample(5:15, num_trees, replace = TRUE),
    canopy_max = sample(20:30, num_trees, replace = TRUE)
  )

  # Add extra trees near the edge in summer
  if (season == "summer") {
    edge_trees <- data.frame(
      x_pos = sample((x_dim * 0.8):x_dim, 20, replace = TRUE),
      y_pos = sample(1:y_dim, 20, replace = TRUE),
      trunk_height = sample(5:15, 20, replace = TRUE),
      canopy_max = sample(20:30, 20, replace = TRUE)
    )
    tree_positions <- dplyr::bind_rows(tree_positions, edge_trees)
  }

  for (i in 1:nrow(tree_positions)) {
    x_pos <- tree_positions$x_pos[i]
    y_pos <- tree_positions$y_pos[i]

    # Determine if this is a shrub
    is_shrub <- (x_pos >= 30 & x_pos <= 50)

    voxel_grid <- add_single_tree(
      voxel_grid,
      x_pos = x_pos,
      y_pos = y_pos,
      trunk_height = tree_positions$trunk_height[i],
      canopy_max = tree_positions$canopy_max[i],
      season = season,
      is_shrub = is_shrub
    )
  }

  # # Add a ground layer (z = 0) with high density
  # ground_layer <- expand.grid(X = 1:x_dim, Y = 1:y_dim, Z = 0) %>%
  #   mutate(density = 1e6)
  #
  # # Combine ground layer with voxel grid
  # voxel_grid <- bind_rows(ground_layer, voxel_grid) %>%
  #   arrange(X, Y, Z)
  return(voxel_grid)
}



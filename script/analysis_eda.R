############################################################
# Position-related neuronal activity and arena boundaries
############################################################

rm(list = ls())

setwd("C:/Users/pietr/OneDrive/Desktop/PhD/Activities/Data Research Camp 2026/Project/Data")

library(ggplot2)
library(patchwork)

all_trials <- readRDS("all_trials.RDS")

str(all_trials, max.level = 3)


############################################################
# 1. Settings and trial mapping
############################################################

percentile_used <- 0.90
number_distance_bins <- 20

trial_key <- data.frame(
  trial = paste0("trial", 1:6),
  day = c(1, 1, 1, 2, 2, 2),
  shape = c(
    "triangle",
    "circle",
    "square",
    "circle",
    "square",
    "triangle"
  ),
  stringsAsFactors = FALSE
)

shape_levels <- c(
  "triangle",
  "circle",
  "square"
)

shape_labels <- c(
  "Triangle",
  "Circle",
  "Square"
)

day_levels <- c(
  "Day 1",
  "Day 2"
)

mouse_names <- names(all_trials)


############################################################
# 2. Basic functions
############################################################

get_trial_name <- function(
    shape,
    day
) {
  
  trial_name <- trial_key$trial[
    trial_key$shape == tolower(shape) &
      trial_key$day == day
  ]
  
  if (length(trial_name) != 1) {
    stop(
      paste(
        "Could not identify the trial for",
        shape,
        "on Day",
        day
      )
    )
  }
  
  trial_name
}


binarize_S <- function(
    S,
    percentile = 0.90
) {
  
  S <- as.matrix(S)
  
  thresholds <- apply(
    S,
    1,
    quantile,
    probs = percentile,
    na.rm = TRUE,
    names = FALSE
  )
  
  binary_S <- sweep(
    S,
    1,
    thresholds,
    FUN = ">"
  )
  
  storage.mode(binary_S) <- "integer"
  binary_S[is.na(binary_S)] <- 0L
  
  binary_S
}


set_condition_factors <- function(data) {
  
  if ("shape" %in% names(data)) {
    
    data$shape <- factor(
      data$shape,
      levels = shape_levels,
      labels = shape_labels
    )
  }
  
  if ("day" %in% names(data)) {
    
    data$day <- factor(
      data$day,
      levels = day_levels
    )
  }
  
  data
}


############################################################
# 3. Collect behavioural positions
############################################################

collect_shape_positions <- function(
    all_trials,
    shape
) {
  
  trials <- trial_key$trial[
    trial_key$shape == shape
  ]
  
  result <- list()
  k <- 1
  
  for (mouse in names(all_trials)) {
    
    for (trial in trials) {
      
      position <-
        all_trials[[mouse]][[trial]]$behav$position
      
      if (is.null(position)) {
        next
      }
      
      position <- as.matrix(position)
      
      valid <- is.finite(position[, 1]) &
        is.finite(position[, 2])
      
      position <- position[
        valid,
        1:2,
        drop = FALSE
      ]
      
      if (nrow(position) == 0) {
        next
      }
      
      result[[k]] <- data.frame(
        x = position[, 1],
        y = position[, 2],
        mouse = mouse,
        trial = trial,
        shape = shape
      )
      
      k <- k + 1
    }
  }
  
  do.call(
    rbind,
    result
  )
}


############################################################
# 4. Estimate arena boundaries
############################################################

estimate_circle_boundary <- function(
    positions,
    boundary_quantile = 0.995
) {
  
  centre <- c(
    x = median(
      positions$x,
      na.rm = TRUE
    ),
    y = median(
      positions$y,
      na.rm = TRUE
    )
  )
  
  radial_distance <- sqrt(
    (positions$x - centre["x"])^2 +
      (positions$y - centre["y"])^2
  )
  
  list(
    type = "circle",
    centre = centre,
    radius = as.numeric(
      quantile(
        radial_distance,
        boundary_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    )
  )
}


estimate_square_boundary <- function(
    positions,
    lower_quantile = 0.005,
    upper_quantile = 0.995
) {
  
  list(
    type = "square",
    
    xmin = as.numeric(
      quantile(
        positions$x,
        lower_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    
    xmax = as.numeric(
      quantile(
        positions$x,
        upper_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    
    ymin = as.numeric(
      quantile(
        positions$y,
        lower_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    
    ymax = as.numeric(
      quantile(
        positions$y,
        upper_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    )
  )
}


triangle_area <- function(
    p1,
    p2,
    p3
) {
  
  abs(
    p1[1] * (p2[2] - p3[2]) +
      p2[1] * (p3[2] - p1[2]) +
      p3[1] * (p1[2] - p2[2])
  ) / 2
}


estimate_triangle_boundary <- function(
    positions,
    max_hull_points = 100
) {
  
  xy <- unique(
    positions[
      ,
      c("x", "y")
    ]
  )
  
  hull <- as.matrix(
    xy[
      chull(
        xy$x,
        xy$y
      ),
      ,
      drop = FALSE
    ]
  )
  
  if (nrow(hull) > max_hull_points) {
    
    selected <- unique(
      round(
        seq(
          1,
          nrow(hull),
          length.out = max_hull_points
        )
      )
    )
    
    hull <- hull[
      selected,
      ,
      drop = FALSE
    ]
  }
  
  combinations <- combn(
    seq_len(nrow(hull)),
    3
  )
  
  areas <- apply(
    combinations,
    2,
    function(index) {
      
      triangle_area(
        hull[index[1], ],
        hull[index[2], ],
        hull[index[3], ]
      )
    }
  )
  
  vertices <- hull[
    combinations[, which.max(areas)],
    ,
    drop = FALSE
  ]
  
  colnames(vertices) <- c(
    "x",
    "y"
  )
  
  list(
    type = "triangle",
    vertices = vertices
  )
}


estimate_arena_definitions <- function(
    all_trials
) {
  
  circle_positions <- collect_shape_positions(
    all_trials,
    "circle"
  )
  
  square_positions <- collect_shape_positions(
    all_trials,
    "square"
  )
  
  triangle_positions <- collect_shape_positions(
    all_trials,
    "triangle"
  )
  
  list(
    circle = estimate_circle_boundary(
      circle_positions
    ),
    
    square = estimate_square_boundary(
      square_positions
    ),
    
    triangle = estimate_triangle_boundary(
      triangle_positions
    )
  )
}


arena_definitions <- estimate_arena_definitions(
  all_trials
)

#print(arena_definitions)


############################################################
# 5. Convert boundaries to plotting data
############################################################

arena_boundary_dataframe <- function(
    arena_definitions,
    number_circle_points = 300
) {
  
  angles <- seq(
    0,
    2 * pi,
    length.out = number_circle_points
  )
  
  circle <- data.frame(
    x = arena_definitions$circle$centre["x"] +
      arena_definitions$circle$radius *
      cos(angles),
    
    y = arena_definitions$circle$centre["y"] +
      arena_definitions$circle$radius *
      sin(angles),
    
    shape = "circle"
  )
  
  square <- data.frame(
    x = c(
      arena_definitions$square$xmin,
      arena_definitions$square$xmax,
      arena_definitions$square$xmax,
      arena_definitions$square$xmin,
      arena_definitions$square$xmin
    ),
    
    y = c(
      arena_definitions$square$ymin,
      arena_definitions$square$ymin,
      arena_definitions$square$ymax,
      arena_definitions$square$ymax,
      arena_definitions$square$ymin
    ),
    
    shape = "square"
  )
  
  triangle_vertices <- rbind(
    arena_definitions$triangle$vertices,
    arena_definitions$triangle$vertices[1, ]
  )
  
  triangle <- data.frame(
    x = triangle_vertices[, 1],
    y = triangle_vertices[, 2],
    shape = "triangle"
  )
  
  set_condition_factors(
    rbind(
      circle,
      square,
      triangle
    )
  )
}


boundary_data <- arena_boundary_dataframe(
  arena_definitions
)


############################################################
# 6. Validate estimated arena boundaries
############################################################

plot_estimated_arenas <- function(
    all_trials,
    arena_definitions
) {
  
  position_data <- do.call(
    rbind,
    lapply(
      shape_levels,
      function(shape) {
        
        collect_shape_positions(
          all_trials,
          shape
        )
      }
    )
  )
  
  position_data <- set_condition_factors(
    position_data
  )
  
  boundary_data <- arena_boundary_dataframe(
    arena_definitions
  )
  
  ggplot() +
    geom_point(
      data = position_data,
      aes(
        x = x,
        y = y
      ),
      size = 0.15,
      alpha = 0.08
    ) +
    geom_path(
      data = boundary_data,
      aes(
        x = x,
        y = y,
        group = shape
      ),
      linewidth = 0.9
    ) +
    facet_wrap(
      vars(shape)
    ) +
    coord_fixed() +
    labs(
      title = "Arena boundaries estimated from mouse positions",
      subtitle = paste(
        "Observed positions and estimated arena boundaries"
      ),
      x = "Position x",
      y = "Position y"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


p_estimated_arenas <- plot_estimated_arenas(
  all_trials,
  arena_definitions
)

#print(p_estimated_arenas)


############################################################
# 7. Align behavioural positions with S frames
############################################################

get_S_frame_times <- function(
    neuron_time,
    number_S_frames,
    behaviour_time
) {
  
  neuron_time <- as.numeric(neuron_time)
  behaviour_time <- as.numeric(behaviour_time)
  
  neuron_time <- neuron_time[
    is.finite(neuron_time)
  ]
  
  behaviour_time <- behaviour_time[
    is.finite(behaviour_time)
  ]
  
  if (length(neuron_time) >= 2) {
    
    if (length(neuron_time) == number_S_frames) {
      return(neuron_time)
    }
    
    return(
      seq(
        min(neuron_time),
        max(neuron_time),
        length.out = number_S_frames
      )
    )
  }
  
  seq(
    min(behaviour_time),
    max(behaviour_time),
    length.out = number_S_frames
  )
}


align_position_to_S <- function(
    behaviour_position,
    behaviour_time,
    S_time
) {
  
  behaviour_position <- as.matrix(
    behaviour_position
  )
  
  behaviour_time <- as.numeric(
    behaviour_time
  )
  
  valid <- is.finite(behaviour_time) &
    is.finite(behaviour_position[, 1]) &
    is.finite(behaviour_position[, 2])
  
  behaviour_time <- behaviour_time[
    valid
  ]
  
  behaviour_position <- behaviour_position[
    valid,
    1:2,
    drop = FALSE
  ]
  
  ordering <- order(
    behaviour_time
  )
  
  behaviour_time <- behaviour_time[
    ordering
  ]
  
  behaviour_position <- behaviour_position[
    ordering,
    ,
    drop = FALSE
  ]
  
  keep <- !duplicated(
    behaviour_time
  )
  
  behaviour_time <- behaviour_time[
    keep
  ]
  
  behaviour_position <- behaviour_position[
    keep,
    ,
    drop = FALSE
  ]
  
  data.frame(
    frame = seq_along(S_time),
    time = S_time,
    
    x = approx(
      behaviour_time,
      behaviour_position[, 1],
      xout = S_time,
      rule = 1
    )$y,
    
    y = approx(
      behaviour_time,
      behaviour_position[, 2],
      xout = S_time,
      rule = 1
    )$y
  )
}


############################################################
# 8. Boundary-distance functions
############################################################

distance_to_segment <- function(
    x,
    y,
    x1,
    y1,
    x2,
    y2
) {
  
  dx <- x2 - x1
  dy <- y2 - y1
  
  segment_length_squared <- dx^2 + dy^2
  
  projection <- (
    (x - x1) * dx +
      (y - y1) * dy
  ) / segment_length_squared
  
  projection <- pmax(
    0,
    pmin(
      1,
      projection
    )
  )
  
  closest_x <- x1 + projection * dx
  closest_y <- y1 + projection * dy
  
  sqrt(
    (x - closest_x)^2 +
      (y - closest_y)^2
  )
}


calculate_boundary_distance <- function(
    x,
    y,
    shape,
    arena_definitions
) {
  
  if (shape == "circle") {
    
    arena <- arena_definitions$circle
    
    radial_distance <- sqrt(
      (x - arena$centre["x"])^2 +
        (y - arena$centre["y"])^2
    )
    
    return(
      abs(
        arena$radius -
          radial_distance
      )
    )
  }
  
  if (shape == "square") {
    
    arena <- arena_definitions$square
    
    return(
      pmin(
        abs(x - arena$xmin),
        abs(arena$xmax - x),
        abs(y - arena$ymin),
        abs(arena$ymax - y)
      )
    )
  }
  
  arena <- arena_definitions$triangle
  vertices <- arena$vertices
  
  d12 <- distance_to_segment(
    x,
    y,
    vertices[1, 1],
    vertices[1, 2],
    vertices[2, 1],
    vertices[2, 2]
  )
  
  d23 <- distance_to_segment(
    x,
    y,
    vertices[2, 1],
    vertices[2, 2],
    vertices[3, 1],
    vertices[3, 2]
  )
  
  d31 <- distance_to_segment(
    x,
    y,
    vertices[3, 1],
    vertices[3, 2],
    vertices[1, 1],
    vertices[1, 2]
  )
  
  pmin(
    d12,
    d23,
    d31
  )
}


############################################################
# 9. Prepare aligned position/activity data
############################################################

prepare_position_activity <- function(
    all_trials,
    mouse,
    shape,
    day,
    percentile,
    arena_definitions
) {
  
  trial_name <- get_trial_name(
    shape,
    day
  )
  
  trial_data <- all_trials[[mouse]][[trial_name]]
  
  binary_S <- binarize_S(
    trial_data$neuron$S,
    percentile
  )
  
  S_time <- get_S_frame_times(
    neuron_time = trial_data$neuron$time,
    number_S_frames = ncol(binary_S),
    behaviour_time = trial_data$behav$time
  )
  
  data <- align_position_to_S(
    behaviour_position = trial_data$behav$position,
    behaviour_time = trial_data$behav$time,
    S_time = S_time
  )
  
  data$boundary_distance <-
    calculate_boundary_distance(
      data$x,
      data$y,
      shape,
      arena_definitions
    )
  
  data$number_active <- colSums(
    binary_S,
    na.rm = TRUE
  )
  
  data$number_neurons <- nrow(
    binary_S
  )
  
  data$proportion_active <-
    data$number_active /
    data$number_neurons
  
  data$mouse <- mouse
  data$trial <- trial_name
  data$shape <- shape
  data$day <- paste(
    "Day",
    day
  )
  
  data
}


prepare_mouse_position_activity <- function(
    all_trials,
    mouse,
    percentile,
    arena_definitions
) {
  
  result <- lapply(
    shape_levels,
    function(shape) {
      
      lapply(
        1:2,
        function(day) {
          
          prepare_position_activity(
            all_trials = all_trials,
            mouse = mouse,
            shape = shape,
            day = day,
            percentile = percentile,
            arena_definitions = arena_definitions
          )
        }
      )
    }
  )
  
  data <- do.call(
    rbind,
    unlist(
      result,
      recursive = FALSE
    )
  )
  
  data <- set_condition_factors(
    data
  )
  
  rownames(data) <- NULL
  
  data
}


prepare_all_mice_position_activity <- function(
    all_trials,
    percentile,
    arena_definitions
) {
  
  data <- do.call(
    rbind,
    lapply(
      names(all_trials),
      function(mouse) {
        
        prepare_mouse_position_activity(
          all_trials = all_trials,
          mouse = mouse,
          percentile = percentile,
          arena_definitions = arena_definitions
        )
      }
    )
  )
  
  rownames(data) <- NULL
  
  data
}


# Prepare the full dataset once
all_position_activity <- prepare_all_mice_position_activity(
  all_trials = all_trials,
  percentile = percentile_used,
  arena_definitions = arena_definitions
)


############################################################
# 10. Activity classes and colours
############################################################

discretize_population_activity <- function(
    proportion_active
) {
  
  cut(
    proportion_active,
    breaks = c(
      -Inf,
      0,
      0.02,
      0.05,
      0.10,
      Inf
    ),
    labels = c(
      "No activity",
      ">0–2%",
      ">2–5%",
      ">5–10%",
      ">10%"
    ),
    include.lowest = TRUE
  )
}


activity_colours <- c(
  "No activity" = "grey85",
  ">0–2%" = "#B3CDE3",
  ">2–5%" = "#74ADD1",
  ">5–10%" = "#FDAE61",
  ">10%" = "#B2182B"
)


add_activity_classes <- function(data) {
  
  data$activity_class <-
    discretize_population_activity(
      data$proportion_active
    )
  
  # Low activity first, high activity last
  data <- data[
    order(
      data$proportion_active,
      data$number_active
    ),
  ]
  
  data
}


############################################################
# 11. Position-related activity for one mouse
############################################################

plot_mouse_position_activity <- function(
    data,
    mouse,
    boundary_data,
    show_title = TRUE
) {
  
  mouse_data <- data[
    data$mouse == mouse &
      is.finite(data$x) &
      is.finite(data$y),
  ]
  
  mouse_data <- add_activity_classes(
    mouse_data
  )
  
  ggplot() +
    
    geom_path(
      data = mouse_data,
      aes(
        x = x,
        y = y,
        group = interaction(
          shape,
          day
        )
      ),
      colour = "grey75",
      linewidth = 0.15,
      alpha = 0.25
    ) +
    
    geom_point(
      data = mouse_data,
      aes(
        x = x,
        y = y,
        colour = activity_class,
        size = number_active
      ),
      alpha = 0.60
    ) +
    
    geom_path(
      data = boundary_data,
      aes(
        x = x,
        y = y,
        group = shape
      ),
      colour = "black",
      linewidth = 0.65
    ) +
    
    facet_grid(
      rows = vars(day),
      cols = vars(shape)
    ) +
    
    scale_colour_manual(
      values = activity_colours,
      limits = names(activity_colours),
      drop = FALSE,
      name = "Active neurons"
    ) +
    
    scale_size_continuous(
      range = c(
        0.15,
        3.2
      ),
      name = "Number active",
      trans = "sqrt"
    ) +
    
    coord_fixed() +
    
    labs(
      title = if (show_title) {
        paste(
          "Mouse",
          mouse
        )
      } else {
        NULL
      },
      
      subtitle = NULL,
      x = "Position x",
      y = "Position y"
    ) +
    
    theme_minimal(base_size = 9) +
    
    theme(
      panel.grid = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 8
      ),
      
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey70",
        linewidth = 0.25
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 8
      ),
      
      axis.text = element_text(
        size = 7
      ),
      
      legend.position = "right",
      
      panel.spacing = grid::unit(
        0.3,
        "lines"
      )
    )
}


############################################################
# 12. Create one activity grid for every mouse
############################################################

position_activity_plots <- lapply(
  mouse_names,
  function(mouse) {
    
    plot_mouse_position_activity(
      data = all_position_activity,
      mouse = mouse,
      boundary_data = boundary_data
    )
  }
)

names(position_activity_plots) <- mouse_names


# Example for one mouse
selected_mouse <- "M3424F"

p_position_activity <-
  position_activity_plots[[selected_mouse]]

print(
  p_position_activity
)


############################################################
# Compact position-activity plot for all mice
#
# Rows: mice
# Columns: shape and day
############################################################

plot_all_mice_position_activity <- function(
    data,
    boundary_data,
    percentile = 0.90
) {
  
  plot_data <- data[
    is.finite(data$x) &
      is.finite(data$y),
  ]
  
  plot_data$activity_class <-
    discretize_population_activity(
      plot_data$proportion_active
    )
  
  # Draw low-activity observations first and
  # high-activity observations last
  plot_data <- plot_data[
    order(
      plot_data$proportion_active,
      plot_data$number_active
    ),
  ]
  
  # Ensure the desired ordering
  plot_data$mouse <- factor(
    plot_data$mouse,
    levels = names(all_trials)
  )
  
  plot_data$shape <- factor(
    plot_data$shape,
    levels = c(
      "Triangle",
      "Circle",
      "Square"
    )
  )
  
  plot_data$day <- factor(
    plot_data$day,
    levels = c(
      "Day 1",
      "Day 2"
    )
  )
  
  ggplot() +
    
    # Mouse trajectory
    geom_path(
      data = plot_data,
      aes(
        x = x,
        y = y,
        group = interaction(
          mouse,
          shape,
          day
        )
      ),
      colour = "grey75",
      linewidth = 0.12,
      alpha = 0.20
    ) +
    
    # Position-related neuronal activity
    geom_point(
      data = plot_data,
      aes(
        x = x,
        y = y,
        colour = activity_class,
        size = number_active
      ),
      alpha = 0.65
    ) +
    
    # Estimated arena boundaries
    geom_path(
      data = boundary_data,
      aes(
        x = x,
        y = y,
        group = shape
      ),
      colour = "black",
      linewidth = 0.55
    ) +
    
    # One row per mouse and six columns:
    # Triangle D1, Triangle D2, Circle D1, ...
    facet_grid(
      rows = vars(mouse),
      cols = vars(shape, day),
      switch = "y"
    ) +
    
    scale_colour_manual(
      values = activity_colours,
      limits = names(activity_colours),
      drop = FALSE,
      name = "Proportion active"
    ) +
    
    scale_size_continuous(
      range = c(
        0.10,
        3.0
      ),
      trans = "sqrt"
    ) +
    guides(
      size = "none"
    ) +
    
    coord_fixed() +
    
    labs(
      title = "Position-related neuronal activity for all mice",
      subtitle = paste0(
        "Colour represents the proportion of active neurons; ",
        "point size represents the number of simultaneously active neurons. ",
        "Activity is defined using each neuron's ",
        percentile * 100,
        "th percentile."
      ),
      x = "Position x",
      y = "Position y"
    ) +
    
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 3,
          alpha = 1
        ),
        nrow = 1
      )
    ) +
    
    theme_minimal(base_size = 9) +
    
    theme(
      panel.grid = element_blank(),
      
      panel.spacing.x = grid::unit(
        0.20,
        "lines"
      ),
      
      panel.spacing.y = grid::unit(
        0.25,
        "lines"
      ),
      
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey70",
        linewidth = 0.25
      ),
      
      strip.text.x = element_text(
        face = "bold",
        size = 8
      ),
      
      strip.text.y.left = element_text(
        face = "bold",
        size = 9,
        angle = 0
      ),
      
      axis.title = element_text(
        size = 9
      ),
      
      axis.text = element_text(
        size = 6
      ),
      
      legend.position = "bottom",
      
      legend.box = "vertical",
      
      legend.title = element_text(
        face = "bold",
        size = 9
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5
      )
    )
}


p_position_activity_all_mice <-
  plot_all_mice_position_activity(
    data = all_position_activity,
    boundary_data = boundary_data,
    percentile = percentile_used
  )

print(
  p_position_activity_all_mice
)


############################################################
# 14. Summarize activity by boundary distance
############################################################

summarize_boundary_activity <- function(
    data,
    number_distance_bins = 20
) {
  
  data <- data[
    is.finite(data$boundary_distance) &
      is.finite(data$proportion_active),
  ]
  
  condition <- interaction(
    data$mouse,
    data$shape,
    data$day,
    drop = TRUE
  )
  
  robust_maximum <- ave(
    data$boundary_distance,
    condition,
    FUN = function(x) {
      
      as.numeric(
        quantile(
          x,
          0.99,
          na.rm = TRUE,
          names = FALSE
        )
      )
    }
  )
  
  robust_maximum[
    !is.finite(robust_maximum) |
      robust_maximum <= 0
  ] <- 1
  
  data$relative_boundary_distance <- pmin(
    data$boundary_distance /
      robust_maximum,
    1
  )
  
  data$distance_bin <- cut(
    data$relative_boundary_distance,
    breaks = seq(
      0,
      1,
      length.out =
        number_distance_bins + 1
    ),
    include.lowest = TRUE,
    labels = FALSE
  )
  
  means <- aggregate(
    cbind(
      relative_boundary_distance,
      proportion_active
    ) ~ mouse + shape + day + distance_bin,
    data = data,
    FUN = mean
  )
  
  frame_counts <- aggregate(
    frame ~ mouse + shape + day + distance_bin,
    data = data,
    FUN = length
  )
  
  names(frame_counts)[
    names(frame_counts) == "frame"
  ] <- "number_frames"
  
  merge(
    means,
    frame_counts,
    by = c(
      "mouse",
      "shape",
      "day",
      "distance_bin"
    )
  )
}


boundary_summary <- summarize_boundary_activity(
  data = all_position_activity,
  number_distance_bins = number_distance_bins
)


############################################################
# 15. Boundary-distance relationship for one mouse
############################################################

plot_mouse_boundary_activity <- function(
    boundary_summary,
    mouse
) {
  
  mouse_data <- boundary_summary[
    boundary_summary$mouse == mouse,
  ]
  
  ggplot(
    mouse_data,
    aes(
      x = relative_boundary_distance,
      y = proportion_active
    )
  ) +
    geom_line(
      linewidth = 0.7
    ) +
    geom_point(
      aes(
        size = number_frames
      ),
      alpha = 0.75
    ) +
    facet_grid(
      rows = vars(shape),
      cols = vars(day)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.25
      )
    ) +
    scale_size_continuous(
      name = "Frames"
    ) +
    labs(
      title = paste(
        "Activity and boundary distance — Mouse",
        mouse
      ),
      subtitle = paste(
        "0 = boundary; 1 = most central observed positions"
      ),
      x = "Relative distance from boundary",
      y = "Mean proportion of active neurons"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


p_boundary_activity <- plot_mouse_boundary_activity(
  boundary_summary,
  selected_mouse
)

print(
  p_boundary_activity
)


############################################################
# 16. Boundary-distance relationship across all mice
############################################################

plot_all_mice_boundary_activity <- function(
    boundary_summary
) {
  
  overall_mean <- aggregate(
    cbind(
      relative_boundary_distance,
      proportion_active
    ) ~ shape + day + distance_bin,
    data = boundary_summary,
    FUN = mean
  )
  
  overall_sd <- aggregate(
    proportion_active ~
      shape + day + distance_bin,
    data = boundary_summary,
    FUN = sd,
    na.rm = TRUE
  )
  
  names(overall_sd)[
    names(overall_sd) == "proportion_active"
  ] <- "activity_sd"
  
  number_mice <- aggregate(
    proportion_active ~
      shape + day + distance_bin,
    data = boundary_summary,
    FUN = function(x) {
      sum(is.finite(x))
    }
  )
  
  names(number_mice)[
    names(number_mice) == "proportion_active"
  ] <- "number_mice"
  
  overall_mean <- merge(
    overall_mean,
    overall_sd,
    by = c(
      "shape",
      "day",
      "distance_bin"
    )
  )
  
  overall_mean <- merge(
    overall_mean,
    number_mice,
    by = c(
      "shape",
      "day",
      "distance_bin"
    )
  )
  
  overall_mean$activity_se <-
    overall_mean$activity_sd /
    sqrt(overall_mean$number_mice)
  
  overall_mean$activity_se[
    !is.finite(overall_mean$activity_se)
  ] <- 0
  
  overall_mean$lower <- pmax(
    0,
    overall_mean$proportion_active -
      overall_mean$activity_se
  )
  
  overall_mean$upper <- pmin(
    1,
    overall_mean$proportion_active +
      overall_mean$activity_se
  )
  
  ggplot() +
    
    geom_line(
      data = boundary_summary,
      aes(
        x = relative_boundary_distance,
        y = proportion_active,
        group = mouse
      ),
      colour = "grey65",
      linewidth = 0.35,
      alpha = 0.45
    ) +
    
    geom_ribbon(
      data = overall_mean,
      aes(
        x = relative_boundary_distance,
        ymin = lower,
        ymax = upper
      ),
      fill = "#74ADD1",
      alpha = 0.25
    ) +
    
    geom_line(
      data = overall_mean,
      aes(
        x = relative_boundary_distance,
        y = proportion_active
      ),
      colour = "#2166AC",
      linewidth = 1
    ) +
    
    facet_grid(
      rows = vars(day),
      cols = vars(shape)
    ) +
    
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.25
      )
    ) +
    
    labs(
      title = paste(
        "Neuronal activity and boundary distance",
        "across all mice"
      ),
      subtitle = paste(
        "Grey lines: individual mice;",
        "blue line: mean; ribbon: ±1 standard error"
      ),
      x = paste(
        "Relative distance from boundary",
        "(0 = boundary, 1 = central position)"
      ),
      y = "Mean proportion of active neurons"
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


p_boundary_activity_all_mice <-
  plot_all_mice_boundary_activity(
    boundary_summary
  )

print(
  p_boundary_activity_all_mice
)


############################################################
# 17. Boundary-distance maps for one mouse
############################################################

distance_colours <- c(
  "0–5%" = "#2166AC",
  ">5–15%" = "#67A9CF",
  ">15–30%" = "#D1E5F0",
  ">30–50%" = "#FDAE61",
  ">50%" = "#B2182B"
)


plot_mouse_boundary_distance_map <- function(
    data,
    mouse,
    boundary_data
) {
  
  mouse_data <- data[
    data$mouse == mouse &
      is.finite(data$x) &
      is.finite(data$y) &
      is.finite(data$boundary_distance),
  ]
  
  condition <- interaction(
    mouse_data$shape,
    mouse_data$day,
    drop = TRUE
  )
  
  maximum_distance <- ave(
    mouse_data$boundary_distance,
    condition,
    FUN = function(x) {
      
      as.numeric(
        quantile(
          x,
          0.99,
          na.rm = TRUE,
          names = FALSE
        )
      )
    }
  )
  
  maximum_distance[
    !is.finite(maximum_distance) |
      maximum_distance <= 0
  ] <- 1
  
  mouse_data$relative_boundary_distance <-
    pmin(
      mouse_data$boundary_distance /
        maximum_distance,
      1
    )
  
  mouse_data$distance_class <- cut(
    mouse_data$relative_boundary_distance,
    breaks = c(
      -Inf,
      0.05,
      0.15,
      0.30,
      0.50,
      Inf
    ),
    labels = names(
      distance_colours
    ),
    include.lowest = TRUE
  )
  
  mouse_data <- mouse_data[
    order(
      mouse_data$relative_boundary_distance
    ),
  ]
  
  ggplot() +
    geom_point(
      data = mouse_data,
      aes(
        x = x,
        y = y,
        colour = distance_class
      ),
      size = 0.7,
      alpha = 0.65
    ) +
    geom_path(
      data = boundary_data,
      aes(
        x = x,
        y = y,
        group = shape
      ),
      colour = "black",
      linewidth = 0.75
    ) +
    facet_grid(
      rows = vars(shape),
      cols = vars(day)
    ) +
    scale_colour_manual(
      values = distance_colours,
      limits = names(distance_colours),
      drop = FALSE,
      name = "Relative distance"
    ) +
    coord_fixed() +
    labs(
      title = paste(
        "Distance from boundary — Mouse",
        mouse
      ),
      subtitle = paste(
        "Blue = close to the boundary;",
        "red = farther from the boundary"
      ),
      x = "Position x",
      y = "Position y"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


p_boundary_distance_map <-
  plot_mouse_boundary_distance_map(
    data = all_position_activity,
    mouse = selected_mouse,
    boundary_data = boundary_data
  )

print(
  p_boundary_distance_map
)

############################################################
# 18. Example raw S time series (informational)
############################################################
plot_example_S_traces <- function(
    all_trials,
    mouse,
    trial,
    number_traces = 5,
    seed = 1
) {
  
  trial_data <- all_trials[[mouse]][[trial]]
  S <- as.matrix(trial_data$neuron$S)
  
  number_traces <- min(number_traces, nrow(S))
  
  set.seed(seed)
  selected_neurons <- sample(
    seq_len(nrow(S)),
    number_traces
  )
  
  time_vector <- as.numeric(trial_data$neuron$time)
  
  if (length(time_vector) != ncol(S)) {
    time_vector <- seq_len(ncol(S))
  }
  
  trace_data <- do.call(
    rbind,
    lapply(
      selected_neurons,
      function(i) {
        data.frame(
          time = time_vector,
          S = S[i, ],
          neuron = paste0("Neuron ", i)
        )
      }
    )
  )
  
  trace_data$neuron <- factor(
    trace_data$neuron,
    levels = paste0("Neuron ", selected_neurons)
  )
  
  ggplot(
    trace_data,
    aes(
      x = time,
      y = S
    )
  ) +
    geom_line(
      linewidth = 0.4,
      colour = "#2166AC"
    ) +
    facet_wrap(
      vars(neuron),
      ncol = 1,
      scales = "free_y"
    ) +
    labs(
      title = "Spike traces",
      subtitle = paste0("Mouse: ",
        mouse,
        ", ",
        trial
      ),
      x = "Time",
      y = "S"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 8
      ),
      plot.title = element_text(
        face = "bold"
      ),
      axis.title.y = element_text(
        angle = 0,
        vjust = 0.5
      )
    )
}

p_S_timeseries_example <- plot_example_S_traces(
  all_trials = all_trials,
  mouse = "M3424F",
  trial = "trial1",
  number_traces = 5
)
print(p_S_timeseries_example)

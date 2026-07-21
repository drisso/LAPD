# Data preparation utilities for the hurdle-GAM analysis.

TRIAL_KEY <- data.frame(
  trial = paste0("trial", 1:6),
  day = c(1L, 1L, 1L, 2L, 2L, 2L),
  shape = c("triangle", "circle", "square", "circle", "square", "triangle"),
  stringsAsFactors = FALSE
)

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required. Install it before continuing.", package),
      call. = FALSE
    )
  }
}

load_trial_data <- function(
    rds_path,
    mouse = 1L,
    trial = 1L,
    cache_path = NULL
) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  if (!file.exists(rds_path)) {
    stop("Could not find the data file: ", rds_path, call. = FALSE)
  }

  all_trials <- readRDS(rds_path)

  mouse_name <- if (is.numeric(mouse)) names(all_trials)[mouse] else mouse
  if (length(mouse_name) != 1L || is.na(mouse_name) ||
      !mouse_name %in% names(all_trials)) {
    stop("The requested mouse was not found in all_trials.", call. = FALSE)
  }

  mouse_trials <- all_trials[[mouse_name]]
  trial_name <- if (is.numeric(trial)) names(mouse_trials)[trial] else trial
  if (length(trial_name) != 1L || is.na(trial_name) ||
      !trial_name %in% names(mouse_trials)) {
    stop("The requested trial was not found for mouse ", mouse_name, ".", call. = FALSE)
  }

  out <- mouse_trials[[trial_name]]
  out$mouse_name <- mouse_name
  out$trial_name <- trial_name
  out$shape <- TRIAL_KEY$shape[match(trial_name, TRIAL_KEY$trial)]

  required <- c("behav", "neuron")
  if (!all(required %in% names(out)) || is.null(out$neuron$S)) {
    stop("The selected trial does not contain behav, neuron, and neuron$S.", call. = FALSE)
  }

  # The source RDS is large. Caching just one trial makes subsequent renders much
  # faster and avoids repeatedly materializing the full nested object.
  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_path, compress = "xz")
  }

  rm(all_trials, mouse_trials)
  invisible(gc())
  out
}

calcium_frame_times <- function(neuron_time, n_frames, behaviour_time) {
  neuron_time <- as.numeric(neuron_time)
  neuron_time <- neuron_time[is.finite(neuron_time)]
  behaviour_time <- as.numeric(behaviour_time)
  behaviour_time <- behaviour_time[is.finite(behaviour_time)]

  if (length(neuron_time) == n_frames) {
    return(neuron_time)
  }

  if (length(neuron_time) >= 2L) {
    # neuron$time follows the original acquisition clock and is approximately
    # twice as long as S. An evenly spaced grid preserves its endpoints without
    # pretending that it is one-to-one with the columns of S.
    return(seq(min(neuron_time), max(neuron_time), length.out = n_frames))
  }

  if (length(behaviour_time) < 2L) {
    stop("At least two finite timestamps are needed to align the data.", call. = FALSE)
  }
  seq(min(behaviour_time), max(behaviour_time), length.out = n_frames)
}

align_behaviour_to_frames <- function(position, behaviour_time, frame_time) {
  position <- as.matrix(position)
  behaviour_time <- as.numeric(behaviour_time)
  if (ncol(position) < 2L || nrow(position) != length(behaviour_time)) {
    stop("Behaviour position must have two columns and one row per timestamp.", call. = FALSE)
  }

  valid <- is.finite(behaviour_time) & is.finite(position[, 1L]) &
    is.finite(position[, 2L])
  position <- position[valid, 1:2, drop = FALSE]
  behaviour_time <- behaviour_time[valid]
  ordering <- order(behaviour_time)
  position <- position[ordering, , drop = FALSE]
  behaviour_time <- behaviour_time[ordering]
  keep <- !duplicated(behaviour_time)
  position <- position[keep, , drop = FALSE]
  behaviour_time <- behaviour_time[keep]

  data.frame(
    frame = seq_along(frame_time),
    time_ms = frame_time,
    time_sec = (frame_time - min(frame_time)) / 1000,
    x = approx(behaviour_time, position[, 1L], xout = frame_time, rule = 2)$y,
    y = approx(behaviour_time, position[, 2L], xout = frame_time, rule = 2)$y
  )
}

triangle_area <- function(p1, p2, p3) {
  abs(
    p1[1L] * (p2[2L] - p3[2L]) +
      p2[1L] * (p3[2L] - p1[2L]) +
      p3[1L] * (p1[2L] - p2[2L])
  ) / 2
}

estimate_arena_boundary <- function(
    positions,
    shape,
    lower_quantile = 0.005,
    upper_quantile = 0.995,
    max_hull_points = 100L
) {
  positions <- as.data.frame(positions)[, c("x", "y"), drop = FALSE]
  positions <- positions[is.finite(positions$x) & is.finite(positions$y), , drop = FALSE]
  if (nrow(positions) < 10L) {
    stop("Too few valid positions to estimate the arena boundary.", call. = FALSE)
  }

  shape <- match.arg(tolower(shape), c("triangle", "circle", "square"))
  x_limits <- stats::quantile(
    positions$x, c(lower_quantile, upper_quantile), na.rm = TRUE, names = FALSE
  )
  y_limits <- stats::quantile(
    positions$y, c(lower_quantile, upper_quantile), na.rm = TRUE, names = FALSE
  )
  trimmed <- positions[
    positions$x >= x_limits[1L] & positions$x <= x_limits[2L] &
      positions$y >= y_limits[1L] & positions$y <= y_limits[2L],
    , drop = FALSE
  ]

  if (shape == "circle") {
    centre <- c(x = median(trimmed$x), y = median(trimmed$y))
    radial <- sqrt((trimmed$x - centre["x"])^2 + (trimmed$y - centre["y"])^2)
    return(list(
      type = shape,
      centre = centre,
      radius = unname(stats::quantile(radial, upper_quantile, names = FALSE))
    ))
  }

  if (shape == "square") {
    return(list(
      type = shape,
      xmin = min(trimmed$x), xmax = max(trimmed$x),
      ymin = min(trimmed$y), ymax = max(trimmed$y)
    ))
  }

  unique_xy <- unique(trimmed)
  hull <- as.matrix(unique_xy[chull(unique_xy$x, unique_xy$y), , drop = FALSE])
  if (nrow(hull) > max_hull_points) {
    idx <- unique(round(seq(1, nrow(hull), length.out = max_hull_points)))
    hull <- hull[idx, , drop = FALSE]
  }
  if (nrow(hull) < 3L) {
    stop("The triangle boundary could not be estimated.", call. = FALSE)
  }
  combinations <- utils::combn(seq_len(nrow(hull)), 3L)
  areas <- apply(combinations, 2L, function(i) {
    triangle_area(hull[i[1L], ], hull[i[2L], ], hull[i[3L], ])
  })
  vertices <- hull[combinations[, which.max(areas)], , drop = FALSE]
  colnames(vertices) <- c("x", "y")
  list(type = shape, vertices = vertices)
}

distance_to_segment <- function(x, y, x1, y1, x2, y2) {
  dx <- x2 - x1
  dy <- y2 - y1
  length_squared <- dx^2 + dy^2
  if (length_squared == 0) {
    return(sqrt((x - x1)^2 + (y - y1)^2))
  }
  projection <- ((x - x1) * dx + (y - y1) * dy) / length_squared
  projection <- pmax(0, pmin(1, projection))
  closest_x <- x1 + projection * dx
  closest_y <- y1 + projection * dy
  sqrt((x - closest_x)^2 + (y - closest_y)^2)
}

distance_to_boundary <- function(x, y, boundary) {
  if (boundary$type == "circle") {
    radial <- sqrt((x - boundary$centre["x"])^2 + (y - boundary$centre["y"])^2)
    return(abs(boundary$radius - radial))
  }

  if (boundary$type == "square") {
    return(pmin(
      abs(x - boundary$xmin), abs(boundary$xmax - x),
      abs(y - boundary$ymin), abs(boundary$ymax - y)
    ))
  }

  vertices <- boundary$vertices
  vertices <- rbind(vertices, vertices[1L, , drop = FALSE])
  edge_distances <- vapply(seq_len(nrow(vertices) - 1L), function(i) {
    distance_to_segment(
      x, y,
      vertices[i, 1L], vertices[i, 2L],
      vertices[i + 1L, 1L], vertices[i + 1L, 2L]
    )
  }, numeric(length(x)))
  if (is.null(dim(edge_distances))) edge_distances else apply(edge_distances, 1L, min)
}

boundary_dataframe <- function(boundary, n_circle = 300L) {
  if (boundary$type == "circle") {
    angle <- seq(0, 2 * pi, length.out = n_circle)
    return(data.frame(
      x = boundary$centre["x"] + boundary$radius * cos(angle),
      y = boundary$centre["y"] + boundary$radius * sin(angle)
    ))
  }
  if (boundary$type == "square") {
    return(data.frame(
      x = c(boundary$xmin, boundary$xmax, boundary$xmax, boundary$xmin, boundary$xmin),
      y = c(boundary$ymin, boundary$ymin, boundary$ymax, boundary$ymax, boundary$ymin)
    ))
  }
  vertices <- rbind(boundary$vertices, boundary$vertices[1L, , drop = FALSE])
  data.frame(x = vertices[, 1L], y = vertices[, 2L])
}

calculate_speed <- function(x, y, time_sec, winsor_quantile = 0.995) {
  dt <- c(NA_real_, diff(time_sec))
  displacement <- c(NA_real_, sqrt(diff(x)^2 + diff(y)^2))
  speed <- displacement / dt
  speed[!is.finite(speed) | speed < 0] <- NA_real_
  speed[1L] <- stats::median(speed, na.rm = TRUE)
  # A short running median suppresses single-frame tracking jumps.
  speed <- stats::runmed(speed, k = 5L, endrule = "median")
  cap <- stats::quantile(speed, winsor_quantile, na.rm = TRUE, names = FALSE)
  pmin(speed, cap)
}

make_frame_covariates <- function(trial_data, holdout_seconds = 30) {
  S <- as.matrix(trial_data$neuron$S)
  n_frames <- ncol(S)
  frame_time <- calcium_frame_times(
    trial_data$neuron$time,
    n_frames,
    trial_data$behav$time
  )
  frames <- align_behaviour_to_frames(
    trial_data$behav$position,
    trial_data$behav$time,
    frame_time
  )
  frames$speed <- calculate_speed(frames$x, frames$y, frames$time_sec)

  holdout_start <- max(frames$time_sec) - holdout_seconds
  frames$set <- ifelse(frames$time_sec > holdout_start, "test", "train")
  if (sum(frames$set == "test") < 10L || sum(frames$set == "train") < 100L) {
    stop("The requested holdout leaves too few training or test frames.", call. = FALSE)
  }

  # Estimate geometry from training positions so the held-out period does not
  # contribute to feature construction.
  shape <- trial_data$shape
  if (length(shape) != 1L || is.na(shape)) {
    shape <- TRIAL_KEY$shape[match(trial_data$trial_name, TRIAL_KEY$trial)]
  }
  boundary <- estimate_arena_boundary(frames[frames$set == "train", ], shape)
  frames$border_distance <- distance_to_boundary(frames$x, frames$y, boundary)

  list(frames = frames, boundary = boundary, S = S)
}

standardize_frame_covariates <- function(
    frames,
    variables = c("x", "y", "border_distance", "speed", "time_sec")
) {
  train <- frames$set == "train"
  scaling <- lapply(variables, function(variable) {
    centre <- mean(frames[[variable]][train], na.rm = TRUE)
    spread <- stats::sd(frames[[variable]][train], na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) spread <- 1
    c(centre = centre, scale = spread)
  })
  names(scaling) <- variables

  z_names <- c(
    x = "x_z", y = "y_z", border_distance = "border_z",
    speed = "speed_z", time_sec = "time_z"
  )
  for (variable in variables) {
    frames[[z_names[[variable]]]] <-
      (frames[[variable]] - scaling[[variable]]["centre"]) /
      scaling[[variable]]["scale"]
  }
  attr(frames, "scaling") <- scaling
  frames
}

event_thresholds <- function(
    S,
    train_frames,
    method = c("three_sd", "positive", "quantile"),
    quantile_probability = 0.90,
    neuron_indices = seq_len(nrow(S))
) {
  method <- match.arg(method)
  train_values <- S[neuron_indices, train_frames, drop = FALSE]
  train_values[!is.finite(train_values) | train_values < 0] <- 0

  if (method == "three_sd") {
    thresholds <- 3 * apply(train_values, 1L, stats::sd, na.rm = TRUE)
  } else if (method == "positive") {
    thresholds <- rep(0, length(neuron_indices))
  } else {
    thresholds <- apply(train_values, 1L, stats::quantile,
      probs = quantile_probability, na.rm = TRUE, names = FALSE
    )
  }
  thresholds[!is.finite(thresholds)] <- Inf
  names(thresholds) <- as.character(neuron_indices)
  thresholds
}

make_long_hurdle_data <- function(
    S,
    frames,
    frame_indices,
    thresholds,
    neuron_indices = seq_len(nrow(S))
) {
  require_namespace("data.table")
  frame_indices <- as.integer(frame_indices)
  neuron_indices <- as.integer(neuron_indices)
  if (!all(as.character(neuron_indices) %in% names(thresholds))) {
    stop("A threshold is required for every selected neuron.", call. = FALSE)
  }

  wide <- data.table::as.data.table(t(S[neuron_indices, frame_indices, drop = FALSE]))
  data.table::setnames(wide, as.character(neuron_indices))
  wide[, frame_row := seq_len(.N)]
  long <- data.table::melt(
    wide,
    id.vars = "frame_row",
    variable.name = "neuron_index",
    value.name = "amplitude",
    variable.factor = FALSE
  )
  rm(wide)
  long[, neuron_index := as.integer(neuron_index)]
  long[!is.finite(amplitude) | amplitude < 0, amplitude := 0]
  long[, threshold := thresholds[as.character(neuron_index)]]
  long[, event := as.integer(amplitude > threshold)]
  long[, neuron := factor(
    paste0("n", neuron_index),
    levels = paste0("n", neuron_indices)
  )]

  frame_subset <- frames[frame_indices, , drop = FALSE]
  covariate_names <- c(
    "frame", "time_sec", "x", "y", "border_distance", "speed",
    "x_z", "y_z", "border_z", "speed_z", "time_z", "set"
  )
  for (variable in covariate_names) {
    long[, (variable) := frame_subset[[variable]][frame_row]]
  }
  long[, frame_row := NULL]
  data.table::setcolorder(long, c(
    "frame", "time_sec", "neuron", "neuron_index", "event", "amplitude",
    "threshold", "x", "y", "border_distance", "speed",
    "x_z", "y_z", "border_z", "speed_z", "time_z", "set"
  ))
  long[]
}

prepare_hurdle_analysis_data <- function(
    trial_data,
    holdout_seconds = 30,
    threshold_method = "three_sd",
    quantile_probability = 0.90,
    neuron_indices = NULL
) {
  prepared <- make_frame_covariates(trial_data, holdout_seconds)
  frames <- standardize_frame_covariates(prepared$frames)
  S <- prepared$S
  if (is.null(neuron_indices)) neuron_indices <- seq_len(nrow(S))
  train_frames <- which(frames$set == "train")
  test_frames <- which(frames$set == "test")
  thresholds <- event_thresholds(
    S, train_frames,
    method = threshold_method,
    quantile_probability = quantile_probability,
    neuron_indices = neuron_indices
  )

  train <- make_long_hurdle_data(S, frames, train_frames, thresholds, neuron_indices)
  test <- make_long_hurdle_data(S, frames, test_frames, thresholds, neuron_indices)

  list(
    train = train,
    test = test,
    frames = frames,
    boundary = prepared$boundary,
    thresholds = thresholds,
    neuron_indices = neuron_indices,
    mouse_name = trial_data$mouse_name,
    trial_name = trial_data$trial_name,
    shape = trial_data$shape
  )
}

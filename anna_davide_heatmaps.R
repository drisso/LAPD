
rm(list=ls())
# working directory

setwd("C:/Users/david/Desktop/DataResearchCamp-Florence2026")

# load data
all_trials <- readRDS("all_trials.RDS")

str(all_trials)


# Extract the matrix
C_raw <- all_trials$M3411$trial1$neuron$C_raw

dim(C_raw)

image(
  x = 1:ncol(C_raw),
  y = 1:nrow(C_raw),
  z = t(C_raw[nrow(C_raw):1, ]),  # flip so neuron 1 is at the top
  col = hcl.colors(100, "Inferno"),
  xlab = "Time",
  ylab = "Neuron",
  useRaster = TRUE
)


trial_key <- data.frame(
  trial = paste0("trial", 1:6),
  day   = c(1, 1, 1, 2, 2, 2),
  shape = c(
    "triangle", "circle", "square",
    "circle", "square", "triangle"
  )
)

trial_key

# Average neighbouring time points
bin_time <- function(mat, bin_size = 10) {
  
  if (bin_size <= 1) {
    return(mat)
  }
  
  n_complete <- floor(ncol(mat) / bin_size) * bin_size
  
  # Remove any incomplete final bin
  mat <- mat[, seq_len(n_complete), drop = FALSE]
  
  bin_id <- rep(
    seq_len(n_complete / bin_size),
    each = bin_size
  )
  
  # Average columns belonging to the same time bin
  binned <- t(
    rowsum(
      t(mat),
      group = bin_id,
      reorder = FALSE
    ) / bin_size
  )
  
  binned
}


binarize_neurons <- function(mat, percentile = 0.90) {
  
  thresholds <- apply(
    mat,
    1,
    quantile,
    probs = percentile,
    na.rm = TRUE
  )
  
  binary_mat <- sweep(
    mat,
    MARGIN = 1,
    STATS = thresholds,
    FUN = ">"
  )
  
  # Convert TRUE/FALSE to 1/0
  binary_mat <- binary_mat * 1L
  
  # Keep missing values from causing plotting issues
  binary_mat[is.na(binary_mat)] <- 0L
  
  binary_mat
}

matrix_to_dataframe <- function(mat, mouse, day, shape) {
  
  data.frame(
    mouse  = mouse,
    day    = paste("Day", day),
    shape  = shape,
    neuron = rep(seq_len(nrow(mat)), times = ncol(mat)),
    time   = rep(seq_len(ncol(mat)), each = nrow(mat)),
    value  = as.vector(mat)
  )
}


get_C_raw <- function(
    all_trials,
    mouse,
    shape,
    day,
    bin_size = 1,
    percentile = 0.90
) {
  
  shape <- match.arg(
    tolower(shape),
    c("triangle", "circle", "square")
  )
  
  trial_name <- trial_key$trial[
    trial_key$shape == shape &
      trial_key$day == day
  ]
  
  if (length(trial_name) != 1) {
    stop("Could not uniquely identify the trial.")
  }
  
  mat <- all_trials[[mouse]][[trial_name]]$neuron$C_raw
  
  if (is.null(mat)) {
    stop(
      paste(
        "C_raw is missing for",
        mouse,
        trial_name
      )
    )
  }
  
  # Binary thresholding should happen before optional time binning
  mat <- binarize_neurons(
    mat,
    percentile = percentile
  )
  
  # Optional: aggregate binary activity over time windows
  if (bin_size > 1) {
    mat <- bin_time(mat, bin_size = bin_size)
  }
  
  mat
}


plot_shape_all_mice <- function(
    all_trials,
    shape,
    percentile = 0.90,
    bin_size = 1
) {
  
  shape <- match.arg(
    tolower(shape),
    c("triangle", "circle", "square")
  )
  
  mice <- names(all_trials)
  
  plot_data <- lapply(mice, function(mouse) {
    
    lapply(1:2, function(day) {
      
      mat <- get_C_raw(
        all_trials = all_trials,
        mouse = mouse,
        shape = shape,
        day = day,
        percentile = percentile,
        bin_size = bin_size
      )
      
      matrix_to_dataframe(
        mat = mat,
        mouse = mouse,
        day = day,
        shape = shape
      )
    })
  })
  
  plot_data <- do.call(
    rbind,
    unlist(plot_data, recursive = FALSE)
  )
  
  plot_data$mouse <- factor(
    plot_data$mouse,
    levels = mice
  )
  
  plot_data$day <- factor(
    plot_data$day,
    levels = c("Day 1", "Day 2")
  )
  
  ggplot(
    plot_data,
    aes(x = time, y = neuron, fill = value)
  ) +
    geom_raster() +
    facet_grid(
      rows = vars(mouse),
      cols = vars(day),
      scales = "free_y"
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 1),
      name = if (bin_size == 1) {
        "Activity"
      } else {
        "Active\nproportion"
      }
    ) +
    scale_y_reverse() +
    labs(
      title = paste(
        tools::toTitleCase(shape),
        "activity: Day 1 versus Day 2"
      ),
      subtitle = paste0(
        "Active when C_raw exceeds the neuron-specific ",
        percentile * 100,
        "th percentile"
      ),
      x = if (bin_size == 1) "Time point" else "Time bin",
      y = "Neuron"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    )
}


p_triangle <- plot_shape_all_mice(
  all_trials,
  shape = "triangle",
  percentile = 0.90,
  bin_size = 10
)

p_circle <- plot_shape_all_mice(
  all_trials,
  shape = "circle",
  percentile = 0.90,
  bin_size = 10
)

p_square <- plot_shape_all_mice(
  all_trials,
  shape = "square",
  percentile = 0.90,
  bin_size = 10
)

p_triangle
p_circle 
p_square

###########

library(patchwork)

plot_mouse_shape <- function(
    all_trials,
    mouse,
    shape,
    percentile = 0.90,
    bin_size = 1
) {
  
  day1 <- get_C_raw(
    all_trials = all_trials,
    mouse = mouse,
    shape = shape,
    day = 1,
    percentile = percentile,
    bin_size = bin_size
  )
  
  day2 <- get_C_raw(
    all_trials = all_trials,
    mouse = mouse,
    shape = shape,
    day = 2,
    percentile = percentile,
    bin_size = bin_size
  )
  
  n_neurons <- min(nrow(day1), nrow(day2))
  
  day1 <- day1[seq_len(n_neurons), , drop = FALSE]
  day2 <- day2[seq_len(n_neurons), , drop = FALSE]
  
  combined <- cbind(day1, day2)
  
  df <- data.frame(
    neuron = rep(
      seq_len(nrow(combined)),
      times = ncol(combined)
    ),
    time = rep(
      seq_len(ncol(combined)),
      each = nrow(combined)
    ),
    value = as.vector(combined)
  )
  
  day_boundary <- ncol(day1) + 0.5
  
  ggplot(
    df,
    aes(x = time, y = neuron, fill = value)
  ) +
    geom_raster() +
    geom_vline(
      xintercept = day_boundary,
      linewidth = 0.6,
      colour = "red"
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 1),
      name = if (bin_size == 1) {
        "Activity"
      } else {
        "Active\nproportion"
      }
    ) +
    scale_y_reverse() +
    labs(
      title = paste(
        mouse,
        "–",
        tools::toTitleCase(shape)
      ),
      subtitle = paste0(
        "C_raw > neuron-specific ",
        percentile * 100,
        "th percentile"
      ),
      x = if (bin_size == 1) "Time point" else "Time bin",
      y = "Neuron"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}


plot_six_mice <- function(
    all_trials,
    shape,
    percentile = 0.90,
    bin_size = 1,
    ncol = 2
) {
  
  mice <- names(all_trials)
  
  plots <- lapply(mice, function(mouse) {
    
    plot_mouse_shape(
      all_trials = all_trials,
      mouse = mouse,
      shape = shape,
      percentile = percentile,
      bin_size = bin_size
    )
  })
  
  wrap_plots(
    plots,
    ncol = ncol,
    guides = "collect"
  ) +
    plot_annotation(
      title = paste(
        tools::toTitleCase(shape),
        "binary activity for all mice"
      ),
      subtitle = paste0(
        "Active when C_raw exceeds each neuron's ",
        percentile * 100,
        "th percentile"
      )
    ) &
    theme(legend.position = "right")
}


triangle_six_mice <- plot_six_mice(
  all_trials,
  shape = "triangle",
  percentile = 0.90,
  bin_size = 10
)

triangle_six_mice

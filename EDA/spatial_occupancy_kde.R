# KERNEL DENSITY ESTIMATION OF SPATIAL OCCUPANCY
#
# Single mouse, six trials
#
# Day 1:
#   Trial 1 = triangular cage
#   Trial 2 = round cage
#   Trial 3 = square cage
#
# Day 2:
#   Trial 4 = round cage
#   Trial 5 = square cage
#   Trial 6 = triangular cage
#
# The KDE is weighted by the time spent at each recorded position.


# 1. Packages --------------------------------------------------

library(tidyverse)
library(purrr)
library(tibble)

if (!requireNamespace("ks", quietly = TRUE)) {
  install.packages("ks")
}


# 2. Load the dataset ------------------------------------------

setwd("/Users/lucacaldera/Documents/Erre/DRC")

all_trials <- readRDS("all_trials.RDS")

mice   <- names(all_trials)
trials <- names(all_trials[[1]])

mice


# 3. User-defined parameters -----------------------------------

mouse_id <- "M3425F"
# "M3411"  "M3412"  "M3421F"  "M3422F"  "M3424F"  "M3425F"

# Kernel bandwidth in spatial-coordinate units.
# Lower values preserve more local detail, while higher values
# produce a smoother density surface.
kde_bandwidth <- 12

# Number of grid points along each axis.
grid_n <- 180

# Hide the weakest tails of the KDE in the final plot.
# A value of 0.005 removes values below 0.5% of each trial's maximum.
tail_threshold <- 0.005

# Cap unusually long time intervals, which may be caused by missing frames.
max_gap_factor <- 3


# 4. Validate the selected mouse -------------------------------

if (!mouse_id %in% names(all_trials)) {
  stop(
    paste0(
      "Mouse not found. Available mice: ",
      paste(names(all_trials), collapse = ", ")
    )
  )
}

trial_names <- names(all_trials[[mouse_id]])

if (length(trial_names) < 6) {
  stop("The selected mouse contains fewer than six trials.")
}

print(trial_names)


# 5. Trial metadata --------------------------------------------

trial_info <- tibble(
  trial_number = 1:6,
  
  trial = paste0("trial", 1:6),
  
  day = factor(
    c(
      "Day 1",
      "Day 1",
      "Day 1",
      "Day 2",
      "Day 2",
      "Day 2"
    ),
    levels = c("Day 1", "Day 2")
  ),
  
  cage = factor(
    c(
      "Triangular",  # Trial 1
      "Round",       # Trial 2
      "Square",      # Trial 3
      "Round",       # Trial 4
      "Square",      # Trial 5
      "Triangular"   # Trial 6
    ),
    levels = c(
      "Triangular",
      "Round",
      "Square"
    )
  )
)

print(trial_info)


# Check that every trial listed above is present in the dataset.
missing_trials <- setdiff(
  trial_info$trial,
  trial_names
)

if (length(missing_trials) > 0) {
  stop(
    paste0(
      "Trials not found: ",
      paste(missing_trials, collapse = ", ")
    )
  )
}


# 6. Extract position data -------------------------------------

position_df <- purrr::map_dfr(
  seq_len(nrow(trial_info)),
  
  function(i) {
    
    trial_id <- trial_info$trial[[i]]
    trial_number_id <- trial_info$trial_number[[i]]
    
    behav <- all_trials[[mouse_id]][[trial_id]]$behav
    
    position <- as.data.frame(behav$position)
    
    if (ncol(position) < 2) {
      stop(
        paste0(
          trial_id,
          ": behav$position must contain at least two columns."
        )
      )
    }
    
    if (length(behav$time) != nrow(position)) {
      stop(
        paste0(
          trial_id,
          ": the length of behav$time does not match ",
          "the number of rows in behav$position."
        )
      )
    }
    
    # Keep invalid coordinates at this stage so that time intervals
    # are calculated from the original sequence of timestamps.
    df <- tibble(
      time = as.numeric(behav$time),
      x = as.numeric(position[[1]]),
      y = as.numeric(position[[2]])
    ) |>
      filter(is.finite(time)) |>
      arrange(time)
    
    if (nrow(df) < 2) {
      stop(
        paste0(
          trial_id,
          ": not enough valid time samples are available."
        )
      )
    }
    
    # Time between consecutive samples, in milliseconds.
    time_difference_ms <- diff(df$time)
    
    typical_interval_ms <- median(
      time_difference_ms[
        is.finite(time_difference_ms) &
          time_difference_ms > 0
      ],
      na.rm = TRUE
    )
    
    if (!is.finite(typical_interval_ms)) {
      stop(
        paste0(
          trial_id,
          ": unable to calculate a valid sampling interval."
        )
      )
    }
    
    max_interval_ms <- typical_interval_ms * max_gap_factor
    
    df |>
      mutate(
        # Assign the median interval to the final sample.
        dt_ms = c(
          time_difference_ms,
          typical_interval_ms
        ),
        
        # Replace zero, negative, or non-finite intervals.
        dt_ms = if_else(
          is.finite(dt_ms) & dt_ms > 0,
          dt_ms,
          typical_interval_ms
        ),
        
        # Prevent a sample after a recording gap from receiving
        # an unrealistically large time weight.
        dt_ms = pmin(
          dt_ms,
          max_interval_ms
        ),
        
        dt_sec = dt_ms / 1000
      ) |>
      filter(
        is.finite(x),
        is.finite(y),
        is.finite(dt_sec),
        dt_sec > 0
      ) |>
      transmute(
        mouse = mouse_id,
        trial_number = trial_number_id,
        trial = trial_id,
        day = trial_info$day[[i]],
        cage = trial_info$cage[[i]],
        time = time,
        x = x,
        y = y,
        dt_sec = dt_sec
      )
  }
)


# 7. Summarize the extracted data ------------------------------

trial_summary <- position_df |>
  group_by(
    mouse,
    day,
    trial_number,
    trial,
    cage
  ) |>
  summarise(
    n_positions = n(),
    duration_sec = sum(dt_sec),
    duration_min = duration_sec / 60,
    min_x = min(x),
    max_x = max(x),
    min_y = min(y),
    max_y = max(y),
    .groups = "drop"
  ) |>
  arrange(trial_number)

print(trial_summary)


# 8. Define common grid limits ---------------------------------
#
# Every trial is estimated on the same grid, making the density
# surfaces directly comparable across panels.

global_x_range <- range(
  position_df$x,
  na.rm = TRUE
)

global_y_range <- range(
  position_df$y,
  na.rm = TRUE
)

x_padding <- diff(global_x_range) * 0.03
y_padding <- diff(global_y_range) * 0.03

kde_min <- c(
  global_x_range[[1]] - x_padding,
  global_y_range[[1]] - y_padding
)

kde_max <- c(
  global_x_range[[2]] + x_padding,
  global_y_range[[2]] + y_padding
)


# 9. Common bandwidth matrix -----------------------------------
#
# Use the same amount of smoothing on both axes and in every trial.

H_common <- diag(
  rep(kde_bandwidth^2, 2)
)

print(H_common)


# 10. Calculate the time-weighted KDE --------------------------

kde_df <- position_df |>
  group_by(
    mouse,
    day,
    cage,
    trial_number,
    trial
  ) |>
  group_split(.keep = TRUE) |>
  purrr::map_dfr(
    
    function(trial_data) {
      
      xy <- trial_data |>
        select(x, y) |>
        as.matrix()
      
      total_time_sec <- sum(
        trial_data$dt_sec,
        na.rm = TRUE
      )
      
      # The KDE function expects weights that sum to one.
      time_weights <- trial_data$dt_sec / total_time_sec
      
      kde_fit <- ks::kde(
        x = xy,
        
        # Apply the same bandwidth to every trial.
        H = H_common,
        
        # Weight each position by the time spent there.
        w = time_weights,
        
        # Estimate all trials over the same spatial grid.
        xmin = kde_min,
        xmax = kde_max,
        
        gridsize = c(
          grid_n,
          grid_n
        ),
        
        bgridsize = c(
          grid_n,
          grid_n
        ),
        
        # Binning improves performance for large datasets.
        binned = TRUE,
        
        # Return the density estimate while avoiding contour calculations.
        density = TRUE,
        
        compute.cont = FALSE
      )
      
      x_grid <- kde_fit$eval.points[[1]]
      y_grid <- kde_fit$eval.points[[2]]
      
      # expand.grid preserves the ordering used by
      # as.vector(kde_fit$estimate).
      kde_grid <- expand.grid(
        x = x_grid,
        y = y_grid,
        KEEP.OUT.ATTRS = FALSE
      ) |>
        as_tibble()
      
      kde_grid$density <- as.vector(
        kde_fit$estimate
      )
      
      # Width and height of one grid cell.
      dx <- mean(diff(x_grid))
      dy <- mean(diff(y_grid))
      
      kde_grid |>
        mutate(
          # KDE values are probability densities per unit area.
          # Multiplying by cell area and trial duration gives the
          # estimated occupancy time within each cell.
          occupancy_sec = density *
            dx *
            dy *
            total_time_sec,
          
          mouse = trial_data$mouse[[1]],
          day = trial_data$day[[1]],
          cage = trial_data$cage[[1]],
          trial_number = trial_data$trial_number[[1]],
          trial = trial_data$trial[[1]]
        )
    }
  )


# 11. Remove weak tails from the plot --------------------------
#
# A Gaussian KDE assigns very small values even far from the cage.
# These values are set to NA for plotting only; the original density
# values are retained for the contour lines.

kde_df <- kde_df |>
  group_by(
    mouse,
    day,
    cage,
    trial_number,
    trial
  ) |>
  mutate(
    panel_max = max(
      occupancy_sec,
      na.rm = TRUE
    ),
    
    occupancy_plot = if_else(
      occupancy_sec >=
        panel_max * tail_threshold,
      occupancy_sec,
      NA_real_
    )
  ) |>
  ungroup()


# 12. Trial labels ---------------------------------------------

label_df <- position_df |>
  distinct(
    mouse,
    day,
    cage,
    trial_number,
    trial
  ) |>
  mutate(
    label = paste0(
      "Trial ",
      trial_number
    ),
    
    x_label = kde_min[[1]] +
      0.025 * diff(c(kde_min[[1]], kde_max[[1]])),
    
    y_label = kde_max[[2]] -
      0.025 * diff(c(kde_min[[2]], kde_max[[2]]))
  )


# 13. KDE plot --------------------------------------------------
#
# Columns:
#   triangular | round | square
#
# Rows:
#   Day 1 | Day 2
#
# This layout allows each cage on Day 1 to be compared directly
# with the same cage on Day 2.

x11()

kde_plot <- ggplot(
  kde_df,
  aes(
    x = x,
    y = y
  )
) +
  
  # KDE surface.
  geom_raster(
    aes(
      fill = occupancy_plot
    ),
    interpolate = TRUE
  ) +
  
  # Density contours.
  geom_contour(
    aes(
      z = occupancy_sec
    ),
    bins = 8,
    colour = "white",
    linewidth = 0.25,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  
  # Trial number shown in the upper-left corner of each panel.
  geom_label(
    data = label_df,
    aes(
      x = x_label,
      y = y_label,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 3.5,
    fill = scales::alpha("white", 0.8),
    colour = "black",
    linewidth = 0,
    label.padding = unit(
      0.18,
      "lines"
    )
  ) +
  
  facet_grid(
    rows = vars(day),
    cols = vars(cage),
    drop = FALSE
  ) +
  
  scale_fill_viridis_c(
    name = "Estimated time\nper cell (s)",
    option = "viridis",
    transform = "sqrt",
    na.value = "white"
  ) +
  
  coord_equal(
    xlim = c(
      kde_min[[1]],
      kde_max[[1]]
    ),
    ylim = c(
      kde_min[[2]],
      kde_max[[2]]
    ),
    expand = FALSE
  ) +
  
  labs(
    title = paste(
      "Kernel density of spatial occupancy — Mouse",
      mouse_id
    ),
    subtitle = paste0(
      "Time-weighted KDE; bandwidth = ",
      kde_bandwidth,
      " spatial units"
    ),
    x = "X position",
    y = "Y position"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey70"
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 10.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right"
  )


# 14. Display the plot -----------------------------------------

print(kde_plot)


# 15. Save the figure ------------------------------------------

ggsave(
  filename = paste0(
    "KDE_spatial_occupancy_",
    mouse_id,
    ".png"
  ),
  plot = kde_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

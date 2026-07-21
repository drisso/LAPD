# Prediction metrics, baselines, and diagnostic summaries for hurdle GAMs.

sampled_concurvity <- function(
    model,
    newdata,
    max_rows = 20000L,
    seed = 2026L,
    exclude_pattern = NULL
) {
  if (!inherits(model, "gam")) {
    stop("model must inherit from class 'gam'.", call. = FALSE)
  }
  if (!nrow(newdata)) stop("newdata has no rows.", call. = FALSE)

  # mgcv::concurvity() constructs the complete fitted-model matrix. For a
  # multi-million-row bam fit, its internal integer dimension calculation can
  # overflow. A sampled linear-predictor matrix estimates the same geometric
  # diagnostic while keeping the dimensions bounded.
  if (nrow(newdata) > max_rows) {
    set.seed(seed)
    row_index <- sort(sample.int(nrow(newdata), max_rows))
    diagnostic_data <- newdata[row_index, , drop = FALSE]
  } else {
    diagnostic_data <- newdata
  }
  X_full <- stats::predict(model, newdata = diagnostic_data, type = "lpmatrix")
  valid <- rowSums(is.na(X_full)) == 0
  X_full <- X_full[valid, , drop = FALSE]

  smooth_labels <- vapply(model$smooth, `[[`, character(1), "label")
  keep_smooth <- rep(TRUE, length(model$smooth))
  if (!is.null(exclude_pattern)) {
    keep_smooth <- !grepl(exclude_pattern, smooth_labels)
  }
  kept <- which(keep_smooth)
  if (!length(kept)) stop("No smooth terms remain for the diagnostic.", call. = FALSE)

  first_smooth_column <- min(vapply(model$smooth, `[[`, numeric(1), "first.para"))
  groups <- list()
  if (first_smooth_column > 1L) {
    groups$para <- seq_len(first_smooth_column - 1L)
  }
  for (i in kept) {
    groups[[smooth_labels[i]]] <-
      model$smooth[[i]]$first.para:model$smooth[[i]]$last.para
  }

  selected_columns <- unlist(groups, use.names = FALSE)
  X <- X_full[, selected_columns, drop = FALSE]
  beta <- stats::coef(model)[selected_columns]
  rm(X_full, diagnostic_data)

  if (nrow(X) <= ncol(X)) {
    stop(
      "The concurvity sample must contain more rows than model columns; ",
      "increase max_rows.",
      call. = FALSE
    )
  }
  X <- qr.R(qr(X, tol = 0, LAPACK = FALSE))
  group_lengths <- lengths(groups)
  stops <- cumsum(group_lengths)
  starts <- stops - group_lengths + 1L
  result <- matrix(
    0,
    nrow = 3L,
    ncol = length(groups),
    dimnames = list(c("worst", "observed", "estimate"), names(groups))
  )

  for (i in seq_along(groups)) {
    target <- starts[i]:stops[i]
    other <- setdiff(seq_len(ncol(X)), target)
    Xi <- X[, other, drop = FALSE]
    Xj <- X[, target, drop = FALSE]
    r <- ncol(Xi)
    R <- qr.R(qr(cbind(Xi, Xj), LAPACK = FALSE, tol = 0))[
      , -(seq_len(r)), drop = FALSE
    ]
    Rt <- qr.R(qr(R, tol = 0))
    target_beta <- beta[target]
    result["worst", i] <-
      svd(forwardsolve(t(Rt), t(R[seq_len(r), , drop = FALSE])))$d[1L]^2
    denominator <- sum((Rt %*% target_beta)^2)
    result["observed", i] <- if (denominator > 0) {
      sum((R[seq_len(r), , drop = FALSE] %*% target_beta)^2) / denominator
    } else {
      NA_real_
    }
    result["estimate", i] <-
      sum(R[seq_len(r), , drop = FALSE]^2) / sum(R^2)
  }

  attr(result, "sample_rows") <- sum(valid)
  attr(result, "excluded_smooths") <- smooth_labels[!keep_smooth]
  result
}

binary_auc <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  n_positive <- sum(observed == 1L)
  n_negative <- sum(observed == 0L)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks <- rank(predicted, ties.method = "average")
  (sum(ranks[observed == 1L]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

binary_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- pmin(pmax(predicted[keep], 1e-8), 1 - 1e-8)
  c(
    AUC = binary_auc(observed, predicted),
    Brier = mean((observed - predicted)^2),
    `log loss` = -mean(observed * log(predicted) + (1 - observed) * log(1 - predicted))
  )
}

continuous_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  c(
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted))
  )
}

calibration_table <- function(observed, predicted, n_bins = 10L) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  ordering <- order(predicted)
  bin <- integer(length(predicted))
  bin[ordering] <- pmin(n_bins, ceiling(seq_along(ordering) * n_bins / length(ordering)))
  out <- aggregate(
    cbind(observed = observed, predicted = predicted, n = rep(1, length(observed))),
    list(bin = bin),
    function(x) c(mean = mean(x), sum = sum(x))
  )
  data.frame(
    bin = out$bin,
    observed = out$observed[, "mean"],
    predicted = out$predicted[, "mean"],
    n = out$n[, "sum"]
  )
}

fit_neuron_baseline <- function(train_data, prior_weight = 100) {
  global_event_rate <- mean(train_data$event)
  global_positive_mean <- mean(train_data$amplitude[train_data$event == 1L])
  if (!is.finite(global_positive_mean)) global_positive_mean <- 0

  counts <- aggregate(event ~ neuron, train_data, function(x) c(events = sum(x), n = length(x)))
  baseline <- data.frame(
    neuron = counts$neuron,
    events = counts$event[, "events"],
    n = counts$event[, "n"]
  )
  baseline$event_probability <-
    (baseline$events + prior_weight * global_event_rate) /
    (baseline$n + prior_weight)

  positive <- train_data[train_data$event == 1L, , drop = FALSE]
  if (nrow(positive)) {
    positive_means <- aggregate(amplitude ~ neuron, positive, mean)
    baseline$positive_mean <- positive_means$amplitude[
      match(as.character(baseline$neuron), as.character(positive_means$neuron))
    ]
  } else {
    baseline$positive_mean <- NA_real_
  }
  baseline$positive_mean[!is.finite(baseline$positive_mean)] <- global_positive_mean
  baseline$expected_amplitude <- baseline$event_probability * baseline$positive_mean
  baseline
}

add_baseline_predictions <- function(data, baseline) {
  idx <- match(as.character(data$neuron), as.character(baseline$neuron))
  if (anyNA(idx)) stop("Baseline is missing one or more neurons.", call. = FALSE)
  data$baseline_event_probability <- baseline$event_probability[idx]
  data$baseline_positive_mean <- baseline$positive_mean[idx]
  data$baseline_expected_amplitude <- baseline$expected_amplitude[idx]
  data
}

evaluate_hurdle_predictions <- function(data, set_label = "Evaluation") {
  required <- c(
    "event", "amplitude", "event_probability", "positive_mean",
    "expected_amplitude", "baseline_event_probability",
    "baseline_positive_mean", "baseline_expected_amplitude"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Prediction data are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  thresholded_amplitude <- ifelse(data$event == 1L, data$amplitude, 0)
  positive <- data$event == 1L

  event_model <- binary_metrics(data$event, data$event_probability)
  event_baseline <- binary_metrics(data$event, data$baseline_event_probability)
  combined_model <- continuous_metrics(thresholded_amplitude, data$expected_amplitude)
  combined_baseline <- continuous_metrics(
    thresholded_amplitude, data$baseline_expected_amplitude
  )
  positive_model <- continuous_metrics(data$amplitude[positive], data$positive_mean[positive])
  positive_baseline <- continuous_metrics(
    data$amplitude[positive], data$baseline_positive_mean[positive]
  )

  make_rows <- function(values, component, predictor) {
    data.frame(
      set = set_label,
      component = component,
      predictor = predictor,
      metric = names(values),
      value = as.numeric(values),
      row.names = NULL
    )
  }
  rbind(
    make_rows(event_model, "Event occurrence", "Hurdle GAM"),
    make_rows(event_baseline, "Event occurrence", "Neuron baseline"),
    make_rows(positive_model, "Positive amplitude", "Hurdle GAM"),
    make_rows(positive_baseline, "Positive amplitude", "Neuron baseline"),
    make_rows(combined_model, "Thresholded amplitude", "Hurdle GAM"),
    make_rows(combined_baseline, "Thresholded amplitude", "Neuron baseline")
  )
}

sample_diagnostic_rows <- function(data, max_rows = 500000L, seed = 2026L) {
  if (nrow(data) <= max_rows) return(data)
  set.seed(seed)
  # A uniform sample preserves prevalence, calibration, and proper scoring rules.
  data[sort(sample.int(nrow(data), max_rows)), ]
}

population_time_summary <- function(data, frame_rate = 15, bin_seconds = 1) {
  data$time_bin <- floor(data$time_sec / bin_seconds) * bin_seconds
  aggregate(
    cbind(
      observed_event_rate = data$event,
      predicted_event_rate = data$event_probability,
      observed_amplitude = ifelse(data$event == 1L, data$amplitude, 0),
      predicted_amplitude = data$expected_amplitude
    ),
    list(time_sec = data$time_bin),
    mean
  )
}

spatial_fit_summary <- function(data, n_bins = 15L) {
  x_breaks <- seq(min(data$x), max(data$x), length.out = n_bins + 1L)
  y_breaks <- seq(min(data$y), max(data$y), length.out = n_bins + 1L)
  data$x_bin <- cut(data$x, x_breaks, include.lowest = TRUE, labels = FALSE)
  data$y_bin <- cut(data$y, y_breaks, include.lowest = TRUE, labels = FALSE)
  out <- aggregate(
    cbind(
      x = data$x,
      y = data$y,
      observed = data$event,
      predicted = data$event_probability,
      n = rep(1, nrow(data))
    ),
    list(x_bin = data$x_bin, y_bin = data$y_bin),
    function(z) c(mean = mean(z), sum = sum(z))
  )
  data.frame(
    x = out$x[, "mean"], y = out$y[, "mean"],
    observed = out$observed[, "mean"],
    predicted = out$predicted[, "mean"],
    n = out$n[, "sum"]
  )
}

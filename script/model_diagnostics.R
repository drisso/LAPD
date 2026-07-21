# Prediction metrics, baselines, and diagnostic summaries for hurdle GAMs.

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

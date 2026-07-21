# Hurdle generalized additive model for deconvolved calcium events.

default_hurdle_formulas <- function(
    position_k = 50L,
    border_k = 8L,
    speed_k = 8L,
    time_k = 12L
) {
  shared <- paste0(
    "s(neuron, bs = 're') + ",
    "s(x_z, y_z, k = ", position_k, ") + ",
    "s(border_z, k = ", border_k, ") + ",
    "s(speed_z, k = ", speed_k, ") + ",
    "s(time_z, k = ", time_k, ")"
  )
  list(
    occurrence = stats::as.formula(paste("event ~", shared)),
    positive = stats::as.formula(paste("amplitude ~", shared))
  )
}

fit_hurdle_gam <- function(
    train_data,
    formulas = default_hurdle_formulas(),
    nthreads = 2L,
    method = "fREML",
    discrete = TRUE,
    select = TRUE
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required to fit the hurdle GAM.", call. = FALSE)
  }
  required <- c(
    "event", "amplitude", "neuron", "x_z", "y_z", "border_z",
    "speed_z", "time_z"
  )
  missing <- setdiff(required, names(train_data))
  if (length(missing)) {
    stop("Training data are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(unique(train_data$event)) < 2L) {
    stop("The occurrence model requires both event and non-event observations.", call. = FALSE)
  }
  positive_data <- train_data[train_data$event == 1L & train_data$amplitude > 0, ]
  if (nrow(positive_data) < 100L) {
    stop(
      "Fewer than 100 positive training events remain after thresholding; ",
      "use a less stringent threshold or a simpler positive-part model.",
      call. = FALSE
    )
  }

  occurrence <- mgcv::bam(
    formula = formulas$occurrence,
    data = train_data,
    family = stats::binomial(link = "logit"),
    method = method,
    discrete = discrete,
    select = select,
    nthreads = nthreads,
    drop.unused.levels = FALSE,
    gc.level = 1L
  )
  positive <- mgcv::bam(
    formula = formulas$positive,
    data = positive_data,
    family = stats::Gamma(link = "log"),
    method = method,
    discrete = discrete,
    select = select,
    nthreads = nthreads,
    drop.unused.levels = FALSE,
    gc.level = 1L
  )

  structure(
    list(
      occurrence = occurrence,
      positive = positive,
      formulas = formulas,
      n_train = nrow(train_data),
      n_positive = nrow(positive_data),
      neuron_levels = levels(train_data$neuron)
    ),
    class = "hurdle_gam"
  )
}

predict.hurdle_gam <- function(object, newdata, ...) {
  newdata$neuron <- factor(newdata$neuron, levels = object$neuron_levels)
  if (anyNA(newdata$neuron)) {
    stop("newdata contains neurons that were not present during training.", call. = FALSE)
  }
  event_probability <- as.numeric(stats::predict(
    object$occurrence, newdata = newdata, type = "response", ...
  ))
  positive_mean <- as.numeric(stats::predict(
    object$positive, newdata = newdata, type = "response", ...
  ))
  event_probability <- pmin(pmax(event_probability, 0), 1)
  positive_mean <- pmax(positive_mean, 0)

  data.frame(
    event_probability = event_probability,
    positive_mean = positive_mean,
    expected_amplitude = event_probability * positive_mean
  )
}

add_hurdle_predictions <- function(model, data) {
  predictions <- predict(model, newdata = data)
  data$event_probability <- predictions$event_probability
  data$positive_mean <- predictions$positive_mean
  data$expected_amplitude <- predictions$expected_amplitude
  data
}

print.hurdle_gam <- function(x, ...) {
  cat("Hurdle GAM\n")
  cat("  training observations:", format(x$n_train, big.mark = ","), "\n")
  cat("  positive observations:", format(x$n_positive, big.mark = ","), "\n")
  cat("\nOccurrence formula:\n  ")
  print(x$formulas$occurrence)
  cat("\nPositive-part formula:\n  ")
  print(x$formulas$positive)
  invisible(x)
}

hurdle_model_summary <- function(model) {
  occurrence_summary <- summary(model$occurrence)
  positive_summary <- summary(model$positive)
  data.frame(
    component = c("Event occurrence", "Positive amplitude"),
    family = c("Binomial-logit", "Gamma-log"),
    observations = c(model$n_train, model$n_positive),
    adjusted_r_squared = c(
      unname(occurrence_summary$r.sq),
      unname(positive_summary$r.sq)
    ),
    deviance_explained = c(
      unname(occurrence_summary$dev.expl),
      unname(positive_summary$dev.expl)
    ),
    stringsAsFactors = FALSE
  )
}

smooth_significance_table <- function(model) {
  extract <- function(fitted_model, component) {
    table <- as.data.frame(summary(fitted_model)$s.table)
    table$term <- rownames(table)
    rownames(table) <- NULL
    table$component <- component
    table[, c("component", "term", setdiff(names(table), c("component", "term")))]
  }
  rbind(
    extract(model$occurrence, "Event occurrence"),
    extract(model$positive, "Positive amplitude")
  )
}

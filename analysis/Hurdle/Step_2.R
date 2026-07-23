# ==============================================================================
# HURDLE MODEL - STEP 2: POSITIVE AMPLITUDE STAGE   S_k(t) | S_k(t) > 0
# ==============================================================================
# Le ampiezze positive sono assunte log-normali:
#     log(S_k(t)) | S_k(t) > 0  ~  N(eta_k(t), sigma_k^2)
# stimate con LASSO gaussiano sui soli frame con evento, riusando dallo Step 1:
#   - la design matrix salvata (fit_t1$prep_data$X);
#   - le colonne ammesse per ciascun neurone (fit_t1$valid_cols_list), che
#     contengono gia' l'esclusione autoreferenziale al lag 0;
#   - gli stessi penalty.factor (distanza fra centroidi, floor 0.05);
#   - la stessa griglia di lambda e gli stessi fold temporali a blocchi,
#     ristretti ai frame positivi.
#
# Questo file contiene SOLO definizioni di funzioni. Caricamento dati,
# esecuzione e salvataggi stanno in main.R. Richiede che Step_1.R sia gia'
# stato caricato (usa prepare_design_matrix() per la previsione out-of-sample).
#
# Contenuto:
#   1. Helper e validazione della struttura Step 1
#   2. fit_positive_amplitude_stage()      stima dello stage positivo
#   3. Previsione hurdle in-sample e out-of-sample
#   4. Combinazione dei due stage (soft / hard / calibrata)
#   5. evaluate_hurdle_predictions()       metriche numeriche
# ==============================================================================

library(glmnet)
library(Matrix)


# ==============================================================================
# 1. HELPER E VALIDAZIONE DELLA STRUTTURA STEP 1
# ==============================================================================

# Lambda usata dallo Step 1 per coefficienti e previsioni.
selected_lambda_stage1 <- function(model) {
  if (inherits(model, "cv.glmnet")) {
    if (!is.finite(model$lambda.1se)) stop("Invalid lambda.1se in a Step 1 model.")
    return(model$lambda.1se)
  }
  model$lambda[mid_lambda_index(model$lambda)]
}

# Lo Step 1 corrente salva le colonne per neurone in valid_cols_list.
# valid_cols resta come fallback per bundle piu' vecchi.
stage1_valid_cols_field <- function(first_fit) {
  if (!is.null(first_fit$valid_cols_list)) return("valid_cols_list")
  if (!is.null(first_fit$valid_cols)) return("valid_cols")
  stop("Step 1 contains neither valid_cols_list nor valid_cols.")
}

stage1_valid_cols_for_neuron <- function(first_fit, neuron_id) {
  field_name <- stage1_valid_cols_field(first_fit)
  valid <- first_fit[[field_name]]
  if (is.list(valid)) {
    if (length(valid) < neuron_id || is.null(valid[[neuron_id]])) {
      stop("Missing Step 1 valid columns for neuron ", neuron_id, ".")
    }
    valid <- valid[[neuron_id]]
  }

  valid <- sort(unique(as.integer(valid)))
  p <- ncol(first_fit$prep_data$X)
  if (length(valid) == 0L || any(!is.finite(valid)) || any(valid < 1L | valid > p)) {
    stop("fit_t1$", field_name, " contains invalid indices for neuron ", neuron_id, ".")
  }
  valid
}

# Numero di basi spline usato nello Step 1: serve per ricostruire la design
# matrix su un trial nuovo. Le covariate fisse sono n_splines^2 + 1 (velocita').
stage1_n_splines <- function(first_fit) {
  n_splines <- first_fit$prep_data$n_splines
  if (is.null(n_splines)) {
    n_splines <- round(sqrt(first_fit$prep_data$n_fixed_cols - 1L))
  }
  as.integer(n_splines)
}

# Recupera la griglia di lambda effettivamente salvata nei modelli Step 1,
# invece di riscriverla a mano nello Step 2.
extract_stage1_lambda_grid <- function(
    first_fit,
    fallback = exp(seq(log(1), log(0.0001), length.out = 20L))) {

  available <- which(!vapply(first_fit$models, is.null, logical(1)))
  if (length(available) == 0L) {
    warning("No Step 1 model is available; using the fallback lambda grid.")
    return(as.numeric(fallback))
  }

  grids <- lapply(first_fit$models[available], function(model) as.numeric(model$lambda))
  reference <- grids[[1L]]
  same_grid <- vapply(
    grids,
    function(grid_i) {
      length(grid_i) == length(reference) &&
        isTRUE(all.equal(grid_i, reference, tolerance = 1e-12, check.attributes = FALSE))
    },
    logical(1)
  )

  if (!all(same_grid)) {
    warning(
      "Not all Step 1 models store exactly the same lambda path. ",
      "Step 2 will use the path from the first available Step 1 model."
    )
  }

  reference
}

# Attiva il parallelo solo se un backend foreach e' gia' registrato.
resolve_parallel_cv <- function(requested = TRUE) {
  if (!isTRUE(requested)) return(FALSE)

  registered <- FALSE
  if (requireNamespace("foreach", quietly = TRUE)) {
    registered <- isTRUE(foreach::getDoParRegistered())
  }

  if (!registered) {
    message(
      "parallel_cv=TRUE, but no foreach backend is registered. ",
      "Step 2 cross-validation will run in serial mode."
    )
  }
  registered
}

validate_stage1_structure <- function(first_fit) {
  required_fit <- c("models", "prep_data", "active_neurons", "foldid")
  missing_fit <- setdiff(required_fit, names(first_fit))
  if (length(missing_fit) > 0L) {
    stop("Missing fields in fit_t1: ", paste(missing_fit, collapse = ", "))
  }
  valid_cols_field <- stage1_valid_cols_field(first_fit)

  prep <- first_fit$prep_data
  required_prep <- c(
    "X", "Z", "K", "T_frames", "dt_ms", "centroids",
    "n_fixed_cols", "max_lag_frames"
  )
  missing_prep <- setdiff(required_prep, names(prep))
  if (length(missing_prep) > 0L) {
    stop("Missing fields in fit_t1$prep_data: ", paste(missing_prep, collapse = ", "))
  }

  if (length(first_fit$foldid) != prep$T_frames) {
    stop("fit_t1$foldid does not have length T_frames.")
  }
  if (any(!is.finite(first_fit$foldid))) {
    stop("fit_t1$foldid contains non-finite values.")
  }

  X <- prep$X
  if (nrow(X) != prep$T_frames) {
    stop("fit_t1$prep_data$X does not have T_frames rows.")
  }
  if (nrow(prep$Z) != prep$K || ncol(prep$Z) != prep$T_frames) {
    stop("fit_t1$prep_data$Z is inconsistent with K and T_frames.")
  }
  if (prep$n_fixed_cols < 1L || prep$n_fixed_cols > ncol(X)) {
    stop("fit_t1$prep_data$n_fixed_cols is invalid.")
  }

  # Design matrix attesa: [fisse | lag 0 | lag 1..L]
  expected_total_cols <- prep$n_fixed_cols + prep$K * (prep$max_lag_frames + 1L)
  if (ncol(X) != expected_total_cols) {
    stop(
      "The Step 1 design matrix is not compatible with lag 0 plus lags 1..L. ",
      "Expected ", expected_total_cols, " columns = ", prep$n_fixed_cols,
      " fixed + ", prep$K, " x (", prep$max_lag_frames, " + 1) network; found ",
      ncol(X), "."
    )
  }

  expected_network_names <- unlist(
    lapply(0:prep$max_lag_frames, function(lag) paste0("N", seq_len(prep$K), "_lag", lag)),
    use.names = FALSE
  )
  observed_network_names <- colnames(X)[(prep$n_fixed_cols + 1L):ncol(X)]
  if (!is.null(observed_network_names) && !identical(observed_network_names, expected_network_names)) {
    warning(
      "The Step 1 network-column names do not exactly match the expected ordering ",
      "lag 0, lag 1, ..., lag L. Column positions from the saved design matrix ",
      "will still be treated as authoritative."
    )
  }

  # Al lag 0, Z_i(t) non deve predire il neurone i stesso.
  for (i in sort(unique(as.integer(first_fit$active_neurons)))) {
    if (i < 1L || i > prep$K || is.null(first_fit$models[[i]])) next
    valid_i <- stage1_valid_cols_for_neuron(first_fit, i)
    self_lag0_col <- prep$n_fixed_cols + i
    if (self_lag0_col %in% valid_i) {
      stop(
        "Invalid Step 1 structure: neuron ", i,
        " retains its own contemporaneous lag-0 predictor."
      )
    }
  }

  message("Step 1 valid-column field detected: fit_t1$", valid_cols_field)
  invisible(TRUE)
}

# Verifica che la matrice di ampiezze S sia quella usata nello Step 1.
check_S_matches_stage1 <- function(S, first_fit) {
  prep <- first_fit$prep_data
  if (!is.matrix(S)) S <- as.matrix(S)

  if (!identical(dim(S), c(as.integer(prep$K), as.integer(prep$T_frames)))) {
    stop(
      "S dimensions do not match Step 1. Expected ",
      prep$K, " x ", prep$T_frames, "; found ", paste(dim(S), collapse = " x "), "."
    )
  }

  Z_from_S <- ifelse(is.finite(S) & S > 0, 1, 0)
  n_mismatch <- sum(Z_from_S != as.matrix(prep$Z), na.rm = TRUE)
  if (n_mismatch > 0L) {
    stop(
      "I(S>0) differs from fit_t1$prep_data$Z in ", n_mismatch,
      " cells. The selected dataset is not the one used in Step 1."
    )
  }

  S
}


# ==============================================================================
# 2. FIT DELLO STAGE POSITIVO (LOG-NORMALE)
# ==============================================================================
# Penalizzazione delle covariate fisse nello Step 2.
#
# Nello Step 1 le covariate fisse sono non penalizzate: 17 parametri liberi su
# ~9000 frame. Nello Step 2 le osservazioni scendono a poche centinaia (i soli
# frame con evento) e le 16 basi spline prodotto sono fortemente collineari:
# lasciarle libere produce coefficienti enormi e quasi cancellanti sul
# sottocampione, che esplodono in exp(eta) quando si predice sui frame senza
# evento. Per questo di default le basi spaziali vengono penalizzate, mentre la
# velocita' (una sola colonna, ben condizionata) resta libera.
build_fixed_penalty <- function(
    fixed_names,
    penalize_spatial_basis = TRUE,
    spatial_penalty_factor = 1,
    unpenalized_fixed = "speed_spline") {

  pen <- rep(0, length(fixed_names))
  if (!isTRUE(penalize_spatial_basis)) return(pen)

  is_spatial <- grepl("^spatial_basis_", fixed_names)
  if (!any(is_spatial)) {
    # Nomi non standard: si penalizza tutto tranne quanto elencato esplicitamente.
    is_spatial <- rep(TRUE, length(fixed_names))
  }
  is_spatial <- is_spatial & !(fixed_names %in% unpenalized_fixed)

  pen[is_spatial] <- spatial_penalty_factor
  pen
}

fit_positive_amplitude_stage <- function(
    first_fit,
    S,
    positive_family = "lognormal",
    min_positive = 15L,
    # lambda_grid = NULL: glmnet calcola lambda_max dai dati di CIASCUN neurone.
    # E' la scelta corretta per lo stage gaussiano, la cui scala di lambda non
    # ha nulla a che vedere con quella della devianza binomiale dello Step 1.
    lambda_grid = NULL,
    nlambda = 100L,
    lambda_min_ratio = NULL,
    penalize_spatial_basis = TRUE,
    spatial_penalty_factor = 1,
    unpenalized_fixed = "speed_spline",
    distance_floor = 0.05,
    maxit = 100000L,
    parallel = FALSE) {

  if (!identical(positive_family, "lognormal")) {
    stop("This version implements positive_family='lognormal'.")
  }

  validate_stage1_structure(first_fit)

  prep <- first_fit$prep_data
  X_full <- as.matrix(prep$X)
  K <- prep$K
  T_frames <- prep$T_frames
  n_fixed <- prep$n_fixed_cols
  max_lag <- prep$max_lag_frames
  lag_values <- 0:max_lag
  n_network_lags <- length(lag_values)
  stage1_foldid <- as.integer(first_fit$foldid)
  nfolds_stage1 <- length(unique(stage1_foldid))

  S <- check_S_matches_stage1(S, first_fit)
  Z_from_S <- ifelse(is.finite(S) & S > 0, 1, 0)

  # La design matrix dello Step 1 e le colonne ammesse per neurone sono
  # autoritative: stessi effetti fissi, stesso lag 0 e stessi lag 1..L, con la
  # colonna lag-0 del neurone target gia' esclusa.
  fixed_cols <- seq_len(n_fixed)

  dist_matrix <- as.matrix(dist(prep$centroids))
  max_dist <- max(dist_matrix)
  if (!is.finite(max_dist) || max_dist <= 0) max_dist <- 1

  models <- vector("list", K)
  valid_cols_by_neuron <- vector("list", K)
  full_coefs_by_neuron <- vector("list", K)
  foldid_by_neuron <- vector("list", K)

  selected_lambda <- rep(NA_real_, K)
  lambda_max_by_neuron <- rep(NA_real_, K)
  lambda_min_by_neuron <- rep(NA_real_, K)
  dispersion <- rep(NA_real_, K)
  n_positive <- rowSums(Z_from_S)
  status <- rep("inactive_stage1", K)
  error_message <- rep(NA_character_, K)
  cv_error_message <- rep(NA_character_, K)

  active_declared <- sort(unique(as.integer(first_fit$active_neurons)))
  active_declared <- active_declared[active_declared >= 1L & active_declared <= K]
  model_available <- which(!vapply(first_fit$models, is.null, logical(1)))
  model_available <- model_available[
    vapply(first_fit$models[model_available], glmnet_fit_is_usable, logical(1))
  ]
  active_stage1 <- intersect(active_declared, model_available)
  missing_stage1_model <- setdiff(active_declared, model_available)

  status[active_stage1] <- "too_few_positive_events"
  status[missing_stage1_model] <- "stage1_fit_missing"

  # Coefficienti di rete dello stage positivo: [sorgente x target x lag 0..L]
  delta_array <- array(
    0,
    dim = c(K, K, n_network_lags),
    dimnames = list(
      source = paste0("N", seq_len(K)),
      target = paste0("N", seq_len(K)),
      lag = paste0("lag", lag_values)
    )
  )

  targets <- intersect(active_stage1, which(n_positive >= min_positive))
  status[targets] <- "pending"
  n_targets <- length(targets)
  fit_start <- Sys.time()

  fixed_names <- colnames(X_full)[fixed_cols]
  if (is.null(fixed_names)) fixed_names <- paste0("X", fixed_cols)

  fixed_penalty <- build_fixed_penalty(
    fixed_names            = fixed_names,
    penalize_spatial_basis = penalize_spatial_basis,
    spatial_penalty_factor = spatial_penalty_factor,
    unpenalized_fixed      = unpenalized_fixed
  )

  data_driven_lambda <- is.null(lambda_grid)

  message(
    sprintf(
      "Starting Step 2: %d neurons, lognormal model, %d Step 1 temporal folds, %s.",
      n_targets, nfolds_stage1,
      if (data_driven_lambda) {
        sprintf("lambda path computed per neuron from the data (nlambda=%d)", nlambda)
      } else {
        sprintf(
          "%d fixed lambda [%.1e, %.1e]",
          length(lambda_grid), max(lambda_grid), min(lambda_grid)
        )
      }
    )
  )
  message(
    "Fixed covariates - unpenalized: ",
    paste(fixed_names[fixed_penalty == 0], collapse = ", "),
    " | penalized (factor ", spatial_penalty_factor, "): ",
    sum(fixed_penalty > 0)
  )
  if (length(missing_stage1_model) > 0L) {
    message(
      "Active neurons without a Step 1 fit, excluded from Step 2: ",
      paste(missing_stage1_model, collapse = ", ")
    )
  }
  flush.console()

  for (target_index in seq_along(targets)) {
    i <- targets[target_index]
    neuron_start <- Sys.time()
    idx_pos <- which(is.finite(S[i, ]) & S[i, ] > 0)
    y_positive <- as.numeric(S[i, idx_pos])

    elapsed_minutes <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
    message(
      sprintf(
        "[%d/%d] Neuron %d - positive events: %d - elapsed: %.1f min",
        target_index, n_targets, i, length(idx_pos), elapsed_minutes
      )
    )
    flush.console()

    if (length(idx_pos) < min_positive) {
      status[i] <- "too_few_positive_events"
      message("    skipped - too few positive events")
      next
    }
    if (length(unique(y_positive)) < 2L) {
      status[i] <- "constant_positive_amplitude"
      message("    skipped - constant positive amplitude")
      next
    }

    # Stessi predittori tenuti dallo Step 1 per questo neurone. La risposta e'
    # ristretta alle ampiezze positive, ma le colonne non vengono ri-selezionate:
    # cosi' Step 1 e Step 2 restano direttamente confrontabili.
    valid_i <- stage1_valid_cols_for_neuron(first_fit, i)
    if (!all(fixed_cols %in% valid_i)) {
      stop("Step 1 fixed covariates are missing from its valid-column set for neuron ", i, ".")
    }
    valid_cols_by_neuron[[i]] <- valid_i

    X_i <- X_full[idx_pos, valid_i, drop = FALSE]
    dimnames(X_i) <- list(NULL, paste0("V", seq_len(ncol(X_i))))
    X_i <- Matrix(X_i, sparse = TRUE)

    # Stesso ordinamento dello Step 1 - effetti fissi, lag 0, lag 1..L - e
    # stessi pesi di distanza sulla parte di rete. Le covariate fisse seguono
    # invece fixed_penalty (vedi build_fixed_penalty).
    raw_pen_lag0 <- dist_matrix[i, ] / max_dist
    pen_lag0 <- pmax(raw_pen_lag0, distance_floor)

    raw_pen_lags <- rep(dist_matrix[i, ] / max_dist, times = max_lag)
    pen_lags <- pmax(raw_pen_lags, distance_floor)

    penalty_full <- c(fixed_penalty, pen_lag0, pen_lags)
    if (length(penalty_full) != ncol(X_full)) {
      stop(
        "Penalty/design mismatch for neuron ", i, ": expected ", ncol(X_full),
        " penalty factors, constructed ", length(penalty_full), "."
      )
    }
    penalty_i <- penalty_full[valid_i]

    y_fit <- log(y_positive)

    # Fold temporali a blocchi dello Step 1 ristretti alle osservazioni positive.
    # La CV gira solo se tutti i fold originali sono rappresentati; altrimenti si
    # usa lo stesso percorso di lambda senza cross-validation.
    fold_i <- stage1_foldid[idx_pos]
    represented_folds <- sort(unique(fold_i))
    all_stage1_folds <- sort(unique(stage1_foldid))
    foldid_by_neuron[[i]] <- fold_i
    can_run_cv <- identical(represented_folds, all_stage1_folds)

    # Argomenti del percorso di lambda: con lambda_grid = NULL si lascia che
    # glmnet ricavi lambda_max dai dati del neurone.
    lambda_args <- if (data_driven_lambda) {
      c(
        list(nlambda = nlambda),
        if (!is.null(lambda_min_ratio)) list(lambda.min.ratio = lambda_min_ratio)
      )
    } else {
      list(lambda = lambda_grid)
    }

    common_args <- c(
      list(
        x = X_i,
        y = y_fit,
        family = "gaussian",
        penalty.factor = penalty_i,
        alpha = 1,
        standardize = TRUE,
        intercept = TRUE,
        maxit = maxit
      ),
      lambda_args
    )

    fit_i <- NULL
    if (can_run_cv) {
      fit_i <- tryCatch(
        do.call(
          cv.glmnet,
          c(
            common_args,
            list(foldid = fold_i, type.measure = "mse", parallel = parallel)
          )
        ),
        error = function(e_cv) {
          cv_error_message[i] <<- conditionMessage(e_cv)
          NULL
        }
      )
    } else {
      cv_error_message[i] <- paste0(
        "Only ", length(represented_folds), " of ", length(all_stage1_folds),
        " Step 1 temporal folds contain positive events."
      )
    }

    if (is.null(fit_i)) {
      fit_i <- tryCatch(
        do.call(glmnet, common_args),
        error = function(e_glmnet) {
          error_message[i] <<- paste(
            "CV:", ifelse(is.na(cv_error_message[i]), "not run", cv_error_message[i]),
            "| fallback:", conditionMessage(e_glmnet)
          )
          NULL
        }
      )
    }

    if (is.null(fit_i)) {
      status[i] <- "fit_failed"
      message("    failed - ", error_message[i])
      flush.console()
      next
    }

    lambda_i <- if (inherits(fit_i, "cv.glmnet")) {
      fit_i$lambda.1se
    } else {
      fit_i$lambda[mid_lambda_index(fit_i$lambda)]
    }

    models[[i]] <- fit_i
    selected_lambda[i] <- lambda_i
    lambda_max_by_neuron[i] <- max(fit_i$lambda)
    lambda_min_by_neuron[i] <- min(fit_i$lambda)
    status[i] <- if (inherits(fit_i, "cv.glmnet")) "cv_glmnet" else "glmnet_fallback"

    # Coefficienti riportati sulla scala della design matrix completa.
    coef_i <- as.matrix(coef(fit_i, s = lambda_i))[, 1]
    full_coef <- numeric(ncol(X_full) + 1L)
    full_coef[1L] <- coef_i[1L]
    full_coef[valid_i + 1L] <- coef_i[-1L]
    full_coefs_by_neuron[[i]] <- full_coef

    delta_vec <- full_coef[(n_fixed + 2L):length(full_coef)]
    delta_array[, i, ] <- matrix(delta_vec, nrow = K, ncol = n_network_lags)

    # Dispersione sigma^2 sulla scala logaritmica.
    eta_i <- as.numeric(predict(fit_i, newx = X_i, s = lambda_i, type = "response"))
    df_i <- max(1L, sum(coef_i != 0) - 1L)
    denom <- max(1L, length(y_positive) - df_i)
    sigma2_i <- sum((y_fit - eta_i)^2) / denom
    dispersion[i] <- if (is.finite(sigma2_i) && sigma2_i >= 0) sigma2_i else NA_real_

    neuron_minutes <- as.numeric(difftime(Sys.time(), neuron_start, units = "mins"))
    message(
      sprintf(
        "    done - status: %s - lambda: %.3e - time: %.1f min",
        status[i], selected_lambda[i], neuron_minutes
      )
    )
    flush.console()
  }

  total_minutes <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
  message(sprintf("Step 2 fit completed in %.1f minutes.", total_minutes))

  # Diagnostica sulla scelta di lambda. Con percorso dipendente dai dati il
  # confronto va fatto contro il percorso di CIASCUN neurone, non contro una
  # griglia globale. Una quota alta al bordo superiore significa che la 1se
  # sceglie il modello nullo (solo covariate non penalizzate).
  cv_ids <- which(status == "cv_glmnet" & is.finite(selected_lambda))
  at_upper <- abs(selected_lambda[cv_ids] - lambda_max_by_neuron[cv_ids]) < 1e-12
  at_lower <- abs(selected_lambda[cv_ids] - lambda_min_by_neuron[cv_ids]) < 1e-12

  upper_boundary_rate <- if (length(cv_ids) > 0L) mean(at_upper) else NA_real_
  lower_boundary_rate <- if (length(cv_ids) > 0L) mean(at_lower) else NA_real_

  lambda_diagnostics <- list(
    selection_rule = "lambda.1se",
    lambda_source = if (data_driven_lambda) {
      "per-neuron path computed by glmnet from the data"
    } else {
      "fixed grid supplied by the caller"
    },
    n_lambda = if (data_driven_lambda) nlambda else length(lambda_grid),
    n_cv_fits = length(cv_ids),
    upper_lambda = if (data_driven_lambda) {
      range(lambda_max_by_neuron[cv_ids])
    } else max(lambda_grid),
    lower_lambda = if (data_driven_lambda) {
      range(lambda_min_by_neuron[cv_ids])
    } else min(lambda_grid),
    upper_boundary_rate = upper_boundary_rate,
    lower_boundary_rate = lower_boundary_rate
  )

  if (length(cv_ids) > 0L) {
    message(
      sprintf(
        "CV lambda diagnostics: %.1f%% at the upper boundary; %.1f%% at the lower boundary.",
        100 * upper_boundary_rate, 100 * lower_boundary_rate
      )
    )
  }

  delta_sum <- apply(delta_array, c(1, 2), sum)

  structure(
    list(
      models = models,
      family = "lognormal",
      selected_lambda = selected_lambda,
      dispersion = dispersion,
      valid_cols = valid_cols_by_neuron,
      full_coefs = full_coefs_by_neuron,
      foldid_by_neuron = foldid_by_neuron,
      stage1_foldid = stage1_foldid,
      cv_scheme = paste(
        "Exact Step 1 interleaved temporal-block folds, restricted to S_k(t)>0;",
        "all original folds required"
      ),
      delta_array = delta_array,
      delta_sum = delta_sum,
      positive_neurons = which(vapply(models, Negate(is.null), logical(1))),
      n_positive = n_positive,
      status = status,
      error_message = error_message,
      cv_error_message = cv_error_message,
      min_positive = min_positive,
      lambda_grid = lambda_grid,
      lambda_max_by_neuron = lambda_max_by_neuron,
      lambda_min_by_neuron = lambda_min_by_neuron,
      fixed_penalty = setNames(fixed_penalty, fixed_names),
      lambda_diagnostics = lambda_diagnostics,
      selection_rule = "lambda.1se; middle grid value for glmnet fallback",
      response_scale = "log-positive-amplitude",
      prep_signature = list(
        K = K,
        T_frames = T_frames,
        n_fixed_cols = n_fixed,
        fixed_covariate_names = fixed_names,
        max_lag_frames = max_lag,
        network_lag_values = lag_values,
        includes_lag0 = TRUE,
        lag0_self_predictor_excluded = TRUE,
        n_splines = stage1_n_splines(first_fit),
        dt_ms = prep$dt_ms,
        centroids = prep$centroids
      ),
      call = match.call()
    ),
    class = "hurdle_positive_fit"
  )
}


# ==============================================================================
# 3. PREVISIONE HURDLE (IN-SAMPLE E OUT-OF-SAMPLE)
# ==============================================================================

# Stage 1: P(S > 0 | X) su una design matrix qualsiasi.
predict_event_stage_from_X <- function(first_fit, X_new = first_fit$prep_data$X) {
  K <- first_fit$prep_data$K
  n <- nrow(X_new)
  pi_hat <- matrix(0, nrow = K, ncol = n)

  for (i in first_fit$active_neurons) {
    model_i <- first_fit$models[[i]]
    if (is.null(model_i)) next
    if (!glmnet_fit_is_usable(model_i)) {
      warning(sprintf("Neuron %d has a malformed Step 1 model: skipped.", i))
      next
    }

    valid_i <- stage1_valid_cols_for_neuron(first_fit, i)
    X_i <- X_new[, valid_i, drop = FALSE]
    dimnames(X_i) <- list(NULL, paste0("V", seq_len(ncol(X_i))))
    X_i <- Matrix(X_i, sparse = TRUE)

    lambda_i <- selected_lambda_stage1(model_i)
    pi_hat[i, ] <- as.numeric(
      predict(model_i, newx = X_i, s = lambda_i, type = "response")
    )
  }
  pi_hat
}

# Stage 2: E[S | S > 0, X] = exp(eta + sigma^2 / 2) su una design matrix qualsiasi.
predict_positive_stage_from_X <- function(second_fit, X_new) {
  K <- second_fit$prep_signature$K
  n <- nrow(X_new)
  m_hat <- matrix(0, nrow = K, ncol = n)

  for (i in second_fit$positive_neurons) {
    model_i <- second_fit$models[[i]]
    cols_i <- second_fit$valid_cols[[i]]
    lambda_i <- second_fit$selected_lambda[i]

    X_i <- X_new[, cols_i, drop = FALSE]
    dimnames(X_i) <- list(NULL, paste0("V", seq_len(ncol(X_i))))
    X_i <- Matrix(X_i, sparse = TRUE)

    eta <- as.numeric(predict(model_i, newx = X_i, s = lambda_i, type = "response"))
    sigma2 <- second_fit$dispersion[i]
    if (!is.finite(sigma2) || sigma2 < 0) sigma2 <- 0
    m_hat[i, ] <- exp(eta + 0.5 * sigma2)
  }
  m_hat
}

# Riusa le probabilita' Step 1 gia' salvate, se compatibili con S (solo in-sample).
extract_saved_stage1_probabilities <- function(saved_predictions, first_fit, S) {
  if (is.null(saved_predictions)) return(NULL)
  required <- c("predicted_probs", "actual_spikes")
  if (!all(required %in% names(saved_predictions))) return(NULL)

  probs <- as.matrix(saved_predictions$predicted_probs)
  actual <- as.matrix(saved_predictions$actual_spikes)
  expected_dim <- c(as.integer(first_fit$prep_data$K), as.integer(first_fit$prep_data$T_frames))

  if (!identical(dim(probs), expected_dim) || !identical(dim(actual), expected_dim)) {
    warning("Saved Step 1 predictions have incompatible dimensions and will be recomputed.")
    return(NULL)
  }

  observed_event <- ifelse(is.finite(S) & S > 0, 1, 0)
  if (any(actual != observed_event, na.rm = TRUE)) {
    warning("Saved actual_spikes does not match I(S>0); event probabilities will be recomputed.")
    return(NULL)
  }

  message("Reusing Step 1 event probabilities from the saved predictions.")
  probs
}

# ------------------------------------------------------------------------------
# Previsione grezza dei due stage, in-sample oppure out-of-sample.
#
#   trial_data = NULL  -> IN-SAMPLE: si usa la design matrix salvata nello Step 1.
#                         S e' la matrice di ampiezze dello stesso trial.
#   trial_data != NULL -> OUT-OF-SAMPLE: la design matrix viene ricostruita con
#                         prepare_design_matrix() (Step_1.R) sugli stessi
#                         max_lag_frames e n_splines dello Step 1. data_w_dist
#                         e' obbligatorio; S, se non fornita, viene presa da
#                         trial_data$neuron$S.
# ------------------------------------------------------------------------------
predict_hurdle_stages <- function(
    first_fit,
    second_fit,
    trial_data = NULL,
    data_w_dist = NULL,
    S = NULL,
    saved_stage1_predictions = NULL) {

  if (is.null(trial_data)) {
    # ---------------------------- IN-SAMPLE ----------------------------------
    if (is.null(S)) stop("In-sample prediction requires the amplitude matrix S.")
    S_obs <- check_S_matches_stage1(S, first_fit)
    X_new <- first_fit$prep_data$X

    pi_hat <- extract_saved_stage1_probabilities(saved_stage1_predictions, first_fit, S_obs)
    event_probability_source <- "saved_stage1_predictions"
    if (is.null(pi_hat)) {
      message("Step 1 predictions unavailable or incompatible; recomputing from the saved models.")
      pi_hat <- predict_event_stage_from_X(first_fit, X_new)
      event_probability_source <- "recomputed_from_saved_stage1_models"
    }
    evaluation_scope <- "in-sample matched trial"

  } else {
    # -------------------------- OUT-OF-SAMPLE --------------------------------
    if (is.null(data_w_dist)) {
      stop("Out-of-sample prediction also requires data_w_dist.")
    }
    if (!exists("prepare_design_matrix", mode = "function")) {
      stop("prepare_design_matrix() not found: source Step_1.R before predicting out-of-sample.")
    }

    # I nodi delle B-spline devono essere quelli della stima, altrimenti
    # spatial_basis_j denota una funzione dello spazio diversa da quella su cui
    # sono stati stimati i coefficienti.
    spline_basis <- first_fit$prep_data$spline_basis
    if (is.null(spline_basis)) {
      warning(
        "fit_t1$prep_data$spline_basis assente: le B-spline verranno rifittate ",
        "sul trial di test e la parte spaziale della previsione out-of-sample ",
        "sara' distorta. Rifitta lo Step 1 per correggere."
      )
    }

    prep_new <- prepare_design_matrix(
      trial_data = trial_data,
      data_w_dist = data_w_dist,
      max_lag_frames = first_fit$prep_data$max_lag_frames,
      n_splines = stage1_n_splines(first_fit),
      spline_basis = spline_basis
    )

    if (prep_new$K != first_fit$prep_data$K) {
      stop(
        "The test trial has ", prep_new$K, " neurons against the ",
        first_fit$prep_data$K, " of Step 1: the fitted models do not apply."
      )
    }
    if (ncol(prep_new$X) != ncol(first_fit$prep_data$X)) {
      stop("The test design matrix does not have the same number of columns as Step 1.")
    }

    X_new <- prep_new$X
    S_obs <- if (is.null(S)) as.matrix(trial_data$neuron$S) else as.matrix(S)
    pi_hat <- predict_event_stage_from_X(first_fit, X_new)
    event_probability_source <- "recomputed_from_saved_stage1_models"
    evaluation_scope <- "out-of-sample trial"
  }

  m_hat <- predict_positive_stage_from_X(second_fit, X_new)

  list(
    design_matrix = X_new,
    event_probability_raw = pi_hat,
    event_probability_source = event_probability_source,
    positive_mean = m_hat,
    observed_amplitude = S_obs,
    observed_event = ifelse(is.finite(S_obs) & S_obs > 0, 1, 0),
    evaluation_scope = evaluation_scope
  )
}


# ==============================================================================
# 4. COMBINAZIONE DEI DUE STAGE
# ==============================================================================
# Metodi disponibili:
#   "soft_raw"              : P(S>0|X) * E[S|S>0,X]   (attesa hurdle canonica)
#   "hard_threshold"        : I(P(S>0|X) >= soglia_i) * E[S|S>0,X]
#   "soft_logit_calibrated" : P_cal(S>0|X) * E[S|S>0,X]
# Soglie e calibrazione sono post-processing: non modificano ne' rifittano i
# modelli glmnet dei due stage.

clip_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

select_event_threshold <- function(y, p, rule = "f1") {
  rule <- match.arg(rule, choices = c("f1", "youden"))

  ok <- is.finite(y) & is.finite(p)
  y <- as.integer(y[ok] > 0)
  p <- as.numeric(p[ok])

  if (length(y) == 0L || sum(y) == 0L) return(Inf)
  if (sum(y) == length(y)) return(-Inf)

  ord <- order(p, decreasing = TRUE)
  p_ord <- p[ord]
  y_ord <- y[ord]

  tp <- cumsum(y_ord == 1L)
  fp <- cumsum(y_ord == 0L)
  n_pos <- sum(y_ord == 1L)
  n_neg <- sum(y_ord == 0L)
  fn <- n_pos - tp

  # Si valuta solo l'ultima occorrenza di ciascuna probabilita' distinta.
  candidate_idx <- which(c(diff(p_ord) != 0, TRUE))

  if (rule == "f1") {
    denominator <- 2 * tp + fp + fn
    score <- ifelse(denominator > 0, 2 * tp / denominator, 0)
  } else {
    sensitivity <- tp / n_pos
    specificity <- 1 - fp / n_neg
    score <- sensitivity + specificity - 1
  }

  candidate_score <- score[candidate_idx]
  best <- candidate_idx[which(candidate_score == max(candidate_score, na.rm = TRUE))]

  # In caso di parita' si sceglie la soglia con prevalenza predetta piu' vicina
  # a quella osservata; se la parita' resta, la soglia piu' alta.
  observed_rate <- mean(y)
  predicted_rate <- best / length(y)
  best <- best[abs(predicted_rate - observed_rate) ==
                 min(abs(predicted_rate - observed_rate))]
  best <- min(best)

  p_ord[best]
}

fit_logit_probability_calibrator <- function(
    y,
    p,
    slope_bounds = c(0.25, 8),
    eps = 1e-6) {

  if (length(slope_bounds) != 2L ||
      !all(is.finite(slope_bounds)) ||
      slope_bounds[1] <= 0 ||
      slope_bounds[1] > slope_bounds[2]) {
    stop("calibration_slope_bounds must contain two positive ordered values.")
  }

  ok <- is.finite(y) & is.finite(p)
  y <- as.integer(y[ok] > 0)
  p <- clip_probability(p[ok], eps)

  if (length(y) == 0L || length(unique(y)) < 2L) {
    return(list(intercept = 0, slope = 1, status = "identity_no_two_classes"))
  }

  x <- qlogis(p)
  calibration_fit <- try(
    suppressWarnings(
      glm.fit(
        x = cbind(`(Intercept)` = 1, raw_logit = x),
        y = y,
        family = binomial()
      )
    ),
    silent = TRUE
  )

  slope <- 1
  fit_status <- "identity_fallback"
  if (!inherits(calibration_fit, "try-error")) {
    coef_fit <- calibration_fit$coefficients
    if (length(coef_fit) >= 2L && is.finite(coef_fit[2L])) {
      slope <- min(max(unname(coef_fit[2L]), slope_bounds[1]), slope_bounds[2])
      fit_status <- if (abs(slope - unname(coef_fit[2L])) < 1e-12) {
        "estimated"
      } else {
        "estimated_slope_clipped"
      }
    }
  }

  # L'intercetta viene ristimata dopo aver vincolato la pendenza, cosi' che
  # mean(P_cal) coincida con il tasso di eventi osservato.
  observed_rate <- mean(y)
  intercept_equation <- function(intercept) {
    mean(plogis(intercept + slope * x)) - observed_rate
  }

  intercept <- tryCatch(
    uniroot(intercept_equation, interval = c(-50, 50), tol = 1e-10)$root,
    error = function(e) NA_real_
  )

  if (!is.finite(intercept)) {
    intercept <- 0
    slope <- 1
    fit_status <- "identity_intercept_failure"
  }

  list(intercept = intercept, slope = slope, status = fit_status)
}

classification_diagnostics <- function(y, predicted_event) {
  y <- as.integer(y > 0)
  predicted_event <- as.integer(predicted_event > 0)

  tp <- sum(predicted_event == 1L & y == 1L)
  fp <- sum(predicted_event == 1L & y == 0L)
  fn <- sum(predicted_event == 0L & y == 1L)
  tn <- sum(predicted_event == 0L & y == 0L)

  precision <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
  sensitivity <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(sensitivity) &&
            precision + sensitivity > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else {
    NA_real_
  }

  c(
    precision = precision,
    sensitivity = sensitivity,
    specificity = specificity,
    f1 = f1
  )
}

# ------------------------------------------------------------------------------
# Costruisce le tre combinazioni dei due stage.
#
# fitted_components = NULL  -> soglie e calibrazione stimate sui dati passati
#                              (uso in-sample).
# fitted_components != NULL -> soglie e calibrazione dell'oggetto fornito
#                              vengono APPLICATE senza ristimarle: e' la forma
#                              corretta per la previsione out-of-sample.
# ------------------------------------------------------------------------------
build_joint_prediction_components <- function(
    base_predictions,
    neuron_ids,
    threshold_rule = "f1",
    calibration_slope_bounds = c(0.25, 8),
    eps = 1e-6,
    fitted_components = NULL) {

  threshold_rule <- match.arg(threshold_rule, choices = c("f1", "youden"))

  pi_raw <- as.matrix(base_predictions$event_probability_raw)
  m_hat <- as.matrix(base_predictions$positive_mean)
  Z <- as.matrix(base_predictions$observed_event)

  if (!identical(dim(pi_raw), dim(m_hat)) || !identical(dim(pi_raw), dim(Z))) {
    stop("The event, amplitude and observed-event matrices have incompatible dimensions.")
  }

  K <- nrow(pi_raw)
  T_frames <- ncol(pi_raw)
  neuron_ids <- sort(unique(as.integer(neuron_ids)))
  neuron_ids <- neuron_ids[neuron_ids >= 1L & neuron_ids <= K]

  reuse <- !is.null(fitted_components)
  if (reuse) {
    threshold_rule <- fitted_components$threshold_rule
    calibration_slope_bounds <- fitted_components$calibration_slope_bounds
    eps <- fitted_components$probability_epsilon
  }

  hard_event <- matrix(0L, nrow = K, ncol = T_frames)
  pi_calibrated <- pi_raw
  thresholds <- rep(NA_real_, K)
  calibration_intercept <- rep(NA_real_, K)
  calibration_slope <- rep(NA_real_, K)
  calibration_status <- rep("not_fitted", K)
  diagnostic_rows <- vector("list", length(neuron_ids))

  for (idx in seq_along(neuron_ids)) {
    i <- neuron_ids[idx]
    y_i <- Z[i, ]
    p_i <- pi_raw[i, ]

    if (reuse) {
      thresholds[i] <- fitted_components$threshold_by_neuron[i]
      calibration_intercept[i] <- fitted_components$calibration_intercept_by_neuron[i]
      calibration_slope[i] <- fitted_components$calibration_slope_by_neuron[i]
      calibration_status[i] <- paste0(
        "applied_from_fit:", fitted_components$calibration_status_by_neuron[i]
      )
      if (!is.finite(calibration_intercept[i]) || !is.finite(calibration_slope[i])) {
        calibration_intercept[i] <- 0
        calibration_slope[i] <- 1
        calibration_status[i] <- "identity_no_fitted_calibration"
      }
    } else {
      thresholds[i] <- select_event_threshold(y = y_i, p = p_i, rule = threshold_rule)

      calibrator_i <- fit_logit_probability_calibrator(
        y = y_i,
        p = p_i,
        slope_bounds = calibration_slope_bounds,
        eps = eps
      )
      calibration_intercept[i] <- calibrator_i$intercept
      calibration_slope[i] <- calibrator_i$slope
      calibration_status[i] <- calibrator_i$status
    }

    hard_event[i, ] <- as.integer(p_i >= thresholds[i])

    p_i_clipped <- clip_probability(p_i, eps)
    pi_calibrated[i, ] <- plogis(
      calibration_intercept[i] + calibration_slope[i] * qlogis(p_i_clipped)
    )

    class_diag <- classification_diagnostics(y_i, hard_event[i, ])

    diagnostic_rows[[idx]] <- data.frame(
      neuron = i,
      observed_event_rate = mean(y_i),
      mean_probability_raw = mean(p_i),
      maximum_probability_raw = max(p_i),
      hard_threshold = thresholds[i],
      hard_predicted_event_rate = mean(hard_event[i, ]),
      hard_precision = unname(class_diag["precision"]),
      hard_sensitivity = unname(class_diag["sensitivity"]),
      hard_specificity = unname(class_diag["specificity"]),
      hard_f1 = unname(class_diag["f1"]),
      calibration_intercept = calibration_intercept[i],
      calibration_slope = calibration_slope[i],
      calibration_status = calibration_status[i],
      mean_probability_calibrated = mean(pi_calibrated[i, ]),
      maximum_probability_calibrated = max(pi_calibrated[i, ]),
      brier_raw = mean((y_i - p_i)^2),
      brier_calibrated = mean((y_i - pi_calibrated[i, ])^2),
      stringsAsFactors = FALSE
    )
  }

  list(
    event_probability_raw = pi_raw,
    event_probability_calibrated = pi_calibrated,
    hard_event = hard_event,
    positive_mean = m_hat,
    observed_amplitude = base_predictions$observed_amplitude,
    observed_event = Z,
    diagnostics = do.call(rbind, diagnostic_rows),
    threshold_by_neuron = thresholds,
    calibration_intercept_by_neuron = calibration_intercept,
    calibration_slope_by_neuron = calibration_slope,
    calibration_status_by_neuron = calibration_status,
    threshold_rule = threshold_rule,
    calibration_slope_bounds = calibration_slope_bounds,
    probability_epsilon = eps,
    available_methods = c(
      "soft_raw",
      "hard_threshold",
      "soft_logit_calibrated"
    ),
    evaluation_scope = base_predictions$evaluation_scope,
    calibration_scope = if (reuse) {
      "thresholds and calibration applied from the in-sample fit"
    } else {
      "in-sample exploratory post-processing"
    }
  )
}

select_joint_prediction <- function(
    components,
    method = c("soft_raw", "hard_threshold", "soft_logit_calibrated")) {

  method <- match.arg(method, choices = components$available_methods)

  event_weight <- switch(
    method,
    soft_raw = components$event_probability_raw,
    hard_threshold = components$hard_event,
    soft_logit_calibrated = components$event_probability_calibrated
  )

  prediction_label <- switch(
    method,
    soft_raw = "Raw probability x positive amplitude",
    hard_threshold = "Thresholded event x positive amplitude",
    soft_logit_calibrated = "Logit-calibrated probability x positive amplitude"
  )

  list(
    prediction_method = method,
    prediction_label = prediction_label,
    event_probability = event_weight,
    event_weight = event_weight,
    event_probability_raw = components$event_probability_raw,
    event_probability_calibrated = components$event_probability_calibrated,
    hard_event = components$hard_event,
    positive_mean = components$positive_mean,
    expected_amplitude = event_weight * components$positive_mean,
    observed_amplitude = components$observed_amplitude,
    observed_event = components$observed_event,
    evaluation_scope = components$evaluation_scope,
    postprocessing_scope = components$calibration_scope
  )
}

# ------------------------------------------------------------------------------
# Interfaccia unica di previsione: dai due fit alla previsione finale.
# Restituisce componenti, tutte le combinazioni disponibili e quella scelta.
# Per l'out-of-sample passare trial_data + data_w_dist e, per non ristimare
# soglie/calibrazione sui dati di test, fitted_components dell'in-sample.
# ------------------------------------------------------------------------------
predict_hurdle <- function(
    first_fit,
    second_fit,
    trial_data = NULL,
    data_w_dist = NULL,
    S = NULL,
    saved_stage1_predictions = NULL,
    method = "soft_raw",
    threshold_rule = "f1",
    calibration_slope_bounds = c(0.25, 8),
    eps = 1e-6,
    fitted_components = NULL,
    neuron_ids = NULL) {

  base_predictions <- predict_hurdle_stages(
    first_fit = first_fit,
    second_fit = second_fit,
    trial_data = trial_data,
    data_w_dist = data_w_dist,
    S = S,
    saved_stage1_predictions = saved_stage1_predictions
  )

  if (is.null(neuron_ids)) neuron_ids <- second_fit$positive_neurons

  components <- build_joint_prediction_components(
    base_predictions = base_predictions,
    neuron_ids = neuron_ids,
    threshold_rule = threshold_rule,
    calibration_slope_bounds = calibration_slope_bounds,
    eps = eps,
    fitted_components = fitted_components
  )

  if (!method %in% components$available_methods) {
    stop(
      "Unknown method: ", method, ". Available methods: ",
      paste(components$available_methods, collapse = ", ")
    )
  }

  options <- setNames(
    lapply(components$available_methods, function(m) select_joint_prediction(components, m)),
    components$available_methods
  )

  list(
    base = base_predictions,
    components = components,
    options = options,
    selected_method = method,
    selected = options[[method]],
    neuron_ids = neuron_ids
  )
}


# ==============================================================================
# 5. VALUTAZIONE NUMERICA DELLE PREVISIONI
# ==============================================================================
evaluate_hurdle_predictions <- function(predictions, neuron_ids = NULL) {
  S <- predictions$observed_amplitude
  Z <- predictions$observed_event
  pi_hat <- predictions$event_probability
  m_hat <- predictions$positive_mean
  expected <- predictions$expected_amplitude

  if (is.null(neuron_ids)) neuron_ids <- seq_len(nrow(S))

  rows <- lapply(neuron_ids, function(i) {
    ok <- is.finite(S[i, ]) & is.finite(expected[i, ])
    pos <- ok & S[i, ] > 0 & m_hat[i, ] > 0

    data.frame(
      neuron = i,
      n_frames = sum(ok),
      n_positive = sum(pos),
      event_brier = if (any(ok)) mean((Z[i, ok] - pi_hat[i, ok])^2) else NA_real_,
      mae_all = if (any(ok)) mean(abs(S[i, ok] - expected[i, ok])) else NA_real_,
      rmse_all = if (any(ok)) sqrt(mean((S[i, ok] - expected[i, ok])^2)) else NA_real_,
      mae_positive = if (any(pos)) mean(abs(S[i, pos] - m_hat[i, pos])) else NA_real_,
      rmse_positive = if (any(pos)) sqrt(mean((S[i, pos] - m_hat[i, pos])^2)) else NA_real_,
      rmsle_positive = if (any(pos)) {
        sqrt(mean((log(S[i, pos]) - log(m_hat[i, pos]))^2))
      } else NA_real_
    )
  })

  per_neuron <- do.call(rbind, rows)
  aggregate <- data.frame(
    metric = c(
      "event_brier", "mae_all", "rmse_all",
      "mae_positive", "rmse_positive", "rmsle_positive"
    ),
    mean_across_neurons = c(
      mean(per_neuron$event_brier, na.rm = TRUE),
      mean(per_neuron$mae_all, na.rm = TRUE),
      mean(per_neuron$rmse_all, na.rm = TRUE),
      mean(per_neuron$mae_positive, na.rm = TRUE),
      mean(per_neuron$rmse_positive, na.rm = TRUE),
      mean(per_neuron$rmsle_positive, na.rm = TRUE)
    )
  )
  list(per_neuron = per_neuron, aggregate = aggregate)
}

# Confronto sintetico fra i tre metodi di combinazione.
compare_joint_prediction_methods <- function(prediction_result, neuron_ids = NULL) {
  if (is.null(neuron_ids)) neuron_ids <- prediction_result$neuron_ids

  evaluations <- lapply(
    prediction_result$options,
    evaluate_hurdle_predictions,
    neuron_ids = neuron_ids
  )

  comparison <- do.call(
    rbind,
    lapply(names(evaluations), function(m) {
      out <- evaluations[[m]]$aggregate
      out$prediction_method <- m
      out[, c("prediction_method", "metric", "mean_across_neurons")]
    })
  )
  row.names(comparison) <- NULL

  list(evaluations = evaluations, comparison = comparison)
}

# ==============================================================================
# HURDLE MODEL - STEP 1: EVENT STAGE   P(S_k(t) > 0 | X)
# ==============================================================================
# Modello logistico penalizzato (LASSO con penalty.factor basato sulla distanza
# fra centroidi) stimato separatamente per ogni neurone.
#
# Design matrix: [covariate fisse (B-spline spaziali + velocita') | lag 0 | lag 1..L]
# Il neurone i esclude la propria colonna al lag 0 (esclusione autoreferenziale).
#
# Questo file contiene SOLO definizioni di funzioni: caricamento dati ed
# esecuzione stanno in main.R.
#
# Contenuto:
#   0. capwords()                       helper di formattazione
#   1. prepare_design_matrix()          costruzione design matrix
#   2. fit_model()                      stima dello stage evento
#   3. predict_model()                  previsione in-sample / out-of-sample
#   4. evaluate_predictions()           ROC / AUC
#   5. helper di estrazione coefficienti
#   6-13. funzioni diagnostiche grafiche
# ==============================================================================

library(dplyr)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(Matrix)
library(pROC)
library(splines)
library(tidyr)


# ==============================================================================
# 0. HELPER GENERALI
# ==============================================================================
# Indice "centrale" di un percorso di lambda, usato quando la CV non e'
# disponibile e si ripiega su glmnet semplice.
#
# ATTENZIONE: round(length(lambda) / 2) vale 0 quando il percorso ha un solo
# valore (R arrotonda 0.5 a 0), e lambda[0] restituisce numeric(0). glmnet
# accorcia il percorso quando il modello satura, quindi il caso si presenta
# davvero e faceva fallire predict() con un errore su dimnames.
mid_lambda_index <- function(lambda) {
  n <- length(lambda)
  if (n < 1L) stop("Percorso di lambda vuoto.")
  max(1L, min(n, round(n / 2)))
}

# Verifica che un oggetto glmnet sia internamente coerente.
#
# Su dati degeneri (risposta separabile, devianza nulla) glmnet puo' restituire
# un oggetto malformato - tipicamente lambda = Inf con length(lambda) = 1 ma
# length(a0) = 20 - su cui coef() e predict() falliscono. Va scartato invece che
# conservato, altrimenti fa crashare la previsione molto piu' a valle.
glmnet_fit_is_usable <- function(fit) {
  if (is.null(fit)) return(FALSE)
  if (inherits(fit, "cv.glmnet")) fit <- fit$glmnet.fit
  n_lambda <- length(fit$lambda)

  n_lambda >= 1L &&
    all(is.finite(fit$lambda)) &&
    !is.null(fit$beta) &&
    ncol(fit$beta) == n_lambda &&
    length(fit$a0) == n_lambda
}

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if (strict) tolower(s) else s},
                           sep = "", collapse = " ")
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


# ==============================================================================
# 1. PREPARAZIONE DESIGN MATRIX (Covariate Fisse + Lag 0 + Lags 1..L)
# ==============================================================================
# spline_basis: se NULL le B-spline vengono fittate sulle posizioni di QUESTO
# trial e gli oggetti bs risultanti vengono restituiti in prep_data$spline_basis.
# Se fornito, i nodi vengono ereditati via predict(): indispensabile per la
# previsione out-of-sample, altrimenti spatial_basis_j denoterebbe una funzione
# dello spazio diversa da quella su cui sono stati stimati i coefficienti.
prepare_design_matrix <- function(trial_data, data_w_dist, max_lag_frames = 5,
                                  n_splines = 4, spline_basis = NULL) {

  # --- 1. Estrarre e Allineare i Tempi Neurali ---
  S <- trial_data$neuron$S                          # [K x T_frames]
  centroids <- trial_data$neuron$centroid           # [K x 2]
  t_raw <- trial_data$neuron$time[, 1]

  K <- nrow(S)
  T_frames <- ncol(S)

  t_neuron <- t_raw[seq(1, length(t_raw), length.out = T_frames)]
  dt_ms <- mean(diff(t_neuron))

  # Matrice binaria degli eventi Z
  Z <- ifelse(S > 0, 1, 0)                          # [K x T_frames]

  # --- 2. Interpolazione Posizione Arena e Velocita' sui Tempi Neurali ---
  t_behav <- trial_data$behav$time[, 1]
  pos_behav <- trial_data$behav$position

  pos_x <- approx(x = t_behav, y = pos_behav[, 1], xout = t_neuron, rule = 2)$y
  pos_y <- approx(x = t_behav, y = pos_behav[, 2], xout = t_neuron, rule = 2)$y

  speed_interp <- approx(x = data_w_dist$time, y = data_w_dist$speed_spline, xout = t_neuron, rule = 2)$y

  # --- 3. Basi Spaziali B-Spline 2D + Covariata Velocita' ---
  if (is.null(spline_basis)) {
    bx <- bs(pos_x, df = n_splines)
    by <- bs(pos_y, df = n_splines)
  } else {
    # predict.bs() riusa knots e Boundary.knots dell'oggetto originale
    bx <- predict(spline_basis$bx, pos_x)
    by <- predict(spline_basis$by, pos_y)
  }

  X_spatial_list <- lapply(1:ncol(bx), function(i) bx[, i] * by)
  X_spatial <- do.call(cbind, X_spatial_list)
  colnames(X_spatial) <- paste0("spatial_basis_", 1:ncol(X_spatial))

  X_fixed <- cbind(X_spatial, speed_spline = speed_interp)
  n_fixed_cols <- ncol(X_fixed)

  # --- 4. Contemporaneous Network Activity (LAG 0) ---
  X_lag0 <- t(Z)                                    # [T_frames x K]
  colnames(X_lag0) <- paste0("N", 1:K, "_lag0")

  # --- 5. Autoregressive Lags (LAG 1..L) ---
  X_lags <- list()
  for (lag in 1:max_lag_frames) {
    Z_lagged <- cbind(matrix(0, nrow = K, ncol = lag), Z[, 1:(T_frames - lag)])
    rownames(Z_lagged) <- paste0("N", 1:K, "_lag", lag)
    X_lags[[lag]] <- Z_lagged
  }
  X_lags <- do.call(rbind, X_lags) %>% t()          # [T_frames x (K * max_lag_frames)]

  # Design Matrix Finale: [X_fixed | X_lag0 | X_lags]
  X_full <- cbind(X_fixed, X_lag0, X_lags)

  return(list(
    X = X_full,
    Z = Z,
    K = K,
    T_frames = T_frames,
    dt_ms = dt_ms,
    centroids = centroids,
    n_fixed_cols = n_fixed_cols,
    max_lag_frames = max_lag_frames,
    n_splines = n_splines,
    # Portatori dei nodi, da riusare in previsione out-of-sample
    spline_basis = list(bx = bx, by = by)
  ))
}


# ==============================================================================
# 2. FIT DEL MODELLO (Esclusione Autoreferenziale al Lag 0)
# ==============================================================================
fit_model <- function(prep_data, family = "binomial", min_spikes = 15, nfolds = 5, block_dur_sec = 2.0) {

  K <- prep_data$K
  X_full <- prep_data$X
  Z <- prep_data$Z
  n_fixed <- prep_data$n_fixed_cols
  max_lag_frames <- prep_data$max_lag_frames
  T_frames <- prep_data$T_frames
  dt_ms <- prep_data$dt_ms

  dist_matrix <- as.matrix(dist(prep_data$centroids))
  max_dist <- max(dist_matrix)

  active_neurons <- which(rowSums(Z) >= min_spikes)

  models <- vector("list", K)
  valid_cols_list <- vector("list", K)
  W_matrix <- matrix(0, nrow = K, ncol = K)

  lambda_grid <- exp(seq(log(1), log(0.0001), length.out = 20))

  # Folds Temporali a Blocchi
  block_size_frames <- max(1, round((block_dur_sec * 1000) / dt_ms))
  n_blocks <- ceiling(T_frames / block_size_frames)
  block_ids <- rep(1:n_blocks, each = block_size_frames, length.out = T_frames)
  foldid <- ((block_ids - 1) %% nfolds) + 1

  col_counts <- colSums(X_full > 0)

  for (i in active_neurons) {
    y <- Z[i, ]

    # 1. Penalita' covariate fisse (0)
    pen_fixed <- rep(0, n_fixed)

    # 2. Penalita' Lag 0 (pesata per distanza dal neurone i)
    raw_pen_lag0 <- dist_matrix[i, ] / max_dist
    pen_lag0 <- pmax(raw_pen_lag0, 0.05)

    # 3. Penalita' Lags 1..L (pesata per distanza)
    raw_pen_lags <- rep(dist_matrix[i, ] / max_dist, times = max_lag_frames)
    pen_lags <- pmax(raw_pen_lags, 0.05)

    penalty_full <- c(pen_fixed, pen_lag0, pen_lags)

    # Indice della colonna Z_i(t) al Lag 0 (da rimuovere per il neurone i)
    col_self_lag0 <- n_fixed + i

    # Seleziona colonne con varianza e rimuovi il neurone i al lag 0
    valid_cols_i <- c(1:n_fixed, which(col_counts >= 2 & (1:ncol(X_full)) > n_fixed))
    valid_cols_i <- setdiff(valid_cols_i, col_self_lag0)

    valid_cols_list[[i]] <- valid_cols_i

    X_sub <- X_full[, valid_cols_i, drop = FALSE]
    penalty_sub <- penalty_full[valid_cols_i]

    dimnames(X_sub) <- list(NULL, paste0("V", 1:ncol(X_sub)))
    X_sparse <- Matrix(X_sub, sparse = TRUE)

    fit_done <- FALSE
    coefs_sub <- NULL

    tryCatch({
      spikes_per_fold <- tapply(y, foldid, sum)
      if (any(spikes_per_fold == 0)) {
        stop("At least one temporal fold contains 0 spikes for this neuron.")
      }

      fit <- cv.glmnet(
        x = X_sparse, y = y, family = family,
        penalty.factor = penalty_sub, foldid = foldid,
        lambda = lambda_grid, type.measure = "deviance",
        alpha = 1, maxit = 100000, parallel = TRUE
      )
      models[[i]] <- fit
      coefs_sub <- as.matrix(coef(fit, s = "lambda.1se"))[-1, 1]
      fit_done <- TRUE
    }, error = function(e) {
      tryCatch({
        fit_fallback <- glmnet(
          x = X_sparse, y = y, family = family,
          penalty.factor = penalty_sub, lambda = lambda_grid,
          alpha = 1, maxit = 100000
        )
        if (!glmnet_fit_is_usable(fit_fallback)) {
          stop("degenerate glmnet fallback (malformed object)")
        }
        # Indice sul percorso EFFETTIVO del modello, non su lambda_grid: glmnet
        # puo' restituirne meno di quanti richiesti.
        mid_idx <- mid_lambda_index(fit_fallback$lambda)
        cs <- as.matrix(coef(fit_fallback))[-1, mid_idx]

        # Si conserva il modello SOLO dopo che i coefficienti sono stati estratti
        # senza errori: un oggetto degenere in models[[i]] farebbe crashare
        # predict() molto piu' a valle. <<- e' necessario perche' qui siamo
        # dentro un gestore d'errore.
        models[[i]] <<- fit_fallback
        coefs_sub <<- cs
        fit_done <<- TRUE
      }, error = function(e2) {
        message(sprintf("Neuron %d CV avoided for numeric instability.", i))
      })
    })

    # Ricostruzione della Matrice di Adiacenza Totale (Lag 0 + Lags 1..L)
    if (fit_done && !is.null(coefs_sub)) {
      full_coefs <- numeric(ncol(X_full))
      full_coefs[valid_cols_i] <- coefs_sub

      coef_lags <- full_coefs[(n_fixed + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = 1 + max_lag_frames)

      W_matrix[, i] <- rowSums(coef_mat)
    }
  }

  return(list(
    models = models,
    adj_matrix = W_matrix,
    prep_data = prep_data,
    active_neurons = active_neurons,
    valid_cols_list = valid_cols_list,
    foldid = foldid
  ))
}


# ==============================================================================
# 3. PREDIZIONE IN-SAMPLE / OUT-OF-SAMPLE
# ==============================================================================
# Per la previsione out-of-sample, data_w_dist deve essere la serie di velocita'
# del trial di TEST, non quella usata in stima.
predict_model <- function(fit_results, test_trial_data, data_w_dist) {

  max_lag <- fit_results$prep_data$max_lag_frames
  n_splines <- fit_results$prep_data$n_splines
  if (is.null(n_splines)) {
    n_splines <- round(sqrt(fit_results$prep_data$n_fixed_cols - 1))
  }

  # Nodi delle B-spline ereditati dalla stima. Assente nei fit salvati prima
  # dell'introduzione del campo: in quel caso le basi vengono rifittate sul
  # trial di test e la parte spaziale risulta distorta.
  spline_basis <- fit_results$prep_data$spline_basis
  if (is.null(spline_basis) && !identical(test_trial_data, NULL)) {
    warning(
      "fit_results$prep_data$spline_basis assente: le B-spline verranno ",
      "rifittate sul trial di test e i coefficienti spaziali non saranno ",
      "confrontabili. Rifitta lo Step 1 per correggere."
    )
  }

  prep_test <- prepare_design_matrix(
    trial_data = test_trial_data,
    data_w_dist = data_w_dist,
    max_lag_frames = max_lag,
    n_splines = n_splines,
    spline_basis = spline_basis
  )

  K <- prep_test$K
  T_test <- prep_test$T_frames
  prob_matrix <- matrix(0, nrow = K, ncol = T_test)

  for (i in fit_results$active_neurons) {
    fit_i <- fit_results$models[[i]]
    valid_cols_i <- fit_results$valid_cols_list[[i]]

    if (!is.null(fit_i) && !glmnet_fit_is_usable(fit_i)) {
      warning(sprintf("Neuron %d has a malformed Step 1 model: skipped.", i))
      next
    }

    if (!is.null(fit_i) && !is.null(valid_cols_i)) {
      X_test_sub <- prep_test$X[, valid_cols_i, drop = FALSE]
      dimnames(X_test_sub) <- list(NULL, paste0("V", 1:ncol(X_test_sub)))
      X_test_sparse <- Matrix(X_test_sub, sparse = TRUE)

      if (inherits(fit_i, "cv.glmnet")) {
        probs <- predict(fit_i, newx = X_test_sparse, s = "lambda.1se", type = "response")
      } else {
        mid_lambda <- fit_i$lambda[mid_lambda_index(fit_i$lambda)]
        probs <- predict(fit_i, newx = X_test_sparse, s = mid_lambda, type = "response")
      }
      prob_matrix[i, ] <- as.vector(probs)
    }
  }

  return(list(
    predicted_probs = prob_matrix,
    actual_spikes   = prep_test$Z,
    time_ms         = prep_test$dt_ms
  ))
}


# ==============================================================================
# 4. VALUTAZIONE PRESTAZIONI (ROC & AUC)
# ==============================================================================
evaluate_predictions <- function(predictions, active_neurons) {

  auc_list <- numeric(length(active_neurons))
  roc_data_list <- list()
  grid_spec <- seq(0, 1, length.out = 200)
  sens_matrix <- matrix(NA, nrow = length(active_neurons), ncol = length(grid_spec))

  for (k in seq_along(active_neurons)) {
    i <- active_neurons[k]
    y_true <- predictions$actual_spikes[i, ]
    y_pred <- predictions$predicted_probs[i, ]

    if (length(unique(y_true)) > 1) {
      roc_i <- roc(y_true, y_pred, quiet = TRUE)
      auc_list[k] <- auc(roc_i)

      ord <- order(roc_i$specificities)
      sens_interp <- approx(
        x = roc_i$specificities[ord],
        y = roc_i$sensitivities[ord],
        xout = grid_spec,
        rule = 2
      )$y
      sens_matrix[k, ] <- sens_interp

      roc_data_list[[k]] <- data.frame(
        Specificity = roc_i$specificities,
        Sensitivity = roc_i$sensitivities,
        Neuron = factor(i)
      )
    }
  }

  df_roc_all <- bind_rows(roc_data_list)
  df_roc_mean <- data.frame(
    Specificity = grid_spec,
    Sensitivity = colMeans(sens_matrix, na.rm = TRUE)
  )

  p_roc <- ggplot() +
    geom_line(data = df_roc_all, aes(x = 1 - Specificity, y = Sensitivity, group = Neuron),
              color = "gray70", alpha = 0.25) +
    geom_line(data = df_roc_mean, aes(x = 1 - Specificity, y = Sensitivity),
              color = "darkred", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = sprintf("ROC Curves (Mean AUC = %.3f)", mean(auc_list, na.rm = TRUE)),
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme_minimal()

  return(list(
    auc_vector = auc_list,
    mean_auc = mean(auc_list, na.rm = TRUE),
    p_roc = p_roc
  ))
}


# ==============================================================================
# 5. HELPER: ESTRAZIONE DEI COEFFICIENTI DI UN NEURONE
# ==============================================================================
# Coefficienti stimati (senza intercetta) riportati sulla scala delle colonne
# della design matrix completa. Restituisce NULL se il neurone non e' stimato.
step1_full_coefs <- function(fit_results, i) {
  fit_i <- fit_results$models[[i]]
  valid_cols_i <- fit_results$valid_cols_list[[i]]
  if (is.null(fit_i) || is.null(valid_cols_i)) return(NULL)

  if (inherits(fit_i, "cv.glmnet")) {
    coefs_sub <- as.matrix(coef(fit_i, s = "lambda.1se"))[-1, 1]
  } else {
    mid_idx <- mid_lambda_index(fit_i$lambda)
    coefs_sub <- as.matrix(coef(fit_i))[-1, mid_idx]
  }

  full_coefs <- numeric(ncol(fit_results$prep_data$X))
  full_coefs[valid_cols_i] <- coefs_sub
  full_coefs
}

# Matrice [K x (1 + max_lag)] dei coefficienti di rete del neurone i.
# La colonna 1 corrisponde al lag 0, la colonna l+1 al lag l.
step1_lag_coef_matrix <- function(fit_results, i) {
  full_coefs <- step1_full_coefs(fit_results, i)
  if (is.null(full_coefs)) return(NULL)

  n_fixed <- fit_results$prep_data$n_fixed_cols
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames

  coef_lags <- full_coefs[(n_fixed + 1):length(full_coefs)]
  matrix(coef_lags, nrow = K, ncol = 1 + max_lag)
}


# ==============================================================================
# 6. DIAGNOSTICA: Heatmap dei Coefficienti Autoregressivi (Lag Singolo)
# ==============================================================================
plot_single_lag <- function(fit_results, lag = 0, only_active = TRUE) {

  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  active_neurons <- fit_results$active_neurons
  dt_ms <- fit_results$prep_data$dt_ms

  # 1. Validation del Lag (accetta da 0 a max_lag)
  if (lag < 0 || lag > max_lag) {
    stop(sprintf("Lag %d non valido! Scegli un valore compreso tra 0 e %d.", lag, max_lag))
  }

  idx_plot <- if (only_active) active_neurons else 1:K
  N_plot <- length(idx_plot)

  W_mat_full <- matrix(0, nrow = K, ncol = K)

  # 2. Estrazione dei coefficienti gamma per ciascun neurone target i
  for (i in active_neurons) {
    coef_mat <- step1_lag_coef_matrix(fit_results, i)
    if (!is.null(coef_mat)) {
      # Seleziona la colonna corrispondente al Lag scelto (lag 0 -> colonna 1)
      W_mat_full[, i] <- coef_mat[, lag + 1]
    }
  }

  # Sotto-matrice per i soli neuroni da visualizzare
  W_sub <- W_mat_full[idx_plot, idx_plot, drop = FALSE]

  # Dataframe per ggplot (struttura a griglia)
  df_single_lag <- data.frame(
    Source = factor(rep(idx_plot, times = N_plot), levels = rev(idx_plot)), # Input (y-axis)
    Target = factor(rep(idx_plot, each = N_plot), levels = idx_plot),       # Target (x-axis)
    Weight = as.vector(W_sub)
  )

  # Titolo dinamico in base al Lag
  title_str <- if (lag == 0) {
    "Adjcency Matrix - Lag 0 (t = 0 ms)"
  } else {
    sprintf("Adjcency Matrix - Lag %d (t - %d ms)", lag, round(lag * dt_ms))
  }

  # 3. Rendering Grafico (Heatmap)
  p <- ggplot(df_single_lag, aes(x = Target, y = Source, fill = Weight)) +
    geom_tile(color = "white", linewidth = 0.1) +

    # Gradiente Divergente: Blu (Inibizione), Bianco (Zero), Rosso (Eccitazione)
    scale_fill_gradient2(
      low = "#2980B9",
      mid = "#FFFFFF",
      high = "#C0392B",
      midpoint = 0,
      name = "Coeff. gamma"
    ) +

    labs(
      x = "Target neuron (i)",
      y = "Input neuron (j)",
      title = title_str
    ) +

    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )

  return(p)
}


# ==============================================================================
# 7. DIAGNOSTICA: Mappa Spaziale Multi-Lag dei Coefficienti Autoregressivi
# ==============================================================================
plot_spatial_multi_lag <- function(fit_results, lags, target_id = NULL) {

  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  active_neurons <- fit_results$active_neurons
  dt_ms <- fit_results$prep_data$dt_ms
  centroids <- fit_results$prep_data$centroids

  # Validation dei lag (accetta anche 0)
  if (any(lags < 0 | lags > max_lag)) {
    stop(sprintf("Uno o piu' lag non sono validi! Scegli valori compresi tra 0 e %d.", max_lag))
  }

  # Etichette dinamiche per i pannelli
  lag_labels <- sapply(lags, function(l) {
    if (l == 0) {
      return("Lag 0 (t = 0 ms)")
    } else {
      return(sprintf("Lag %d (t - %d ms)", l, round(l * dt_ms)))
    }
  })

  all_edges_list <- list()
  all_self_list  <- list()

  # --- 1. Estrazione dei coefficienti gamma_ij e gamma_ii per ogni Lag ---
  for (idx in seq_along(lags)) {
    lag <- lags[idx]
    lag_label <- lag_labels[idx]

    W_mat_lag <- matrix(0, nrow = K, ncol = K)

    for (i in active_neurons) {
      coef_mat <- step1_lag_coef_matrix(fit_results, i)
      if (!is.null(coef_mat)) {
        # lag 0 -> Colonna 1, lag 1 -> Colonna 2, etc.
        W_mat_lag[, i] <- coef_mat[, lag + 1]
      }
    }

    # 1A. Estrazione Autoconnessioni (gamma_ii sulla diagonale - 0 al Lag 0)
    gamma_ii <- diag(W_mat_lag)
    self_ids <- which(gamma_ii != 0)

    if (length(self_ids) > 0) {
      all_self_list[[length(all_self_list) + 1]] <- data.frame(
        Neuron = self_ids,
        X = centroids[self_ids, 1],
        Y = centroids[self_ids, 2],
        Weight = gamma_ii[self_ids],
        Lag = lag_label
      )
    }

    # 1B. Rimuoviamo la diagonale per le frecce inter-neurone
    diag(W_mat_lag) <- 0

    for (j in 1:K) {
      for (i in 1:K) {
        w <- W_mat_lag[j, i]
        if (w != 0) {
          all_edges_list[[length(all_edges_list) + 1]] <- data.frame(
            Source = j, Target = i,
            x_src = centroids[j, 1], y_src = centroids[j, 2],
            x_tgt = centroids[i, 1], y_tgt = centroids[i, 2],
            Weight = w,
            Lag = lag_label
          )
        }
      }
    }
  }

  df_edges <- if (length(all_edges_list) > 0) bind_rows(all_edges_list) else data.frame()
  df_self  <- if (length(all_self_list) > 0) bind_rows(all_self_list) else data.frame()

  # Filtro per eventuale neurone target isolato
  if (!is.null(target_id) && nrow(df_edges) > 0) {
    df_edges <- df_edges %>% filter(Source == target_id | Target == target_id)
  }

  if (nrow(df_edges) > 0) df_edges$Lag <- factor(df_edges$Lag, levels = lag_labels)
  if (nrow(df_self) > 0)  df_self$Lag  <- factor(df_self$Lag, levels = lag_labels)

  # Sfondo nodi dell'ippocampo replicati per ciascun pannello
  df_nodes <- bind_rows(lapply(lag_labels, function(lbl) {
    data.frame(
      Neuron = 1:K,
      X = centroids[, 1],
      Y = centroids[, 2],
      Lag = factor(lbl, levels = lag_labels)
    )
  }))

  # Nodi collegati (inter-neurone)
  if (nrow(df_edges) > 0) {
    df_connected <- df_edges %>%
      group_by(Lag) %>%
      do({
        conn_ids <- unique(c(.$Source, .$Target))
        data.frame(
          Neuron = conn_ids,
          X = centroids[conn_ids, 1],
          Y = centroids[conn_ids, 2]
        )
      }) %>%
      ungroup()
  } else {
    df_connected <- data.frame()
  }

  # --- 2. Rendering Grafico ggplot2 ---
  p <- ggplot() +
    # Layer 1: Nodi dell'ippocampo
    geom_point(data = df_nodes, aes(x = X, y = Y),
               color = "#E0E0E0", size = 1.2, alpha = 0.6) +

    # Layer 2: Frecce orientate per i collegamenti gamma_ij (j -> i)
    {
      if (nrow(df_edges) > 0) {
        geom_segment(
          data = df_edges,
          aes(x = x_src, y = y_src, xend = x_tgt, yend = y_tgt,
              color = Weight, alpha = Weight),
          arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
          linewidth = 0.55
        )
      }
    } +

    # Layer 3: Nodi coinvolti in connessioni inter-neurone
    {
      if (nrow(df_connected) > 0) {
        geom_point(data = df_connected, aes(x = X, y = Y),
                   color = "#2C3E50", size = 1.8, alpha = 0.9)
      }
    } +

    # Layer 4: NEURONI CON AUTO-CONNESSIONE (gamma_ii != 0)
    {
      if (nrow(df_self) > 0) {
        geom_point(data = df_self, aes(x = X, y = Y),
                   color = "#39B5FF", size = 1.8, alpha = 0.9)
      }
    } +

    # Disposizione dei grafici su 2 righe
    facet_wrap(~ Lag, nrow = 2) +

    # GRADIENTE DIVERGENTE CON ZERO = BIANCO PERFETTO
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "darkred",
      midpoint = 0,
      name = "Coeff. gamma"
    ) +

    scale_alpha_continuous(range = c(0.45, 0.95), guide = "none") +
    coord_fixed() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "gray30"),
      strip.background = element_rect(fill = "#EAECEE", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid = element_line(color = "gray95"),
      legend.position = "bottom"
    )

  return(p)
}


# ==============================================================================
# 8. DIAGNOSTICA: Goodness of Fit (Mappa di Devianza Spaziale)
# ==============================================================================
plot_neuron_deviance_map <- function(fit_results) {

  K <- fit_results$prep_data$K
  centroids <- fit_results$prep_data$centroids
  active_neurons <- fit_results$active_neurons

  dev_1se_vec <- rep(NA, K)

  for (i in active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i) && inherits(fit_i, "cv.glmnet")) {
      idx_1se <- fit_i$index["1se"]
      if (is.na(idx_1se) || is.null(idx_1se)) {
        idx_1se <- which(fit_i$lambda == fit_i$lambda.1se)
      }
      dev_1se_vec[i] <- fit_i$cvm[idx_1se]
    }
  }

  df_map <- data.frame(
    Neuron = 1:K,
    X = centroids[, 1],
    Y = centroids[, 2],
    Deviance_1se = dev_1se_vec
  )

  p <- ggplot(df_map, aes(x = X, y = Y)) +
    geom_point(data = dplyr::filter(df_map, is.na(Deviance_1se)), color = "#E0E0E0", size = 1.5, alpha = 0.5) +
    geom_point(data = dplyr::filter(df_map, !is.na(Deviance_1se)), aes(color = Deviance_1se), size = 2.8, alpha = 0.95) +
    scale_color_viridis_c(option = "plasma", name = "Deviance (1se)", direction = -1) +
    coord_fixed() +
    labs(
      title = "Spatial Map of Neuronal Deviance",
      subtitle = "Binomial Deviance at lambda.1se across Hippocampal Centroids",
      x = "X Centroid Position", y = "Y Centroid Position"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "bottom"
    )

  return(p)
}


# ==============================================================================
# 9. DIAGNOSTICA: Heatmap Sovrapposta Osservato vs Predetto
# ==============================================================================
plot_spike_overlay_heatmap <- function(predictions,
                                       active_neurons = NULL,
                                       frame_range = c(1, 1000),
                                       active_only = TRUE,
                                       threshold = 0.3) { # Threshold parameter chosen by user

  Z_obs  <- predictions$actual_spikes    # Observed matrix [K x T] (0/1)
  Z_pred <- predictions$predicted_probs  # Estimated matrix [K x T] (probabilities)
  dt_ms  <- predictions$time_ms

  K <- nrow(Z_obs)
  T_total <- ncol(Z_obs)

  # 1. Time window management (frames)
  if (is.null(frame_range)) {
    frame_range <- c(1, T_total)
  }

  start_frame <- max(1, frame_range[1])
  end_frame   <- min(T_total, frame_range[2])
  frames_idx  <- start_frame:end_frame

  # Time in seconds for the X-axis
  time_sec <- (frames_idx * dt_ms) / 1000
  time_step <- mean(diff(time_sec))

  # 2. Neuron selection (all or active only)
  if (active_only && !is.null(active_neurons)) {
    neuron_idx <- active_neurons
  } else {
    neuron_idx <- 1:K
  }

  # Subsetting matrices
  Z_obs_sub  <- Z_obs[neuron_idx, frames_idx, drop = FALSE]
  Z_pred_sub <- Z_pred[neuron_idx, frames_idx, drop = FALSE]

  # Binarizing predictions based on THRESHOLD
  Z_pred_bin <- ifelse(Z_pred_sub >= threshold, 1, 0)

  # 3. Data Frame Preparation (Filtering only Value == 1 for performance and overlay)
  neuron_factor <- factor(rep(neuron_idx, times = length(frames_idx)), levels = rev(neuron_idx))
  time_vector   <- rep(time_sec, each = length(neuron_idx))

  df_obs <- data.frame(
    Neuron = neuron_factor,
    Time   = time_vector,
    Value  = as.vector(Z_obs_sub)
  ) %>% dplyr::filter(Value == 1)

  df_pred <- data.frame(
    Neuron = neuron_factor,
    Time   = time_vector,
    Value  = as.vector(Z_pred_bin)
  ) %>% dplyr::filter(Value == 1)

  # 4. Overlaid Heatmap
  p <- ggplot() +
    # Layer 1: OBSERVED Spikes (Black - Full height)
    geom_tile(
      data = df_obs,
      aes(x = Time, y = Neuron, fill = "Observed (Z = 1)"),
      width = time_step, height = 0.85, alpha = 0.9
    ) +
    # Layer 2: PREDICTED Spikes (Red - Inner stripe to highlight overlaps)
    geom_tile(
      data = df_pred,
      aes(x = Time, y = Neuron, fill = "Predicted (P >= Threshold)"),
      width = time_step, height = 0.40, alpha = 0.85
    ) +
    scale_fill_manual(
      name = "Class",
      values = c(
        "Observed (Z = 1)" = "black",
        "Predicted (P >= Threshold)" = "red"
      )
    ) +
    labs(
      title = paste0("Observed vs. Predicted Spikes Overlay (Threshold = ", threshold, ")"),
      x = "Time (seconds)",
      y = "Neurons"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "gray30", hjust = 0),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top"
    )

  return(p)
}


# ==============================================================================
# 10. DIAGNOSTICA: Residual ACF e PACF
# ==============================================================================
plot_residual_acf_pacf <- function(fit_results, predictions, max_lag_plot = 20, residual_type = "pearson") {

  active_neurons <- fit_results$active_neurons
  Z <- predictions$actual_spikes
  P <- predictions$predicted_probs
  T_frames <- ncol(Z)

  # Matrices to store ACF and PACF values across active neurons [n_active x max_lag_plot]
  acf_mat  <- matrix(NA, nrow = length(active_neurons), ncol = max_lag_plot)
  pacf_mat <- matrix(NA, nrow = length(active_neurons), ncol = max_lag_plot)

  for (idx in seq_along(active_neurons)) {
    i <- active_neurons[idx]
    y <- Z[i, ]
    p <- P[i, ]

    # 1. Calculate GLM Residuals
    if (residual_type == "pearson") {
      eps <- 1e-5
      res <- (y - p) / sqrt(pmax(p * (1 - p), eps))
    } else if (residual_type == "raw") {
      res <- y - p
    }

    # 2. Compute ACF & PACF up to max_lag_plot
    res_acf  <- acf(res, lag.max = max_lag_plot, plot = FALSE)
    acf_mat[idx, ] <- res_acf$acf[2:(max_lag_plot + 1)] # Exclude lag 0 (always 1)

    res_pacf <- pacf(res, lag.max = max_lag_plot, plot = FALSE)
    pacf_mat[idx, ] <- res_pacf$acf[1:max_lag_plot]
  }

  # 95% Confidence threshold for Gaussian white noise
  ci_bound <- 1.96 / sqrt(T_frames)

  # Data frames for plotting
  df_acf <- data.frame(
    Lag    = 1:max_lag_plot,
    Mean   = colMeans(acf_mat, na.rm = TRUE),
    Median = apply(acf_mat, 2, median, na.rm = TRUE)
  )

  df_pacf <- data.frame(
    Lag    = 1:max_lag_plot,
    Mean   = colMeans(pacf_mat, na.rm = TRUE),
    Median = apply(pacf_mat, 2, median, na.rm = TRUE)
  )

  # 3. ACF Plot
  p_acf <- ggplot(df_acf, aes(x = Lag, y = Mean)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "#2980B9", linewidth = 0.8) +
    geom_point(color = "#2980B9", size = 2) +
    geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "#C0392B", linewidth = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      title = "Average Residual ACF Across Active Neurons",
      subtitle = paste0(capwords(residual_type), " residuals"),
      x = "Lag (frames)", y = "Autocorrelation (ACF)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "gray30")
    )

  # 4. PACF Plot
  p_pacf <- ggplot(df_pacf, aes(x = Lag, y = Mean)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "#8E44AD", linewidth = 0.8) +
    geom_point(color = "#8E44AD", size = 2) +
    geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "#C0392B", linewidth = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      title = "Average Residual PACF Across Active Neurons",
      subtitle = paste0(capwords(residual_type), " residuals"),
      x = "Lag (frames)", y = "Partial Autocorrelation (PACF)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "gray30")
    )

  grid.arrange(p_acf, p_pacf, ncol = 1)
}


# ==============================================================================
# 11. DIAGNOSTICA: Mappa Spaziale dei Coefficienti della Velocita'
# ==============================================================================
plot_neuron_speed_beta_map <- function(fit_results) {

  K <- fit_results$prep_data$K
  centroids <- fit_results$prep_data$centroids
  active_neurons <- fit_results$active_neurons
  n_fixed <- fit_results$prep_data$n_fixed_cols

  beta_speed_vec <- rep(NA, K)

  for (i in active_neurons) {
    # speed_spline e' l'ultima covariata fissa: colonna n_fixed della design
    # matrix completa. Si estrae da full_coefs perche' coef() restituisce nomi
    # posizionali V1..Vp (le colonne vengono rinominate dentro fit_model).
    full_coefs <- step1_full_coefs(fit_results, i)
    if (!is.null(full_coefs)) {
      beta_speed_vec[i] <- full_coefs[n_fixed]
    }
  }

  df_map <- data.frame(
    Neuron = 1:K,
    X = centroids[, 1],
    Y = centroids[, 2],
    Beta_Speed = beta_speed_vec
  )

  p <- ggplot(df_map, aes(x = X, y = Y)) +
    # Inactive / NA neurons in light grey
    geom_point(data = dplyr::filter(df_map, is.na(Beta_Speed)),
               color = "#E0E0E0", size = 1.5, alpha = 0.5) +
    # Active neurons colored by Speed Beta
    geom_point(data = dplyr::filter(df_map, !is.na(Beta_Speed)),
               aes(color = Beta_Speed), size = 2.8, alpha = 0.95) +
    # Diverging color palette centered at zero
    scale_color_gradient2(
      low = "#2980B9",
      mid = "#F7F9F9",
      high = "#C0392B",
      midpoint = 0,
      name = expression(beta["speed"])
    ) +
    coord_fixed() +
    labs(
      title = "Spatial Map of Velocity Coefficients",
      x = "X Centroid Position",
      y = "Y Centroid Position"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "right"
    )

  return(p)
}


# ==============================================================================
# 12. DIAGNOSTICA: Parametro di Penalizzazione (lambda) per Neurone
# ==============================================================================
plot_neuron_lambda_map <- function(fit_results) {

  K <- fit_results$prep_data$K
  active_neurons <- fit_results$active_neurons

  lambda_vec <- rep(NA, K)
  fit_type_vec <- rep(NA, K)

  for (i in active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i)) {
      if (inherits(fit_i, "cv.glmnet")) {
        lambda_vec[i] <- fit_i$lambda.1se
        fit_type_vec[i] <- "cv.glmnet (1se)"
      } else {
        mid_idx <- mid_lambda_index(fit_i$lambda)
        lambda_vec[i] <- fit_i$lambda[mid_idx]
        fit_type_vec[i] <- "glmnet (fallback)"
      }
    }
  }

  df_lambda <- data.frame(
    Neuron = 1:K,
    Lambda = lambda_vec,
    FitType = fit_type_vec
  ) %>% dplyr::filter(!is.na(Lambda))

  p <- ggplot(df_lambda, aes(x = Neuron, y = Lambda)) +
    # Vertical stems for a lollipop/pin style
    geom_segment(aes(x = Neuron, xend = Neuron, y = min(Lambda) * 0.8, yend = Lambda),
                 color = "gray80", linewidth = 0.4) +
    # Points for each neuron
    geom_point(aes(color = FitType), size = 2.4, alpha = 0.85) +
    # Logarithmic scale for lambda
    scale_y_log10() +
    scale_color_manual(
      values = c("cv.glmnet (1se)" = "#2E86C1", "glmnet (fallback)" = "#E67E22"),
      name = "Fit Method"
    ) +
    labs(
      title = expression("Choice of the Lasso Regularization Parameter"),
      x = "Neuron ID",
      y = expression(lambda * " (log scale)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank()
    ) +
    theme(
      legend.position = "none"
    )

  return(p)
}


# ==============================================================================
# 13. DIAGNOSTICA: Place Field del Neurone sulla Mappa dell'Arena
# ==============================================================================
plot_neuron_place_field <- function(fit_results, trial_data, neuron_id, grid_size = 100) {

  # --- 1. Allineamento Coordinate e Recupero Basi Spaziali ---
  t_behav <- trial_data$behav$time[, 1]
  pos_behav <- trial_data$behav$position

  t_raw <- trial_data$neuron$time[, 1]
  T_frames <- ncol(trial_data$neuron$S)
  t_neuron <- t_raw[seq(1, length(t_raw), length.out = T_frames)]

  pos_x <- approx(t_behav, pos_behav[, 1], xout = t_neuron, rule = 2)$y
  pos_y <- approx(t_behav, pos_behav[, 2], xout = t_neuron, rule = 2)$y

  # X_fixed = [X_spatial (n_spatial cols) | speed_spline (1 col)]
  n_fixed <- fit_results$prep_data$n_fixed_cols
  n_spatial <- n_fixed - 1
  n_splines <- round(sqrt(n_spatial))

  # Fit originale delle B-splines sulle traiettorie reali per ereditare i knots
  bx_fit <- bs(pos_x, df = n_splines)
  by_fit <- bs(pos_y, df = n_splines)

  # --- 2. Creazione della Griglia 2D dell'Arena ---
  x_grid <- seq(min(pos_x, na.rm = TRUE), max(pos_x, na.rm = TRUE), length.out = grid_size)
  y_grid <- seq(min(pos_y, na.rm = TRUE), max(pos_y, na.rm = TRUE), length.out = grid_size)
  grid_2d <- expand.grid(X = x_grid, Y = y_grid)

  # Valutazione delle B-splines sulla nuova griglia
  bx_grid <- predict(bx_fit, grid_2d$X)
  by_grid <- predict(by_fit, grid_2d$Y)

  # Prodotto tensoriale per generare la matrice di predittori della griglia
  X_spatial_grid <- do.call(cbind, lapply(1:ncol(bx_grid), function(i) bx_grid[, i] * by_grid))

  # --- 3. Estrazione dei Coefficienti Spaziali per il Neurone ---
  fit_i <- fit_results$models[[neuron_id]]
  valid_cols_i <- fit_results$valid_cols_list[[neuron_id]]

  if (is.null(fit_i) || is.null(valid_cols_i)) {
    stop(sprintf("Il neurone %d non e' tra quelli attivi stimati dal modello.", neuron_id))
  }

  if (inherits(fit_i, "cv.glmnet")) {
    raw_coefs <- as.matrix(coef(fit_i, s = "lambda.1se"))
  } else {
    mid_idx <- mid_lambda_index(fit_i$lambda)
    raw_coefs <- as.matrix(coef(fit_i))[, mid_idx, drop = FALSE]
  }

  # Intercetta beta_0
  beta_0 <- raw_coefs[1, 1]

  # Ripristino del vettore completo dei coefficienti su X_full
  full_coefs <- numeric(ncol(fit_results$prep_data$X))
  full_coefs[valid_cols_i] <- raw_coefs[-1, 1]

  # Selezione esplicita delle sole prime n_spatial colonne B-Spline
  beta_spatial <- full_coefs[1:n_spatial]

  # --- 4. Calcolo Probabilita' Baseline Spaziale P(Spike | X, Y) ---
  eta <- beta_0 + X_spatial_grid %*% beta_spatial
  prob_grid <- 1 / (1 + exp(-eta))                     # Inverse Logit
  prob_matrix <- matrix(prob_grid, nrow = grid_size, ncol = grid_size)

  soft_palette <- colorRampPalette(c(
    "grey90",
    "#5499C7",
    "#A9DFBF",
    "#F9E79F",
    "#F5B041",
    "#E67E22"
  ))

  # --- 5. Rendering Grafico (Mappa di Calore + Traiettoria Topo) ---
  filled.contour(
    x = x_grid,
    y = y_grid,
    z = prob_matrix,
    color.palette = soft_palette,
    xlab = "Position X",
    ylab = "Position Y",
    main = sprintf("Neuron %d", neuron_id),
    plot.axes = {
      axis(1)
      axis(2)
      # Sovrapponi la traiettoria reale del topo in nero semitrasparente
      lines(pos_x, pos_y, col = rgb(0, 0, 0, 0.25), lwd = 0.5)
    }
  )
}

rm(list = ls())

# ==============================================================================
# LIBRERIE E CARICAMENTO DATI
# ==============================================================================
library(glmnet)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pROC)
library(gridExtra)

all_trials <- readRDS("~/Downloads/all_trials.RDS")
load("~/Downloads/mouse_trial.RData")
data_w_dist <- df.neuron[, c('speed_spline', 'dist_to_hull', 'time')]

experiment <- all_trials$M3424F$trial2


# ==============================================================================
# FASE 1: Preparazione Design Matrix (Covariate Fisse + Lags Autoregressivi)
# ==============================================================================
prepare_design_matrix <- function(trial_data, data_w_dist, max_lag_frames = 3) {
  
  # --- 1. Estrarre e Allineare i Tempi Neurali ---
  S <- trial_data$neuron$S                         # [K x T_frames]
  centroids <- trial_data$neuron$centroid           # [K x 2]
  t_raw <- trial_data$neuron$time[, 1]              
  
  K <- nrow(S)
  T_frames <- ncol(S)
  
  # Downsampling dei tempi del calcio per matchare i T_frames di S
  t_neuron <- t_raw[seq(1, length(t_raw), length.out = T_frames)]
  dt_ms <- mean(diff(t_neuron))
  
  # Matrice binaria degli eventi Z
  Z <- ifelse(S > 0, 1, 0)                         # [K x T_frames]
  
  # --- 2. Interpolazione Covariate Fisse (data_w_dist) sui Tempi Neurali ---
  speed_interp <- approx(x = data_w_dist$time, y = data_w_dist$speed_spline, xout = t_neuron, rule = 2)$y
  dist_interp  <- approx(x = data_w_dist$time, y = data_w_dist$dist_to_hull,  xout = t_neuron, rule = 2)$y
  
  X_fixed <- cbind(
    speed_spline = speed_interp, 
    dist_to_hull = dist_interp
  )
  n_fixed_cols <- ncol(X_fixed) # 2 covariate fisse
  
  # --- 3. Storico Autoregressivo di Rete (Lags) ---
  X_lags <- list()
  for (lag in 1:max_lag_frames) {
    Z_lagged <- cbind(matrix(0, nrow = K, ncol = lag), Z[, 1:(T_frames - lag)])
    rownames(Z_lagged) <- paste0("N", 1:K, "_lag", lag)
    X_lags[[lag]] <- Z_lagged
  }
  X_lags <- do.call(rbind, X_lags) %>% t()        # [T_frames x (K * max_lag_frames)]
  
  # Final Design Matrix: [X_fixed | X_lags]
  X_full <- cbind(X_fixed, X_lags)
  
  return(list(
    X = X_full,
    Z = Z,
    K = K,
    T_frames = T_frames,
    dt_ms = dt_ms,
    centroids = centroids,
    n_fixed_cols = n_fixed_cols,
    max_lag_frames = max_lag_frames
  ))
}


# ==============================================================================
# FASE 2: Fit del Modello con CV Temporale e Lasso Pesato Spazialmente
# ==============================================================================
fit_model <- function(prep_data, family = "binomial", min_spikes = 15, nfolds = 5, block_dur_sec = 2.0) {
  
  K <- prep_data$K
  X_full <- prep_data$X
  Z <- prep_data$Z
  n_fixed <- prep_data$n_fixed_cols
  max_lag_frames <- prep_data$max_lag_frames
  T_frames <- prep_data$T_frames
  dt_ms <- prep_data$dt_ms
  
  # Matrice distanze FOV
  dist_matrix <- as.matrix(dist(prep_data$centroids))
  max_dist <- max(dist_matrix)
  
  # Consideriamo solo neuroni con un numero sufficiente di spike
  active_neurons <- which(rowSums(Z) >= min_spikes)
  
  models <- vector("list", K)
  W_matrix <- matrix(0, nrow = K, ncol = K)
  
  # Griglia Lambda Fissa
  lambda_grid <- exp(seq(log(0.1), log(0.0001), length.out = 50))
  
  # --- Cross-Validation Temporale a Blocchi (Interleaved) ---
  block_size_frames <- max(1, round((block_dur_sec * 1000) / dt_ms))
  n_blocks <- ceiling(T_frames / block_size_frames)
  block_ids <- rep(1:n_blocks, each = block_size_frames, length.out = T_frames)
  foldid <- ((block_ids - 1) %% nfolds) + 1
  
  for (i in active_neurons) {
    y <- Z[i, ]
    
    # Covariate fisse non penalizzate (0), Lags penalizzati in base alla distanza
    pen_fixed <- rep(0, n_fixed)
    raw_pen_lags <- rep(dist_matrix[i, ] / max_dist, times = max_lag_frames)
    pen_lags <- pmax(raw_pen_lags, 0.05)
    penalty_full <- c(pen_fixed, pen_lags)
    
    # Filtraggio colonne a varianza zero
    col_counts <- colSums(X_full > 0)
    valid_cols <- c(1:n_fixed, which(col_counts >= 2 & (1:ncol(X_full)) > n_fixed))
    
    X_sub <- X_full[, valid_cols, drop = FALSE]
    penalty_sub <- penalty_full[valid_cols]
    
    dimnames(X_sub) <- list(NULL, paste0("V", 1:ncol(X_sub)))
    X_sparse <- Matrix(X_sub, sparse = TRUE)
    
    fit_done <- FALSE
    coefs_sub <- NULL
    
    # Fit con CV
    tryCatch({
      spikes_per_fold <- tapply(y, foldid, sum)
      if (any(spikes_per_fold == 0)) {
        stop("At least one temporal fold contains 0 spikes.")
      }
      
      fit <- cv.glmnet(
        x = X_sparse, y = y,
        family = family,
        penalty.factor = penalty_sub,
        foldid = foldid,
        lambda = lambda_grid,                     
        type.measure = "deviance",
        alpha = 1, maxit = 100000, parallel = TRUE
      )
      models[[i]] <- fit
      coefs_sub <- as.matrix(coef(fit, s = "lambda.1se"))[-1, 1]
      fit_done <- TRUE
    }, error = function(e) {
      # Fallback su glmnet senza CV
      tryCatch({
        fit_fallback <- glmnet(
          x = X_sparse, y = y,
          family = family,
          penalty.factor = penalty_sub,
          lambda = lambda_grid,
          alpha = 1, maxit = 100000
        )
        models[[i]] <- fit_fallback
        mid_idx <- round(length(lambda_grid) / 2)
        coefs_sub <<- as.matrix(coef(fit_fallback))[-1, mid_idx]
        fit_done <<- TRUE
      }, error = function(e2) {
        message(sprintf("Neuron %d avoided for numeric instability.", i))
      })
    })
    
    # Ricostruzione Matrice di Adiacenza W
    if (fit_done && !is.null(coefs_sub)) {
      full_coefs <- numeric(ncol(X_full))
      full_coefs[valid_cols] <- coefs_sub
      
      coef_lags <- full_coefs[(n_fixed + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag_frames)
      
      W_matrix[, i] <- rowSums(coef_mat)
    }
  }
  
  return(list(
    models = models,
    adj_matrix = W_matrix,
    prep_data = prep_data,
    active_neurons = active_neurons,
    valid_cols = valid_cols,
    foldid = foldid
  ))
}


# ==============================================================================
# FASE 3: Predizione In-Sample / Out-of-Sample
# ==============================================================================
predict_model <- function(fit_results, test_trial_data, data_w_dist) {
  
  max_lag <- fit_results$prep_data$max_lag_frames
  valid_cols <- fit_results$valid_cols
  
  prep_test <- prepare_design_matrix(
    test_trial_data, 
    data_w_dist = data_w_dist, 
    max_lag_frames = max_lag
  )
  
  X_test_sub <- prep_test$X[, valid_cols, drop = FALSE]
  dimnames(X_test_sub) <- list(NULL, paste0("V", 1:ncol(X_test_sub)))
  X_test_sparse <- Matrix(X_test_sub, sparse = TRUE)
  
  K <- prep_test$K
  T_test <- prep_test$T_frames
  
  prob_matrix <- matrix(0, nrow = K, ncol = T_test)
  
  for (i in fit_results$active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i)) {
      if (inherits(fit_i, "cv.glmnet")) {
        probs <- predict(fit_i, newx = X_test_sparse, s = "lambda.1se", type = "response")
      } else {
        mid_lambda <- fit_i$lambda[round(length(fit_i$lambda)/2)]
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
# FASE 4: Valutazione Prestazioni (ROC & AUC)
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
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    theme_minimal()
  
  return(list(
    auc_vector = auc_list,
    mean_auc = mean(auc_list, na.rm = TRUE),
    p_roc = p_roc
  ))
}


# ==============================================================================
# FASE 5: Visualizzazione Effetti Covariate Fisse per Singolo Neurone
# ==============================================================================
plot_neuron_covariate_effects <- function(fit_results, data_w_dist, neuron_id, grid_size = 80) {
  
  fit_i <- fit_results$models[[neuron_id]]
  if (is.null(fit_i)) {
    stop(sprintf("Neuron %d is not active in the fitted model.", neuron_id))
  }
  
  if (inherits(fit_i, "cv.glmnet")) {
    coefs <- as.matrix(coef(fit_i, s = "lambda.1se"))
  } else {
    mid_idx <- round(length(fit_i$lambda) / 2)
    coefs <- as.matrix(coef(fit_i))[, mid_idx, drop = FALSE]
  }
  
  beta_0 <- coefs[1, 1]
  beta_speed <- coefs[2, 1]
  beta_dist  <- coefs[3, 1]
  
  speed_seq <- seq(min(data_w_dist$speed_spline), max(data_w_dist$speed_spline), length.out = grid_size)
  dist_seq  <- seq(min(data_w_dist$dist_to_hull), max(data_w_dist$dist_to_hull), length.out = grid_size)
  
  grid_2d <- expand.grid(Speed = speed_seq, Dist = dist_seq)
  grid_2d$Eta <- beta_0 + beta_speed * grid_2d$Speed + beta_dist * grid_2d$Dist
  grid_2d$Prob <- 1 / (1 + exp(-grid_2d$Eta))
  
  p <- ggplot(grid_2d, aes(x = Speed, y = Dist, fill = Prob)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", name = "P(Spike)") +
    labs(
      title = sprintf("Covariate Response Surface - Neuron %d", neuron_id),
      subtitle = sprintf("Beta Speed: %.3f | Beta Dist to Hull: %.3f", beta_speed, beta_dist),
      x = "Speed (spline)",
      y = "Distance to Hull"
    ) +
    theme_minimal()
  
  return(p)
}


# ==============================================================================
# FASE 6: Heatmap dei Coefficienti Autoregressivi (Lag Singolo)
# ==============================================================================
plot_single_lag <- function(fit_results, lag = 1, only_active = TRUE) {
  
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  n_fixed <- fit_results$prep_data$n_fixed_cols
  active_neurons <- fit_results$active_neurons
  valid_cols <- fit_results$valid_cols
  
  if (lag < 1 || lag > max_lag) {
    stop(sprintf("Lag %d is invalid! Must be between 1 and %d.", lag, max_lag))
  }
  
  idx_plot <- if (only_active) active_neurons else 1:K
  N_plot <- length(idx_plot)
  
  W_mat_full <- matrix(0, nrow = K, ncol = K)
  
  for (i in active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i)) {
      if (inherits(fit_i, "cv.glmnet")) {
        coefs_sub <- as.matrix(coef(fit_i, s = "lambda.1se"))[-1, 1]
      } else {
        mid_idx <- round(length(fit_i$lambda) / 2)
        coefs_sub <- as.matrix(coef(fit_i))[-1, mid_idx]
      }
      
      full_coefs <- numeric(ncol(fit_results$prep_data$X))
      full_coefs[valid_cols] <- coefs_sub
      
      coef_lags <- full_coefs[(n_fixed + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag)
      
      W_mat_full[, i] <- coef_mat[, lag]
    }
  }
  
  W_sub <- W_mat_full[idx_plot, idx_plot, drop = FALSE]
  
  df_single_lag <- data.frame(
    Source = factor(rep(idx_plot, times = N_plot), levels = rev(idx_plot)),
    Target = factor(rep(idx_plot, each = N_plot), levels = idx_plot),
    Weight = as.vector(W_sub)
  )
  
  p <- ggplot(df_single_lag, aes(x = Target, y = Source, fill = Weight)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#F8F9F9", "#CBD5E1", "#F5B041", "#E74C3C", "#78281F"),
      name = "Coeff. gamma"
    ) +
    labs(x = "Target neuron", y = "Input neuron", title = sprintf("Adjacency Matrix - Lag %d", lag)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  
  return(p)
}


# ==============================================================================
# FASE 7: Mappa Spaziale Multi-Lag dei Coefficienti Autoregressivi
# ==============================================================================
plot_spatial_multi_lag <- function(fit_results, lags = c(1, 2, 3), target_id = NULL) {
  
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  n_fixed <- fit_results$prep_data$n_fixed_cols
  active_neurons <- fit_results$active_neurons
  valid_cols <- fit_results$valid_cols
  dt_ms <- fit_results$prep_data$dt_ms
  centroids <- fit_results$prep_data$centroids
  
  if (any(lags < 1 | lags > max_lag)) {
    stop(sprintf("Uno o più lag non sono validi! Scegli valori compresi tra 1 e %d.", max_lag))
  }
  
  lag_labels <- sprintf("Lag %d (t - %d ms)", lags, round(lags * dt_ms))
  
  all_edges_list <- list()
  all_self_list  <- list()
  
  # --- 1. Estrazione dei coefficienti gamma_ij e gamma_ii per ogni Lag ---
  for (idx in seq_along(lags)) {
    lag <- lags[idx]
    lag_label <- lag_labels[idx]
    
    W_mat_lag <- matrix(0, nrow = K, ncol = K)
    
    for (i in active_neurons) {
      fit_i <- fit_results$models[[i]]
      if (!is.null(fit_i)) {
        if (inherits(fit_i, "cv.glmnet")) {
          coefs_sub <- as.matrix(coef(fit_i, s = "lambda.1se"))[-1, 1]
        } else {
          mid_idx <- round(length(fit_i$lambda) / 2)
          coefs_sub <- as.matrix(coef(fit_i))[-1, mid_idx]
        }
        
        full_coefs <- numeric(ncol(fit_results$prep_data$X))
        full_coefs[valid_cols] <- coefs_sub
        
        coef_lags <- full_coefs[(n_fixed + 1):length(full_coefs)]
        coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag)
        
        W_mat_lag[, i] <- coef_mat[, lag]
      }
    }
    
    # 1A. Estrazione Autoconnessioni (gamma_ii sulla diagonale)
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
    
    # 1B. Rimuoviamo la diagonale per le frecce di connessione inter-neurone
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
    
    facet_wrap(~ Lag, nrow = 1) +
    
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
      legend.position = "right"
    )
  
  return(p)
}




# ==============================================================================
# FASE 8: Goodness of Fit (Mappa di Devianza Spaziale)
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
      legend.position = "right"
    )
  
  return(p)
}


# ==============================================================================
# PIPELINE DI ESECUZIONE COMPLETA
# ==============================================================================

# 1. Preparazione Design Matrix
prep_t1 <- prepare_design_matrix(
  trial_data = experiment, 
  data_w_dist = data_w_dist, 
  max_lag_frames = 3
)

# 2. Fit del Modello
fit_t1 <- fit_model(
  prep_data = prep_t1, 
  family = "binomial", 
  min_spikes = 15, 
  nfolds = 5, 
  block_dur_sec = 2.0
)

# 3. Predizione In-Sample
predictions_t2 <- predict_model(
  fit_results = fit_t1, 
  test_trial_data = experiment, 
  data_w_dist = data_w_dist
)

# 4. Valutazione ROC
eval_res <- evaluate_predictions(predictions_t2, fit_t1$active_neurons)
print(eval_res$p_roc)

# 5. Effetto Covariate su un Neurone
neuron_da_mostrare <- fit_t1$active_neurons[1]
p_cov_effect <- plot_neuron_covariate_effects(fit_t1, data_w_dist, neuron_id = neuron_da_mostrare)
print(p_cov_effect)

# 6. Heatmaps dei Lags
p_lag1 <- plot_single_lag(fit_t1, lag = 1, only_active = TRUE)
print(p_lag1)

# 7. Grafico Spaziale Multi-Lag
p_multi <- plot_spatial_multi_lag(fit_t1, lags = c(1, 2, 3))
print(p_multi)

# 8. Mappa di Devianza Spaziale
p_deviance_map <- plot_neuron_deviance_map(fit_t1)
print(p_deviance_map)

# 9. Heatmap Sovrapposta
p_overlay <- plot_spike_overlay_heatmap(
  predictions = predictions_t2,
  active_neurons = fit_t1$active_neurons,
  frame_range = c(1, 120),
  active_only = TRUE,
  threshold = 0.15
)
print(p_overlay)

# 10. Salvataggio Risultati
save(
  fit_t1, 
  predictions_t2, 
  prep_t1, 
  file = "~/Downloads/fit_t1_results.rda"
)

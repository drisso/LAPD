# CV MODIFICATA E AGGIUNTO EFFETTO FISSO DELLA VELOCITà

rm(list=ls())

library(dplyr)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(Matrix)
library(pROC)
library(splines
library(tidyr)

all_trials <- readRDS("~/Downloads/all_trials.RDS")
load("~/Downloads/mouse_trial.RData")
data_w_dist <- df.neuron[, c('speed_spline', 'time')]

experiment <- all_trials$M3424F$trial2


# ==============================================================================
# FASE 1: Preparazione Dati (B-Splines Spaziali + Velocità + Lags) ####
# ==============================================================================
prepare_design_matrix <- function(trial_data, data_w_dist, max_lag_frames = 5, n_splines = 4) {
  
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
  
  # --- 2. Interpolazione Posizione Arena e Velocità sui Tempi Neurali ---
  t_behav <- trial_data$behav$time[, 1]
  pos_behav <- trial_data$behav$position
  
  pos_x <- approx(x = t_behav, y = pos_behav[, 1], xout = t_neuron, rule = 2)$y
  pos_y <- approx(x = t_behav, y = pos_behav[, 2], xout = t_neuron, rule = 2)$y
  
  # Interpolazione della velocità da data_w_dist
  speed_interp <- approx(x = data_w_dist$time, y = data_w_dist$speed_spline, xout = t_neuron, rule = 2)$y
  
  # --- 3. Basi Spaziali B-Spline 2D + Covariata Velocità ---
  bx <- bs(pos_x, df = n_splines)
  by <- bs(pos_y, df = n_splines)
  
  X_spatial_list <- lapply(1:ncol(bx), function(i) bx[, i] * by)
  X_spatial <- do.call(cbind, X_spatial_list)      # [T_frames x (n_splines^2)]
  colnames(X_spatial) <- paste0("spatial_basis_", 1:ncol(X_spatial))
  
  # Matrice delle covariate fisse (Spazio 2D + Velocità)
  X_fixed <- cbind(X_spatial, speed_spline = speed_interp)
  n_fixed_cols <- ncol(X_fixed)
  
  # --- 4. Storico Autoregressivo di Rete (Lags) ---
  X_lags <- list()
  for (lag in 1:max_lag_frames) {
    Z_lagged <- cbind(matrix(0, nrow = K, ncol = lag), Z[, 1:(T_frames - lag)])
    rownames(Z_lagged) <- paste0("N", 1:K, "_lag", lag)
    X_lags[[lag]] <- Z_lagged
  }
  X_lags <- do.call(rbind, X_lags) %>% t()        # [T_frames x (K * max_lag_frames)]
  
  # Design Matrix Finale: [X_fixed | X_lags]
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
# FASE 2: Fit del Modello ####
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
  lambda_grid <- exp(seq(log(1), log(0.0001), length.out = 20))
  
  # Cross-Validation Temporale a Blocchi
  block_size_frames <- max(1, round((block_dur_sec * 1000) / dt_ms))
  n_blocks <- ceiling(T_frames / block_size_frames)
  block_ids <- rep(1:n_blocks, each = block_size_frames, length.out = T_frames)
  foldid <- ((block_ids - 1) %% nfolds) + 1
  
  for (i in active_neurons) {
    y <- Z[i, ]
    
    # Penalità: 0 per covariate fisse (spazio + velocità), pesata per i lags
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
        stop("At least one temporal fold contains 0 spikes for this neuron.")
      }
      
      fit <- cv.glmnet(
        x = X_sparse,
        y = y,
        family = family,
        penalty.factor = penalty_sub,
        foldid = foldid,
        lambda = lambda_grid,                     
        type.measure = "deviance",
        alpha = 1,
        maxit = 100000, 
        parallel = TRUE
      )
      models[[i]] <- fit
      coefs_sub <- as.matrix(coef(fit, s = "lambda.1se"))[-1, 1]
      fit_done <- TRUE
    }, error = function(e) {
      # Fallback su glmnet senza CV
      tryCatch({
        fit_fallback <- glmnet(
          x = X_sparse,
          y = y,
          family = family,
          penalty.factor = penalty_sub,
          lambda = lambda_grid,
          alpha = 1,
          maxit = 100000
        )
        models[[i]] <- fit_fallback
        
        mid_idx <- round(length(lambda_grid) / 2)
        coefs_sub <<- as.matrix(coef(fit_fallback))[-1, mid_idx]
        fit_done <<- TRUE
      }, error = function(e2) {
        message(sprintf("Neuron %d CV avoided for numeric instability.", i))
      })
    })
    
    # Ricostruzione matrice W
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
# FASE 3: Predizione In-Sample / Out-of-Sample ####
# ==============================================================================
predict_model <- function(fit_results, test_trial_data, data_w_dist) {
  
  max_lag <- fit_results$prep_data$max_lag_frames
  valid_cols <- fit_results$valid_cols
  
  # Prepara la design matrix del test incorporando data_w_dist
  prep_test <- prepare_design_matrix(
    trial_data = test_trial_data, 
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
# FASE 4: Valutazione Prestazioni (ROC & AUC) ####
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
# FASE 5: Visualizzazione Effetti Covariate Fisse per Singolo Neurone ####
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
# FASE 6: Heatmap dei Coefficienti Autoregressivi (Lag Singolo) ####
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


# FASE 7: Mappa Spaziale Multi-Lag dei Coefficienti Autoregressivi ####

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
    
    # Disposizione su 2 righe
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
# FASE 8: Goodness of Fit (Mappa di Devianza Spaziale) ####
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
# FASE 9: Heatmap Sovrapposta ####
# ==============================================================================

plot_spike_overlay_heatmap <- function(predictions, 
                                       active_neurons = NULL, 
                                       frame_range = c(1, 1000), 
                                       active_only = TRUE,
                                       threshold = 0.3) { # Parametro soglia scelto dall'utente
  
  Z_obs  <- predictions$actual_spikes    # Matrice [K x T] osservata (0/1)
  Z_pred <- predictions$predicted_probs # Matrice [K x T] stimata (probabilità)
  dt_ms  <- predictions$time_ms
  
  K <- nrow(Z_obs)
  T_total <- ncol(Z_obs)
  
  # 1. Gestione finestra temporale (frames)
  if (is.null(frame_range)) {
    frame_range <- c(1, T_total)
  }
  
  start_frame <- max(1, frame_range[1])
  end_frame   <- min(T_total, frame_range[2])
  frames_idx  <- start_frame:end_frame
  
  # Tempo in secondi per l'asse X
  time_sec <- (frames_idx * dt_ms) / 1000
  time_step <- mean(diff(time_sec))
  
  # 2. Selezione dei neuroni (tutti o solo gli attivi)
  if (active_only && !is.null(active_neurons)) {
    neuron_idx <- active_neurons
  } else {
    neuron_idx <- 1:K
  }
  
  # Subsetting delle matrici
  Z_obs_sub  <- Z_obs[neuron_idx, frames_idx, drop = FALSE]
  Z_pred_sub <- Z_pred[neuron_idx, frames_idx, drop = FALSE]
  
  # Binarizzazione delle predizioni in base alla SOGLIA
  Z_pred_bin <- ifelse(Z_pred_sub >= threshold, 1, 0)
  
  # 3. Preparazione Data Frame (Filtriamo solo per Value == 1 per velocizzare e sovrapporre)
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
  
  # 4. Heatmap Sovrapposta
  p <- ggplot() +
    # Layer 1: Spikes OSSERVATI (Nero - Altezza completa)
    geom_tile(
      data = df_obs, 
      aes(x = Time, y = Neuron, fill = "Osservato (Z = 1)"), 
      width = time_step, height = 0.85, alpha = 0.9
    ) +
    # Layer 2: Spikes PREDETTI (Rosso - Striscia centrale per evidenziare sovrapposizioni)
    geom_tile(
      data = df_pred, 
      aes(x = Time, y = Neuron, fill = "Predetto (P >= Soglia)"), 
      width = time_step, height = 0.40, alpha = 0.85
    ) +
    scale_fill_manual(
      name = "Classe",
      values = c(
        "Osservato (Z = 1)" = "black", 
        "Predetto (P >= Soglia)" = "red"
      )
    ) +
    labs(
      title = paste0("Sovrapposizione Spikes Osservati vs. Predetti (Soglia = ", threshold, ")"),
      subtitle = "Nero = Spike Reale | Rosso = Spike Predetto | Rosso dentro Nero = Match (TP)",
      x = "Tempo (secondi)", 
      y = "Neuroni"
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
# FASE 10: Residual ACF e PACF ####
# ==============================================================================

plot_residual_acf_pacf <- function(fit_results, predictions, max_lag_plot = 20, residual_type = "pearson") {
  
  active_neurons <- fit_results$active_neurons
  Z <- predictions$actual_spikes
  P <- predictions$predicted_probs
  T_frames <- ncol(Z)
  fitted_max_lag <- fit_results$prep_data$max_lag_frames
  
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
    geom_vline(xintercept = fitted_max_lag + 0.5, linetype = "dotted", color = "#27AE60", linewidth = 1) +
    annotate("text", x = fitted_max_lag + 0.8, y = max(df_acf$Mean, ci_bound * 1.5), 
             label = paste0("Fitted Lag = ", fitted_max_lag), color = "#27AE60", hjust = 0, size = 3.5, fontface = "bold") +
    labs(
      title = "Average Residual ACF Across Active Neurons",
      subtitle = paste0("Residuals: ", toupper(residual_type), " | Dashed Red: 95% CI White Noise Limit"),
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
    geom_vline(xintercept = fitted_max_lag + 0.5, linetype = "dotted", color = "#27AE60", linewidth = 1) +
    annotate("text", x = fitted_max_lag + 0.8, y = max(df_pacf$Mean, ci_bound * 1.5), 
             label = paste0("Fitted Lag = ", fitted_max_lag), color = "#27AE60", hjust = 0, size = 3.5, fontface = "bold") +
    labs(
      title = "Average Residual PACF Across Active Neurons",
      subtitle = paste0("Residuals: ", toupper(residual_type), " | Dashed Red: 95% CI White Noise Limit"),
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
# FASE 11: Spatial Map of Speed Beta Coefficients ####
# ==============================================================================
plot_neuron_speed_beta_map <- function(fit_results) {
  
  K <- fit_results$prep_data$K
  centroids <- fit_results$prep_data$centroids
  active_neurons <- fit_results$active_neurons
  
  beta_speed_vec <- rep(NA, K)
  
  for (i in active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i)) {
      if (inherits(fit_i, "cv.glmnet")) {
        coefs <- as.matrix(coef(fit_i, s = "lambda.1se"))
      } else {
        mid_idx <- round(length(fit_i$lambda) / 2)
        coefs <- as.matrix(coef(fit_i))[, mid_idx, drop = FALSE]
      }
      
      # Extract speed_spline coefficient (1st predictor after intercept)
      if ("speed_spline" %in% rownames(coefs)) {
        beta_speed_vec[i] <- coefs["speed_spline", 1]
      } else {
        beta_speed_vec[i] <- coefs[2, 1] # fallback index
      }
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
      title = "Spatial Map of Speed Beta Coefficients",
      subtitle = expression("Speed Spline Coefficient (" * beta["speed"] * ") at " * lambda["1se"] * " across Centroids"),
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
# FASE 12: Penalization Parameter (lambda) per Neuron  ####
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
        mid_idx <- round(length(fit_i$lambda) / 2)
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
      title = expression("Penalization Parameter (" * lambda * ") per Neuron"),
      subtitle = expression("Optimal L1 Regularization Penalty (" * lambda["1se"] * ") across Neurons"),
      x = "Neuron ID",
      y = expression(lambda * " (log scale)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

# ==============================================================================
# ESECUZIONE ####
# ==============================================================================

# 1. Preparazione Dati
prep_t1 <- prepare_design_matrix(
  trial_data = experiment, 
  data_w_dist = data_w_dist, 
  max_lag_frames = 5
)

# 2. Fit del Modello
fit_t1 <- fit_model(
  prep_data = prep_t1, 
  family = "binomial", 
  min_spikes = 15, 
  nfolds = 5
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
p_multi <- plot_spatial_multi_lag(fit_t1, lags = c(1, 2, 3, 4, 5))
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

# 10. Run residual ACF/PACF diagnostic up to 20 lags
plot_residual_acf_pacf(
  fit_results = fit_t1, 
  predictions = predictions_t2, 
  max_lag_plot = 20, 
  residual_type = "pearson"
)

# 11. Plot Spatial Speed Beta Map
p_speed_beta <- plot_neuron_speed_beta_map(fit_t1)
print(p_speed_beta)

# 12. Plot Lambda parameters per neuron
p_lambda <- plot_neuron_lambda_map(fit_t1)
print(p_lambda)

# FINAL. Salvataggio Risultati
save(
  fit_t1, 
  predictions_t2, 
  prep_t1, 
  file = "~/Downloads/fit_t1_results versione 3.rda"
)

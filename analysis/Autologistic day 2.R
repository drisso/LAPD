rm(list=ls())

all_trials <- readRDS("~/Downloads/all_trials.RDS")

library(glmnet)
library(splines)
library(Matrix)
library(dplyr)
library(splines)
library(ggplot2)
library(tidyr)
library(pROC)

# ==============================================================================
# FASE 1: Preparazione e Allineamento Dati (Baseline Spaziale + Lags)
# ==============================================================================
prepare_design_matrix <- function(trial_data, max_lag_frames = 3, n_splines = 4) {
  
  # --- 1. Estrarre e Allineare i Tempi Neurali ---
  S <- trial_data$neuron$S                         # [K x T_frames] (es. 462 x 10870)
  centroids <- trial_data$neuron$centroid           # [K x 2]
  t_raw <- trial_data$neuron$time[, 1]              # [21740]
  
  K <- nrow(S)
  T_frames <- ncol(S)
  
  # Downsampling dei tempi del calcio per matchare i 10870 frame di S
  t_neuron <- t_raw[seq(1, length(t_raw), length.out = T_frames)]
  dt_ms <- mean(diff(t_neuron))
  
  # Matrice binaria degli eventi Z
  Z <- ifelse(S > 0, 1, 0)                         # [K x T_frames]
  
  # --- 2. Interpolazione Posizione Arena (behav) sui Tempi del Calcio ---
  t_behav <- trial_data$behav$time[, 1]
  pos_behav <- trial_data$behav$position
  
  pos_x <- approx(x = t_behav, y = pos_behav[, 1], xout = t_neuron, rule = 2)$y
  pos_y <- approx(x = t_behav, y = pos_behav[, 2], xout = t_neuron, rule = 2)$y
  
  # --- 3. Basi Spaziali B-Spline 2D per mu_k(X,Y) [Place Fields] ---
  bx <- bs(pos_x, df = n_splines)
  by <- bs(pos_y, df = n_splines)
  
  # Prodotto tensoriale per generare la griglia di basi spaziali 2D
  X_spatial_list <- lapply(1:ncol(bx), function(i) bx[, i] * by)
  X_spatial <- do.call(cbind, X_spatial_list)      # [T_frames x (n_splines^2)]
  colnames(X_spatial) <- paste0("spatial_basis_", 1:ncol(X_spatial))
  
  # --- 4. Storico Autoregressivo di Rete (Lags) ---
  X_lags <- list()
  for (lag in 1:max_lag_frames) {
    Z_lagged <- cbind(matrix(0, nrow = K, ncol = lag), Z[, 1:(T_frames - lag)])
    rownames(Z_lagged) <- paste0("N", 1:K, "_lag", lag)
    X_lags[[lag]] <- Z_lagged
  }
  X_lags <- do.call(rbind, X_lags) %>% t()        # [T_frames x (K * max_lag_frames)]
  
  # Design Matrix Finale: [X_spaziale | X_lags]
  X_full <- cbind(X_spatial, X_lags)
  
  return(list(
    X = X_full,
    Z = Z,
    K = K,
    T_frames = T_frames,
    dt_ms = dt_ms,
    centroids = centroids,
    n_spatial_cols = ncol(X_spatial),
    max_lag_frames = max_lag_frames
  ))
}


# ==============================================================================
# FASE 2: Stima dei Parametri tramite Poisson/Binomial GLM Regularizzato
# ==============================================================================
fit_model <- function(prep_data, family = "binomial", min_spikes = 15, nfolds = 5) {
  
  K <- prep_data$K #number of neurons
  X_full <- prep_data$X
  Z <- prep_data$Z
  n_spatial <- prep_data$n_spatial_cols
  max_lag_frames <- prep_data$max_lag_frames
  
  # Matrice distanze FOV
  dist_matrix <- as.matrix(dist(prep_data$centroids))
  max_dist <- max(dist_matrix)
  
  # Consideriamo solo neuroni con un numero sufficiente di spike
  active_neurons <- which(rowSums(Z) >= min_spikes)
  
  models <- vector("list", K)
  W_matrix <- matrix(0, nrow = K, ncol = K)
  
  # 1. GRIGLIA LAMBDA FISSA: Forza tutti i fold ad avere la stessa identica dimensione
  lambda_grid <- exp(seq(log(0.1), log(0.0001), length.out = 50))
  
  for (i in active_neurons) {
    y <- Z[i, ]
    
    # Penalità spaziale e lags
    pen_spatial <- rep(0, n_spatial)
    raw_pen_lags <- rep(dist_matrix[i, ] / max_dist, times = max_lag_frames)
    pen_lags <- pmax(raw_pen_lags, 0.05)
    penalty_full <- c(pen_spatial, pen_lags)
    
    # Filtraggio colonne a varianza zero
    col_counts <- colSums(X_full > 0)
    valid_cols <- c(1:n_spatial, which(col_counts >= 2 & (1:ncol(X_full)) > n_spatial))
    
    X_sub <- X_full[, valid_cols, drop = FALSE]
    penalty_sub <- penalty_full[valid_cols]
    
    # 2. PULIZIA TOTALE DIMNAMES (Risolve il bug di allineamento C++)
    dimnames(X_sub) <- list(NULL, paste0("V", 1:ncol(X_sub)))
    X_sparse <- Matrix(X_sub, sparse = TRUE)
    
    # Stratificazione dei fold
    pos_idx <- which(y == 1)
    neg_idx <- which(y == 0)
    
    foldid <- numeric(length(y))
    foldid[pos_idx] <- sample(rep(1:nfolds, length.out = length(pos_idx)))
    foldid[neg_idx] <- sample(rep(1:nfolds, length.out = length(neg_idx)))
    
    fit_done <- FALSE
    coefs_sub <- NULL
    
    # 3. TENTATIVO CON CROSS-VALIDATION
    tryCatch({
      fit <- cv.glmnet(
        x = X_sparse,
        y = y,
        family = family,
        penalty.factor = penalty_sub,
        foldid = foldid,
        lambda = lambda_grid,                     
        type.measure = "mse",
        alpha = 1,
        maxit = 100000, parallel = TRUE
      )
      models[[i]] <- fit
      coefs_sub <- as.matrix(coef(fit, s = "lambda.1se"))[-1, 1]
      fit_done <- TRUE
    }, error = function(e) {
      
      # FALLBACK SICURO: Se la CV fallisce sul neurone i, usa glmnet semplice
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
        
        # Estraiamo i coefficienti a metà della griglia lambda
        mid_idx <- round(length(lambda_grid) / 2)
        coefs_sub <<- as.matrix(coef(fit_fallback))[-1, mid_idx]
        fit_done <<- TRUE
      }, error = function(e2) {
        message(sprintf("Neuron %d CV avoided for numeric instability.", i))
      })
    })
    
    # Ricostruzione della matrice W
    if (fit_done && !is.null(coefs_sub)) {
      full_coefs <- numeric(ncol(X_full))
      full_coefs[valid_cols] <- coefs_sub
      
      coef_lags <- full_coefs[(n_spatial + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag_frames)
      
      W_matrix[, i] <- rowSums(coef_mat)
    }
  }
  
  return(list(
    models = models,
    adj_matrix = W_matrix,
    prep_data = prep_data,
    active_neurons = active_neurons,
    valid_cols = valid_cols
  ))
}

# 1. Preparazione Dati per Trial 1 (Train)
prep_t1 <- prepare_design_matrix(all_trials$M3424F$trial1, max_lag_frames = 3)

# 2. Fit del Modello
fit_t1 <- fit_model(prep_t1, family = "binomial", min_spikes = 15, nfolds = 5)


# ==============================================================================
# FASE 3: Funzione di Previsione (Forecasting su Nuovi Trial / Test Set)
# ==============================================================================
predict_model <- function(fit_results, test_trial_data) {
  
  max_lag <- fit_results$prep_data$max_lag_frames
  valid_cols <- fit_results$valid_cols # Recupera le colonne valide del train
  
  # 1. Prepara la matrice grezza del nuovo trial di test
  prep_test <- prepare_design_matrix(
    test_trial_data, 
    max_lag_frames = max_lag
  )
  
  # 2. FILTRA LA MATRICE DI TEST SULLE STESSE COLONNE DEL TRAIN (1117 colonne)
  X_test_sub <- prep_test$X[, valid_cols, drop = FALSE]
  dimnames(X_test_sub) <- list(NULL, paste0("V", 1:ncol(X_test_sub)))
  X_test_sparse <- Matrix(X_test_sub, sparse = TRUE)
  
  K <- prep_test$K
  T_test <- prep_test$T_frames
  
  prob_matrix <- matrix(0, nrow = K, ncol = T_test)
  
  # 3. Predizione neurone per neurone
  for (i in fit_results$active_neurons) {
    fit_i <- fit_results$models[[i]]
    if (!is.null(fit_i)) {
      if (inherits(fit_i, "cv.glmnet")) {
        probs <- predict(fit_i, newx = X_test_sparse, s = "lambda.1se", type = "response")
      } else {
        # Gestione modello di fallback
        mid_lambda <- fit_i$lambda[round(length(fit_i$lambda)/2)]
        probs <- predict(fit_i, newx = X_test_sparse, s = mid_lambda, type = "response")
      }
      prob_matrix[i, ] <- as.vector(probs)
    }
  }
  
  return(list(
    predicted_probs = prob_matrix,  # [K x T_test] Probabilità predette
    actual_spikes   = prep_test$Z,  # [K x T_test] Activity reale
    time_ms         = prep_test$dt_ms
  ))
}

# 3. Previsione sul Trial 2 (Test)
predictions_t2 <- predict_model(fit_t1, test_trial_data = all_trials$M3412$trial1)

# 1. Preparazione parametri
active_ids <- fit_t1$active_neurons
grid_spec <- seq(0, 1, length.out = 200) # Griglia di campionamento per la media
sens_matrix <- matrix(NA, nrow = length(active_ids), ncol = length(grid_spec))
auc_vector <- c()

first_valid <- TRUE

# 2. Ciclo su tutti i neuroni attivi per disegnare le curve individuali
for (k in seq_along(active_ids)) {
  i <- active_ids[k]
  y_true <- predictions_t2$actual_spikes[i, ]
  y_pred <- predictions_t2$predicted_probs[i, ]
  
  # Verifichiamo che il neurone abbia sparato almeno una volta nel test set
  if (length(unique(y_true)) > 1) {
    roc_i <- roc(y_true, y_pred, quiet = TRUE)
    auc_vector <- c(auc_vector, auc(roc_i))
    
    # Inizializza il grafico con la prima curva, poi aggiungi le successive (add = TRUE)
    if (first_valid) {
      plot(roc_i, col = rgb(0.5, 0.5, 0.5, 0.15), 
           main = sprintf("ROC curve (Active Neurons %d)", length(active_ids)),
           lwd = 1)
      first_valid <- FALSE
    } else {
      plot(roc_i, add = TRUE, col = rgb(0.5, 0.5, 0.5, 0.15), lwd = 1)
    }
    
    # Interpolazione delle sensibilità sulla griglia di specificità comune
    ord <- order(roc_i$specificities)
    sens_interp <- approx(
      x = roc_i$specificities[ord], 
      y = roc_i$sensitivities[ord], 
      xout = grid_spec, 
      rule = 2
    )$y
    
    sens_matrix[k, ] <- sens_interp
  }
}

# 3. Calcolo e sovrapposizione della curva ROC MEDIA (in ROSSO)
mean_sens <- colMeans(sens_matrix, na.rm = TRUE)
lines(x = grid_spec, y = mean_sens, col = "red", lwd = 3.5)

mean_auc <- mean(auc_vector, na.rm = TRUE)



#############################################################
plot_neuron_place_field <- function(fit_results, trial_data, neuron_id, grid_size = 100) {
  
  # --- 1. Allineamento Coordinate e Recupero Basi Spaziali ---
  t_behav <- trial_data$behav$time[, 1]
  pos_behav <- trial_data$behav$position
  
  t_raw <- trial_data$neuron$time[, 1]
  T_frames <- ncol(trial_data$neuron$S)
  t_neuron <- t_raw[seq(1, length(t_raw), length.out = T_frames)]
  
  pos_x <- approx(t_behav, pos_behav[, 1], xout = t_neuron, rule = 2)$y
  pos_y <- approx(t_behav, pos_behav[, 2], xout = t_neuron, rule = 2)$y
  
  n_spatial <- fit_results$prep_data$n_spatial_cols
  n_splines <- sqrt(n_spatial)
  
  # Fit originale delle B-splines sulle traiettorie reali (serve per ereditare i knots)
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
  
  # --- 3. Estrazione Coefficienti per il Neurone ---
  fit_i <- fit_results$models[[neuron_id]]
  if (is.null(fit_i)) {
    stop(sprintf("Il neurone %d non è tra quelli attivi stimati dal modello.", neuron_id))
  }
  
  if (inherits(fit_i, "cv.glmnet")) {
    coefs <- as.matrix(coef(fit_i, s = "lambda.1se"))
  } else {
    mid_idx <- round(length(fit_i$lambda) / 2)
    coefs <- as.matrix(coef(fit_i))[, mid_idx, drop = FALSE]
  }
  
  # L'intercetta è alla riga 1, i coefficienti spaziali occupano le righe da 2 a (n_spatial + 1)
  beta_0 <- coefs[1, 1]
  beta_spatial <- coefs[2:(n_spatial + 1), 1]
  
  # --- 4. Calcolo Probabilità Baseline P(Spike | X, Y) ---
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
  
  # --- 5. Grafico (Mappa di Calore + Traiettoria Topo) ---
  filled.contour(
    x = x_grid, 
    y = y_grid, 
    z = prob_matrix,
    color.palette = soft_palette,
    xlab = "Position X", 
    ylab = "Position Y",
    plot.axes = {
      axis(1)
      axis(2)
      # Sovrapponi la traiettoria reale del topo in nero semitrasparente
      lines(pos_x, pos_y, col = rgb(0, 0, 0, 0.3), lwd = 0.5)
    }
  )
}

neuron_da_mostrare <- fit_t1$active_neurons[18]

plot_neuron_place_field(
  fit_results = fit_t1, 
  trial_data = all_trials$M3412$trial1, 
  neuron_id = neuron_da_mostrare,
  grid_size = 120
)


################################################################
# Heatmap autoregressive coefficient

plot_single_lag <- function(fit_results, lag = 1, only_active = TRUE) {
  
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  n_spatial <- fit_results$prep_data$n_spatial_cols
  active_neurons <- fit_results$active_neurons
  valid_cols <- fit_results$valid_cols
  dt_ms <- fit_results$prep_data$dt_ms
  
  if (lag < 1 || lag > max_lag) {
    stop(sprintf("Il lag %d non è valido! Deve essere un numero compreso tra 1 e %d.", lag, max_lag))
  }
  
  idx_plot <- if (only_active) active_neurons else 1:K
  N_plot <- length(idx_plot)
  
  # 1. Matrice numerica [K x K] completa per il lag selezionato
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
      
      coef_lags <- full_coefs[(n_spatial + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag)
      
      # Salviamo la colonna del lag di interesse
      W_mat_full[, i] <- coef_mat[, lag]
    }
  }
  
  # Subsetting sui soli neuroni richiesti
  W_sub <- W_mat_full[idx_plot, idx_plot, drop = FALSE]
  
  # 2. VETTORIZZAZIONE: Creazione diretta del data.frame (senza cicli for)
  df_single_lag <- data.frame(
    Source = factor(rep(idx_plot, times = N_plot), levels = rev(idx_plot)),
    Target = factor(rep(idx_plot, each = N_plot), levels = idx_plot),
    Weight = as.vector(W_sub)
  )
  
  lag_ms <- round(lag * dt_ms)
  
  # 3. Plotting con ggplot2
  p <- ggplot(df_single_lag, aes(x = Target, y = Source, fill = Weight)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#F8F9F9", "#CBD5E1", "#F5B041", "#E74C3C", "#78281F"),
      name = "Coeff gamma"
    ) +
    labs(
      x = "Target neuron",
      y = "Input neuron"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
  
  return(p)
}

# Heatmap solo per il Lag 1 
p_lag1 <- plot_single_lag(fit_t1, lag = 1, only_active = T)
print(p_lag1)

# Heatmap solo per il Lag 2 
p_lag2 <- plot_single_lag(fit_t1, lag = 2, only_active = T)
print(p_lag2)

# Heatmap solo per il Lag 3 
p_lag3 <- plot_single_lag(fit_t1, lag = 3, only_active = TRUE)
print(p_lag3)



plot_spatial_from_lag_heatmap <- function(fit_results, lag = 1, target_id = NULL, min_weight = 0) {
  
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  n_spatial <- fit_results$prep_data$n_spatial_cols
  active_neurons <- fit_results$active_neurons
  valid_cols <- fit_results$valid_cols
  dt_ms <- fit_results$prep_data$dt_ms
  centroids <- fit_results$prep_data$centroids
  
  if (lag < 1 || lag > max_lag) {
    stop(sprintf("Lag non valido! Scegli un valore compreso tra 1 e %d.", max_lag))
  }
  
  # --- 1. Estrazione Matrice W esatta per il solo LAG selezionato ---
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
      
      coef_lags <- full_coefs[(n_spatial + 1):length(full_coefs)]
      coef_mat <- matrix(coef_lags, nrow = K, ncol = max_lag)
      
      W_mat_lag[, i] <- coef_mat[, lag]
    }
  }
  diag(W_mat_lag) <- 0  # Rimuoviamo le autoconnessioni
  
  # --- 2. Costruzione della lista dei collegamenti non nulli per questo Lag ---
  edges_list <- list()
  for (j in 1:K) {
    for (i in 1:K) {
      w <- W_mat_lag[j, i]
      if (w > min_weight) {
        edges_list[[length(edges_list) + 1]] <- data.frame(
          Source = j, Target = i,
          x_src = centroids[j, 1], y_src = centroids[j, 2],
          x_tgt = centroids[i, 1], y_tgt = centroids[i, 2],
          Weight = w
        )
      }
    }
  }
  
  df_edges <- if (length(edges_list) > 0) bind_rows(edges_list) else data.frame()
  
  # Filtro opzionale se si vuole isolare un singolo neurone target
  if (!is.null(target_id) && nrow(df_edges) > 0) {
    df_edges <- df_edges %>% filter(Source == target_id | Target == target_id)
  }
  
  df_nodes <- data.frame(
    Neuron = 1:K,
    X = centroids[, 1],
    Y = centroids[, 2]
  )
  
  lag_ms <- round(lag * dt_ms)
  
  # --- 3. Rendering Grafico Spaziale ---
  p <- ggplot() +
    # Sfondo: Centroidi di tutti i neuroni
    geom_point(data = df_nodes, aes(x = X, y = Y), 
               color = "#E0E0E0", size = 1.2, alpha = 0.6) +
    
    # Frecce orientate per i collegamenti definiti nella heatmap
    {
      if (nrow(df_edges) > 0) {
        geom_segment(
          data = df_edges,
          aes(x = x_src, y = y_src, xend = x_tgt, yend = y_tgt, 
              color = Weight, alpha = Weight),
          arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
          linewidth = 0.6
        )
      }
    } +
    
    # Evidenziazione dei neuroni coinvolti nei collegamenti
    {
      if (nrow(df_edges) > 0) {
        connected_ids <- unique(c(df_edges$Source, df_edges$Target))
        df_connected <- df_nodes %>% filter(Neuron %in% connected_ids)
        
        geom_point(data = df_connected, aes(x = X, y = Y), 
                   color = "#2C3E50", size = 1.8)
      }
    } +
    
    # Evidenziazione del neurone target se specificato
    {
      if (!is.null(target_id)) {
        df_target <- df_nodes %>% filter(Neuron == target_id)
        geom_point(data = df_target, aes(x = X, y = Y), 
                   color = "black", fill = "#00FF66", shape = 21, size = 3, stroke = 1.2)
      }
    } +
    
    scale_color_gradientn(
      colors = c("#F5B041", "#E74C3C", "#78281F"),
      name = "Coeff. gamma"
    ) +
    scale_alpha_continuous(range = c(0.45, 0.95), guide = "none") +
    coord_fixed() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      panel.grid = element_line(color = "gray95")
    )
  
  return(p)
}

p_net_lag1 <- plot_spatial_from_lag_heatmap(fit_t1, lag = 1)
print(p_net_lag1)

p_net_lag2 <- plot_spatial_from_lag_heatmap(fit_t1, lag = 2)
print(p_net_lag2)

p_net_lag3 <- plot_spatial_from_lag_heatmap(fit_t1, lag = 3)
print(p_net_lag3)


library(ggplot2)
library(dplyr)

plot_spatial_multi_lag <- function(fit_results, lags = c(1, 2, 3), target_id = NULL, min_weight = 0) {
  
  K <- fit_results$prep_data$K
  max_lag <- fit_results$prep_data$max_lag_frames
  n_spatial <- fit_results$prep_data$n_spatial_cols
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
        
        coef_hawkes <- full_coefs[(n_spatial + 1):length(full_coefs)]
        coef_mat <- matrix(coef_hawkes, nrow = K, ncol = max_lag)
        
        W_mat_lag[, i] <- coef_mat[, lag]
      }
    }
    
    # 1A. Estrazione Autoconnessioni (gamma_ii sulla diagonale)
    gamma_ii <- diag(W_mat_lag)
    self_ids <- which(gamma_ii > min_weight)
    
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
        if (w > min_weight) {
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
    # Layer 1: Sfondo nodi dell'ippocampo (grigio chiaro)
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
    
    # Layer 3: Nodi coinvolti in connessioni inter-neurone (blu scuro)
    {
      if (nrow(df_connected) > 0) {
        geom_point(data = df_connected, aes(x = X, y = Y), 
                   color = "#2C3E50", size = 1.6, alpha = 0.9)
      }
    } +
    
    # Layer 4: NEURONI CON AUTO-CONNESSIONE (gamma_ii > 0) IN ROSSO
    {
      if (nrow(df_self) > 0) {
        geom_point(data = df_self, aes(x = X, y = Y), 
                   color = "#78281F", size = 1.6, alpha = 0.9)
      }
    } +
    
    # Layer 5: Neurone target evidenziato (se specificato)
    {
      if (!is.null(target_id)) {
        df_target <- df_nodes %>% filter(Neuron == target_id)
        geom_point(data = df_target, aes(x = X, y = Y), 
                   color = "black", fill = "#00FF66", shape = 21, size = 3.5, stroke = 1.2)
      }
    } +
    
    facet_wrap(~ Lag, nrow = 1) +
    
    scale_color_gradientn(
      colors = c("#F5B041", "#E74C3C", "#78281F"),
      name = "Coeff. gamma_ij"
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

# Disegna i 3 lag affiancati in una sola figura
p_multi <- plot_spatial_multi_lag(fit_t1, lags = c(1, 2, 3))
print(p_multi)



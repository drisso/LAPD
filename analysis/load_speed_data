# ==============================================================================
# SCRIPT DI IMPORTAZIONE GENERICO DATI MOUSE
# ==============================================================================

load_mouse_speed_data <- function(mouse_id = "M3424F", 
                                  file_prefix = "mouse", 
                                  base_dir = "~/Downloads", 
                                  n_trials = 6,
                                  target_cols = c("speed_spline", "dist_to_hull", "time")) {
  
  speed_list  <- vector("list", n_trials)
  trial_names <- paste0("trial", 1:n_trials)
  names(speed_list) <- trial_names
  
  for (i in 1:n_trials) {
    # Costruzione dinamica del nome file (es. mouse_trial1.RData o M3424F_trial1.RData)
    file_name <- sprintf("%s_trial%d.RData", file_prefix, i)
    file_path <- file.path(base_dir, file_name)
    
    if (file.exists(file_path)) {
      temp_env <- new.env()
      load(file_path, envir = temp_env)
      
      if ("df.neuron" %in% ls(temp_env)) {
        df_temp <- temp_env$df.neuron
        
        missing_cols <- setdiff(target_cols, colnames(df_temp))
        if (length(missing_cols) == 0) {
          speed_list[[trial_names[i]]] <- df_temp[, target_cols]
        } else {
          warning(sprintf("[%s] Trial %d: Colonne mancanti in df.neuron: %s", 
                          mouse_id, i, paste(missing_cols, collapse = ", ")))
        }
      } else {
        warning(sprintf("[%s] Trial %d: Oggetto 'df.neuron' non trovato in %s", 
                        mouse_id, i, file_path))
      }
    } else {
      warning(sprintf("[%s] File non trovato: %s", mouse_id, file_path))
    }
  }
  
  return(speed_list)
}

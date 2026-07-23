# ==============================================================================
# HURDLE MODEL - MAIN
# ==============================================================================
# Modello hurdle in due step per la serie di ampiezze neurali S_k(t), nulla
# quasi ovunque e log-normale sui valori positivi:
#
#   Step 1 (evento)    P(S_k(t) > 0 | X)              -> logistica LASSO
#   Step 2 (ampiezza)  log(S_k(t)) | S_k(t) > 0       -> gaussiana LASSO
#   Previsione hurdle  E[S_k(t) | X] = P(S>0|X) * E[S|S>0,X]
#
# Questo file carica i dati, esegue il fit e produce le previsioni. Le funzioni
# stanno in Step_1.R e Step_2.R. Nulla viene scritto su disco: i risultati
# restano oggetti nella sessione R.
#
# Due modalita' d'uso, controllate da FIT_STEP1:
#   FIT_STEP1 = TRUE  -> stima entrambi gli step da zero
#   FIT_STEP1 = FALSE -> carica fit_t1 da un .rda gia' salvato e stima solo lo Step 2
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Percorsi e impostazioni
# ------------------------------------------------------------------------------
PROJECT_DIR <- "C:/Users/pietr/OneDrive/Desktop/PhD/Activities/Data Research Camp 2026/Project"
HURDLE_DIR  <- file.path(PROJECT_DIR, "Hurdle")
DATA_DIR    <- file.path(PROJECT_DIR, "Data")

ALL_TRIALS_FILE  <- file.path(DATA_DIR, "all_trials.RDS")

# Trial su cui viene stimato il modello
MOUSE_ID <- "M3424F"
TRIAL_ID <- "trial1"

# --- Sorgente della velocita' ---------------------------------------------
# "per_trial"   : load_mouse_speed_data() restituisce una lista per trial. E' la
#                 modalita' corretta, perche' ogni trial ha la propria serie di
#                 velocita' e la previsione out-of-sample deve usare quella del
#                 trial di TEST, non quella di training.
# "single_file" : un unico data.frame da mouse_trial.RData (df.neuron). Da usare
#                 solo se load_speed_data.R non e' disponibile; in quel caso la
#                 previsione out-of-sample usa una velocita' non allineata al
#                 trial di test, quindi va considerata inaffidabile.
SPEED_SOURCE <- "per_trial"

# Usati se SPEED_SOURCE = "per_trial"
LOAD_SPEED_SCRIPT <- file.path(HURDLE_DIR, "load_speed_data.R")
SPEED_BASE_DIR    <- DATA_DIR
SPEED_N_TRIALS    <- 6L

# ATTENZIONE: in load_mouse_speed_data() i nomi file dipendono da file_prefix,
# non da mouse_id (che compare solo nei messaggi di warning). Quindi il prefisso
# va mappato esplicitamente sul topo, altrimenti train e test di topi diversi
# leggerebbero gli stessi file.
SPEED_FILE_PREFIX <- c(M3424F = "mouse")

# Usato se SPEED_SOURCE = "single_file"
MOUSE_TRIAL_FILE <- file.path(DATA_DIR, "mouse_trial.RData")

# Step 1: stimarlo (TRUE) oppure caricarlo da un .rda gia' salvato (FALSE)
FIT_STEP1        <- FALSE
STEP1_FIT_FILE   <- file.path(DATA_DIR, "fit_M3424F trial 1_results.rda")   # usato solo se FIT_STEP1 = FALSE

# Iperparametri Step 1
MAX_LAG_FRAMES <- 5      # numero di lag temporali passati (oltre al lag 0)
N_SPLINES      <- 4      # gradi di liberta' delle B-spline spaziali per asse
MIN_SPIKES     <- 15     # neuroni considerati attivi
NFOLDS         <- 5      # blocchi di cross-validation temporale
BLOCK_DUR_SEC  <- 2.0    # durata in secondi di ogni blocco temporale

# Iperparametri Step 2
POSITIVE_FAMILY <- "lognormal"
MIN_POSITIVE    <- 15L   # eventi positivi minimi per stimare un neurone
PARALLEL_CV     <- TRUE  # effettivo solo con un backend foreach registrato

# Percorso di lambda dello Step 2. NULL = glmnet calcola lambda_max dai dati di
# ciascun neurone: la scala di lambda dello stage gaussiano su log(S) non ha
# nulla a che vedere con quella della devianza binomiale dello Step 1, quindi
# ereditarne la griglia porta la 1se a scegliere il modello nullo.
# Per tornare al comportamento precedente: extract_stage1_lambda_grid(fit_t1).
STEP2_LAMBDA_GRID <- NULL
STEP2_NLAMBDA     <- 100L

# Penalizzazione delle covariate fisse nello Step 2. Le 16 basi spline prodotto
# sono collineari e verrebbero stimate su poche centinaia di frame positivi:
# lasciarle libere fa esplodere exp(eta) fuori dal supporto di stima.
PENALIZE_SPATIAL_BASIS <- TRUE
SPATIAL_PENALTY_FACTOR <- 1
UNPENALIZED_FIXED      <- "speed_spline"

# Combinazione dei due stage nella previsione finale:
#   "soft_raw" (attesa hurdle canonica) | "hard_threshold" | "soft_logit_calibrated"
JOINT_PREDICTION_METHOD <- "soft_raw"
HARD_THRESHOLD_RULE     <- "f1"

# Previsione out-of-sample: trial di test (NULL per saltarla)
TEST_MOUSE_ID <- "M3424F"
TEST_TRIAL_ID <- "trial6"


# ------------------------------------------------------------------------------
# 1. Funzioni
# ------------------------------------------------------------------------------
source(file.path(HURDLE_DIR, "Step_1.R"))
source(file.path(HURDLE_DIR, "Step_2.R"))

# Carica fit_t1 (ed eventuali previsioni Step 1) da un .rda salvato in precedenza
load_stage1_bundle <- function(path) {
  if (!file.exists(path)) stop("Step 1 file not found: ", path)
  if (!grepl("\\.(rda|RData)$", path, ignore.case = TRUE)) {
    warning("The Step 1 file does not have an .rda/.RData extension: ", path)
  }

  env <- new.env(parent = globalenv())
  loaded_names <- load(path, envir = env)

  if (!"fit_t1" %in% loaded_names) {
    stop(
      "fit_t1 was not found in the file. Available objects: ",
      paste(loaded_names, collapse = ", ")
    )
  }

  list(
    fit_t1 = env$fit_t1,
    predictions_t1 = if ("predictions_t2" %in% loaded_names) env$predictions_t2 else NULL,
    prep_t1 = if ("prep_t1" %in% loaded_names) env$prep_t1 else NULL,
    loaded_names = loaded_names,
    source_file = normalizePath(path, mustWork = TRUE)
  )
}


# ------------------------------------------------------------------------------
# 2. Dati
# ------------------------------------------------------------------------------
all_trials <- readRDS(ALL_TRIALS_FILE)

# Restituisce la serie di velocita' del trial richiesto.
get_speed_data <- function(mouse_id, trial_id) {
  if (identical(SPEED_SOURCE, "per_trial")) {
    prefix <- SPEED_FILE_PREFIX[[mouse_id]]
    if (is.null(prefix) || is.na(prefix)) {
      stop(
        "Nessun prefisso file di velocita' definito per il topo ", mouse_id,
        ". Aggiungilo a SPEED_FILE_PREFIX."
      )
    }

    mouse_speed <- load_mouse_speed_data(
      mouse_id    = mouse_id,
      file_prefix = prefix,
      base_dir    = SPEED_BASE_DIR,
      n_trials    = SPEED_N_TRIALS
    )
    speed <- mouse_speed[[trial_id]]
    if (is.null(speed)) {
      stop("Speed data not found for ", mouse_id, " - ", trial_id, ".")
    }
    return(speed)
  }

  # single_file: unica serie condivisa da tutti i trial
  env <- new.env()
  load(MOUSE_TRIAL_FILE, envir = env)
  env$df.neuron[, c("speed_spline", "time")]
}

if (identical(SPEED_SOURCE, "per_trial")) {
  if (!file.exists(path.expand(LOAD_SPEED_SCRIPT))) {
    stop(
      "LOAD_SPEED_SCRIPT non trovato: ", LOAD_SPEED_SCRIPT,
      ". Correggi il percorso oppure imposta SPEED_SOURCE <- \"single_file\"."
    )
  }
  source(LOAD_SPEED_SCRIPT)
}

experiment <- all_trials[[MOUSE_ID]][[TRIAL_ID]]
if (is.null(experiment)) {
  stop("Trial not found in all_trials: ", MOUSE_ID, " - ", TRIAL_ID, ".")
}
S_train <- as.matrix(experiment$neuron$S)

# Velocita' del trial di stima
data_w_dist <- get_speed_data(MOUSE_ID, TRIAL_ID)

message("Trial: ", MOUSE_ID, " - ", TRIAL_ID,
        " | neuroni: ", nrow(S_train), " | frame: ", ncol(S_train))


# ------------------------------------------------------------------------------
# 3. STEP 1 - stage evento
# ------------------------------------------------------------------------------
if (FIT_STEP1) {

  prep_t1 <- prepare_design_matrix(
    trial_data     = experiment,
    data_w_dist    = data_w_dist,
    max_lag_frames = MAX_LAG_FRAMES,
    n_splines      = N_SPLINES
  )

  fit_t1 <- fit_model(
    prep_data     = prep_t1,
    family        = "binomial",
    min_spikes    = MIN_SPIKES,
    nfolds        = NFOLDS,
    block_dur_sec = BLOCK_DUR_SEC
  )

  # Previsioni Step 1 in-sample (probabilita' di evento)
  predictions_t1 <- predict_model(
    fit_results     = fit_t1,
    test_trial_data = experiment,
    data_w_dist     = data_w_dist
  )

  stage1_source <- "fitted in this session"

} else {

  stage1_bundle  <- load_stage1_bundle(STEP1_FIT_FILE)
  fit_t1         <- stage1_bundle$fit_t1
  predictions_t1 <- stage1_bundle$predictions_t1
  prep_t1        <- stage1_bundle$prep_t1
  stage1_source  <- stage1_bundle$source_file

  message("Step 1 caricato da: ", stage1_source)
  message("Oggetti nel file: ", paste(stage1_bundle$loaded_names, collapse = ", "))
}

validate_stage1_structure(fit_t1)

# Verifica che il trial caricato sia esattamente quello dello Step 1
S_train <- check_S_matches_stage1(S_train, fit_t1)

message(
  "Setup Step 1: n_fixed=", fit_t1$prep_data$n_fixed_cols,
  "; lag di rete=0..", fit_t1$prep_data$max_lag_frames,
  "; fold temporali=", length(unique(fit_t1$foldid)),
  "; neuroni attivi=", length(fit_t1$active_neurons), "."
)


# ------------------------------------------------------------------------------
# 4. STEP 2 - stage ampiezza positiva (log-normale)
# ------------------------------------------------------------------------------
use_parallel <- resolve_parallel_cv(PARALLEL_CV)

fit_t2 <- fit_positive_amplitude_stage(
  first_fit              = fit_t1,
  S                      = S_train,
  positive_family        = POSITIVE_FAMILY,
  min_positive           = MIN_POSITIVE,
  lambda_grid            = STEP2_LAMBDA_GRID,
  nlambda                = STEP2_NLAMBDA,
  penalize_spatial_basis = PENALIZE_SPATIAL_BASIS,
  spatial_penalty_factor = SPATIAL_PENALTY_FACTOR,
  unpenalized_fixed      = UNPENALIZED_FIXED,
  parallel               = use_parallel
)

print(table(fit_t2$status))
print(fit_t2$lambda_diagnostics)


# ------------------------------------------------------------------------------
# 5. Modello hurdle completo
# ------------------------------------------------------------------------------
hurdle_model <- structure(
  list(
    stage1_event    = fit_t1,
    stage2_positive = fit_t2,
    mouse           = MOUSE_ID,
    trial           = TRIAL_ID,
    fixed_covariates = colnames(fit_t1$prep_data$X)[seq_len(fit_t1$prep_data$n_fixed_cols)],
    stage1_source   = stage1_source,
    lambda_grid     = fit_t2$lambda_grid,
    fixed_penalty   = fit_t2$fixed_penalty,
    response        = "neuron$S",
    fitted_on       = Sys.time()
  ),
  class = "neural_hurdle_model"
)


# ------------------------------------------------------------------------------
# 6. Previsione IN-SAMPLE
# ------------------------------------------------------------------------------
hurdle_in <- predict_hurdle(
  first_fit                = fit_t1,
  second_fit               = fit_t2,
  S                        = S_train,
  saved_stage1_predictions = predictions_t1,
  method                   = JOINT_PREDICTION_METHOD,
  threshold_rule           = HARD_THRESHOLD_RULE
)

hurdle_predictions_in <- hurdle_in$selected
hurdle_evaluation_in  <- evaluate_hurdle_predictions(
  hurdle_predictions_in,
  neuron_ids = fit_t2$positive_neurons
)

message("Metodo di combinazione selezionato: ", hurdle_in$selected_method)
print(hurdle_evaluation_in$aggregate)

# Confronto fra i tre metodi di combinazione (nessun rifit)
hurdle_comparison_in <- compare_joint_prediction_methods(hurdle_in)
print(hurdle_comparison_in$comparison)


# ------------------------------------------------------------------------------
# 7. Previsione OUT-OF-SAMPLE (opzionale)
# ------------------------------------------------------------------------------
# La design matrix del trial di test viene ricostruita con gli stessi
# max_lag_frames e n_splines dello Step 1. Soglie e calibrazione NON vengono
# ristimate sui dati di test: si applicano quelle dell'in-sample.
if (!is.null(TEST_MOUSE_ID) && !is.null(TEST_TRIAL_ID)) {

  experiment_test <- all_trials[[TEST_MOUSE_ID]][[TEST_TRIAL_ID]]
  if (is.null(experiment_test)) {
    stop("Test trial not found in all_trials: ", TEST_MOUSE_ID, " - ", TEST_TRIAL_ID, ".")
  }

  # Velocita' del trial di TEST, non quella di training.
  data_w_dist_test <- get_speed_data(TEST_MOUSE_ID, TEST_TRIAL_ID)

  hurdle_out <- predict_hurdle(
    first_fit         = fit_t1,
    second_fit        = fit_t2,
    trial_data        = experiment_test,
    data_w_dist       = data_w_dist_test,
    method            = JOINT_PREDICTION_METHOD,
    fitted_components = hurdle_in$components
  )

  hurdle_predictions_out <- hurdle_out$selected
  hurdle_evaluation_out  <- evaluate_hurdle_predictions(
    hurdle_predictions_out,
    neuron_ids = fit_t2$positive_neurons
  )

  message("Previsione out-of-sample: ", TEST_MOUSE_ID, " - ", TEST_TRIAL_ID)
  print(hurdle_evaluation_out$aggregate)
}

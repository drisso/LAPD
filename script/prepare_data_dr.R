library(tidyverse)
source("script/data_preparation.R")

project_root <- "."

trial_data <- load_trial_data(
  rds_path = file.path(project_root, "all_trials.RDS"),
  mouse = 2L,
  trial = 1L,
  cache_path = NULL
)

prepared <- make_frame_covariates(trial_data)
S <- as.data.frame(t(prepared$S))
colnames(S) <- paste0("n", seq_len(ncol(S)))
wide <- cbind(prepared$frames, S)
long <- pivot_longer(wide, names_to = "neuron", values_to = "intensity", cols=grep("^n", colnames(wide)))
long <- mutate(long, activation = ifelse(intensity>0, 1, 0))

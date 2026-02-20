# ============================================================================
# script4_bigRegression_FUNCTIONS.R
# Helper functions for the regression pipeline (Step 4). Provides data
# preparation, column derivation, external predictor joining, and model
# fitting routines for both logistic and linear models.
#
# Inputs: haploDF data frames (passed by caller)
# Outputs: prepared data frames, fitted model objects
# Dependencies: functions_pick_strains_functions.R, plot_functions.R
#               (sourced by caller)
# ============================================================================

############################################################
# Regression pipeline: detected (binomial) + delay (lm)
# Tidied + streamlined
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(MMWRweek)
  library(car)
  library(broom)
  library(ggplot2)
  library(patchwork)
})

# ---- user-supplied helper functions ----
# must provide these:
# - summarise_entropy_by_group()
# - assert_cumcount_matches()
# - summarise_detection_delay()
# - scale_columns()  (ensure it does NOT scale detected; best to only scale vars you pass)
# - log_reg_summary_table()
# - glm_summary_table()
# - forestPlot()          (for logistic)
# - forest_plot_lm()      (for linear)

# source(".../plot_functions.R")
# source(".../functions_pick_strains_functions.R")


# -----------------------------
# Helpers
# -----------------------------

ensure_dirs <- function(paths) {
  walk(paths, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))
  invisible(TRUE)
}

parse_subtype_from_filename <- function(x) {
  # returns H1 / H3
  str_match(x, "_(H\\d)_")[,2]
}

parse_threshold_from_filename <- function(x) {
  # matches: haploDF_noMich_H1_0.005_2021_v2.csv
  # returns: 0.005
  str_match(x, "^haploDF_noMich_H\\d_([^_]+)_\\d{4}")[,2]
}

add_derived_columns <- function(df) {
  df %>%
    mutate(
      season    = USseason_dominant,
      SEASON    = USseason_dominant,
      source_file = as.character(source_file),
      subtype   = parse_subtype_from_filename(source_file),
      threshold = parse_threshold_from_filename(source_file),
      haplo_length = nchar(Haplotype),
      haplo_days_to_extinction = as.numeric(last_haplo - firstUSDate),
      detected = if_else(!is.na(delay), 1L, 0L),
      WEEK = MMWRweek(firstUSDate)$MMWRweek
    )
}

drop_problem_rows <- function(df) {
  # you explicitly dropped H1 2021 + rows where firstUSDate is NA
  df %>%
    filter(!(subtype == "H1" & season == "2021")) %>%
    filter(!is.na(firstUSDate))
}

join_entropy <- function(df) {
  res <- summarise_entropy_by_group(df)
  
  # sanity checks you already do
  if (any(!res$cumcount_equals_sumcount, na.rm = TRUE)) {
    message("WARNING: cumcount_equals_sumcount has FALSE rows. Printing a few:")
    print(res %>% filter(!cumcount_equals_sumcount) %>% head(20))
  }
  assert_cumcount_matches(res)
  
  df %>%
    left_join(res %>% select(season, subtype, threshold, entropy),
              by = c("season", "subtype", "threshold"))
}

join_external_predictors <- function(df,
                                     episize_path = "../episize_stats.csv",
                                     state_metrics_path = "../state_metrics.csv",
                                     state_seq_path = "../state_seq_numbers.csv") {
  
  episize_stats <- read.csv(episize_path)
  state_metrics <- read.csv(state_metrics_path)
  state_seq     <- read.csv(state_seq_path)
  
  df %>%
    left_join(episize_stats,  by = c("SEASON", "WEEK")) %>%
    left_join(state_metrics,  by = c("firstState" = "state")) %>%
    left_join(state_seq,      by = c("season", "firstState")) %>%
    mutate(
      pop_size = pop_size / 100000,
      dist_to_Mich = dist_to_Mich / 1000
    )
}

fit_detected_models <- function(df_scaled, threshold,
                                out_dir_models, out_dir_figs,
                                predictors,
                                ref_state = "Illinois") {
  
  # Drop lines where first detection is Michigan (your intent)
  df_scaled <- df_scaled %>% filter(firstState != "Michigan")
  
  # keep only modeled columns + complete cases
  model_df <- df_scaled %>%
    select(all_of(c("detected", predictors))) %>%
    filter(complete.cases(.))
  
  # Base “check” model for alias/vif (optional)
  check_model <- glm(detected ~ ., data = model_df, family = binomial())
  # alias(check_model); vif(check_model)  # keep if you want console output
  
  # Example rule you used: drop FLUWEEK due to correlation with EPIDEMIC_PROP
  if ("FLUWEEK" %in% colnames(model_df)) {
    model_df <- model_df %>% select(-FLUWEEK)
  }
  
  final_modelB1 <- glm(detected ~ ., data = model_df, family = binomial())
  
  tab_B1 <- log_reg_summary_table(final_modelB1, value = "z") %>%
    mutate(threshold = threshold, model = "DETECTED_final_modelB1")
  
  write.csv(
    tab_B1,
    file = file.path(out_dir_models, paste0("detected_TF_modelB1_threshold_", threshold, ".csv")),
    row.names = FALSE
  )
  saveRDS(final_modelB1, file = file.path(out_dir_models, paste0("final_modelB1_t", threshold, ".rds")))
  
  gB1 <- forestPlot(final_modelB1)
  ggsave(
    filename = file.path(out_dir_figs, paste0("forest_detected_modelB1_t", threshold, ".jpg")),
    plot = gB1,
    width = 10, height = 6
  )
  
  # State-only model (effect of state)
  state_df <- df_scaled %>%
    select(detected, firstState) %>%
    mutate(firstState = relevel(factor(firstState), ref = ref_state)) %>%
    filter(complete.cases(.))
  
  final_modelB4 <- glm(detected ~ firstState, data = state_df, family = binomial())
  
  tab_B4 <- log_reg_summary_table(final_modelB4, value = "z") %>%
    mutate(threshold = threshold, model = "DETECTED_state_only")
  
  write.csv(
    tab_B4,
    file = file.path(out_dir_models, paste0("detected_TF_stateOnly_threshold_", threshold, ".csv")),
    row.names = FALSE
  )
  saveRDS(final_modelB4, file = file.path(out_dir_models, paste0("final_modelB4_stateOnly_t", threshold, ".rds")))
  
  gB4 <- forestPlot(final_modelB4)
  ggsave(
    filename = file.path(out_dir_figs, paste0("forest_detected_stateOnly_t", threshold, ".jpg")),
    plot = gB4,
    width = 10, height = 6
  )
  
  list(
    modelB1 = final_modelB1,
    modelB4 = final_modelB4,
    plotB1 = gB1,
    plotB4 = gB4,
    tabB1 = tab_B1,
    tabB4 = tab_B4
  )
}

fit_delay_models <- function(df_scaled, threshold,
                             out_dir_models, out_dir_figs,
                             predictors) {
  
  # Only those with observed delay
  model_df <- df_scaled %>%
    filter(!is.na(delay)) %>%
    filter(delay > 0) %>%                 # you had this
    select(all_of(c("delay", predictors))) %>%
    filter(complete.cases(.))
  
  if (nrow(model_df) < 10) {
    message("WARNING: very few rows for delay model at threshold=", threshold,
            " (n=", nrow(model_df), "). Skipping.")
    return(NULL)
  }
  
  # Optional collinearity checks
  check_model <- glm(delay ~ ., data = model_df)
  # alias(check_model); vif(check_model)
  
  # Your earlier logic: remove frequency_order + FLUWEEK if present
  drop_vars <- intersect(c("frequency_order", "FLUWEEK"), colnames(model_df))
  if (length(drop_vars) > 0) model_df2 <- model_df %>% select(-all_of(drop_vars)) else model_df2 <- model_df
  
  final_model1C <- glm(delay ~ ., data = model_df2)
  
  tab_1C <- glm_summary_table(final_model1C) %>%
    mutate(threshold = threshold, model = "DELAY_final_model1C")
  
  write.csv(
    tab_1C,
    file = file.path(out_dir_models, paste0("delay_model1C_threshold_", threshold, ".csv")),
    row.names = FALSE
  )
  saveRDS(final_model1C, file = file.path(out_dir_models, paste0("final_model1C_t", threshold, ".rds")))
  
  gC1 <- forest_plot_lm(final_model1C)
  ggsave(
    filename = file.path(out_dir_figs, paste0("forest_delay_model1C_t", threshold, ".jpg")),
    plot = gC1,
    width = 10, height = 6
  )
  
  # State-only delay model
  state_df <- df_scaled %>%
    filter(!is.na(delay), delay > 0) %>%
    select(delay, firstState) %>%
    filter(complete.cases(.))
  
  final_model3C <- glm(delay ~ firstState, data = state_df)
  
  tab_3C <- glm_summary_table(final_model3C) %>%
    mutate(threshold = threshold, model = "DELAY_state_only")
  
  write.csv(
    tab_3C,
    file = file.path(out_dir_models, paste0("delay_stateOnly_threshold_", threshold, ".csv")),
    row.names = FALSE
  )
  saveRDS(final_model3C, file = file.path(out_dir_models, paste0("final_model3C_stateOnly_t", threshold, ".rds")))
  
  gC3 <- forest_plot_lm(final_model3C)
  ggsave(
    filename = file.path(out_dir_figs, paste0("forest_delay_stateOnly_t", threshold, ".jpg")),
    plot = gC3,
    width = 10, height = 6
  )
  
  list(
    model1C = final_model1C,
    model3C = final_model3C,
    plot1C = gC1,
    plot3C = gC3,
    tab1C = tab_1C,
    tab3C = tab_3C
  )
}


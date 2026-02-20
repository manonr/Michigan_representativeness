# ============================================================================
# script4_bigRegression_EXECUTE.R
# Step 4: Regression analysis of haplotype detection. Reads haplotype data
# from Step 3, joins external predictors, and fits logistic (detection
# yes/no) and linear (delay) regression models per threshold.
#
# Inputs: haploDF_noMich_*_v2.csv
# Outputs: all_delays.csv, model CSVs/RDS, forest plots
# Dependencies: functions_pick_strains_functions.R, plot_functions.R,
#               regression_FUNCTIONS.R
# ============================================================================

library(dplyr)
library(readr)
library(purrr)
library(stringr)

# -----------------------------
# EXECUTION
# -----------------------------
rootDir <- "/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/haplo_stats3"

# source your helper functions
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/plot_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/regression_FUNCTIONS.R")


pattern_main = "haploDF_noMich"
pattern_version = "_v2"
out_dir_figs   = "../figures"
out_dir_models = "../detectionDelay_regModels"

setwd(rootDir)
ensure_dirs(c(out_dir_models, out_dir_figs))

# ---- read input files ----
files <- dir()
myFiles <- files[str_detect(files, fixed(pattern_main)) & str_detect(files, fixed(pattern_version))]
if (length(myFiles) == 0) stop("No files matched patterns in rootDir.")

all_data <- map_df(myFiles, ~ read_csv(.x, show_col_types = FALSE) %>%
                     mutate(source_file = basename(.x)))

# ---- derived cols + filtering ----
all_data <- all_data %>%
  add_derived_columns() %>%
  drop_problem_rows()

cat("[script4] after filtering | total haplotypes:", nrow(all_data),
    "| H1:", sum(all_data$subtype == "H1"), "| H3:", sum(all_data$subtype == "H3"), "\n")
message("Basic table(firstUS_season x subtype x threshold):")
print(table(all_data$firstUS_season, all_data$subtype, all_data$threshold))

# ---- entropy ----
all_data <- join_entropy(all_data)

# ---- external predictors ----
all_data <- join_external_predictors(all_data)

# save “big table”
write.csv(all_data, file = "all_delays.csv", row.names = FALSE)

# ---- per-threshold models ----
thresholds <- sort(unique(all_data$threshold))
thresholds <- thresholds[!is.na(thresholds)]
thresholds <- as.character(thresholds)

# predictors (you can edit here)
predictors_detected <- c(
  "FLUWEEK", "SEASON", "EPIDEMIC_PROP", "haplo_length", "frequency_order",
  "haplo_prop", "haplo_days_to_extinction", "pop_size", "dist_to_Mich",
  "entropy", "nSeqSubmitted"
)

predictors_delay <- predictors_detected

all_results <- list()

for (thr in thresholds) {
  message("Threshold: ", thr)
  
  df_thr <- all_data %>% filter(threshold == thr)
  
  # Summaries 
  out_sum <- summarise_detection_delay(df_thr)
  write.csv(out_sum,
            file = file.path(out_dir_models, paste0("summary_detection_delay_threshold_", thr, ".csv")),
            row.names = FALSE)
  
  # scale selected columns (IMPORTANT: do not include detected/delay in the scaling list)
  df_thr_scaled <- scale_columns(df_thr, predictors_detected)
  
  # Fit detected models
  det_res <- fit_detected_models(
    df_scaled = df_thr_scaled,
    threshold = thr,
    out_dir_models = out_dir_models,
    out_dir_figs   = out_dir_figs,
    predictors = predictors_detected
  )
  
  # Fit delay models
  del_res <- fit_delay_models(
    df_scaled = df_thr_scaled,
    threshold = thr,
    out_dir_models = out_dir_models,
    out_dir_figs   = out_dir_figs,
    predictors = predictors_delay
  )
  
  thr_key <- str_replace_all(as.character(thr), "\\.", "p")  # "0.005" -> "0p005"
  all_results[[thr_key]] <- list(
    summary = out_sum,
    detected = det_res,
    delay = del_res
  )
  print("threshold done")
}



combine_model_csvs <- function(
    folder = "/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/detectionDelay_regModels",
    out_file = file.path(folder, "ALL_detected_delay_models_combined.csv")
) {
  
  files <- list.files(
    folder,
    pattern = "^(detected_|delay_).*\\.csv$",
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop("No detected_ or delay_ CSVs found in: ", folder)
  }
  
  message("Found ", length(files), " model CSVs. Combining...")
  
  combined <- map_df(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    
    df %>%
      mutate(
        source_file = basename(f),
        model_type = if_else(
          str_starts(source_file, "detected_"),
          "detected",
          "delay"
        ),
        threshold_from_name = str_extract(source_file, "(?<=threshold_)[0-9.]+"),
        model_from_name = str_remove(source_file, "\\.csv$")
      )
  })
  
  write_csv(combined, out_file)
  message("Wrote combined CSV: ", out_file)
  
  invisible(combined)
}

# run it
combined_df <- combine_model_csvs()



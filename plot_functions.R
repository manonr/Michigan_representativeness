# ============================================================================
# plot_functions.R
# Plotting and table-formatting functions for regression results.
# Provides forest plots, faceted scatter plots, and model summary table
# extractors for both logistic and linear regression.
#
# Inputs: model objects and data frames from regression pipeline
# Outputs: ggplot objects, formatted summary tables
# Dependencies: ggplot2, dplyr, tidyr
#
# KNOWN ISSUES:
#   - plotFunction_2Y() hardcodes season factor levels (2021-2024) which may
#     be fragile if seasons change
# ============================================================================

# plot functions

# plot delays
plot_delay_vs_x <- function(df, x_var, x_label, ylim) {
  ggplot(df, aes(x = .data[[x_var]], y = delay_plot, color = label, shape = subtype)) +
    geom_point() +
    scale_y_continuous(
      name = "Delay in detection (days)",
      limits = c(0, ylim),
      breaks = c(seq(0, ylim, 50)),
      labels = c(seq(0, ylim-50, 50), "NA")
    ) +
    scale_color_manual(
      values = c(
        "Observed in UOM, same season" = "blue",
        "Observed in UOM, next season" = "darkgreen",
        "Never observed in UOM" = "red"
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    labs(title = "", color = "", x = x_label)
}

library(tidyverse)


# Plot
plotFunction_1Y <- function(df, xAxis, yAxis, yAxisLabel){
  ggplot(df, aes(x = {{xAxis}}, y = {{yAxis}}, color = threshold, group = threshold, shape = subtype)) +
    geom_point(size = 2, 
               position = position_jitter(width = 0.05, height = 0)) +
    facet_wrap(~subtype, ncol = 1,scales = "free_y") +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(breaks = scales::pretty_breaks()) +  # round y-axis ticks
    theme_bw(base_size = 14) +
    labs(
      x = "Season",
      y = yAxisLabel,
      color = "Mutation threshold",
      shape = "Subtype"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}





plotFunction_2Y <- function(df, xVar, yVar1, yVar2, yAxisLabel) {
  # Pivot to long format
  df_long <- df %>%
    select(all_of(c("subtype", "threshold", xVar, yVar1, yVar2))) %>%
    pivot_longer(cols = c(all_of(yVar1), all_of(yVar2)),
                 names_to = "measure", values_to = "value")
  
  # Make sure x-axis and threshold are factors as needed
  df_long <- df_long %>%
    mutate(
      !!xVar := as.character(.data[[xVar]]),                     # convert to character
      !!xVar := ifelse(.data[[xVar]] == "2021.2022.2023.2024",
                       "2021-2024",
                       .data[[xVar]]),                           # replace
      !!xVar := factor(.data[[xVar]],
                       levels = c("2021", "2022", "2023", "2024", "2021-2024")),
      threshold = factor(threshold),
      measure = factor(measure)
    )
  
  ggplot(df_long, aes(x = .data[[xVar]], y = value, color = threshold, shape = measure)) +
    geom_point(size = 2) +
    facet_wrap(~subtype, ncol = 1) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 14) +
    labs(
      x = "Season",
      y = yAxisLabel,
      color = "Threshold",
      shape = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

log_reg_summary_table <- function(final_model, value = c("t","z")){
  print("only makes sense for logistic regression")
  model_summary1 <- summary(final_model)
  #p_values <- coef(model_summary1)[, "Pr(>|z|)"]
  p_values <- coef(model_summary1)[, paste0("Pr(>|", value, "|)")]
  odds_ratios <- exp(coef(final_model))
  #print(odds_ratios)
  conf_int <- confint(final_model)  # Profile likelihood CIs
  OR_ci <- exp(conf_int)
  #print(OR_ci)
  summary_table <- data.frame(
    Predictor = names(odds_ratios),
    OR = odds_ratios,
    CI_lower = OR_ci[, 1],
    CI_upper = OR_ci[, 2],
    p_value = p_values
  )
  
  return(summary_table)
}

glm_summary_table <- function(final_model){
  if (!inherits(final_model, "lm")) {
    warning("This function is designed for linear regression (lm) models.")
  }
  
  # Extract summary
  model_summary <- summary(final_model)
  coef_table <- model_summary$coefficients
  
  # Confidence intervals
  ci <- confint(final_model)
  
  # Build nice table
  summary_table <- data.frame(
    Predictor = rownames(coef_table),
    Estimate  = coef_table[, "Estimate"],
    CI_lower  = ci[, 1],
    CI_upper  = ci[, 2],
    p_value   = coef_table[, "Pr(>|t|)"],
    row.names = NULL
  )
  
  return(summary_table)
}

forestPlot <- function(final_model){
  # Tidy model output with odds ratios and 95% CIs
  tidy_model <- broom::tidy(final_model, conf.int = TRUE, exponentiate = TRUE)
  
  # Remove intercept
  tidy_model <- tidy_model[tidy_model$term != "(Intercept)", ]
  
  # only show significant terms on plot
  tidy_model <- tidy_model %>%filter(p.value<0.01)
  
  # Forest plot with log-scale y-axis
  g1 <- ggplot(tidy_model, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_y_log10() +
    coord_flip() +
    labs(
      title = "",
      x = "Predictor",
      y = "Odds Ratio (log scale)"
    ) +
    theme_minimal(base_size = 14)
  g1
}

forest_plot_lm <- function(model,
                           drop_intercept = TRUE,
                           drop_small = FALSE,
                           small_thresh = 0.1,
                           drop_bigP = T,
                           xlab = "Estimate",
                           ylab = "Predictor",
                           title = NULL) {
  # summary & coefficients
  sm <- summary(model)
  coef_tab <- sm$coefficients   # matrix: Estimate, Std. Error, t value, Pr(>|t|)
  
  # confidence intervals
  ci <- confint(model)
  
  # build data frame
  df <- data.frame(
    term     = rownames(coef_tab),
    estimate = coef_tab[, "Estimate"],
    se       = coef_tab[, "Std. Error"],
    p_value  = coef_tab[, "Pr(>|t|)"],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    row.names = NULL
  )
  
  # drop intercept if requested
  if (drop_intercept) {
    df <- df %>% filter(term != "(Intercept)")
  }
  
  # drop when p> 0.05
  if (drop_bigP) {
    df <- df %>%
      filter(!(p_value > 0.05))
  }
  
  # optionally drop very small coefficients (near zero)
  if (drop_small) {
    df <- df %>%
      filter(!(estimate >= -small_thresh & estimate <= small_thresh))
  }
  
  # order predictors by estimate (or abs(estimate) if preferred)
  df <- df %>%
    mutate(term = factor(term, levels = rev(term[order(estimate)])))
  
  # build forest plot
  p <- ggplot(df, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
    geom_point(size = 2) +
    labs(
      x = xlab,
      y = ylab,
      title = title
    ) +
    theme_minimal(base_size = 14)
  
  return(p)
}



make_lm_coef_df <- function(model, drop_intercept = TRUE) {
  sm <- summary(model)
  coef_tab <- sm$coefficients
  ci <- confint(model)
  
  # Detect correct p-value column (works for lm and glm)
  p_col_name <- grep("^Pr\\(>\\|", colnames(coef_tab), value = TRUE)
  
  df <- data.frame(
    term     = rownames(coef_tab),
    estimate = coef_tab[, "Estimate"],
    se       = coef_tab[, "Std. Error"],
    p_value  = coef_tab[, p_col_name],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    row.names = NULL
  )
  
  if (drop_intercept) {
    df <- df[df$term != "(Intercept)", ]
  }
  
  df
}


library(ggplot2)
library(dplyr)

forest_plot_signedlog <- function(df,
                                  drop_small   = FALSE,
                                  small_thresh = 0.5,
                                  drop_nonsig  = TRUE,
                                  title = NULL) {
  
  # Signed log transformation
  signed_log <- function(x) sign(x) * log10(1 + abs(x))
  
  df_plot <- df
  
  # 1) Drop tiny effects around 0 if requested
  if (drop_small) {
    df_plot <- df_plot %>%
      filter(!(estimate >= -small_thresh & estimate <= small_thresh))
  }
  
  # 2) Drop CIs that overlap 0 if requested
  if (drop_nonsig) {
    df_plot <- df_plot %>%
      filter(!(CI_lower <= 0 & CI_upper >= 0))
  }
  
  # Signed-log transform
  df_plot$est_t  <- signed_log(df_plot$estimate)
  df_plot$low_t  <- signed_log(df_plot$CI_lower)
  df_plot$high_t <- signed_log(df_plot$CI_upper)
  
  # Order by effect size magnitude
  df_plot <- df_plot %>%
    mutate(abs_est = abs(estimate)) %>%
    arrange(abs_est) %>%
    mutate(term = factor(term, levels = unique(term)))
  
  # Tick marks (in original units)
  break_vals <- c(
    -6000, -1000, -100, -10, -1, 
    0,
    1, 10, 100
  )
  break_trans <- signed_log(break_vals)
  
  p <- ggplot(df_plot, aes(x = est_t, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = low_t, xmax = high_t), height = 0.2) +
    geom_point(size = 2.5) +
   scale_x_continuous(
      breaks = break_trans,
      labels = break_vals
    ) +
    labs(
      x = "Beta estimate (signed log scale)",
      y = NULL,
      title = NULL,
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  p
}

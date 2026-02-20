# =============================================================================
# script7b_stateAUC_analysis.R
#
# Computes area-under-the-curve (AUC) of cumulative haplotype coverage over
# detection delay for each US state, fits regression models predicting AUC
# from log(nSeq), season, and subtype, and compares each state's observed AUC
# to a null distribution from Michigan-sized replicates.
#
# Outputs:
#   - plot_state_AUCs.csv       (AUC per state with replicate summaries)
#   - plot_state_AUCs_CIs.csv   (CI data for plotting)
#   - regression_table.csv      (formatted regression coefficients)
#
# Companion plotting script: plot_script7b_stateAUC.R
#
# KNOWN ISSUES:
#   - [FIXED] df_auc was merged with state_season_counts twice (duplicate removed)
#
# M Ragonnet-Cronin, 2026
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(pracma)
library(stringr)
library(lme4)

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/state_analysis_functions.R")


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/")
setwd(rootDir)

all_data<- dir()[grep("haploStats_stateProxyAnalysis", dir())]
all_data <- read.csv(paste0(all_data[which.max(as.Date(sub("haploStats_stateProxyAnalysis_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(all_data))))]))
cat("[script8] state proxy data | rows:", nrow(all_data),
    "| unique states:", length(unique(all_data$thisState)),
    "| H1:", sum(all_data$subtype == "H1"), "| H3:", sum(all_data$subtype == "H3"), "\n")

rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/reps")
setwd(rootDir)

all_data2<- dir()[grep("haploStats_stateProxyAnalysis_MichReps_", dir())]
all_data2 <- read.csv(paste0(all_data2[which.max(as.Date(sub("haploStats_stateProxyAnalysis_MichReps_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(all_data2))))]))
cat("[script8] Michigan replicate data | rows:", nrow(all_data2),
    "| H1:", sum(all_data2$subtype == "H1"), "| H3:", sum(all_data2$subtype == "H3"), "\n")

# ---- 1) Prep: remove NAs in delay, compute cumulative per state within each combo ----
df_prepped <-all_data %>%
  filter(!is.na(delay), !is.na(haplo_prop)     ) %>%
  mutate(state = factor(thisState, levels = state_levels)) %>%
  arrange(subtype, threshold, season, state, delay) %>%
  group_by(subtype, threshold, season, state) %>%
  mutate(cum_prop = cumsum(haplo_prop)) %>%
  ungroup()

df_prepped2 <-all_data2%>%
  filter(!is.na(delay), !is.na(haplo_prop) ) %>%
  mutate(state = factor(thisState, levels = state_levels)) %>%
  arrange(subtype, threshold, season, state, delay, rep) %>%
  group_by(subtype, threshold, season, state, rep) %>%
  mutate(cum_prop = cumsum(haplo_prop)) %>%
  ungroup()

# ---- 2) Compute AUC per state ----

df_auc <- df_prepped %>%
  filter(delay <= 365) %>%
  group_by(state, subtype, threshold, season, delay) %>%
  summarise(cum_prop = max(cum_prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(state, subtype, threshold, season, delay) %>%
  group_by(state, subtype, threshold, season) %>%
  summarise(
    n_points      = n(),
    x_min         = min(delay),
    x_max         = max(delay),
    auc           = if (n() >= 2) trapz(delay, cum_prop) else NA_real_,
    max_cum_prop  = max(cum_prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(auc_norm_365= auc / 365)


df_auc <-  merge(df_auc, state_season_counts, by=c("state", "season", "subtype"))

# ---- 3) Regression: AUC ~ log(nSeq) + season + subtype ----

cor.test(df_auc$auc_norm_365, df_auc$nSeq)
cor.test(df_auc$auc_norm_365, log(df_auc$nSeq))

df_auc$state <- relevel(df_auc$state, ref = "Michigan")
#
state_metrics <- read.csv("../../../state_metrics.csv")
df_auc2 <- merge(df_auc, state_metrics, by="state", all.x=T)
df_auc2 <- df_auc2[-which(df_auc2$season=="2021.2022.2023.2024"),]

df_auc2 <- df_auc2 |>
  mutate(
    state   = factor(state),
    season  = factor(season),
    subtype = factor(subtype),
    nSeq_sc = scale(nSeq),
    lognseq = log(nSeq)
  ) %>%
  filter(season !="2021.2022.2023.2024")
df_auc2$season <- relevel(df_auc2$season, ref = "2024")




# final model
mod_auc <- glm(
  auc_norm_365 ~  log(nSeq) + season + subtype,
  data = df_auc2,
  family = gaussian()
)
performance::r2(mod_auc)
summary(mod_auc)
car::Anova(mod_auc, type = "III")


# models not used, less good fit

#XXXX


#mod_auc <- glm(
#  auc_norm_365 ~  log(nSeq)+state+ season + subtype,
#  data = df_auc2,
#  family = gaussian()
#)
#performance::r2(mod_auc)



mod_final <- lm(auc_norm_365 ~ log(nSeq) + season + subtype, data = df_auc2)

tab <- format_regression_table(mod_final)
print(tab, row.names = FALSE)

#mod_lm <- lm(
#  auc_norm_365 ~ log(nSeq)*state + season + subtype,
#  data = df_auc2
#)
effectsize::eta_squared(car::Anova(mod_final, type = "III"), partial = TRUE)

#mod_final <- lm(auc_norm_365 ~ log(nSeq) + as.numeric(season) + subtype, data = df_auc2)
#summary(mod_final)
#performance::r2(mod_final)
#car::Anova(mod_final, type = "III")

tab <- format_regression_table(mod_final)
print(tab, row.names = FALSE)
write.csv(tab, "regression_table.csv", row.names = FALSE)

mod_final <- lm(auc_norm_365 ~ log(nSeq) + season + subtype, data = df_auc2)
summary(mod_final)
confint(mod_final)

mod_final <- lm(auc_norm_365 ~ log(nSeq) + season + subtype, data = df_auc2)
summary(mod_final)$adj.r.squared

# ---- 4) Compute AUC for Michigan-sized replicates ----

df_auc2 <- df_prepped2 %>%
  filter(delay <= 365) %>%
  group_by(state, subtype, threshold, season, delay, rep) %>%
  summarise(cum_prop = max(cum_prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(state, subtype, threshold, season, rep, delay) %>%
  group_by(state, subtype, threshold, season, rep) %>%
  summarise(
    n_points      = n(),
    x_min         = min(delay),
    x_max         = max(delay),
    auc           = if (n() >= 2) trapz(delay, cum_prop) else NA_real_,
    max_cum_prop  = max(cum_prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(auc_norm_365= auc / 365) %>%
  group_by(state, subtype, threshold, season) %>%
  summarise(
    n_rep         = n(),
    all_aucs      = paste(auc, collapse=","),
    x_min         = mean(x_min, na.rm = TRUE),
    x_max         = mean(x_max, na.rm = TRUE),
    auc_mean      = mean(auc, na.rm = TRUE),
    auc_sd        = sd(auc, na.rm = TRUE),
    auc_se        = auc_sd / sqrt(n_rep),
    t_crit        = qt(0.975, df = n_rep - 1),
    auc_CI_lower  = auc_mean - t_crit * auc_se,
    auc_CI_upper  = auc_mean + t_crit * auc_se,
    auc365_mean      = mean(auc_norm_365, na.rm = TRUE),
    auc365_sd        = sd(auc_norm_365, na.rm = TRUE),
    auc365_se        = auc365_sd / sqrt(n_rep),
    auc365_CI_lower  = auc365_mean - t_crit * auc365_se,
    auc365_CI_upper  = auc365_mean + t_crit * auc365_se,
    max_cum_prop  = mean(max_cum_prop, na.rm = TRUE),
    auc_norm_365  = mean(auc_norm_365, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(-t_crit)


# df_auc already merged with state_season_counts above (line ~83) #claude-fixed-error: duplicate merge removed
df_auc2 <-  merge(df_auc2, state_season_counts, by=c("state", "season", "subtype"))
names(df_auc2)[5:length(names(df_auc2))] <- paste0(names(df_auc2)[5:length(names(df_auc2))], "_MichReps")


auc_df <- merge(df_auc, df_auc2, by=c( "state", "season", "subtype", "threshold"))

auc_df$nSeq.y <- NULL
auc_df$nSeq<- auc_df$nSeq.x
auc_df$nSeq.x <- NULL

# ---- 5) Empirical p-values: state AUC vs Michigan replicate null ----

auc_pvals <- auc_df %>%
  dplyr::mutate(
    null_values_raw = str_split(all_aucs_MichReps, ",\\s*"),
    null_values = lapply(null_values_raw, as.numeric)
  ) %>%
  rowwise() %>%
  mutate(
    null_mean = mean(unlist(null_values), na.rm = TRUE),
    n_null = sum(!is.na(unlist(null_values))),
    p_empirical = (sum(abs(unlist(null_values) - null_mean) >= abs(auc - null_mean)) + 1) / (n_null + 1)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_empirical, method = "BH"))

auc_pvals <- auc_pvals %>%
  select(state, season, subtype, threshold, auc, null_mean, p_empirical, p_adj)%>%
  arrange(p_adj)

# ---- 6) Export CSVs for plotting ----

write.csv(auc_df, "plot_state_AUCs.csv", row.names = F)

plot_long <- auc_df %>%
  select(state, season, subtype, threshold, nSeq,
         auc_norm_365,auc365_mean_MichReps, auc365_CI_lower_MichReps,auc365_CI_upper_MichReps) %>%
  pivot_longer(
    cols = c(auc_norm_365, auc365_mean_MichReps),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = dplyr::case_when(
      metric == "auc_norm_365" ~ "State AUC",
      metric == "auc365_mean_MichReps" ~ "Michigan replicate AUC",
      TRUE ~ metric
    ),
    state = factor(state, levels = state_levels)
)

ci_df <- auc_df %>%
  select(state, season, subtype, threshold, nSeq,
         auc365_mean_MichReps, auc365_CI_lower_MichReps,auc365_CI_upper_MichReps) %>%
  mutate(
    metric = "Michigan replicate AUC",
    state  = factor(state, levels = state_levels)
  ) %>%
  dplyr::rename(value = auc365_mean_MichReps,
         ymin  = auc365_CI_lower_MichReps,
         ymax  = auc365_CI_upper_MichReps)

write.csv(ci_df, "plot_state_AUCs_CIs.csv", row.names = F)

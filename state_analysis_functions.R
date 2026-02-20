# ============================================================================
# state_analysis_functions.R
# Functions and setup data for state-level haplotype detection analysis.
# On load, reads H1/H3 mutation RDS files to compute per-state sequence
# counts. Provides functions for subsampling Michigan sequences to match
# other states, processing per-state haplotype files, and multi-panel
# state comparison plots.
#
# Inputs: H1_ALLmuts.rds, H3_ALLmuts.rds
# Outputs: per-state sequence counts, subsampled data, comparison plots
# Dependencies: functions_pick_strains_functions.R
#
# KNOWN ISSUES:
#   - Top-level code runs on source() (reads RDS files, computes state_season_counts)
#     -- will fail if RDS files are not available at the hardcoded paths
#   - Commented-out alternative color palette definitions (lines ~36-49)
# ============================================================================

# state analysis functions

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")


# state colours
H1df <-  readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/H1_ALLmuts.rds")) %>% select(location, date)%>% mutate(subtype="H1")
H3df <-  readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/H3_ALLmuts.rds")) %>% select(location, date)%>% mutate(subtype="H3")

stateCounts1 <- rbind(H1df, H3df)
stateCounts1$season <- unlist(lapply(stateCounts1$date, season_from_single_date))

stateCounts_df <- as.data.frame(table(stateCounts1$location))  %>%
  #filter(n >= 20) %>%                # exclude <20 sequences
  arrange(desc(Freq))%>%
  dplyr::rename(state = "Var1")
state_levels <- stateCounts_df$state

state_season_counts <- as.data.frame(table(stateCounts1$location, stateCounts1$season,stateCounts1$subtype))
names(state_season_counts) <- c("state", "season", "subtype", "nSeq")

# Expand and fill
state_season_counts2 <- state_season_counts%>%
  group_by(state, subtype) %>%
  summarise(nSeq = sum(nSeq, na.rm = TRUE))%>%
  mutate(season = "2021.2022.2023.2024")

state_season_counts <- rbind(state_season_counts, state_season_counts2)
# Color mapping: continuous gradient by counts, fixed once here

state_colors <- c("#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"
                      ,"#4eb3d3","#2b8cbe","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b","#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#006d2c"
                      ,"#00441b","#ffffcc","#ffeda0","#fed976","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#fc9272","#fb6a4a",
                      "#ef3b2c","#cb181d","#a50f15","#67000d", "#d9d9d9","#bdbdbd"  , "#969696","#737373","#525252","#252525","#000000")

#state_colors <- rev(c("#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"
 #                     ,"#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b","#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#006d2c"
  #                    ,"#00441b","#ffffcc","#ffeda0","#fed976","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#fc9272","#fb6a4a"
   #                   ,"#ef3b2c","#cb181d","#a50f15","#67000d"))


#user_palette <- rev(c(
#  "#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177",
#  "#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d",
#  "#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b",
#  "#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#006d2c","#00441b",
#  "#ffffcc","#ffeda0","#fed976","#fdd0a2","#fdae6b","#fd8d3c","#f16913",
#  "#d94801","#a63603","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"
#))

names(state_colors) <- stateCounts_df$state


# downsampling
keepMich_n <- function(mutDF, mySamplingFreq){
  #allUMich <- length(which(mutDF$UMich=="UMich"))
  #print(paste("all mich seq = ", allUMich))
  #sampleUMich_drop <- round((1-sample_prop)*allUMich,0)
  
  mySamplingFreq<- mySamplingFreq %>%
    mutate(mySeasons = as.integer(mySeasons)) %>%   # ensure integer
    complete(mySeasons = 2021:2024, fill = list(nStateSeq = 0))
  
  safe_sample <- function(available_idx, n_available, n_request) {
    # Handle empty case
    if (n_available == 0) {
      warning("No available items to sample from.")
      return(integer(0))
    }
    
    # Cap sample size to available count
    sample_size <- min(n_available, n_request)
    
    sample(available_idx, size = sample_size, replace = FALSE)
  }
  
  
  michSample2021 <- safe_sample(available_idx = which(mutDF$UMich == "UMich" & mutDF$season == 2021), 
                                n_available = length(which(mutDF$UMich == "UMich" & mutDF$season == 2021)), 
                                n_request = mySamplingFreq$nStateSeq[mySamplingFreq$mySeasons == 2021])
  michSample2022 <- safe_sample(available_idx = which(mutDF$UMich == "UMich" & mutDF$season == 2022), 
                                n_available = length(which(mutDF$UMich == "UMich" & mutDF$season == 2022)), 
                                n_request = mySamplingFreq$nStateSeq[mySamplingFreq$mySeasons == 2022])
  michSample2023 <- safe_sample(available_idx = which(mutDF$UMich == "UMich" & mutDF$season == 2023), 
                                n_available = length(which(mutDF$UMich == "UMich" & mutDF$season == 2023)), 
                                n_request = mySamplingFreq$nStateSeq[mySamplingFreq$mySeasons == 2023]) 
  michSample2024 <- safe_sample(available_idx = which(mutDF$UMich == "UMich" & mutDF$season == 2024), 
                                n_available = length(which(mutDF$UMich == "UMich" & mutDF$season == 2024)), 
                                n_request = mySamplingFreq$nStateSeq[mySamplingFreq$mySeasons == 2024])
  
  
  
  nonMich <- which(mutDF$UMich!="UMich")
  newmutDF <- mutDF[c(nonMich, michSample2021,michSample2022,michSample2023,michSample2024 ),]
  print(paste("kept n mich = ",length(which(newmutDF$UMich=="UMich"))))
  return(newmutDF)
}





get_state <- function(path) {
  bn <- basename(path)
  # capture the part after haploDF_H1_ or haploDF_H3_ and before _stateProxyAnalysis
  sub("haploDF_[^_]+_([^_]+)_stateProxyAnalysis.*", "\\1", bn)
}


# ---------------------------------------
# 2) Page builder: one page per subtype
# ---------------------------------------
plot_panel_for_subtype <- function(dat, st, state_levels, state_colors) {
  # Show only states present on this page, but keep global count-based order
  present_levels <- intersect(state_levels, as.character(unique(dat$thisState)))
  
  ggplot(dat, aes(x = delay, y = cum_prop, color = thisState, group = thisState)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    facet_grid(rows = vars(threshold), cols = vars(season), labeller = label_both) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))) +
    scale_color_manual(
      values = state_colors,
      breaks = present_levels,  # legend order by counts among states present
      guide  = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))
    ) +
    labs(
      title = paste0("Cumulative haplo_prop vs delay — subtype = ", st),
      x = "Delay (weeks)",
      y = "Cumulative haplo_prop",
      color = "State (by sequence count)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.spacing = unit(0.7, "lines"),
      strip.text = element_text(face = "bold")
    )
}

format_regression_table <- function(model, digits = 2) {
  coefs <- summary(model)$coefficients
  ci <- confint(model)
  
  est <- round(coefs[, "Estimate"], digits)
  lo  <- round(ci[, 1], digits)
  hi  <- round(ci[, 2], digits)
  pv  <- coefs[, "Pr(>|t|)"]
  
  # Format p-values: use "<0.001" for very small values
  p_fmt <- ifelse(pv < 0.001, "<0.001", format(round(pv, 3), nsmall = 3))
  
  data.frame(
    Variable           = rownames(coefs),
    `Estimate (95% CI)` = paste0(format(est, nsmall = digits), " (",
                                 format(lo, nsmall = digits), ", ",
                                 format(hi, nsmall = digits), ")"),
    `P value`          = p_fmt,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


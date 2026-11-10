# ============================================================================
# scriptX_empirical_mut_freq.R
# Reviewer response: quantify and plot EMPIRICAL mutation frequencies, to show
# where the predetermined thresholds (script1) fall in the real distribution.
#
# For each subtype x season, at every variable HA amino-acid position:
#   (a) sequence-weighted non-dominant frequency  = 1 - freq(most common AA)
#       (all sequences counted; the empirical frequency of variation)
#   (b) unique-haplotype-weighted non-dominant frequency = what getMutList()
#       in script1 actually thresholds (computed on !duplicated(aa) sequences,
#       positions with any NA dropped, seasons with <20 sequences skipped)
#   (c) frequency of every residue (dominant and non-dominant) at each site
#
# Inputs : haplo_stats3/H1_ALLmuts.rds, H3_ALLmuts.rds  (from script1)
# Outputs: haplo_stats3/empirical_mut_freq/
#            site_freqs_*.csv, mutation_freqs_*.csv, sites_above_threshold_*.csv,
#            fig1-4_*.pdf
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

dataDir <- "/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/"
outDir  <- paste0(dataDir, "empirical_mut_freq/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

subtypes   <- c("H1", "H3")
seasons    <- 2021:2024
thresholds <- c(0.2, 0.15, 0.1, 0.05, 0.01, 0.005)   # same as script1
badAA      <- c("X", "*", "-", "")                     # ambiguous / stop / gap

# ---- helpers ---------------------------------------------------------------

# Site-level table for one season's sequences
site_freqs_one <- function(df) {
  aaposCols <- grep("^aapos", names(df), value = TRUE)
  uniqDF    <- df[!duplicated(df$aa), ]   # unique sequences, as in getMutList()

  site_rows <- lapply(aaposCols, function(col) {
    x <- df[[col]]
    x <- x[!is.na(x) & !x %in% badAA]
    if (length(x) == 0) return(NULL)
    tb <- table(x)
    # unique-haplotype measure, replicating getMutList()
    u <- uniqDF[[col]]
    nondom_unique <- if (nrow(df) >= 20 && !any(is.na(u))) {
      tu <- table(u)
      1 - max(tu) / sum(tu)
    } else NA_real_
    tibble(site = as.integer(sub("aapos", "", col)),
           n_seq = length(x),
           n_aa = length(tb),
           dominant_aa = names(tb)[which.max(tb)],
           nondom_freq_seq = 1 - max(tb) / length(x),
           nondom_freq_unique = nondom_unique)
  })
  bind_rows(site_rows)
}

# Residue-level table: every residue (dominant and non-dominant) at each site
mutation_freqs_one <- function(df) {
  aaposCols <- grep("^aapos", names(df), value = TRUE)
  rows <- lapply(aaposCols, function(col) {
    x <- df[[col]]
    x <- x[!is.na(x) & !x %in% badAA]
    if (length(x) == 0) return(NULL)
    tb <- sort(table(x), decreasing = TRUE)
    site <- as.integer(sub("aapos", "", col))
    tibble(site = site,
           aa = names(tb),
           is_dominant = seq_along(tb) == 1,
           count = as.integer(tb),
           n_seq = length(x),
           freq = as.integer(tb) / length(x),
           residue = paste0(site, names(tb)))
  })
  bind_rows(rows)
}

# ---- compute ---------------------------------------------------------------

siteList <- list(); mutList <- list()
for (subtype in subtypes) {
  bigMutDF <- readRDS(paste0(dataDir, subtype, "_ALLmuts.rds"))
  for (s in seasons) {
    d <- bigMutDF %>% filter(season == s)
    if (nrow(d) == 0) next
    cat("[scriptX]", subtype, s, "| n =", nrow(d), "\n")
    siteList[[paste(subtype, s)]] <- site_freqs_one(d) %>%
      mutate(subtype = subtype, season = s, n_total_seq = nrow(d))
    mutList[[paste(subtype, s)]] <- mutation_freqs_one(d) %>%
      mutate(subtype = subtype, season = s)
  }
}
siteDF <- bind_rows(siteList)
mutDF  <- bind_rows(mutList)
write.csv(siteDF, paste0(outDir, "site_freqs_", Sys.Date(), ".csv"), row.names = FALSE)
write.csv(mutDF,  paste0(outDir, "mutation_freqs_", Sys.Date(), ".csv"), row.names = FALSE)

# ---- how many sites does each threshold capture? ---------------------------
# per season, and the union across seasons (script1's "2021-2024" definition)
count_above <- function(measure) {
  per_season <- bind_rows(lapply(thresholds, function(th) {
    siteDF %>%
      filter(!is.na(.data[[measure]])) %>%
      group_by(subtype, season) %>%
      summarise(n_sites = sum(.data[[measure]] > th), .groups = "drop") %>%
      mutate(threshold = th, season = as.character(season))
  }))
  union_df <- bind_rows(lapply(thresholds, function(th) {
    siteDF %>%
      filter(!is.na(.data[[measure]]), .data[[measure]] > th) %>%
      group_by(subtype) %>%
      summarise(n_sites = n_distinct(site), .groups = "drop") %>%
      mutate(threshold = th, season = "2021-2024 (union)")
  }))
  bind_rows(per_season, union_df) %>% mutate(measure = measure)
}
aboveDF <- bind_rows(count_above("nondom_freq_seq"), count_above("nondom_freq_unique")) %>%
  arrange(measure, subtype, season, desc(threshold))
write.csv(aboveDF, paste0(outDir, "sites_above_threshold_", Sys.Date(), ".csv"), row.names = FALSE)
print(aboveDF %>% pivot_wider(names_from = threshold, values_from = n_sites), n = Inf)

# ---- plots -----------------------------------------------------------------
measure_labels <- c(nondom_freq_seq    = "Sequence-weighted",
                    nondom_freq_unique = "Unique-haplotype-weighted (used in script1)")
siteLong <- siteDF %>%
  select(subtype, season, site, nondom_freq_seq, nondom_freq_unique) %>%
  pivot_longer(starts_with("nondom"), names_to = "measure", values_to = "freq") %>%
  filter(!is.na(freq), freq > 0) %>%
  mutate(measure = factor(measure_labels[measure], levels = measure_labels),
         season = factor(season))

thr_lines <- geom_vline(xintercept = thresholds, linetype = "dashed",
                        colour = "grey40", linewidth = 0.3)

# Residue-level data (all residues at variable sites), linear scale plots
mutPlot <- mutDF %>%
  group_by(subtype, season, site) %>%
  filter(n() > 1) %>%    # drop invariant sites (a single residue at frequency 1)
  ungroup() %>%
  mutate(season = factor(season),
         residue_type = ifelse(is_dominant, "Dominant", "Non-dominant"))

# 1. Distribution of residue frequencies (linear scale, to see any discontinuity)
p1 <- ggplot(filter(mutPlot, !is_dominant), aes(freq)) +
  geom_histogram(binwidth = 0.005, boundary = 0, fill = "firebrick",
                 colour = "white", linewidth = 0.1) +
  coord_cartesian(xlim = c(0, max(mutPlot$freq[!mutPlot$is_dominant]))) +  # x ends at max observed frequency
  facet_grid(subtype ~ season) +   # fixed axes: same y range in every panel, ending at the max bin count
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Non-dominant residue frequency", y = "Number of residues (site x amino acid)",
       title = "Figure 1. Empirical distribution of non-dominant residue frequencies",
       subtitle = "Variable sites only (sites with a single residue excluded); bin width 0.005") +
  theme_bw()
ggsave(paste0(outDir, "fig1_residue_freq_distribution.pdf"), p1, width = 11, height = 6)

# 1b. Table behind Figure 1: number of non-dominant residues per 0.005 frequency bin
#     (empty bins included, so gaps in the distribution are visible)
binBreaks <- round(seq(0, 1, by = 0.005), 3)   # all thresholds fall on bin edges
nondomDF <- filter(mutPlot, !is_dominant)
binTable <- nondomDF %>%
  mutate(bin = cut(freq, breaks = binBreaks, include.lowest = TRUE)) %>%   # (lower, upper]
  count(subtype, season, bin, .drop = FALSE) %>%
  mutate(bin_idx = as.integer(bin),
         bin_lower = binBreaks[bin_idx], bin_upper = binBreaks[bin_idx + 1]) %>%
  rename(n_residues = n) %>%
  group_by(subtype, season) %>%
  mutate(cum_n_residues_from_top = rev(cumsum(rev(n_residues)))) %>%   # residues with freq > bin_lower
  ungroup() %>%
  left_join(distinct(siteDF, subtype, season, n_total_seq) %>% mutate(season = factor(season)),
            by = c("subtype", "season")) %>%
  select(subtype, season, n_total_seq, bin_lower, bin_upper, n_residues, cum_n_residues_from_top)
write.csv(binTable, paste0(outDir, "fig1_binned_counts_", Sys.Date(), ".csv"), row.names = FALSE)
print(binTable %>% filter(n_residues > 0), n = 60)

# 1c. Same data, log y-axis, so the sparse higher-frequency bins are visible
p1b <- binTable %>%
  filter(n_residues > 0) %>%
  ggplot(aes(xmin = bin_lower, xmax = bin_upper, ymin = 0.5, ymax = n_residues)) +
  geom_rect(fill = "firebrick", colour = "white", linewidth = 0.1) +
  thr_lines +   # grey dashed lines at the script1 thresholds
  coord_cartesian(xlim = c(0, max(nondomDF$freq))) +
  scale_y_log10() +
  facet_grid(subtype ~ season) +
  labs(x = "Non-dominant residue frequency (bin midpoint)",
       y = "Number of residues (log scale)",
       title = "Figure 1b. As Figure 1, log y-axis; empty bins not shown",
       subtitle = "Variable sites only; bin width 0.005") +
  theme_bw()
ggsave(paste0(outDir, "fig1b_residue_freq_distribution_logy.pdf"), p1b, width = 11, height = 6)

# 2. Number of sites exceeding each frequency (reverse cumulative)
revcum <- siteLong %>%
  group_by(subtype, season, measure) %>%
  arrange(desc(freq)) %>%
  mutate(n_sites_ge = row_number()) %>%
  ungroup()
p2 <- ggplot(revcum, aes(freq, n_sites_ge, colour = season)) +
  geom_step() +
  thr_lines +
  scale_x_log10() +
  facet_grid(measure ~ subtype, scales = "free_y") +
  labs(x = "Non-dominant amino acid frequency (log scale)",
       y = "Number of sites with frequency ≥ x",
       title = "Figure 2. Sites captured as a function of threshold",
       subtitle = "Dashed lines = thresholds used in script1") +
  theme_bw()
ggsave(paste0(outDir, "fig2_sites_vs_threshold.pdf"), p2, width = 10, height = 7)

# 3. Residue frequencies along HA (all residues at variable sites)
p3 <-ggplot(mutPlot, aes(site, freq, colour = residue_type)) +
  geom_hline(yintercept = thresholds, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 1) +
  geom_text(data = filter(mutPlot, !is_dominant, freq >= 0.1), aes(label = residue),
            size = 2, vjust = -0.6, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = c(Dominant = "grey40", `Non-dominant` = "firebrick")) +
  facet_grid(season ~ subtype) +
  labs(x = "HA amino acid position", y = "Residue frequency",
       colour = NULL,
       title = "Figure 3. Empirical frequency of every residue at variable sites",
       subtitle = "Dashed lines = thresholds used in script1; non-dominant residues labelled if ≥ 10%") +
  theme_bw()
ggsave(paste0(outDir, "fig3_residue_freq_by_position.pdf"), p3, width = 12, height = 10)

# 4. Sequence-weighted vs unique-haplotype-weighted (what script1 thresholds)
p4 <- siteDF %>%
  filter(!is.na(nondom_freq_unique)) %>%
  ggplot(aes(nondom_freq_seq, nondom_freq_unique, colour = factor(season))) +
  geom_abline(linetype = "dotted") +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~subtype) +
  labs(x = "Sequence-weighted non-dominant frequency",
       y = "Unique-haplotype-weighted non-dominant frequency",
       colour = "Season",
       title = "Figure 4. Sequence- vs unique-haplotype-weighted site frequency") +
  theme_bw()
ggsave(paste0(outDir, "fig4_seq_vs_unique_freq.pdf"), p4, width = 9, height = 4.5)

cat("[scriptX] done. Outputs in", outDir, "\n")

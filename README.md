# Influenza HA Haplotype Representativeness Analysis

R scripts for analysing the representativeness of influenza hemagglutinin (HA) haplotype surveillance at the University of Michigan (UOM) relative to the broader US population. This project tracks how amino-acid haplotypes at high-frequency mutation sites are detected over time at UOM compared with the rest of the US, across H1 and H3 subtypes and multiple flu seasons (2021–2024).

## Repository structure

```
.
├── functions_pick_strains_functions.R   # Core utility functions (sequence I/O, haplotypes, entropy, delays, downsampling)
├── state_analysis_functions.R           # Functions and setup data for state-level haplotype detection analysis
├── regression_FUNCTIONS.R               # Helper functions for regression data prep, column derivation, and model fitting
├── plot_functions.R                     # Shared plotting helpers (forest plots, scatter plots, model summary tables)
├── predictors_for_regression.R          # Prepare external predictor variables (ILINet, Census, distance, sequencing counts)
├── prep_nexstrain_michigan_v2.R         # Prepare Michigan sequences for Nextstrain builds (GISAID + GenBank merging)
│
│   # --- Analysis pipeline (run in order) ---
├── script1_getMuts_michigan_US_subsampling_FINAL.R   # Step 1: Identify mutation sites & define haplotypes
├── script2_detecting_haplos_later_FINAL.R            # Step 2: Track cross-season haplotype detection
├── script3_haplotype_stats_FINAL.R                   # Step 3: Compute coverage metrics & detection delays
├── script4_bigRegression_FINAL.R                     # Step 4: Logistic & linear regression on haplotype detection
│
│   # --- Sensitivity & extended analyses ---
├── script5_downsamplingAnalysis.R    # Downsampling Michigan sequences to test robustness
├── script6_epitopeAnalysis.R         # Repeat analysis restricted to antigenic epitope sites
├── script7_otherStateProxies.R       # Extend analysis to all US states as Michigan proxies
├── script8_stateAUC_analysis.R       # State-level AUC computation & regression
└── script9_rarefaction_analysis.R    # Rarefaction / iNEXT coverage analysis
```

## Pipeline overview

1. **Sequence processing** (`script1`): Reads influenza HA sequences, identifies high-frequency mutation sites, and extracts amino-acid haplotypes for UOM and the rest of the US. Computes Shannon entropy and builds ranked haplotype frequency tables per subtype, threshold, and season.

2. **Cross-season detection** (`script2`): Applies mutation site definitions from Step 1 across all seasons to track when each haplotype is first observed in the US vs. UOM.

3. **Haplotype statistics** (`script3`): Updates first-detection dates across seasons and computes coverage metrics (p50, p90) and counts of haplotypes missing from Michigan.

4. **Regression analysis** (`script4` + `regression_FUNCTIONS.R` + `predictors_for_regression.R`): Fits logistic (detection yes/no) and linear (delay) regression models to identify predictors of haplotype detection and delay (e.g., haplotype frequency, epidemic timing, season, subtype). Produces forest plots of model coefficients.

5. **Sensitivity analyses**:
   - `script5`: Downsamples Michigan sequences at varying proportions to assess robustness of results.
   - `script6`: Repeats the full pipeline restricted to antigenic epitope sites on HA.
   - `script7`: Extends the haplotype representativeness analysis to all US states, comparing each state's detection performance.
   - `script8`: Computes AUC of cumulative haplotype coverage over detection delay per state, fits regression models, and compares observed vs. null distributions.

6. **Rarefaction** (`script9`): Uses iNEXT to estimate haplotype diversity coverage and determine the number of sequences needed for 95% coverage.

## Helper scripts

- **`functions_pick_strains_functions.R`**: Core functions used across the pipeline — sequence translation, haplotype extraction, entropy calculation, detection delay tracking, and downsampling utilities.
- **`state_analysis_functions.R`**: Per-state sequence counts, Michigan downsampling functions, and state comparison plotting utilities.
- **`regression_FUNCTIONS.R`**: Data preparation, predictor joining, and model fitting routines shared by the regression scripts.
- **`plot_functions.R`**: Shared plotting functions for forest plots, faceted scatter plots, and regression summary tables.
- **`predictors_for_regression.R`**: Processes ILINet surveillance data, US Census population estimates, geographic distances to Michigan, and per-state sequencing counts into output CSVs (`episize_stats.csv`, `state_metrics.csv`, `state_seq_numbers.csv`).
- **`prep_nexstrain_michigan_v2.R`**: Combines GISAID and GenBank sequences, deduplicates, assigns US state labels and Michigan flags, and produces Nextstrain input files (`epi.tsv` and FASTA).

## Dependencies

Key R packages required:

- **Bioinformatics**: `bioseq`, `Biostrings`, `ape`
- **Data wrangling**: `dplyr`, `tidyr`, `purrr`, `stringr`, `readr`, `reshape2`, `tibble`, `lubridate`
- **Plotting**: `ggplot2`, `ggpubr`, `ggforce`, `viridisLite`, `scales`, `patchwork`
- **Statistics**: `lme4`, `car`, `performance`, `effectsize`, `broom`, `pracma`
- **Diversity**: `iNEXT`

## Data

Input data (FASTA files and processed CSVs) are not included in this repository. Sequences were obtained from GISAID and GenBank. The analysis expects data files in a directory structure under `Representativeness/haplo_stats3/`.

## Author

M Ragonnet-Cronin, 2024

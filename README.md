# Influenza HA Haplotype Representativeness Analysis

R scripts for analysing the representativeness of influenza hemagglutinin (HA) haplotype surveillance at the University of Michigan (UOM) relative to the broader US population. This project tracks how amino-acid haplotypes at high-frequency mutation sites are detected over time at UOM compared with the rest of the US, across H1 and H3 subtypes and multiple flu seasons (2021-2024).

## Repository structure

```
.
├── functions_pick_strains_functions.R   # Core utility functions (sequence I/O, haplotypes, entropy, delays)
├── state_analysis_functions.R           # Functions for state-level proxy analysis
├── plotting/
│   └── plot_functions.R                 # Shared plotting helper functions
│
│   # --- Analysis pipeline (run in order) ---
├── script1_getMuts_michigan_US_v5_subsampling.R   # Step 1: Read sequences, identify mutations & haplotypes
├── script2_detecting_haplos_later_v2.R            # Step 2: Track cross-season haplotype detection
├── script3_haplotype_stats_v2.R                   # Step 3: Compute haplotype statistics (coverage, delay, entropy)
├── script4_bigRegression_FUNCTIONS.R              # Step 4a: Regression model functions
├── script4_bigRegression_EXECUTE.R                # Step 4b: Execute logistic/linear regression analysis
├── predictors_for_regression.R                    # Step 4c: Prepare predictor variables for regression
│
│   # --- Sensitivity & extended analyses ---
├── script5_downsamplingAnalysis.R                 # Downsampling Michigan sequences to test robustness
├── script6_epitopeAnalysis.R                      # Repeat analysis restricted to antigenic epitope sites
├── script7_otherStateProxies.R                    # State-level proxy analysis (data prep & haplotype tracking)
├── script7b_stateAUC_analysis.R                   # State-level AUC computation & regression
├── rarefaction_analysis.R                         # Rarefaction / iNEXT coverage analysis
│
│   # --- Data preparation ---
├── prep_nexstrain_michigan_v2.R                   # Prepare Michigan sequences for Nextstrain builds
│
│   # --- Manuscript numbers & figures ---
├── numbers_for_paper.R                            # Summary statistics cited in the manuscript
├── ms_fig4_Forest_plots.R                         # Figure 4: forest plots of regression coefficients
├── plotting/
│   ├── plot_fig1_script1_general_plots.R          # Figure 1: sequence counts, haplotype frequencies
│   ├── ms_fig2_plot.R                             # Figure 2: representativeness summary
│   ├── ms_fig3_MichDelays_entropy.R               # Figure 3: Michigan detection delays & entropy
│   ├── FigS1_detectionDelayMich.R                 # Fig S1: Michigan-specific delay scatter plots
│   ├── plot_script2_representativeness_combined.R # Supp: representativeness metrics (main + epitope)
│   ├── plot_script4_plot_effects_downsampling.R   # Supp: downsampling effect plots
│   ├── plot_script5_plot_delays_combined.R        # Supp: detection delay plots (main + epitope)
│   ├── plot_script7b_stateAUC.R                   # Supp: state AUC comparison plots
│   ├── plot_state_comparisons.R                   # Supp: cumulative coverage curves per state
│   ├── plot_state_proxies.R                       # Supp: state proxy scatter/heatmap plots
│   ├── plot_sequence_source_mich_nonMich.R        # Supp: sequence source bar charts
│   └── SupTables.R                                # Supplementary tables
│
│   # --- Archived ---
├── old_scripts/                                   # Superseded or exploratory scripts (not part of final analysis)
```

## Pipeline overview

1. **Sequence processing** (`script1`): Reads influenza HA FASTA files, translates DNA to amino acids, identifies high-frequency mutation sites, and extracts amino-acid haplotypes for UOM and the rest of the US.

2. **Cross-season detection** (`script2`): Tracks when each haplotype is first observed in the US vs. UOM, computing detection delays across seasons.

3. **Haplotype statistics** (`script3`): Computes coverage metrics (proportion of US haplotypes observed at UOM), Shannon entropy of haplotype distributions, and summary delay statistics.

4. **Regression analysis** (`script4` + `predictors_for_regression`): Fits logistic and linear regression models to identify predictors of haplotype detection and delay (e.g., haplotype frequency, epidemic timing, season, subtype).

5. **Sensitivity analyses**:
   - `script5`: Downsamples Michigan sequences to assess robustness
   - `script6`: Repeats the analysis restricted to antigenic epitope sites
   - `script7`/`script7b`: Extends analysis to other US states as Michigan proxies, computes AUC of cumulative coverage, and tests whether states differ from Michigan-sized random samples

6. **Rarefaction** (`rarefaction_analysis`): Uses iNEXT to estimate haplotype diversity coverage.

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

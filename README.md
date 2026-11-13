# R Statistical Mutation Analyzer

An R script to process mutation summary files across peptide screening subdirectories, perform Fisher's Exact Tests, Chi-Square tests, permutation tests, Bonferroni-adjusted pairwise comparisons, and generate ggplot visualizations.

## Features

- Recursively scans project subdirectories for mutation summary Excel files.
- Categorizes amino acids into physicochemical groups.
- Executes simulated Fisher's Exact tests, Chi-Square tests, and non-parametric permutation tests (`coin` package).
- Exports high-resolution visualization plots (`300 DPI`).

## Installation & Usage

Make sure R is installed with required packages (`readxl`, `dplyr`, `coin`, `ggplot2`), then run:
```R
Rscript statistical_analysis.R
```

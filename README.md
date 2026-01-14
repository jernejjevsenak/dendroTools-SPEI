# dendroTools–SPEI climate–growth correlations (R scripts)

This repository contains **R scripts for calculating climate–growth correlations** using tree-ring chronologies and **SPEI-based water balance metrics**, for both:

- **Daily climate data** (daily PET → daily water deficit → daily SPEI across many day-scales)
- **Monthly climate data** (daily-to-monthly aggregation → monthly PET / water deficit → monthly SPEI across many month-scales)

The workflow is designed to help identify **seasonal windows (time scales and timing)** where climate variability (via SPEI) is most strongly correlated with tree growth.

## What the scripts do

- Read a tree-ring width chronology (`.crn`) and climate time series (temperature + precipitation).
- Compute **PET** (Hargreaves) and **water deficit** (Prec − PET).
- Compute **SPEI** for a range of aggregation scales (days or months).
- Correlate SPEI values with the chronology (Pearson correlation) over moving seasonal windows.
- Visualize results as a **heatmap** (season length × day-of-year/month).

## Requirements

R packages used across scripts include (not exhaustive):

`dplyr`, `lubridate`, `reshape2`, `dendroTools`, `SPEI`, `zoo`, `dplR`, `lmomco`, `ggplot2`

## Inputs

Typical input files:
- Tree-ring chronology: `my_chron.crn`
- Climate CSVs (daily): `my_chron_Tmean.csv`, `my_chron_Tmax.csv`, `my_chron_Tmin.csv`, `my_chron_Psum.csv`

Expected climate columns (daily scripts): `Y`, `M`, `D`, and the variable column (e.g., `Tmean`, `Tmax`, `Tmin`, `Prec/Psum`).

## Outputs

- Correlation matrix (scale × time position)
- Heatmap figure saved as `.png` (example names in scripts):
  - `SPEI_example2.png` (daily)
  - `SPEI_example_monthly.png` (monthly)

## Notes

- Daily scripts may include an option to append **previous-year** climate predictors.
- Monthly scripts aggregate daily data to monthly means/sums before SPEI calculation.




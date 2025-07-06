# HTM-CRB

[Script1](Script1_FilteringAndVisualization.R) : Filtering and Visualization  


## Overview

This script performs the core data processing and exploratory analysis of protein group LFQ (Label-Free Quantification) intensity data from mass spectrometry experiments. It prepares the data and generates key visualizations to assess variability and condition differences.

## Requirements

- R (recommended version 4.0 or higher)  
- R packages: `dplyr`, `tidyr`, `ggplot2`, `car`, `vegan`, `pheatmap`, `cowplot`, `showtext`  
  *(the script automatically installs missing packages if needed)*

## Installation and Loading of Libraries

The script starts by checking for and installing required packages, then loading them.  
It also sets up the "Open Sans" font for plotting via the `showtext` package and defines global plot settings including enlarged font sizes and high resolution for better clarity.

## Main Steps

1. **Data Loading and Preparation**
   - Selects only LFQ intensity columns from the protein groups data.
   - Converts data to a numeric matrix, replacing zeros with `NA` to avoid errors in log transformation.
   - Applies a log2 transformation.
   - Imputes missing or zero values with a small value below the global minimum observed intensity to allow visualization and statistical analyses.

2. **Heatmap Generation**
   - Creates a clustered heatmap of log2 LFQ intensities across samples.
   - Includes sample condition annotations.
   - Skips heatmap if insufficient variance or data quality.

3. **NMDS and ANOSIM Analysis**
   - Performs Non-metric Multidimensional Scaling (NMDS) on Euclidean distances between samples.
   - Runs ANOSIM to test for differences between condition groups.
   - Produces an NMDS plot annotated with ANOSIM statistics.
   - Handles cases with insufficient samples, groups, or variance gracefully.

4. **Batch Processing for Multiple Works**
   - Iterates over multiple defined experimental comparisons ("works"), each with specific condition patterns.
   - Calls the main processing function with error handling to skip failed analyses.
   - Outputs variability metrics and saves all plots and results to disk.

---

## Usage

- Set the base path and condition patterns for each work.
- Run the script to process all defined works sequentially.
- Inspect generated heatmaps, NMDS plots, and ANOSIM results for exploratory insights into proteomic variability and condition effects.

---

## Notes

- The script is designed to be robust to missing data and zero intensities.
- Imputation uses a small value relative to the minimum detected intensity to avoid biases.
- Works with multiple experimental comparisons via a flexible batch execution loop.
- Subsequent scripts in the project continue downstream analyses and result integration.

---


# HTM-CRB

[1. Filtering and Visualization](Script1_FilteringAndVisualization.R)  

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

[2. Hypothesis Testing Pipeline for Proteomics Differential Analysis](Script2_HypothesisTesting.R)

## Overview
The script performs hypothesis testing for differential protein expression using multiple statistical methods. It generates and saves diagnostic plots to evaluate the distribution and ranking of adjusted p-values across methods.

## Features

- Supports multiple statistical tests:
  - t-Student
  - t-Welch
  - limma
  - DEqMS
  - MSstats
- Density plots of `-log10(adjusted p-value)` across all methods
- Ranked line plots of `-log10(adjusted p-value)` for method comparison
- Automatic detection of protein ID columns and adjusted p-value columns
- Handles NULL or missing data gracefully
- Saves plots in high-resolution PNG format
- Parallelization with `future` package (reset to sequential at the end)

## Outputs

### 1. Density Plot
- **File**: `density_adj_pvalue_all_<condition1>_vs_<condition2>.png`
- **Description**: Compares the density distributions of `-log10(adjusted p-value)` for each method.
- **Highlight**: A vertical red dashed line indicates the significance threshold in `log10` scale.

### 2. Ranked Plot
- **File**: `ranked_adj_pvalue_lines_<condition1>_vs_<condition2>.png`
- **Description**: Plots proteins ranked by their significance (`-log10(adjusted p-value)`), comparing the statistical power and ranking behavior of each method.

## Requirements

- R (recommended version 4.0 or higher)  
- R packages: `ggplot2`, `dplyr`, `readr`, `stringr`, `future`, `janitor`, `tidyr`    
  
## Usage 
This script is part of a larger analysis pipeline. Make sure the following inputs are available before running:

* Data frames for each statistical method with:
  *    Adjusted p-values column (e.g., `adj_p_value`, `adj_p_val`, `adj_pvalue`)
  *    Protein identifiers (`protein_i_ds`, `majority_protein_i_ds`, or `protein_id` depending on method)
 
## Notes 
The script dynamically detects the correct protein ID and adjusted p-value columns.  
Warnings and debug messages guide the user if expected columns are missing.  
Parallel processing is disabled at the end of the analysis for clean exit.  
  
  
----
  
[3. Fold change analysis and protein significance overlap](Script3_BiologicalRelevance.R)

## Overview

This script performs an in-depth analysis of significant proteins derived from multiple statistical tests. It focuses on evaluating the consistency of protein significance across methods by examining intersections and differences, particularly regarding fold change (FC) filtering and Bayesian filtering. The goal is to assess the robustness and biological relevance of detected differential proteins.

## Features

- Processes results from various statistical tests (`tstudent`, `twelch`, `msstats`, `bayes`, etc.).
- Evaluates and classifies proteins based on their inclusion/exclusion after FC and Bayes filtering.
- Computes descriptive statistics for absolute log2 Fold Changes by group.
- Generates informative visualizations:
  - **UpSet plots** showing overlap and discarded protein intersections.
  - **Boxplots** comparing FC distributions across intersection groups.
- Performs statistical testing (ANOVA, Tukey HSD) to assess group differences in fold change.

## Outputs

The script generates several outputs for each statistical test analyzed:

- **`<W#>_BioRel_<test_name>_comparison.txt`**
A detailed text summary including:
  - FC statistics per group.
  - Results from ANOVA and Tukey HSD tests.
  - Protein classification counts.

- **`UpSet_discarded_significant_<test_name>.png`**
UpSet plot visualizing overlaps among discarded significant proteins due to:
  - Fold Change filtering.
  - Bayesian filtering.

- **`Boxplot_absFC_by_Group_<test_name>.png`**
Boxplot showing distribution of absolute log2 fold changes by intersection group:
  - Significant proteins retained by both filters.
  - Discarded proteins.
  - Other intersection-based classifications.

 - **`Discarded_Significant_Counts_<test_name>.txt`**
Table summarizing the number of significant proteins discarded by FC and Bayes filters.

## Requirements

- **R** version ≥ 4.0
- R packages: `dplyr`, `ggplot2`, `UpSetR`, `tibble`, `readr` (optional but recommended for file reading).  
  
Ensure your working directory is properly set up with:
- Filtered and unfiltered result files from multiple statistical tests.
- Biological relevance filtered files (`*_FC.txt`, `*_Bayes.txt`).
- Total protein list (`*_proteingroups_filtered.txt`).

## Usage

1. **Prepare Input Files**:
   Ensure result files are named with the following format:
   - `W<work_number>_res_<test_name>.txt`
   - `W<work_number>_biolrel_<test_name>_FC.txt`
   - `W<work_number>_biolrel_<test_name>_Bayes.txt`
   - `W<work_number>_proteingroups_filtered.txt`

2. **Define Parameters**:
   Modify `work_number`, `biol_rel_folder_path_current`, `biol_comparison_folder_path_current`, and `ht_results_folder_path_current` as needed in your script.

3. **Run the Script**:
   Execute the script in R or RStudio:
   
       source("your_script_name.R")

4. Review Outputs:
Navigate to your output folder to explore the result .txt summaries, UpSet plots, and boxplots.

## Notes
Proteins are grouped into intersection-based categories such as:

 * Sig_FC_Bayes: Significant in original results and retained by both filters.
 * Sig_Discarded_FC_Bayes: Significant in original but discarded by both FC and Bayes.
 * Others include combinations of retention and exclusion.
 * Statistical analysis requires at least 2 proteins per group. If group sizes are too small, ANOVA and Tukey HSD are skipped with a warning.
 * UpSet plots are generated only if at least 2 non-empty sets of discarded proteins are available.
 * The script handles missing or malformed files gracefully with informative messages to aid debugging.

----

[3. ClueGO File Preparation and Intersection Tables Generation](Script4_ClueGO_preparation.R)  

## Overview

This R script automates the preparation of input files for ClueGO analysis by separating UP and DOWN regulated proteins from multiple statistical method results. It also generates intersection tables of protein sets across these methods to assist downstream functional enrichment analysis.


## Input Files

For each dataset/work number W#, the script expects the following files organized under the base directory:

- `W#/W#_BioRel/` folder containing:
  - Files named as `W#_biolrel_<method>_(FC|Bayes).txt`
    - `<method>` can be `tstudent`, `twelch`, `limma`, `deqms`, or `msstats`.
  - `W#_resbiolrel_bayes.txt` (copied Bayes results)
  
- Optional backup Bayes results located at:
  - `W#/W#_H_testing/W#_res_bayes.txt`

Each input file must contain at least the columns:
- `Protein.IDs`
- `log2FC`


## Output Files

The script creates and writes output files under:

- `W#/W#_ClueGO/`

Output files include:

- UP and DOWN protein lists for each input file, named:
  - `W#_<method>_<criteria>_ClueGO_up.txt`
  - `W#_<method>_<criteria>_ClueGO_down.txt`

- Intersection tables for UP and DOWN protein sets by criterion (`FC` or `Bayes`):
  - `Intersection_UP_FC.txt`
  - `Intersection_UP_Bayes.txt`
  - `Intersection_DOWN_FC.txt`
  - `Intersection_DOWN_Bayes.txt`


## Workflow

1. **Setup Paths:** Defines paths for the current work number dataset.
2. **File Collection:** Gathers all relevant input files matching expected patterns.
3. **UP/DOWN Separation:** For each input file, proteins are separated based on `log2FC` value into UP (`>0`) and DOWN (`<0`) sets.
4. **Save ClueGO Files:** The separated protein lists are saved as new files in the ClueGO output directory.
5. **Intersection Tables:** Generates intersection binary tables comparing protein presence across methods for both UP and DOWN proteins, separately for FC and Bayes criteria.


## Usage

1. Set the base directory path:

       base_analysis_path <- "path/to/base_directory"

2. Select the work numbers (datasets) to process:
     
       work_numbers_to_process <- c(1, 2, 3, 4, 5)
     
4. Source or run the script:

       source("script4_cluego_processing.R")

## Requirements

- **R** version ≥ 4.0
- R packages: `purr`, `stringr`
  
## Notes  
Input files must have consistent naming conventions and contain the required columns.

The script assumes previous steps have prepared these input files correctly.

Duplicate or redundant files are internally handled via unique filtering.

Intersection tables generated provide a binary presence/absence matrix useful for comparative analyses.

Proper directory structure and file organization are essential for smooth execution.







# HTM-CRB

## Table of Contents

1. [Filtering and Visualization](#1-filtering-and-visualization)  
2. [Hypothesis Testing Pipeline for Differential Proteomics Analysis](#2-hypothesis-testing-pipeline-for-differential-proteomics-analysis)  
3. [Fold Change Analysis and Overlap of Significant Proteins](#3-fold-change-analysis-and-overlap-of-significant-proteins)  
4. [Preparation of Files for ClueGO and Intersection Tables Generation](#4-preparation-of-files-for-cluego-and-intersection-tables-generation)  
5. [ClueGO Analysis Pipeline](#5-cluego-analysis-pipeline)  
6. [Meta-Analysis of ClueGO Enrichment Results](#6-meta-analysis-of-cluego-enrichment-results)

---

## 1. Filtering and Visualization

### Summary  
This [script](Script1_FilteringAndVisualization.R) processes proteomics LFQ intensity data, generating visualizations to assess variability and differences between conditions.

### Requirements  
- R ≥ 4.0  
- Packages: `dplyr`, `tidyr`, `ggplot2`, `car`, `vegan`, `pheatmap`, `cowplot`, `showtext`  
*(The script installs missing packages automatically.)*

### Main Workflow  
- Load and prepare data (select LFQ columns, log2-transform, impute zeros).  
- Generate clustered heatmaps with condition annotations.  
- Perform NMDS and ANOSIM to evaluate group differences.  
- Process multiple comparisons in batch mode with error handling.  
- Save results and plots.

### Usage  
- Set the base path and condition patterns for each work.
- Run the script to process all defined works sequentially.
- Inspect generated heatmaps, NMDS plots, and ANOSIM results for exploratory insights into proteomic variability and condition effects.

### Notes

- The script is designed to be robust to missing data and zero intensities.
- Imputation uses a small value relative to the minimum detected intensity to avoid biases.
- Works with multiple experimental comparisons via a flexible batch execution loop.
- Subsequent scripts in the project continue downstream analyses and result integration.

---

## 2. Hypothesis Testing Pipeline for Differential Proteomics Analysis

### Summary  
The [script](Script2_HypothesisTesting.R) applies various statistical tests to detect differential proteins and generates diagnostic plots of adjusted p-values.

### Included Methods  
- Student’s t-test, Welch’s t-test, limma, DEqMS, MSstats

### Main Outputs  
- Density and ranking plots of -log10(adjusted p-values) per method.  
- High-resolution PNG files.

### Requirements  
- R ≥ 4.0  
- Packages: `ggplot2`, `dplyr`, `readr`, `stringr`, `future`, `janitor`, `tidyr`

### Usage 
This script is part of a larger analysis pipeline. Make sure the following inputs are available before running:

* Data frames for each statistical method with:
  *    Adjusted p-values column (e.g., `adj_p_value`, `adj_p_val`, `adj_pvalue`)
  *    Protein identifiers (`protein_i_ds`, `majority_protein_i_ds`, or `protein_id` depending on method)
 
### Notes 
- The script dynamically detects the correct protein ID and adjusted p-value columns.
- Warnings and debug messages guide the user if expected columns are missing.
- Parallel processing is disabled at the end of the analysis for clean exit. 


---

## 3. Fold Change Analysis and Overlap of Significant Proteins

### Summary  
This [script](Script3_BiologicalRelevance.R) analyzes the consistency of significant proteins across methods and filters (Fold Change and Bayesian), producing statistics and visualizations (UpSet plots, boxplots).

### Outputs  
- Statistical summary per method (text files).  
- UpSet plots and boxplots per method.  
- Tables counting discarded proteins.

### Requirements  
- R ≥ 4.0  
- Packages: `dplyr`, `ggplot2`, `UpSetR`, `tibble`, `readr`  

### Usage

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

4. **Review Outputs:**
Navigate to your output folder to explore the result .txt summaries, UpSet plots, and boxplots.

### Notes
Proteins are grouped into intersection-based categories such as:

 * Sig_FC_Bayes: Significant in original results and retained by both filters.
 * Sig_Discarded_FC_Bayes: Significant in original but discarded by both FC and Bayes.
 * Others include combinations of retention and exclusion.
 * Statistical analysis requires at least 2 proteins per group. If group sizes are too small, ANOVA and Tukey HSD are skipped with a warning.
 * UpSet plots are generated only if at least 2 non-empty sets of discarded proteins are available.
 * The script handles missing or malformed files gracefully with informative messages to aid debugging.

---

## 4. Preparation of Files for ClueGO and Intersection Tables Generation

### Summary  
The [script](Script4_ClueGO_preparation.R) separates up- and down-regulated proteins based on statistical results and generates intersection tables for downstream functional analysis.

### Inputs  
- Files with biologically relevant results (`*_FC.txt`, `*_Bayes.txt`)  
- Filtered protein files

### Outputs  
- UP/DOWN protein lists for ClueGO  
- Binary intersection tables by method and criterion

### Requirements  
- R ≥ 4.0  
- Packages: `purrr`, `stringr`

### Usage

1. Set the base directory path:

       base_analysis_path <- "path/to/base_directory"

2. Select the work numbers (datasets) to process:
     
       work_numbers_to_process <- c(1, 2, 3, 4, 5)
     
4. Source or run the script:

       source("script4_cluego_processing.R")

## Notes  
- Input files must have consistent naming conventions and contain the required columns.
- The script assumes previous steps have prepared these input files correctly.
- Duplicate or redundant files are internally handled via unique filtering.
- Intersection tables generated provide a binary presence/absence matrix useful for comparative analyses.

Proper directory structure and file organization are essential for smooth execution.

---

## 5. ClueGO Analysis Pipeline

### Summary  
This [script](Script5_EnrichmentComparisons.R) consolidates ClueGO results, performs statistical tests (Kruskal-Wallis and Dunn), and produces annotated visualizations.

### Requirements  
- R ≥ 4.0  
- Packages: `tidyverse`, `ggpubr`, `dunn.test`

### Outputs  
- Consolidated TSV files  
- Statistical results and PNG plots per job and overall

### Usage
- Set the required input and output folder paths.
- Define the lists of works, methods, and metric names as required.
- Run the script. It will:  
  * Consolidate `.xls` files into TSVs, filling missing data if necessary.
  * Process each consolidated TSV to compute metrics and transformations.
  * Aggregate results globally for combined statistical testing and visualization.
- Check the output folders for results and plots.

### Notes

- Ensure all input folders and files are correctly named and accessible.  
- The script is designed to be robust to missing or incomplete data, but consistent file naming is critical.  
- Adjust the lists of methods, works, and metrics as necessary to fit your data.  
- This script is the final stage of a multi-step ClueGO analysis pipeline, integrating results and providing comprehensive summary statistics and plots.

---

## 6. Meta-Analysis of ClueGO Enrichment Results

### Summary  
The [script](Script6_Metaanalysis.R) performs a comprehensive meta-analysis on consolidated ClueGO enrichment results, including global and ontology-specific tests, sensitivity analysis, and final visualizations.

### Requirements  
- R ≥ 4.0  
- Packages: `tidyverse`, `ggplot2`, `ggpubr`, `dunn.test`, `patchwork`, `readr`

### Outputs  
- Consolidated TSV files  
- Global and ontology-specific statistics  
- High-quality PNG figures  
- Sensitivity analysis results

### Usage
- Ensure all required raw ClueGO files are placed under `base_analysis_path`.
- Set all necessary configuration variables and lists prior to running the script.
- Run the script in an R environment with required packages installed (`tidyverse`, `ggplot2`, `ggpubr`, `dunn.test`, `patchwork`, `readr`).
- Check console logs for progress messages and warnings.
- Review outputs in the designated result folders.

---

## Contact

For questions or support, please contact [Your Name or Team] at [Alfonso Olaya](AlfonsoOA).

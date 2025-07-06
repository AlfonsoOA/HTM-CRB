# --- GLOBAL INITIAL CONFIGURATION AND FOLDER PREPARATION ---
# Load necessary libraries
# Make sure you have all these libraries installed. You can install them with:
# install.packages(c("dplyr", "tidyr", "ggplot2", "car", "limma", "DEqMS", "MSstats",
# "VennDiagram", "purrr", "future", "future.apply", "rstanarm", "posterior",
# "UpSetR", "matrixStats", "parallel", "janitor"))
# Some libraries like BiocGenerics, Biostrings, XVector, IRanges, S4Vectors, grid, futile.logger
# are dependencies of the above and are usually loaded automatically.

library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(limma)
library(DEqMS)
library(MSstats)
library(VennDiagram)
library(purrr)
library(future)
library(future.apply)
library(rstanarm)
library(posterior)
library(UpSetR)
library(parallel)
library(matrixStats)
library(janitor) # Added for clean_names()

# --- 1. MAIN CONFIGURATION (MODIFY FOR EACH WORK!) ---
# Define the base path for all projects (Work_1, Work_2, etc.)
# MAKE SURE THIS PATH EXISTS ON YOUR SYSTEM AND USE FORWARD SLASHES / OR DOUBLE BACKSLASHES \\!
base_analysis_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/7_Especial_IJMS/analisis27062025"

# Define the current work number and the names of the two conditions to compare.
# These names will be used for file naming and condition identification in results.
work_number <- 5
condition1_name <- "PM" # Example: "Cold", "PM", "TreatmentA"
condition2_name <- "C"   # Example: "RT", "Control", "TreatmentB"

# --- 2. REPLICATE DEFINITION FOR MSstats (MODIFY ACCORDING TO EXPERIMENTAL DESIGN!) ---
# Define the number of biological replicates for each condition.
num_bio_replicates_per_condition <- 4 # For Work 1: 3 biological replicates (Cold1, Cold2, Cold3)

# Define the number of technical replicates for each biological replicate.
num_tech_replicates_per_bio_replicate <- 1 # For Work 1: 2 technical replicates per biological replicate (e.g., Cold1_1, Cold1_2)

# --- END OF MAIN CONFIGURATION SECTION ---

# --- PATH DEFINITION AND FOLDER CREATION (ALL HERE, AT THE BEGINNING) ---

# Path of the main folder for the current work (e.g., W1, W5)
# Based on the new structure: .../analisis27062025/W1/
current_main_folder_path <- file.path(base_analysis_path, paste0("W", work_number))

# Path of the folder where the filtered file is located (from phase 1 of the pipeline)
# Based on the new structure: .../W1/1_ProteinGroups_Analysis_W1/
current_filtered_folder <- file.path(current_main_folder_path, paste0("1_ProteinGroups_Analysis_W", work_number))

# Paths of the main input files
# Filtered file: 1_proteingroups_filtered_W1.txt
current_input_filtered_file <- file.path(current_filtered_folder, paste0("1_proteingroups_filtered_W", work_number, ".txt"))

# For MSstats, we need the original proteinGroups.txt and evidence.txt files from MaxQuant
# Assuming they are directly inside the main work folder (e.g., W1)
current_original_protein_groups_file <- file.path(current_main_folder_path, paste0("W", work_number, "_proteinGroups.txt"))
current_original_evidence_file <- file.path(current_main_folder_path, paste0("W", work_number, "_evidence.txt"))

# Create the folder for hypothesis testing results (e.g., W1_H_testing)
current_results_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_H_testing"))
if (!dir.exists(current_results_folder_path)) {
  dir.create(current_results_folder_path, recursive = TRUE)
  cat(paste("Results folder for Work", work_number, "created at:", current_results_folder_path, "\n"))
} else {
  cat(paste("The results folder for Work", work_number, "already exists at:", current_results_folder_path, "\n"))
}

# Output file paths for each analysis method
current_output_tStudent_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_tstudent.txt"))
current_output_tWelch_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_twelch.txt"))
current_output_limma_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_limma.txt"))
current_output_DEqMS_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_deqms.txt"))
current_output_MSStats_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_msstats.txt"))
current_output_Bayes_file <- file.path(current_results_folder_path, paste0("W", work_number, "_res_bayes.txt")) # Not yet implemented in this script

# --- GENERATING AND WRITING THE ANNOTATION FILE FOR MSstats ---

# Define simplified condition names for MSstats 'Condition' column and 'Run' column.
# These will be 'u' and 'y', matching the LFQ.intensity column prefixes in your MaxQuant output.
# based on your MaxQuant LFQ column names (e.g., LFQ.intensity.u1).
# Example: If condition1_name is "u", msstats_run_prefix1 will be "u".
msstats_run_prefix1 <- condition1_name 
msstats_run_prefix2 <- condition2_name 

# Annotation file specific for MSstats
current_msstats_annotation_file <- file.path(current_main_folder_path, paste0("annotation_work", work_number, ".csv"))

cat(paste("\n### Generating and writing annotation file for MSstats (Work", work_number, ":", msstats_run_prefix1, "vs", msstats_run_prefix2, ") ###\n"))

# Calculate total number of 'runs' per condition
# These numbers are based on the global num_bio_replicates_per_condition and num_tech_replicates_per_bio_replicate.
# Example: If num_bio_replicates_per_condition = 4 and num_tech_replicates_per_bio_replicate = 1, num_runs_cond1 will be 4.
num_runs_cond1 <- num_bio_replicates_per_condition * num_tech_replicates_per_bio_replicate
num_runs_cond2 <- num_bio_replicates_per_condition * num_tech_replicates_per_bio_replicate


# --- 1. Create RAW file names ('Run' column) ---
# IMPORTANT! These must EXACTLY match the experiment names inferred by MaxQuant/MSstats
# from your LFQ.intensity columns (e.g., "u1", "u2", "u3", "u4" for your current Work 2).
# Modification: Removed the "_1" suffix from the original paste0.
# Example: "u" + "1" becomes "u1".
runs_cond1 <- paste0(msstats_run_prefix1, 1:num_bio_replicates_per_condition)
runs_cond2 <- paste0(msstats_run_prefix2, 1:num_bio_replicates_per_condition)
all_runs <- c(runs_cond1, runs_cond2)


# --- 2. Define conditions ('Condition' column) ---
# Use the simplified condition names (e.g., "u", "y") for this column as expected by MSstats
# for defining contrasts.
# Modification: Changed from using 'msstats_condition1_name' (e.g., "u_Cond") to 'msstats_run_prefix1' (e.g., "u").
conditions_for_runs <- c(rep(msstats_run_prefix1, num_runs_cond1),
                         rep(msstats_run_prefix2, num_runs_cond2))

# --- 3. Define biological replicate IDs ('BioReplicate' column) ---
# THESE MUST BE UNIQUE ACROSS THE ENTIRE EXPERIMENT!
# Using the simplified prefix for consistency in BioReplicate naming.
# Example: "u_BioRep1", "u_BioRep2", etc.
# Modification: Changed from 'condition1_name' (e.g., "u_Cond") to 'msstats_run_prefix1' (e.g., "u").
bioreplicates_cond1 <- paste0(msstats_run_prefix1, "_BioRep", 1:num_bio_replicates_per_condition)
bioreplicates_cond2 <- paste0(msstats_run_prefix2, "_BioRep", 1:num_bio_replicates_per_condition)
all_bioreplicates <- c(bioreplicates_cond1, bioreplicates_cond2)


# --- 4. Define technical replicate IDs ('TechReplicate' column) ---
# If each BioReplicate has more than one 'run', TechReplicate must differentiate them.
# For Work 2, num_tech_replicates_per_bio_replicate is 1, so all values in this column will be 1.
all_techreplicates <- rep(1:num_tech_replicates_per_bio_replicate, times = num_bio_replicates_per_condition)


# --- 5. Create the final annotation dataframe ---
# Column names are CRUCIAL and must be 'Run', 'Condition', 'BioReplicate', 'TechReplicate'!
msstats_annotation_df <- data.frame(
  Run = all_runs,
  Condition = conditions_for_runs,
  BioReplicate = all_bioreplicates,
  TechReplicate = all_techreplicates,
  stringsAsFactors = FALSE # Important to avoid issues with text strings
)

# --- 6. Optional: Print the dataframe to verify its content ---
cat("\nDEBUG: Content of 'msstats_annotation_df' dataframe before writing:\n")
print(msstats_annotation_df)

# --- 7. Write the dataframe to the annotation file ---
write.table(msstats_annotation_df,
            file = current_msstats_annotation_file,
            sep = ",",         # Use comma as separator since the file is .csv
            row.names = FALSE, # DO NOT include R row numbers
            quote = FALSE,     # DO NOT put quotes around text values
            col.names = TRUE)  # Include column names as header

cat(paste("Annotation file for MSstats (Work", work_number, ") created at:", current_msstats_annotation_file, "\n"))

# --- Global Parallelization Setup (for rstanarm, not used here, but good practice) ---
num_cores_to_use <- parallel::detectCores() - 1
if (num_cores_to_use < 1) num_cores_to_use <- 1 # Ensure at least 1 core
plan(multisession, workers = num_cores_to_use)
cat(paste("Parallelization configured with", num_cores_to_use, "cores for Bayesian analysis.\n"))

# --- 1. Read the filtered input file from the previous phase ---
cat(paste("Attempting to read filtered file for Work", work_number, "at:\n", current_input_filtered_file, "\n"))
if (!file.exists(current_input_filtered_file)) {
  stop(paste("Error: Filtered protein file not found for Work", work_number, "at:", current_input_filtered_file))
}
protein_groups_filtered_current <- read.delim(current_input_filtered_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cat("\nDEBUG: Full column names of protein_groups_filtered_current:\n")
print(colnames(protein_groups_filtered_current))
cat("First rows of filtered protein file for Work", work_number, ":\n")
print(head(protein_groups_filtered_current[, 1:min(ncol(protein_groups_filtered_current), 10)]))

# --- Auxiliary Functions ---
# Function to expand rows with multiple IDs (useful for various analyses)
expand_rows_ids <- function(df, id_col, other_cols) {
  if (nrow(df) == 0) {
    all_expected_cols <- unique(c("Protein.ID", other_cols))
    return(data.frame(matrix(ncol = length(all_expected_cols), nrow = 0, dimnames = list(NULL, all_expected_cols)),
                      stringsAsFactors = FALSE))
  }
  
  df[[id_col]] <- as.character(df[[id_col]])
  
  split_ids <- strsplit(df[[id_col]], ";")
  n_ids <- sapply(split_ids, length)
  
  if (length(other_cols) == 0 || all(other_cols %in% id_col)) {
    expanded_df <- data.frame(
      Protein.ID = unlist(split_ids),
      stringsAsFactors = FALSE
    )
  } else {
    expanded_df <- data.frame(
      Protein.ID = unlist(split_ids),
      df[rep(1:nrow(df), n_ids), other_cols, drop = FALSE],
      stringsAsFactors = FALSE
    )
  }
  return(expanded_df)
}

# Function to perform t-test (Student or Welch) per protein WITH LOG2
run_ttest_log2 <- function(row, cond1_cols, cond2_cols, var.equal = TRUE, alpha = 0.05) {
  cond1_intensities <- as.numeric(row[cond1_cols])
  cond2_intensities <- as.numeric(row[cond2_cols])
  
  # Log2 transform intensities. Zeros or NAs will become -Inf or NA, respectively.
  # t.test handles NAs by removing them. -Inf will also be removed by is.finite.
  cond1_valid_log2 <- log2(cond1_intensities)
  cond2_valid_log2 <- log2(cond2_intensities)
  
  # Filter NAs and infinities that may result from non-numeric values or log2(0)
  cond1_valid_log2 <- cond1_valid_log2[is.finite(cond1_valid_log2) & !is.na(cond1_valid_log2)]
  cond2_valid_log2 <- cond2_valid_log2[is.finite(cond2_valid_log2) & !is.na(cond2_valid_log2)]
  
  result <- data.frame(p.value = NA, log2FC = NA, significant = FALSE)
  
  # Perform t-test only if there are enough valid data in both groups
  if (length(cond1_valid_log2) >= 2 && length(cond2_valid_log2) >= 2) {
    ttest_result <- t.test(cond1_valid_log2, cond2_valid_log2, var.equal = var.equal)
    mean_cond1_log2 <- mean(cond1_valid_log2, na.rm = TRUE)
    mean_cond2_log2 <- mean(cond2_valid_log2, na.rm = TRUE)
    log2fc <- mean_cond1_log2 - mean_cond2_log2 # Log2 FC (Cond1 / Cond2)
    result$p.value <- ttest_result$p.value
    result$log2FC <- log2fc
    result$significant <- result$p.value < alpha
  }
  return(result)
}

# --- Dynamic identification of LFQ intensity columns ---
# Extraer el prefijo de una sola letra de los nombres de las condiciones (ej., 'u' de 'u_Cond')
cond1_prefix <- sub("_Cond$", "", condition1_name) # Esto convertirá "u_Cond" en "u"
cond2_prefix <- sub("_Cond$", "", condition2_name) # Esto convertirá "y_Cond" en "y"

# REGEX: Modificado para buscar "LFQ.intensity." seguido del prefijo de la condición y uno o más dígitos.
lfq_cols_cond1 <- grep(paste0("^LFQ\\.intensity\\.", cond1_prefix, "\\d+$"), colnames(protein_groups_filtered_current), value = TRUE)
lfq_cols_cond2 <- grep(paste0("^LFQ\\.intensity\\.", cond2_prefix, "\\d+$"), colnames(protein_groups_filtered_current), value = TRUE)

if (length(lfq_cols_cond1) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition1_name, "'. Check regex or column names."))
if (length(lfq_cols_cond2) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition2_name, "'. Check regex or column names."))

cat(paste0("\nLFQ intensity columns identified for condition '", condition1_name, "': ", paste(lfq_cols_cond1, collapse = ", "), "\n"))
cat(paste0("LFQ intensity columns identified for condition '", condition2_name, "': ", paste(lfq_cols_cond2, collapse = ", "), "\n"))

# --- ANALYSIS WITH T-TEST (Student and Welch) ---
cat(paste("\n### Analysis with t-test (Work", work_number, ":", condition1_name, "vs", condition2_name, ") ###\n"))

# Apply Student's t-test
cat(paste0("\nApplying Student's t-test (with log2 transformation, ", condition1_name, " vs ", condition2_name, ")...\n"))
student_results_list_log2 <- apply(protein_groups_filtered_current, 1, run_ttest_log2,
                                   cond1_cols = lfq_cols_cond1, cond2_cols = lfq_cols_cond2,
                                   var.equal = TRUE)
student_results_df_log2 <- do.call(rbind, student_results_list_log2)

# Join with the original filtered dataframe to keep all columns
student_output_full <- protein_groups_filtered_current %>%
  select(Protein.IDs, everything()) %>%
  bind_cols(student_results_df_log2)
student_output_full$adj.p.value <- p.adjust(student_output_full$p.value, method = "BH")

# Reorder columns (statistics first, then the rest)
cols_to_move_student <- c("p.value", "log2FC", "significant", "adj.p.value")
first_cols_student <- c("Protein.IDs", cols_to_move_student)
other_cols_student <- setdiff(colnames(student_output_full), first_cols_student)
student_output_reordered <- student_output_full[, c(first_cols_student, other_cols_student)]

# Filter by non-NA adj.p.value and then by significance
sig_student_df <- student_output_reordered %>%
  filter(!is.na(adj.p.value) & adj.p.value < 0.05)

if (dir.exists(current_results_folder_path)) {
  write.table(sig_student_df, file = current_output_tStudent_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significantly different proteins (Student's t-test, log2, adj.p < 0.05) saved to:", current_output_tStudent_file, "\n"))
  print(head(sig_student_df[, 1:min(ncol(sig_student_df), 10)]))
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

# Apply Welch's t-test
cat(paste0("\nApplying Welch's t-test (with log2 transformation, ", condition1_name, " vs ", condition2_name, ")...\n"))
welch_results_list_log2 <- apply(protein_groups_filtered_current, 1, run_ttest_log2,
                                 cond1_cols = lfq_cols_cond1, cond2_cols = lfq_cols_cond2,
                                 var.equal = FALSE)
welch_results_df_log2 <- do.call(rbind, welch_results_list_log2)

# Join with the original filtered dataframe
welch_output_full <- protein_groups_filtered_current %>%
  select(Protein.IDs, everything()) %>%
  bind_cols(welch_results_df_log2)
welch_output_full$adj.p.value <- p.adjust(welch_output_full$p.value, method = "BH")

# Reorder columns
cols_to_move_welch <- c("p.value", "log2FC", "significant", "adj.p.value")
first_cols_welch <- c("Protein.IDs", cols_to_move_welch)
other_cols_welch <- setdiff(colnames(welch_output_full), first_cols_welch)
welch_output_reordered <- welch_output_full[, c(first_cols_welch, other_cols_welch)]

# Filter by non-NA adj.p.value and then by significance
sig_welch_df <- welch_output_reordered %>%
  filter(!is.na(adj.p.value) & adj.p.value < 0.05)

if (dir.exists(current_results_folder_path)) {
  write.table(sig_welch_df, file = current_output_tWelch_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significantly different proteins (Welch's t-test, log2, adj.p < 0.05) saved to:", current_output_tWelch_file, "\n"))
  print(head(sig_welch_df[, 1:min(ncol(sig_welch_df), 10)]))
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

# --- ANALYSIS WITH limma ---
cat(paste("\n### Analysis with limma (Work", work_number, ":", condition1_name, "vs", condition2_name, ") ###\n"))

lfq_cols <- c(lfq_cols_cond1, lfq_cols_cond2)
if (length(lfq_cols_cond1) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition1_name, "' for limma. Check regex or column names."))
if (length(lfq_cols_cond2) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition2_name, "' for limma. Check regex or column names."))
if (length(lfq_cols_cond1) < 2 || length(lfq_cols_cond2) < 2) stop("Error: Not enough replicates in one or both groups for limma. At least 2 are needed.")

# Create the expression matrix from LFQ intensity values transformed to log2
expression_matrix <- as.matrix(protein_groups_filtered_current[, lfq_cols])
# Handle zeros or values <= 0: log2(0) produces -Inf. lmFit will treat -Inf as NA.
# No epsilon imputation.
expression_matrix_log2 <- log2(expression_matrix)

# Set row names of the expression matrix with Protein.IDs
rownames(expression_matrix_log2) <- protein_groups_filtered_current$Protein.IDs

# Create the design matrix
group <- factor(c(rep(condition1_name, length(lfq_cols_cond1)), rep(condition2_name, length(lfq_cols_cond2))))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit the linear model
fit <- lmFit(expression_matrix_log2, design)

# Define the contrast of interest
contrast.matrix <- makeContrasts(paste0(condition1_name, " - ", condition2_name), levels = design)

# Apply the contrast and Bayes correction
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
limma_results <- topTable(fit2, coef = paste0(condition1_name, " - ", condition2_name), number = Inf)

# Filter by adj.P.Val and then merge with original data
limma_results_filtered_na <- limma_results[!is.na(limma_results$adj.P.Val), ]
sig_limma_df <- limma_results_filtered_na[limma_results_filtered_na$adj.P.Val < 0.05, ]
significant_protein_ids_limma <- rownames(sig_limma_df)

# Merge limma results with the original filtered dataframe columns
original_data_for_merge_limma <- protein_groups_filtered_current %>%
  select(Protein.IDs, all_of(lfq_cols_cond1), all_of(lfq_cols_cond2), everything()) # Select important columns first

sig_limma_output <- merge(
  data.frame(Protein.IDs = significant_protein_ids_limma, sig_limma_df),
  original_data_for_merge_limma,
  by = "Protein.IDs",
  all.x = TRUE # Keep only proteins that were significant in limma
)

# Reorder columns for better visualization
limma_specific_cols <- c("logFC", "P.Value", "adj.P.Val", "AveExpr", "t", "B")
# Get the rest of the original columns excluding those already in limma_specific_cols or are IDs
other_original_cols_limma <- setdiff(colnames(protein_groups_filtered_current), c("Protein.IDs", lfq_cols_cond1, lfq_cols_cond2))
final_limma_cols_order <- c("Protein.IDs", limma_specific_cols, lfq_cols_cond1, lfq_cols_cond2, other_original_cols_limma)
# Ensure that only columns that actually exist in sig_limma_output are selected
sig_limma_output <- sig_limma_output[, intersect(final_limma_cols_order, colnames(sig_limma_output))]

if (dir.exists(current_results_folder_path)) {
  write.table(sig_limma_output, file = current_output_limma_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significantly different proteins (limma, adj.P.Val < 0.05) saved to:", current_output_limma_file, "\n"))
  print(head(sig_limma_output[, 1:min(ncol(sig_limma_output), 10)]))
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

# --- ANALYSIS WITH DEqMS ---
cat(paste("\n### Analysis with DEqMS (Work", work_number, ":", condition1_name, "vs", condition2_name, ") ###\n"))

cat("First rows of filtered protein file for DEqMS (Work", work_number, ":", condition1_name, "vs", condition2_name, "):\n")
head(protein_groups_filtered_current[, 1:min(ncol(protein_groups_filtered_current), 10)])

lfq_cols_deqms <- c(lfq_cols_cond1, lfq_cols_cond2)
if (length(lfq_cols_cond1) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition1_name, "' for DEqMS. Check regex or column names."))
if (length(lfq_cols_cond2) == 0) stop(paste0("Error: No LFQ intensity columns found for condition '", condition2_name, "' for DEqMS. Check regex or column names."))

# Prepare the intensity matrix
df.LFQ <- protein_groups_filtered_current[, lfq_cols_deqms, drop = FALSE]
rownames(df.LFQ) <- protein_groups_filtered_current$Protein.IDs

# Filter proteins that have too many NAs/zeros in both groups
na_or_zero_count_cond1 <- apply(df.LFQ[, lfq_cols_cond1, drop = FALSE], 1, function(x) sum(is.na(x) | x == 0))
na_or_zero_count_cond2 <- apply(df.LFQ[, lfq_cols_cond2, drop = FALSE], 1, function(x) sum(is.na(x) | x == 0))

# At least 2 non-NA/zero values per condition are required for a t-test (internal to limma/DEqMS)
# Adjust this number according to your data if necessary (e.g., >0 to have at least one)
df.LFQ.filter <- df.LFQ[na_or_zero_count_cond1 <= (length(lfq_cols_cond1) - 2) &
                          na_or_zero_count_cond2 <= (length(lfq_cols_cond2) - 2), , drop = FALSE]

if (nrow(df.LFQ.filter) == 0) stop("Error: No proteins remaining after filtering for NAs/zeros for DEqMS. Consider relaxing the filter.")

# Prepare the peptide count table
cat("\n### Preparing peptide count table for DEqMS ###\n")

# Attempt to find condition-specific peptide count columns first.
# Example column names: Razor...unique.peptides.Cold1_1, Razor...unique.peptides.RT3_2
# The corrected regex now includes the three escaped dots (\\.\\.\\.)
# and the technical replicate part (\\d+_\\d+).
peptide_count_cols_pattern <- paste0("^Razor\\.\\.\\.unique\\.peptides\\.(", condition1_name, "|", condition2_name, ")\\d+_\\d+$")
peptide_count_cols <- grep(peptide_count_cols_pattern, colnames(protein_groups_filtered_current), value = TRUE)

if (length(peptide_count_cols) == 0) {
  # If no condition-specific columns are found, look for the general "Razor.unique.peptides" column.
  # Note: Based on your input, condition-specific columns are expected,
  # but this is a common fallback in some scripts.
  general_peptide_count_col <- grep("^Razor\\.\\.\\.unique\\.peptides$", colnames(protein_groups_filtered_current), value = TRUE)
  if (length(general_peptide_count_col) == 1) {
    peptide_count_cols <- general_peptide_count_col
    cat("Warning: No condition-specific peptide count columns found for DEqMS. Using general 'Razor...unique.peptides' column.\n")
  } else {
    stop("Error: No peptide count columns found (neither condition-specific nor general). Check regex or column names.")
  }
} else {
  cat(paste0("Found ", length(peptide_count_cols), " condition-specific peptide count columns.\n"))
}

# If peptide count columns are found, proceed.
if (length(peptide_count_cols) > 0) {
  # If a general column was found (should be a single column name)
  if (length(peptide_count_cols) == 1 && grepl("^Razor\\.\\.\\.unique\\.peptides$", peptide_count_cols)) {
    pep.count.table <- data.frame(count = protein_groups_filtered_current[[peptide_count_cols]],
                                  row.names = protein_groups_filtered_current$Protein.IDs)
  } else {
    # If condition-specific columns are found, use the minimum of those.
    # Ensure they are numeric for rowMins.
    peptide_count_matrix <- as.matrix(protein_groups_filtered_current[, peptide_count_cols, drop = FALSE])
    mode(peptide_count_matrix) <- 'numeric' # Ensure numeric type for calculations
    pep.count.table <- data.frame(count = rowMins(peptide_count_matrix, na.rm = TRUE), # Use the minimum across conditions
                                  row.names = protein_groups_filtered_current$Protein.IDs)
  }
  cat("Peptide count table successfully created.\n")
  print(head(pep.count.table))
} else {
  stop("Error: No peptide count columns found after processing logic. This should not happen if the previous checks passed.")
}


# DEqMS requires peptide count to be > 0, so 1 is added
pep.count.table$count <- pep.count.table$count + 1

# Ensure that only proteins present in both tables are included
common_proteins_deqms <- intersect(rownames(df.LFQ.filter), rownames(pep.count.table))
if (length(common_proteins_deqms) == 0) stop("Error: No common proteins between filtered LFQ dataframe and peptide count table for DEqMS.")

# Prepare expression matrix (log2) and filtered peptide count table
# Log2 transformation will result in -Inf for 0 values; limma handles these as NA.
protein.matrix <- log2(as.matrix(df.LFQ.filter[common_proteins_deqms, , drop = FALSE]))
pep.count.table_filtered <- pep.count.table[common_proteins_deqms, , drop = FALSE]

# Create design matrix for limma
num_cond1 <- length(lfq_cols_cond1)
num_cond2 <- length(lfq_cols_cond2)
class <- as.factor(c(rep(condition1_name, num_cond1), rep(condition2_name, num_cond2)))
design <- model.matrix(~0 + class)
colnames(design) <- levels(class)

# Fit the limma model
fit1_deqms <- lmFit(protein.matrix, design = design)
cont_deqms <- makeContrasts(paste0(condition1_name, " - ", condition2_name), levels = design)
fit2_deqats <- contrasts.fit(fit1_deqms, contrasts = cont_deqms)
fit3_deqms <- eBayes(fit2_deqats)

# Add peptide count to the limma object
fit3_deqms$count <- pep.count.table_filtered$count
# Remove possible NAs in peptide count (even though +1 was added, for safety)
fit3_deqms <- fit3_deqms[!is.na(fit3_deqms$count), ]
if (nrow(fit3_deqms) == 0) stop("Error: No proteins remaining after adding peptide count for DEqMS.")
fit4_deqms <- spectraCounteBayes(fit3_deqms)
DEqMS.results <- outputResult(fit4_deqms, coef_col = 1)
# Merge DEqMS results with the original filtered dataframe columns
original_data_for_merge_deqms <- protein_groups_filtered_current %>%
  select(Protein.IDs, all_of(lfq_cols_cond1), all_of(lfq_cols_cond2), everything())
deqms_output <- merge(
  data.frame(Protein.IDs = rownames(DEqMS.results), DEqMS.results),
  original_data_for_merge_deqms,
  by = "Protein.IDs",
  all.x = TRUE
)
deqms_specific_cols <- c("logFC", "P.Value", "adj.P.Val", "meanSd", "t", "df", "B", "count")
other_original_cols_deqms <- setdiff(colnames(protein_groups_filtered_current), c("Protein.IDs", lfq_cols_cond1, lfq_cols_cond2))
final_deqms_cols_order <- c("Protein.IDs", deqms_specific_cols, lfq_cols_cond1, lfq_cols_cond2, other_original_cols_deqms)
deqms_output <- deqms_output[, intersect(final_deqms_cols_order, colnames(deqms_output))]
sig_deqms_df <- deqms_output %>%
  filter(!is.na(adj.P.Val) & adj.P.Val <= 0.05)
if (dir.exists(current_results_folder_path)) {
  write.table(sig_deqms_df, file = current_output_DEqMS_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significant DEqMS results (log2 FC, adj.P.Val <= 0.05) saved to:", current_output_DEqMS_file, "\n"))
  head(sig_deqms_df[, 1:min(ncol(sig_deqms_df), 10)])
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

# --- ANALYSIS WITH MSstats ---
cat(paste("\n### Analysis with MSstats (Work", work_number, ":", condition1_name, "vs", condition2_name, ") ###\n"))

# NOTE: The generation and writing of 'msstats_annotation_df'
# and the variables 'msstats_condition1_name' and 'msstats_condition2_name'
# should be performed BEFORE this block, as we have reorganized.
# Therefore, they are not included here.

# MaxQuant input files (the original ones, not the filtered ones, as required by MSstats)
if (!file.exists(current_original_protein_groups_file)) {
  stop(paste("Error: Original proteinGroups.txt file not found at:", current_original_protein_groups_file))
}
if (!file.exists(current_original_evidence_file)) {
  stop(paste("Error: Original evidence.txt file not found at:", current_original_evidence_file))
}

protein_groups_df_original <- read.delim(current_original_protein_groups_file, sep = "\t", header = TRUE)
evidence_df_original <- read.delim(current_original_evidence_file, sep = "\t", header = TRUE)

# --- KEY TO THE SOLUTION: Use the 'msstats_annotation_df' dataframe directly ---
# Pass the 'msstats_annotation_df' dataframe directly to MaxQtoMSstatsFormat
# This assumes 'msstats_annotation_df' is defined and is the correct dataframe.
msstats_raw_data <- MaxQtoMSstatsFormat(
  evidence = evidence_df_original,
  annotation = msstats_annotation_df, # <--- PASS THE DATAFRAME DIRECTLY!
  proteinGroups = protein_groups_df_original,
  mergeFeatures = FALSE,
  featureColumns = "Sequence" # MSstats will typically look for "Intensity" or "LFQ intensity" columns by default
)

msstats_processed_data <- dataProcess(msstats_raw_data, censoredInt = 'NA')

# Define the contrast matrix
# Use the modified condition variables for MSstats (e.g., PM_Cond vs C_Cond)
# Define the contrast matrix for MSstats groupComparison
# The column names of this matrix MUST exactly match the condition levels
# as processed internally by MSstats (i.e., "u" and "y").
comparison_msstats_current <- matrix(c(1, -1), nrow = 1, byrow = TRUE,
                                     dimnames = list(paste0(msstats_run_prefix1, "_vs_", msstats_run_prefix2),
                                                     c(msstats_run_prefix1, msstats_run_prefix2))) 
msstats_results_data <- groupComparison(contrast.matrix = comparison_msstats_current, data = msstats_processed_data)

all_msstats_results <- msstats_results_data$ComparisonResult
msstats_expanded <- expand_rows_ids(all_msstats_results, "Protein", setdiff(colnames(all_msstats_results), "Protein"))
msstats_clean <- msstats_expanded[!grepl("^CON__", msstats_expanded$Protein.ID), ]

cat("\n--- MSstats Diagnosis: Before significance filter ---\n")
cat("Total number of proteins in msstats_clean (after expansion and contaminant removal):", nrow(msstats_clean), "\n")
cat("Column names in msstats_clean:\n")
print(colnames(msstats_clean))

if ("log2FC" %in% colnames(msstats_clean) && "pvalue" %in% colnames(msstats_clean) && "adj.pvalue" %in% colnames(msstats_clean)) {
  cat("\nFirst 10 rows of msstats_clean (Protein.ID, log2FC, pvalue, adj.pvalue):\n")
  print(head(msstats_clean[, c("Protein.ID", "log2FC", "pvalue", "adj.pvalue")], 10))
  
  na_logFC_count <- sum(is.na(msstats_clean$log2FC))
  cat(paste0("\nNumber of proteins with log2FC = NA: ", na_logFC_count, " (", round(na_logFC_count/nrow(msstats_clean)*100, 2), "%)\n"))
  inf_logFC_count <- sum(is.infinite(msstats_clean$log2FC))
  cat(paste0("Number of proteins with log2FC = Inf/-Inf: ", inf_logFC_count, " (", round(inf_logFC_count/nrow(msstats_clean)*100, 2), "%)\n"))
  zero_adj_pvalue_count <- sum(msstats_clean$adj.pvalue == 0, na.rm = TRUE)
  cat(paste0("Number of proteins with adj.pvalue = 0: ", zero_adj_pvalue_count, "\n"))
  sig_na_logFC_count <- sum(msstats_clean$adj.pvalue < 0.05 & is.na(msstats_clean$log2FC), na.rm = TRUE)
  cat(paste0("Number of proteins with adj.pvalue < 0.05 AND log2FC = NA: ", sig_na_logFC_count, "\n"))
  sig_inf_logFC_count <- sum(msstats_clean$adj.pvalue < 0.05 & is.infinite(msstats_clean$log2FC), na.rm = TRUE)
  cat(paste0("Number of proteins with adj.pvalue < 0.05 AND log2FC = Inf/-Inf: ", sig_inf_logFC_count, "\n"))
  current_filter_count <- sum(msstats_clean$adj.pvalue < 0.05 & !is.na(msstats_clean$log2FC), na.rm = TRUE)
  cat(paste0("Number of proteins that would meet the current filter (adj.p < 0.05 & !is.na(log2FC)): ", current_filter_count, "\n"))
} else {
  cat("Warning: Not all expected columns (Protein.ID, log2FC, pvalue, adj.pvalue) are present in msstats_clean. Detailed diagnosis cannot be performed.\n")
  cat("Current columns in msstats_clean: ", paste(colnames(msstats_clean), collapse = ", "), "\n")
}
cat("----------------------------------------------------------------\n")

if (nrow(msstats_clean) > 0 && "log2FC" %in% colnames(msstats_clean) && "adj.pvalue" %in% colnames(msstats_clean)) {
  sig_msstats_df_only_stats <- msstats_clean %>%
    filter(adj.pvalue < 0.05 & !is.na(log2FC))
} else {
  cat("Warning: Not enough proteins found in msstats_clean or crucial columns ('log2FC', 'adj.pvalue') are missing. An empty MSstats results file will be generated.\n")
  sig_msstats_df_only_stats <- data.frame(
    Protein.ID = character(0), Label = character(0), log2FC = numeric(0), se = numeric(0),
    Tvalue = numeric(0), DF = numeric(0), pvalue = numeric(0), adj.pvalue = numeric(0),
    stringsAsFactors = FALSE
  )
}

# Prepare the original dataframe with desired columns for merging
# Assume protein_groups_filtered_current is already defined and available.
original_data_for_merge_msstats <- protein_groups_filtered_current %>%
  select(Protein.IDs, all_of(lfq_cols_cond1), all_of(lfq_cols_cond2), everything())

# Expand protein_groups_filtered_current
original_df_other_cols_msstats <- setdiff(colnames(protein_groups_filtered_current), "Protein.IDs")
protein_groups_expanded_for_merge <- expand_rows_ids(
  protein_groups_filtered_current,
  id_col = "Protein.IDs",
  other_cols = original_df_other_cols_msstats
)
if ("Protein.ID" %in% colnames(protein_groups_expanded_for_merge) && !"Protein.IDs" %in% colnames(protein_groups_expanded_for_merge)) {
  colnames(protein_groups_expanded_for_merge)[colnames(protein_groups_expanded_for_merge) == "Protein.ID"] <- "Protein.IDs"
}

# Perform the merge of MSstats results with the EXPANDED original dataframe
msstats_output_merged <- merge(
  sig_msstats_df_only_stats,
  protein_groups_expanded_for_merge,
  by.x = "Protein.ID",
  by.y = "Protein.IDs",
  all = FALSE
)

msstats_specific_cols <- c("Protein.ID", "Label", "log2FC", "se", "Tvalue", "DF", "pvalue", "adj.pvalue")
# Adjusted to find LFQ.intensity columns for merging
lfq_cols_merged <- grep(paste0("^LFQ\\.intensity\\.", condition1_name, "[0-9]+$|^LFQ\\.intensity\\.", condition2_name, "[0-9]+$"), colnames(msstats_output_merged), value = TRUE)
all_lfq_cols <- unique(lfq_cols_merged)

other_original_cols_msstats_final <- setdiff(colnames(msstats_output_merged), c(msstats_specific_cols, all_lfq_cols))
final_msstats_cols_order <- c(msstats_specific_cols, all_lfq_cols, other_original_cols_msstats_final)

sig_msstats_df_final <- msstats_output_merged[, intersect(final_msstats_cols_order, colnames(msstats_output_merged))]

if (dir.exists(current_results_folder_path)) {
  write.table(sig_msstats_df_final, file = current_output_MSStats_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significant, clean, and expanded MSstats results (adj.p < 0.05) saved to:", current_output_MSStats_file, "\n"))
  cat(paste("Number of significant proteins (MSstats, clean, expanded, adj.p < 0.05) AND present in filtered file:", nrow(sig_msstats_df_final), "\n"))
  print(head(sig_msstats_df_final[, 1:min(10, ncol(sig_msstats_df_final))]))
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

# --- BAYESIAN ANALYSIS WITH rstanarm ---
# Find the overall minimum non-zero LFQ intensity for imputation
# This value will be used to replace 0 or NA intensities before log2 transformation.
all_lfq_values <- as.numeric(unlist(protein_groups_filtered_current[, c(lfq_cols_cond1, lfq_cols_cond2), drop = FALSE]))
min_non_zero_lfq <- min(all_lfq_values[all_lfq_values > 0 & is.finite(all_lfq_values)], na.rm = TRUE)
imputation_value <- min_non_zero_lfq / 2 # Common practice: impute with half of the minimum non-zero value

cat(paste("\n### Starting Bayesian Analysis with rstanarm (Work", work_number, ":", condition1_name, "vs", condition2_name, ") ###\n"))

# Modify pre-filter to require at least 2 valid (finite and > 0) LFQ intensity values per condition
# This implements a relaxed "50% rule" for N=3 or N=4 replicates, ensuring some data per group.
min_reps_for_bayes <- 3 # Relaxed from 3 to 2

cat(paste("\nApplying pre-filter for Bayesian analysis: Requiring at least", min_reps_for_bayes, "valid LFQ intensity values (finite and > 0) per condition...\n"))

count_valid_lfq_bayes <- function(row_data, lfq_cols) {
  sum(is.finite(as.numeric(row_data[lfq_cols])) & as.numeric(row_data[lfq_cols]) > 0)
}

# Apply the function to count valid replicates for each protein and condition
num_valid_cond1_bayes <- apply(protein_groups_filtered_current, 1, count_valid_lfq_bayes, lfq_cols = lfq_cols_cond1)
num_valid_cond2_bayes <- apply(protein_groups_filtered_current, 1, count_valid_lfq_bayes, lfq_cols = lfq_cols_cond2)

# Filter the dataframe for Bayesian analysis
protein_groups_for_bayes <- protein_groups_filtered_current %>%
  filter(num_valid_cond1_bayes >= min_reps_for_bayes & num_valid_cond2_bayes >= min_reps_for_bayes)

cat(paste("Number of proteins before Bayesian pre-filter:", nrow(protein_groups_filtered_current), "\n"))
cat(paste("Number of proteins after Bayesian pre-filter (>= ", min_reps_for_bayes, " valid replicates per group):", nrow(protein_groups_for_bayes), "\n"))

if (nrow(protein_groups_for_bayes) == 0) {
  stop("Error: No proteins remaining after pre-filtering with the criterion of ", min_reps_for_bayes, " valid replicates per group for Bayesian analysis. Consider relaxing this filter if necessary.")
}

run_bayesian_analysis_single_protein_v5 <- function(row, cond1_cols, cond2_cols, fold_change_threshold_log2 = 1, imputation_value_arg) {
  protein_id <- row["Protein.IDs"]
  
  lfq_cond1 <- as.numeric(row[cond1_cols])
  lfq_cond2 <- as.numeric(row[cond2_cols])
  
  # Impute zero/NA values before log2 transformation
  # Values that are NA or 0 are replaced with the imputation_value_arg
  lfq_cond1_imputed <- ifelse(is.na(lfq_cond1) | lfq_cond1 == 0, imputation_value_arg, lfq_cond1)
  lfq_cond2_imputed <- ifelse(is.na(lfq_cond2) | lfq_cond2 == 0, imputation_value_arg, lfq_cond2)
  
  result <- list(prob_biol_relevance = NA_real_, log2FC_mean = NA_real_, max_rhat = NA_real_, min_ess = NA_real_, protein_id = protein_id, error_message = "")
  
  # The check for minimum valid replicates is now performed *before* calling this function
  # when constructing protein_groups_for_bayes. So, we no longer need an internal check here.
  data_bayes <- data.frame(
    intensity_log2 = c(log2(lfq_cond1_imputed), log2(lfq_cond2_imputed)), # Use imputed values
    condition = factor(c(rep(condition1_name, length(lfq_cond1_imputed)),
                         rep(condition2_name, length(lfq_cond2_imputed))),
                       levels = c(condition2_name, condition1_name)) # Cond2 as reference
  )
  
  current_fit_warnings <- character(0)
  fit_error <- NULL
  
  fit_bayes <- tryCatch({
    withCallingHandlers({
      stan_glm(
        intensity_log2 ~ condition,
        data = data_bayes,
        family = gaussian(),
        chains = 4,
        iter = 4000,
        warmup = 2000,
        adapt_delta = 0.99,
        refresh = 0,
        prior_intercept = normal(0, 5),
        prior = normal(0, 2)
      )
    },
    warning = function(w) {
      current_fit_warnings <<- c(current_fit_warnings, w$message)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    fit_error <<- e$message
    NULL
  })
  
  if (is.null(fit_bayes) || !inherits(fit_bayes, "stanreg")) {
    if (!is.null(fit_error)) {
      result$error_message <- paste("ERROR in stan_glm fit:", fit_error, sep=" | ")
    } else {
      result$error_message <- "ERROR: stan_glm returned NULL/invalid object after fitting."
    }
    if (length(current_fit_warnings) > 0) {
      result$error_message <- paste(result$error_message, " | Raw warnings from fit:", paste(current_fit_warnings, collapse=" | "), sep="")
    }
    return(result)
  }
  
  if (length(current_fit_warnings) > 0) {
    if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
    result$error_message <- paste(result$error_message, "WARNINGS from stan_glm fit:", paste(current_fit_warnings, collapse=" | "), sep="")
  }
  
  posterior_matrix <- tryCatch({
    # The coefficient name will be `condition` followed by the NON-reference condition name
    as.matrix(fit_bayes)
  }, error = function(e) {
    if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
    result$error_message <- paste(result$error_message, "ERROR extracting posterior (as.matrix):", e$message, sep = "")
    NULL
  })
  
  if (!is.null(posterior_matrix)) {
    contrast_coef_name <- paste0("condition", condition1_name) # e.g., conditionPM
    if (contrast_coef_name %in% colnames(posterior_matrix)) {
      posterior_diff <- posterior_matrix[, contrast_coef_name]
      result$prob_biol_relevance <- max(prob_gt_0, prob_lt_0)
      result$log2FC_mean <- mean(posterior_diff)
    } else {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, paste0("ERROR: '", contrast_coef_name, "' parameter not found in posterior matrix. Check condition levels."), sep = "")
    }
  } else {
    if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
    result$error_message <- paste(result$error_message, "Posterior matrix extraction failed (returned NULL).", sep = "")
  }
  
  mcmc_summary_result <- tryCatch({
    s = summary(fit_bayes)
    if ("Rhat" %in% colnames(s) && "ESS" %in% colnames(s)) {
      return(s[, c("Rhat", "ESS"), drop = FALSE])
    } else {
      stop("Rhat or ESS columns not found in summary output.")
    }
  }, error = function(e) {
    if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
    result$error_message <- paste(result$error_message, "ERROR extracting summary (Rhat/ESS):", e$message, sep = " | ")
    NULL
  })
  
  if (!is.null(mcmc_summary_result) && "Rhat" %in% colnames(mcmc_summary_result) && "ESS" %in% colnames(mcmc_summary_result)) {
    result$max_rhat <- max(mcmc_summary_result[, "Rhat"], na.rm = TRUE)
    result$min_ess <- min(mcmc_summary_result[, "ESS"], na.rm = TRUE)
  } else {
    if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
    result$error_message <- paste(result$error_message, "WARNING: Rhat/ESS summary failed or returned unexpected format.", sep = " | ")
  }
  
  return(result)
}

cat(paste("Executing Bayesian analysis for Work", work_number, "in parallel with", num_cores_to_use, "cores...\n"))

bayesian_results_list <- future_lapply(
  1:nrow(protein_groups_for_bayes),
  function(i) {
    run_bayesian_analysis_single_protein_v5(
      row = protein_groups_for_bayes[i, ],
      cond1_cols = lfq_cols_cond1,
      cond2_cols = lfq_cols_cond2,
      fold_change_threshold_log2 = 1,
      imputation_value_arg = imputation_value # Pass the imputation value
    )
  },
  future.seed = TRUE
)

bayes_output_df_full <- purrr::map_dfr(bayesian_results_list, ~data.frame(
  prob_biol_relevance = .x$prob_biol_relevance,
  log2FC_mean = .x$log2FC_mean,
  max_rhat = .x$max_rhat,
  min_ess = .x$min_ess,
  error_message = .x$error_message,
  stringsAsFactors = FALSE
), .id = NULL)

bayes_output_df_full$Protein.IDs <- protein_groups_for_bayes$Protein.IDs

cat("\n--- Error and Diagnostic Summary (Bayesian Analysis) ---\n")
num_errors <- sum(nchar(bayes_output_df_full$error_message) > 0)
cat(paste("Total proteins with error/info messages:", num_errors, "\n"))

if (num_errors > 0) {
  cat("\nCount of each message type:\n")
  error_messages_present <- bayes_output_df_full$error_message[nchar(bayes_output_df_full$error_message) > 0]
  print(table(error_messages_present))
  
  problematic_proteins_details <- bayes_output_df_full %>%
    filter(grepl("ERROR:|WARNING \\(Divergence\\):|INFO:", error_message, ignore.case = TRUE) & nchar(error_message) > 0) %>%
    select(Protein.IDs, error_message, prob_biol_relevance, log2FC_mean, max_rhat, min_ess)
  
  if (nrow(problematic_proteins_details) > 0) {
    cat("\nDetails for some problematic proteins (first 20 rows):\n")
    print(head(problematic_proteins_details, 20))
  }
} else {
  cat("\nNo explicit error or info messages found in the 'error_message' column.\n")
}
cat("\n--- End of Diagnostic Summary (Bayesian Analysis) ---\n")

# Merge Bayesian results with the original pre-filtered dataframe for Bayes
protein_groups_with_bayes <- protein_groups_for_bayes %>%
  left_join(bayes_output_df_full, by = "Protein.IDs")

cat(paste("\nTotal proteins after Bayesian analysis (before significance filter):", nrow(protein_groups_with_bayes), "\n"))

threshold_prob_biol_relevance <- 0.8
max_rhat_threshold <- 1.01 # These thresholds are kept but not applied in final filter as per user request
min_ess_threshold <- 200    # These thresholds are kept but not applied in final filter as per user request

cat("\n--- GENERATING MAIN BAYESIAN FILE --- \n(Filtering by prob_biol_relevance and excluding fatal rstanarm errors)\n")

sig_bayes_df_final <- protein_groups_with_bayes %>%
  filter(
    !is.na(prob_biol_relevance) & prob_biol_relevance >= threshold_prob_biol_relevance,
    (nchar(error_message) == 0 | !grepl("ERROR:|Divergent transitions|stan_glm returned NULL|Posterior matrix extraction failed|parameter not found", error_message, ignore.case = TRUE))
  ) %>%
  mutate(significant_bayes = TRUE)

cols_to_move_bayes <- c("prob_biol_relevance", "log2FC_mean", "max_rhat", "min_ess", "significant_bayes", "error_message")
first_cols_bayes <- c("Protein.IDs", cols_to_move_bayes)
other_cols_bayes <- setdiff(colnames(sig_bayes_df_final), first_cols_bayes)
cols_to_select <- c(first_cols_bayes, other_cols_bayes)
cols_to_select <- cols_to_select[cols_to_select %in% colnames(sig_bayes_df_final)]
sig_bayes_output_reordered_final <- sig_bayes_df_final[, cols_to_select]

if (dir.exists(current_results_folder_path)) {
  write.table(sig_bayes_output_reordered_final, file = current_output_Bayes_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Significant Bayesian results (prob_biol_relevance >", threshold_prob_biol_relevance, ") saved to:", current_output_Bayes_file, "\n"))
  cat(paste("Number of significant Bayesian proteins:", nrow(sig_bayes_output_reordered_final), "\n"))
} else {
  cat(paste("Error: Results folder does not exist:", current_results_folder_path, "\n"))
}

cat("\n--- First 20 non-significant proteins (based on final filter) ---\n")
non_sig_proteins_details_final <- protein_groups_with_bayes %>%
  anti_join(sig_bayes_df_final, by = "Protein.IDs") %>%
  select(Protein.IDs, prob_biol_relevance, log2FC_mean, max_rhat, min_ess, error_message)

if (nrow(non_sig_proteins_details_final) > 0) {
  print(head(non_sig_proteins_details_final, 20))
} else {
  cat("All proteins were filtered or no non-significant proteins to show.\n")
}
cat("\n--- End of Non-Significant Protein Summary (final) ---\n")

# --- COMPARISONS AND PLOTS (Dynamic for Work_X) ---
cat(paste0("\n### Generating comparisons and plots (Work_", work_number, ": ", condition1_name, " vs ", condition2_name, ") ###\n"))

# --- 1) Bar plot of the number of significant proteins ---
# Read significant results from each test from the results folder dynamically
# Apply janitor::clean_names() right after reading to standardize column names.

# t-Student
sig_student_df_comp <- tryCatch(
  read.delim(current_output_tStudent_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read t-Student file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)
# t-Welch
sig_welch_df_comp <- tryCatch(
  read.delim(current_output_tWelch_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read t-Welch file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)
# limma
sig_limma_df_comp <- tryCatch(
  read.delim(current_output_limma_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read limma file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)
# DEqMS
sig_deqms_df_comp <- tryCatch(
  read.delim(current_output_DEqMS_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read DEqMS file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)
# MSstats
sig_msstats_df_comp <- tryCatch(
  read.delim(current_output_MSStats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read MSstats file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)
# Bayesian
sig_bayes_df_comp <- tryCatch(
  read.delim(current_output_Bayes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    janitor::clean_names(), # Clean column names
  error = function(e) {
    cat(paste("WARNING: Could not read Bayesian file for comparisons -", e$message, "\n"))
    return(NULL)
  }
)

# Auxiliary function to get the number of rows of a dataframe or 0 if NULL
get_nrow <- function(df) {
  if (is.null(df)) {
    return(0)
  } else {
    return(nrow(df))
  }
}

n_sig <- data.frame(
  Test = c("t-Student", "t-Welch", "limma", "DEqMS", "MSstats", "Bayesian"),
  Significant = c(get_nrow(sig_student_df_comp), get_nrow(sig_welch_df_comp), get_nrow(sig_limma_df_comp),
                  get_nrow(sig_deqms_df_comp), get_nrow(sig_msstats_df_comp), get_nrow(sig_bayes_df_comp))
)

barplot_n_sig <- ggplot(n_sig, aes(x = Test, y = Significant, fill = Test)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Significant), vjust = -0.3) +
  labs(title = paste0("Number of significant proteins (adj.p < 0.05 or Prob. > ", threshold_prob_biol_relevance, ") - ", condition1_name, " vs ", condition2_name),
       x = "Test", y = "Number of proteins") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), # Eliminar líneas de cuadrícula mayores
        panel.grid.minor = element_blank()) # Eliminar líneas de cuadrícula menores

print(barplot_n_sig)
ggsave(file.path(current_results_folder_path, paste0("barplot_significant_proteins_", condition1_name, "_vs_", condition2_name, ".png")), plot = barplot_n_sig, width = 10, height = 6)
cat(paste0("Bar plot of significant proteins saved to: ", file.path(current_results_folder_path, paste0("barplot_significant_proteins_", condition1_name, "_vs_", condition2_name, ".png")), "\n"))

# ---
### Visualization of intersections with UpSet Plot
### ---
cat("\n--- Preparing data for UpSet Plot and intersection table ---\n")

# Create the list of Protein.IDs vectors for the UpSet plot
# Column names are already cleaned by janitor::clean_names()
list_sig_proteins_venn <- list()

if (!is.null(sig_student_df_comp) && "protein_i_ds" %in% colnames(sig_student_df_comp)) {
  list_sig_proteins_venn$Student <- sig_student_df_comp$protein_i_ds
}
if (!is.null(sig_welch_df_comp) && "protein_i_ds" %in% colnames(sig_welch_df_comp)) {
  list_sig_proteins_venn$Welch <- sig_welch_df_comp$protein_i_ds
}
if (!is.null(sig_limma_df_comp) && "protein_i_ds" %in% colnames(sig_limma_df_comp)) {
  list_sig_proteins_venn$Limma <- sig_limma_df_comp$protein_i_ds
}
if (!is.null(sig_deqms_df_comp) && "protein_i_ds" %in% colnames(sig_deqms_df_comp)) {
  list_sig_proteins_venn$DEqMS <- sig_deqms_df_comp$protein_i_ds
}
if (!is.null(sig_msstats_df_comp) && "protein_id" %in% colnames(sig_msstats_df_comp)) { # MSstats uses 'protein_id'
  list_sig_proteins_venn$MSstats <- sig_msstats_df_comp$protein_id
}
if (!is.null(sig_bayes_df_comp) && "protein_i_ds" %in% colnames(sig_bayes_df_comp)) {
  list_sig_proteins_venn$Bayesian <- sig_bayes_df_comp$protein_i_ds
}

# Remove empty or NULL elements from the list for the UpSet plot
list_sig_proteins_venn_cleaned <- list_sig_proteins_venn[sapply(list_sig_proteins_venn, function(x) length(x) > 0)]

# Ensure IDs are character type to avoid issues in UpSetR
list_sig_proteins_venn_cleaned_char <- lapply(list_sig_proteins_venn_cleaned, as.character)

# Check if there are enough sets for a meaningful UpSet Plot
if (length(list_sig_proteins_venn_cleaned_char) >= 2) {
  cat("\n--- Generating UpSet Plot of significant protein overlap ---\n")
  
  # Generate the UpSet Plot
  upset_plot <- UpSetR::upset(
    UpSetR::fromList(list_sig_proteins_venn_cleaned_char),
    nsets = length(list_sig_proteins_venn_cleaned_char),
    nintersects = NA,
    order.by = "freq",
    decreasing = TRUE,
    point.size = 3.5,
    line.size = 2,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
  )
  
  # Save the UpSet Plot to a PNG file
  png(file.path(current_results_folder_path, paste0("UpSet_diagram_", condition1_name, "_vs_", condition2_name, ".png")),
      width = 2500, height = 2000, res = 300)
  print(upset_plot)
  dev.off()
  
  cat(paste0("UpSet Plot saved to: ", file.path(current_results_folder_path, paste0("UpSet_diagram_", condition1_name, "_vs_", condition2_name, ".png")), "\n"))
  
} else {
  cat("\nCannot generate UpSet Plot: At least 2 sets of significant proteins are needed for comparison (after cleaning empty sets).\n")
}


# --- To get intersections in tabular form ---
cat(paste0("\n--- Calculating intersections of significant proteins (", condition1_name, " vs ", condition2_name, ") ---\n"))

n_sets_cleaned <- length(list_sig_proteins_venn_cleaned)
all_intersections <- list()

if (n_sets_cleaned > 0) {
  for (i in 1:n_sets_cleaned) {
    combinations <- combn(names(list_sig_proteins_venn_cleaned), m = i, simplify = FALSE)
    for (combo in combinations) {
      intersection <- list_sig_proteins_venn_cleaned[[combo[1]]]
      if (length(combo) > 1) {
        for (j in 2:length(combo)) {
          intersection <- intersect(intersection, list_sig_proteins_venn_cleaned[[combo[j]]])
        }
      }
      combo_name <- paste(combo, collapse = "&")
      all_intersections[[combo_name]] <- length(intersection)
    }
  }
  
  cat(paste0("\n--- Intersections of significant proteins (", condition1_name, " vs ", condition2_name, ") ---\n"))
  print(all_intersections)
  
  capture.output(print(all_intersections), file = file.path(current_results_folder_path, paste0("Intersections_", condition1_name, "_vs_", condition2_name, ".txt")))
  cat(paste0("Intersection details saved to ", paste0("Intersections_", condition1_name, "_vs_", condition2_name, ".txt"), " in: ", current_results_folder_path, "\n"))
} else {
  cat("\nNo sets of significant proteins to calculate intersections.\n")
}


# ---
### Distribution of adjusted p-values
# ---
cat("\n--- Generating distribution of -log10(adjusted p-values) ---\n")

# Read full results of methods from files in W_X_H_testing
# janitor::clean_names() is also applied here
student_output_df_all <- tryCatch(
  read.delim(current_output_tStudent_file, sep = "\t", stringsAsFactors = FALSE) %>% janitor::clean_names(),
  error = function(e) {
    cat(paste("WARNING: Could not read", current_output_tStudent_file, "for density plot -", e$message, "\n"))
    return(NULL)
  }
)
welch_output_df_all <- tryCatch(
  read.delim(current_output_tWelch_file, sep = "\t", stringsAsFactors = FALSE) %>% janitor::clean_names(),
  error = function(e) {
    cat(paste("WARNING: Could not read", current_output_tWelch_file, "for density plot -", e$message, "\n"))
    return(NULL)
  }
)
limma_results_all <- tryCatch(
  read.delim(current_output_limma_file, sep = "\t", stringsAsFactors = FALSE) %>% janitor::clean_names(),
  error = function(e) {
    cat(paste("WARNING: Could not read", current_output_limma_file, "for density plot -", e$message, "\n"))
    return(NULL)
  }
)
deqms_output_all <- tryCatch(
  read.delim(current_output_DEqMS_file, sep = "\t", stringsAsFactors = FALSE) %>% janitor::clean_names(),
  error = function(e) {
    cat(paste("WARNING: Could not read", current_output_DEqMS_file, "for density plot -", e$message, "\n"))
    return(NULL)
  }
)
msstats_output_all <- tryCatch(
  read.delim(current_output_MSStats_file, sep = "\t", stringsAsFactors = FALSE) %>% janitor::clean_names(),
  error = function(e) {
    cat(paste("WARNING: Could not read", current_output_MSStats_file, "for density plot -", e$message, "\n"))
    return(NULL)
  }
)

# Calculate the intercept value for adjusted p-value
cutoff_log10p <- -log10(0.05)
cat(paste0("Adjusted p-value intercept value (Work_", work_number, "): ", cutoff_log10p, "\n"))

# Prepare transformed data (using all adjusted p-values)
dfs_for_density_plot <- list()

# Limma
if (!is.null(limma_results_all) && "adj_p_val" %in% colnames(limma_results_all) && nrow(limma_results_all) > 0) {
  dfs_for_density_plot[["Limma"]] <- data.frame(transformed_p = -log10(limma_results_all$adj_p_val), method = "Limma")
} else {
  cat("Limma data not available for density plot, missing 'adj_p_val', or no rows of data.\n")
}

# Student's t-test
if (!is.null(student_output_df_all) && "adj_p_value" %in% colnames(student_output_df_all) && nrow(student_output_df_all) > 0) {
  dfs_for_density_plot[["Student"]] <- data.frame(transformed_p = -log10(student_output_df_all$adj_p_value), method = "Student")
} else {
  cat("Student data not available for density plot, missing 'adj_p_value', or no rows of data.\n")
}

# Welch's t-test (assuming it also uses 'adj_p_value' column)
if (!is.null(welch_output_df_all) && "adj_p_value" %in% colnames(welch_output_df_all) && nrow(welch_output_df_all) > 0) {
  dfs_for_density_plot[["Welch"]] <- data.frame(transformed_p = -log10(welch_output_df_all$adj_p_value), method = "Welch")
} else {
  cat("Welch data not available for density plot, missing 'adj_p_value', or no rows of data.\n")
}

# DEqMS (assuming it uses 'adj_p_val' column)
if (!is.null(deqms_output_all) && "adj_p_val" %in% colnames(deqms_output_all) && nrow(deqms_output_all) > 0) {
  dfs_for_density_plot[["DEqMS"]] <- data.frame(transformed_p = -log10(deqms_output_all$adj_p_val), method = "DEqMS")
} else {
  cat("DEqMS data not available for density plot, missing 'adj_p_val', or no rows of data.\n")
}

# MSstats (assuming it uses 'adj_p_value' column based on your previous logs)
if (!is.null(msstats_output_all) && nrow(msstats_output_all) > 0) {
  # DEBUG: Print column names to verify what clean_names() produced
  cat("DEBUG: Columns in msstats_output_all after clean_names(): ")
  print(colnames(msstats_output_all))
  
  # CORRECCIÓN: Usar "adj_pvalue" según la salida de depuración
  if ("adj_pvalue" %in% colnames(msstats_output_all)) {
    # Filter out non-finite or zero p-values before transformation
    valid_p_values <- msstats_output_all$adj_pvalue[is.finite(msstats_output_all$adj_pvalue) & msstats_output_all$adj_pvalue > 0]
    if (length(valid_p_values) > 0) {
      dfs_for_density_plot[["MSstats"]] <- data.frame(transformed_p = -log10(valid_p_values), method = "MSstats")
      cat("DEBUG: MSstats data successfully added for density plot.\n")
    } else {
      cat("WARNING: MSstats data has no finite, positive adjusted p-values after clean_names(). Not added to density plot.\n")
    }
  } else {
    # This warning should now only appear if 'adj_pvalue' itself is truly missing
    cat("WARNING: MSstats data does not contain 'adj_pvalue' column after clean_names(). Not added to density plot.\n")
  }
} else {
  cat("WARNING: msstats_output_all is NULL or has no rows. Not added to density plot.\n")
}

# Combine all dataframes into one for plotting
if (length(dfs_for_density_plot) > 0) {
  combined_density_df <- do.call(rbind, dfs_for_density_plot)
  
  # Check if there's actual data after combining before plotting
  if (nrow(combined_density_df) > 0) {
    # Create the density plot of -log10(adjusted p-value) for all results with a cutoff line
    transformed_density_plot_cutoff_all <- ggplot(combined_density_df, aes(x = transformed_p, color = method)) + 
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = cutoff_log10p, linetype = "dashed", color = "red") +
      labs(title = paste0("Distribution of -log10 (adjusted p-value) - ", condition1_name, " vs ", condition2_name),
           x = "-log10(adjusted p-value)",
           y = "Density",
           color = "Method") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
    
    # Save the plot
    filename_logp_cutoff_all <- file.path(current_results_folder_path, paste0("density_adj_pvalue_all_", condition1_name, "_vs_", condition2_name, ".png"))
    ggsave(filename = filename_logp_cutoff_all, plot = transformed_density_plot_cutoff_all, width = 10, height = 6, dpi = 300)
    cat(paste0("Density plot of -log10(adjusted p-value) (all with cutoff) saved to: ", filename_logp_cutoff_all, "\n"))
    
    print(transformed_density_plot_cutoff_all)
  } else {
    cat("\nNo combined data with rows available to generate the adjusted p-value density plot.\n")
  }
} else {
  cat("\nNo data available from any method for the density plot.\n")
}


# ---
### Ranked -log10(adjusted p-value) line plots
# ---
cat(paste0("\n--- Generating ranked -log10(adjusted p-value) line plots (Work_", work_number, ": ", condition1_name, " vs ", condition2_name, ") ---\n"))

dfs_for_ranked_plot <- list()

# --- t-Student ---
if (!is.null(student_output_df_all)) {
  # Force protein ID columns to character immediately to prevent type issues
  # This targets the original dataframes directly
  if ("majority_protein_i_ds" %in% colnames(student_output_df_all)) {
    student_output_df_all$majority_protein_i_ds <- as.character(student_output_df_all$majority_protein_i_ds)
  }
  if ("protein_i_ds" %in% colnames(student_output_df_all)) {
    student_output_df_all$protein_i_ds <- as.character(student_output_df_all$protein_i_ds)
  }
  
  current_cols <- colnames(student_output_df_all) # Re-evaluate columns after potential modification
  adj_p_value_col_name <- "adj_p_value"
  
  # Define possible protein ID column names after janitor::clean_names()
  possible_protein_id_cols <- c("majority_protein_i_ds", "protein_i_ds")
  found_protein_id_col <- NULL
  
  for (col_name_check in possible_protein_id_cols) {
    if (col_name_check %in% current_cols) {
      found_protein_id_col <- col_name_check
      break
    }
  }
  
  if (!is.null(found_protein_id_col) && adj_p_value_col_name %in% current_cols) {
    cat(paste0("DEBUG: Adding t-Student data to ranked plot. Found protein ID column: '", found_protein_id_col, "'.\n"))
    dfs_for_ranked_plot[["Student"]] <- student_output_df_all %>%
      # Now that original column is forced to character, select and rename
      select(p_adj = all_of(adj_p_value_col_name), protein = all_of(found_protein_id_col)) %>%
      mutate(test = "t-Student")
  } else {
    cat(paste0("WARNING: t-Student data not available or required columns ('", adj_p_value_col_name, "', and a protein ID column from ", paste(possible_protein_id_cols, collapse = ", "), ") missing for ranked plot.\n"))
    cat("DEBUG: Current columns of student_output_df_all: ")
    print(current_cols)
  }
} else {
  cat("WARNING: student_output_df_all is NULL.\n")
}

# --- t-Welch ---
if (!is.null(welch_output_df_all)) {
  # Force protein ID columns to character immediately to prevent type issues
  if ("majority_protein_i_ds" %in% colnames(welch_output_df_all)) {
    welch_output_df_all$majority_protein_i_ds <- as.character(welch_output_df_all$majority_protein_i_ds)
  }
  if ("protein_i_ds" %in% colnames(welch_output_df_all)) {
    welch_output_df_all$protein_i_ds <- as.character(welch_output_df_all$protein_i_ds)
  }
  
  current_cols <- colnames(welch_output_df_all) # Re-evaluate columns after potential modification
  adj_p_value_col_name <- "adj_p_value"
  
  possible_protein_id_cols <- c("majority_protein_i_ds", "protein_i_ds")
  found_protein_id_col <- NULL
  
  for (col_name_check in possible_protein_id_cols) {
    if (col_name_check %in% current_cols) {
      found_protein_id_col <- col_name_check
      break
    }
  }
  
  if (!is.null(found_protein_id_col) && adj_p_value_col_name %in% current_cols) {
    cat(paste0("DEBUG: Adding t-Welch data to ranked plot. Found protein ID column: '", found_protein_id_col, "'.\n"))
    dfs_for_ranked_plot[["Welch"]] <- welch_output_df_all %>%
      # Now that original column is forced to character, select and rename
      select(p_adj = all_of(adj_p_value_col_name), protein = all_of(found_protein_id_col)) %>%
      mutate(test = "t-Welch")
  } else {
    cat(paste0("WARNING: t-Welch data not available or required columns ('", adj_p_value_col_name, "', and a protein ID column from ", paste(possible_protein_id_cols, collapse = ", "), ") missing for ranked plot.\n"))
    cat("DEBUG: Current columns of welch_output_df_all: ")
    print(current_cols)
  }
} else {
  cat("WARNING: welch_output_df_all is NULL.\n")
}

# --- limma ---
if (!is.null(limma_results_all)) {
  # Force protein ID columns to character immediately to prevent type issues
  if ("majority_protein_i_ds" %in% colnames(limma_results_all)) {
    limma_results_all$majority_protein_i_ds <- as.character(limma_results_all$majority_protein_i_ds)
  }
  if ("protein_i_ds" %in% colnames(limma_results_all)) {
    limma_results_all$protein_i_ds <- as.character(limma_results_all$protein_i_ds)
  }
  
  current_cols <- colnames(limma_results_all) # Re-evaluate columns after potential modification
  adj_p_val_col_name <- "adj_p_val"
  
  possible_protein_id_cols <- c("majority_protein_i_ds", "protein_i_ds")
  found_protein_id_col <- NULL
  
  for (col_name_check in possible_protein_id_cols) {
    if (col_name_check %in% current_cols) {
      found_protein_id_col <- col_name_check
      break
    }
  }
  
  if (!is.null(found_protein_id_col) && adj_p_val_col_name %in% current_cols) {
    cat(paste0("DEBUG: Adding limma data to ranked plot. Found protein ID column: '", found_protein_id_col, "'.\n"))
    dfs_for_ranked_plot[["Limma"]] <- limma_results_all %>%
      # Now that original column is forced to character, select and rename
      select(p_adj = all_of(adj_p_val_col_name), protein = all_of(found_protein_id_col)) %>%
      mutate(test = "Limma")
  } else {
    cat(paste0("WARNING: limma data not available or required columns ('", adj_p_val_col_name, "', and a protein ID column from ", paste(possible_protein_id_cols, collapse = ", "), ") missing for ranked plot.\n"))
    cat("DEBUG: Current columns of limma_results_all: ")
    print(current_cols)
  }
} else {
  cat("WARNING: limma_results_all is NULL.\n")
}

# --- DEqMS ---
if (!is.null(deqms_output_all)) {
  # Force protein ID columns to character immediately to prevent type issues
  if ("majority_protein_i_ds" %in% colnames(deqms_output_all)) {
    deqms_output_all$majority_protein_i_ds <- as.character(deqms_output_all$majority_protein_i_ds)
  }
  if ("protein_i_ds" %in% colnames(deqms_output_all)) {
    deqms_output_all$protein_i_ds <- as.character(deqms_output_all$protein_i_ds)
  }
  
  current_cols <- colnames(deqms_output_all) # Re-evaluate columns after potential modification
  adj_p_val_col_name <- "adj_p_val"
  
  possible_protein_id_cols <- c("majority_protein_i_ds", "protein_i_ds")
  found_protein_id_col <- NULL
  
  for (col_name_check in possible_protein_id_cols) {
    if (col_name_check %in% current_cols) {
      found_protein_id_col <- col_name_check
      break
    }
  }
  
  if (!is.null(found_protein_id_col) && adj_p_val_col_name %in% current_cols) {
    cat(paste0("DEBUG: Adding DEqMS data to ranked plot. Found protein ID column: '", found_protein_id_col, "'.\n"))
    dfs_for_ranked_plot[["DEqMS"]] <- deqms_output_all %>%
      # Now that original column is forced to character, select and rename
      select(p_adj = all_of(adj_p_val_col_name), protein = all_of(found_protein_id_col)) %>%
      mutate(test = "DEqMS")
  } else {
    cat(paste0("WARNING: DEqMS data not available or required columns ('", adj_p_val_col_name, "', and a protein ID column from ", paste(possible_protein_id_cols, collapse = ", "), ") missing for ranked plot.\n"))
    cat("DEBUG: Current columns of deqms_output_all: ")
    print(current_cols)
  }
} else {
  cat("WARNING: deqms_output_all is NULL.\n")
}

# --- MSstats ---
if (!is.null(msstats_output_all)) {
  # Force protein ID columns to character immediately to prevent type issues
  if ("protein_id" %in% colnames(msstats_output_all)) {
    msstats_output_all$protein_id <- as.character(msstats_output_all$protein_id)
  }
  
  current_cols <- colnames(msstats_output_all) # Re-evaluate columns after potential modification
  adj_pvalue_col_name <- "adj_pvalue"
  # MSstats specifically uses 'protein_id'
  possible_protein_id_cols <- c("protein_id") 
  found_protein_id_col <- NULL
  
  for (col_name_check in possible_protein_id_cols) {
    if (col_name_check %in% current_cols) {
      found_protein_id_col <- col_name_check
      break
    }
  }
  
  if (adj_pvalue_col_name %in% current_cols && !is.null(found_protein_id_col)) {
    cat(paste0("DEBUG: Adding MSstats data to ranked plot. Found protein ID column: '", found_protein_id_col, "'.\n"))
    dfs_for_ranked_plot[["MSstats"]] <- msstats_output_all %>%
      # Now that original column is forced to character, select and rename
      select(p_adj = all_of(adj_pvalue_col_name), protein = all_of(found_protein_id_col)) %>%
      mutate(test = "MSstats")
  } else {
    cat(paste0("WARNING: MSstats data not available or required columns ('", adj_pvalue_col_name, "', and a protein ID column from ", paste(possible_protein_id_cols, collapse = ", "), ") missing for ranked plot.\n"))
    cat("DEBUG: Current columns of msstats_output_all: ")
    print(current_cols)
  }
} else {
  cat("WARNING: msstats_output_all is NULL.\n")
}


# Combine all dataframes and prepare for ranking if data exists
if (length(dfs_for_ranked_plot) > 0) {
  all_sig_proteins_ranked <- bind_rows(dfs_for_ranked_plot) %>%
    filter(is.finite(p_adj) & p_adj > 0) %>% # Ensure p_adj is finite and > 0 for log10
    mutate(neg_log10_p_adj = -log10(p_adj)) %>%
    group_by(test) %>%
    arrange(desc(neg_log10_p_adj)) %>% # Order by -log10(adjusted p-value) descending for ascending rank
    mutate(rank = row_number()) %>% # Assign rank
    ungroup()
  
  if (nrow(all_sig_proteins_ranked) > 0) {
    # Create the ranked line plot
    ranked_plot <- ggplot(all_sig_proteins_ranked, aes(x = rank, y = neg_log10_p_adj, color = test)) +
      geom_line(linewidth = 1) + # Modificado: 'size' a 'linewidth'
      labs(title = paste0("Rank of -log10(adjusted p-value) by method - ", condition1_name, " vs ", condition2_name),
           x = "Rank (Proteins ordered by -log10(adjusted p-value) from high to low)",
           y = "-log10(adjusted p-value)",
           color = "Method") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), # Eliminar líneas de cuadrícula mayores
            panel.grid.minor = element_blank()) # Eliminar líneas de cuadrícula menores
    
    # Save the plot
    filename_ranked_plot <- file.path(current_results_folder_path, paste0("ranked_adj_pvalue_lines_", condition1_name, "_vs_", condition2_name, ".png"))
    ggsave(filename = filename_ranked_plot, plot = ranked_plot, width = 10, height = 6, dpi = 300)
    cat(paste0("Ranked -log10(adjusted p-value) line plot saved to: ", filename_ranked_plot, "\n"))
    
    print(ranked_plot)
  } else {
    cat("\nNot enough data to generate the ranked line plot.\n")
  }
} else {
  cat("\nNo data available to generate the ranked adjusted p-value line plot.\n")
}

# --- FINALIZING PARALLEL SESSION ---
plan(sequential) # Deactivates parallelization and cleans up workers.
cat("Hypothesis testing analysis completed. Parallel workers shut down.\n")


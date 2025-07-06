# Load necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(VennDiagram) # For Venn diagrams
library(UpSetR)      # For UpSet plots
library(future)      # For parallel processing
library(future.apply) # For parallel apply functions
library(gridExtra)   # For grid.arrange if needed, though grid.draw for Venn is direct.
library(grid)        # For grid.draw with VennDiagram

cat("Libraries loaded successfully.\n")

# --- Configuration ---
# Set the base parent directory where W1, W2, etc., folders are located
base_data_parent_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/7_Especial_IJMS/analisis27062025"
cat(paste0("Base data parent path set to: ", base_data_parent_path, "\n"))

# --- Parallel Processing Configuration ---
# To force sequential processing for debugging, use: plan(sequential)
# To enable parallel processing, use: plan(multisession, workers = X) where X is the number of cores.
# Or: plan(multisession, workers = availableCores()) to use all available cores.

# *** CHOOSE ONE OF THE FOLLOWING TWO LINES, AND COMMENT THE OTHER ***

# For DEBUGGING (Sequential processing):
# plan(sequential)
# cat("Forcing sequential processing for debugging purposes.\n")

# For PRODUCTION (Parallel processing):
plan(multisession, workers = 7) # Using 7 cores. Adjust 'workers' as needed.
cat(paste0("Enabling parallel processing with ", future::nbrOfWorkers(), " workers.\n"))

# --- END Parallel Processing Configuration ---


# --- AUXILIARY FUNCTIONS DEFINITION ---

# Auxiliary function to get protein ID column name based on the method
get_id_col_name <- function(method) {
  if (method == "msstats") {
    return("Protein.ID")
  } else {
    return("Protein.IDs") # Default for MaxQuant outputs
  }
}

# Function for Bayesian analysis for a single protein (IMPROVED Version v7)
run_bayesian_analysis_single_protein_v7 <- function(protein_row_data,
                                                    cond1_ibaq_cols,
                                                    cond2_ibaq_cols,
                                                    cond1_label, # This is the label for the model (e.g., "Cold", "u", "NaCl")
                                                    cond2_label, # This is the label for the model (e.g., "RT", "y", "Control")
                                                    fold_change_threshold_log2 = 1) {
  
  protein_id_col_name <- ifelse("Protein.IDs" %in% names(protein_row_data), "Protein.IDs", "Protein.ID")
  protein_id <- protein_row_data[[protein_id_col_name]]
  
  # Ensure all column names are treated as characters for safe indexing
  protein_row_data_list <- as.list(protein_row_data) # Convert to list to preserve names for as.numeric
  
  ibaq_cond1_values <- as.numeric(protein_row_data_list[cond1_ibaq_cols])
  ibaq_cond2_values <- as.numeric(protein_row_data_list[cond2_ibaq_cols])
  
  ibaq_cond1_filtered <- ibaq_cond1_values[is.finite(ibaq_cond1_values) & ibaq_cond1_values > 0]
  ibaq_cond2_filtered <- ibaq_cond2_values[is.finite(ibaq_cond2_values) & ibaq_cond2_values > 0]
  
  result <- list(prob_biol_relevance = NA_real_,
                 log2FC_mean = NA_real_,
                 max_rhat = NA_real_,
                 min_ess = NA_real_,
                 protein_id = protein_id,
                 error_message = "")
  
  if (length(ibaq_cond1_filtered) >= 3 && length(ibaq_cond2_filtered) >= 3) {
    data_bayes <- data.frame(
      intensity_log2 = c(log2(ibaq_cond1_filtered), log2(ibaq_cond2_filtered)),
      condition = factor(c(rep(cond1_label, length(ibaq_cond1_filtered)),
                           rep(cond2_label, length(ibaq_cond2_filtered))),
                         levels = c(cond2_label, cond1_label)) # Ensure correct contrast for log2FC interpretation (cond1 vs cond2)
    )
    
    current_fit_warnings <- character(0)
    fit_error <- NULL
    
    cat(paste0("   DEBUG: Attempting stan_glm fit for protein: ", protein_id, ". Data points: ", nrow(data_bayes), "\n"))
    
    fit_bayes <- tryCatch({
      withCallingHandlers({
        rstanarm::stan_glm(
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
      cat(paste0("   DEBUG: stan_glm FAILED for protein ", protein_id, ". Error/Warnings: ", result$error_message, "\n"))
      return(result)
    }
    
    cat(paste0("   DEBUG: stan_glm fit SUCCESSFUL for protein: ", protein_id, ". Extracting posterior and diagnostics.\n"))
    
    if (length(current_fit_warnings) > 0) {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, "WARNINGS from stan_glm fit:", paste(current_fit_warnings, collapse=" | "), sep="")
    }
    
    posterior_matrix <- tryCatch({
      as.matrix(fit_bayes)
    }, error = function(e) {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, "ERROR extracting posterior (as.matrix):", e$message, sep = "")
      NULL
    })
    
    if (!is.null(posterior_matrix)) {
      coeff_name <- paste0("condition", cond1_label) 
      
      if (coeff_name %in% colnames(posterior_matrix)) {
        posterior_diff <- posterior_matrix[, coeff_name]
        result$prob_biol_relevance <- mean(abs(posterior_diff) > fold_change_threshold_log2)
        result$log2FC_mean <- mean(posterior_diff)
      } else {
        if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
        result$error_message <- paste(result$error_message,
                                      paste0("ERROR: '", coeff_name, "' parameter not found in posterior matrix. Check factor levels in stan_glm."),
                                      sep = "")
      }
    } else {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, "Posterior matrix extraction failed (returned NULL).", sep = "")
    }
    
    mcmc_summary_data <- tryCatch({
      s_summary <- summary(fit_bayes, pars = c("condition", "(Intercept)"), digits = 4, probs = c(0.025, 0.975))
      if ("Rhat" %in% colnames(s_summary$coefficients) && "ESS" %in% colnames(s_summary$coefficients)) {
        return(s_summary$coefficients[, c("Rhat", "ESS"), drop = FALSE])
      } else {
        stop("Rhat or ESS columns not found in summary(fit_bayes)$coefficients.")
      }
    }, error = function(e) {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, "ERROR extracting summary (Rhat/ESS):", e$message, sep = " | ")
      NULL
    })
    
    if (!is.null(mcmc_summary_data) && "Rhat" %in% colnames(mcmc_summary_data) && "ESS" %in% colnames(mcmc_summary_data)) {
      result$max_rhat <- max(mcmc_summary_data[, "Rhat"], na.rm = TRUE)
      result$min_ess <- min(mcmc_summary_data[, "ESS"], na.rm = TRUE)
    } else {
      if (result$error_message != "") result$error_message <- paste(result$error_message, " | ", sep="")
      result$error_message <- paste(result$error_message, "WARNING: Rhat/ESS summary failed or returned unexpected format (mcmc_summary_data is NULL or missing columns).", sep = " | ")
    }
    
    cat(paste0("   INFO: Protein ", protein_id, " - Prob. Relevance: ", round(result$prob_biol_relevance, 4),
               ", Log2FC Mean: ", round(result$log2FC_mean, 4),
               ", Max Rhat: ", round(result$max_rhat, 4),
               ", Min ESS: ", round(result$min_ess, 0),
               ", Error Message: '", result$error_message, "'\n"))
    
    return(result)
    
  } else {
    result$error_message <- "INFO: Not enough valid replicates (less than 3 per group)."
    cat(paste0("   INFO: Protein ", protein_id, " - Skipped (Not enough replicates). Error Message: '", result$error_message, "'\n"))
    return(result)
  }
}

# --- MAIN FUNCTION FOR BIOLOGICAL RELEVANCE FILTERS AND OVERLAP PLOTS ---

process_overlap_analysis <- function(work_number, condition1_name, condition2_name, base_parent_path) {
  
  cat(paste0("\n--- STARTING OVERLAP AND FC ANALYSIS FOR WORK_", work_number, " (", condition1_name, " vs ", condition2_name, ") ---\n"))
  
  current_main_folder_path <- file.path(base_parent_path, paste0("W", work_number))
  current_ht_results_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_H_testing"))
  current_biol_rel_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_BioRel"))
  current_biol_comparison_folder_path <- file.path(current_biol_rel_folder_path, paste0("W", work_number, "_BioRel_comparison"))
  
  if (!dir.exists(current_biol_rel_folder_path)) {
    dir.create(current_biol_rel_folder_path, recursive = TRUE)
    cat(paste("Biological relevance results folder created at:", current_biol_rel_folder_path, "\n"))
  } else {
    cat(paste("Biological relevance results folder already exists at:", current_biol_rel_folder_path, "\n"))
  }
  
  if (!dir.exists(current_biol_comparison_folder_path)) {
    dir.create(current_biol_comparison_folder_path, recursive = TRUE)
    cat(paste("Biological relevance comparison folder created at:", current_biol_comparison_folder_path, "\n"))
  } else {
    cat(paste("Biological relevance comparison folder already exists at:", current_biol_comparison_folder_path, "\n"))
  }
  
  fc_cutoff <- 2
  log2_fc_cutoff <- log2(fc_cutoff)
  threshold_prob_biol_relevance <- 0.95
  
  cat(paste("\n--- START: Applying Fold Change filter (absolute log2FC >= ", round(log2_fc_cutoff, 2), ") to Work ", work_number, " ---\n"))
  
  test_names_fc <- c("tstudent", "twelch", "limma", "deqms", "msstats")
  files_to_process_fc <- setNames(as.list(test_names_fc), paste0("W", work_number, "_res_", test_names_fc, ".txt"))
  proteins_passing_fc_filter <- list()
  
  for (file_name in names(files_to_process_fc)) {
    method_name <- files_to_process_fc[[file_name]]
    file_path <- file.path(current_ht_results_folder_path, file_name)
    
    cat(paste("\nProcessing file for FC filter:", file_name, " (Method:", method_name, ")\n"))
    
    if (!file.exists(file_path)) {
      cat(paste("WARNING: File not found:", file_path, ". Skipping this method for FC filter.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_FC.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(data.frame(), file = output_file_path)
      cat(paste("Empty FC output file created (due to missing input) at:", output_file_path, "\n"))
      proteins_passing_fc_filter[[method_name]] <- character(0)
      next
    }
    
    res_df <- tryCatch({
      read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      cat(paste("ERROR reading file for FC filter:", file_path, ". Message:", e$message, ". Skipping this method.\n"))
      return(NULL)
    })
    
    if (is.null(res_df) || nrow(res_df) == 0) {
      cat(paste("INFO: File", file_name, "is empty or could not be read. No proteins to filter by FC.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_FC.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(res_df[0, ], file = output_file_path)
      cat(paste("Empty FC output file created with header at:", output_file_path, "\n"))
      proteins_passing_fc_filter[[method_name]] <- character(0)
      next
    }
    
    current_log2fc_col <- if (method_name %in% c("tstudent", "twelch", "msstats")) "log2FC" else "logFC"
    
    if (!current_log2fc_col %in% colnames(res_df)) {
      cat(paste("WARNING: Column '", current_log2fc_col, "' not found in file",
                file_name, ". Skipping FC filter for this method.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_FC.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(data.frame(), file = output_file_path)
      proteins_passing_fc_filter[[method_name]] <- character(0)
      next
    }
    
    proteins_fc_filtered <- res_df %>%
      filter(abs(.[[current_log2fc_col]]) >= log2_fc_cutoff)
    
    protein_id_col <- get_id_col_name(method_name)
    
    if (!(protein_id_col %in% colnames(proteins_fc_filtered))) {
      cat(paste("WARNING: Protein ID column '", protein_id_col, "' not found in filtered data for", method_name, ".\n"))
      proteins_passing_fc_filter[[method_name]] <- character(0)
    } else if (nrow(proteins_fc_filtered) > 0) {
      proteins_passing_fc_filter[[method_name]] <- unique(proteins_fc_filtered[[protein_id_col]])
      cat(paste("Method", method_name, ":", length(proteins_passing_fc_filter[[method_name]]), "proteins passed FC filter.\n"))
    } else {
      proteins_passing_fc_filter[[method_name]] <- character(0)
      cat(paste("Method", method_name, ": No proteins passed FC filter.\n"))
    }
    
    output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_FC.txt")
    output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
    write_tsv(proteins_fc_filtered, file = output_file_path)
    cat(paste("FC filtered results written to:", output_file_path, "\n"))
  }
  
  cat(paste("--- END: Fold Change filter for Work ", work_number, " ---\n"))
  
  # --- 4. Bayesian Filtering ---
  
  cat(paste("\n--- START: Applying Bayesian filter (probability of relevance >= ", threshold_prob_biol_relevance, ") to Work ", work_number, " ---\n"))
  
  test_names_bayes_filter <- c("tstudent", "twelch", "limma", "deqms", "msstats")
  files_to_process_bayes_filter <- setNames(as.list(test_names_bayes_filter), paste0("W", work_number, "_res_", test_names_bayes_filter, ".txt"))
  proteins_passing_bayes_filter <- list()
  
  for (file_name in names(files_to_process_bayes_filter)) {
    method_name <- files_to_process_bayes_filter[[file_name]]
    file_path <- file.path(current_ht_results_folder_path, file_name)
    
    cat(paste("\nProcessing file for Bayesian filter:", file_name, " (Method:", method_name, ")\n"))
    
    if (!file.exists(file_path)) {
      cat(paste("WARNING: File not found:", file_path, ". Skipping this method for Bayesian filter.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_Bayes.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(data.frame(), file = output_file_path)
      cat(paste("Empty Bayesian output file created (due to missing input) at:", output_file_path, "\n"))
      proteins_passing_bayes_filter[[method_name]] <- character(0)
      next
    }
    
    res_df_freq <- tryCatch({
      read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      cat(paste("ERROR reading file for Bayesian filter:", file_path, ". Message:", e$message, ". Skipping this method.\n"))
      return(NULL)
    })
    
    if (is.null(res_df_freq) || nrow(res_df_freq) == 0) {
      cat(paste("INFO: File", file_name, "is empty or could not be read. No proteins for Bayesian analysis.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_Bayes.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(res_df_freq[0, ], file = output_file_path)
      cat(paste("Empty Bayesian output file created with header at:", output_file_path, "\n"))
      proteins_passing_bayes_filter[[method_name]] <- character(0)
      next
    }
    
    cat("DEBUG: Loaded data for Bayesian analysis head and structure:\n")
    print(head(res_df_freq))
    print(str(res_df_freq))
    
    cat(paste("DEBUG: Attempting to identify LFQ.intensity columns for conditions: '", condition1_name, "' and '", condition2_name, "'.\n", sep=""))
    cat(paste("DEBUG: Full column names in current file:", paste(colnames(res_df_freq), collapse=", "), "\n"))
    
    # The condition1_name and condition2_name are now *expected* to be the exact prefixes in the LFQ.intensity column names.
    escaped_cond1_name <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", condition1_name)
    escaped_cond2_name <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", condition2_name)
    
    # This regex now handles both "LFQ.intensity.Cold1_1" (matches "_1") and "LFQ.intensity.u1" (no "_1")
    # It looks for "LFQ.intensity.YourPrefix" followed by one or more digits,
    # optionally followed by an underscore and more digits (for technical replicates).
    regex_cond1_pattern <- paste0("^LFQ\\.intensity\\.", escaped_cond1_name, "[0-9]+(_[0-9]+)?$")
    regex_cond2_pattern <- paste0("^LFQ\\.intensity\\.", escaped_cond2_name, "[0-9]+(_[0-9]+)?$")
    
    ibaq_cols_cond1_current_file <- colnames(res_df_freq)[grepl(regex_cond1_pattern, colnames(res_df_freq))]
    ibaq_cols_cond2_current_file <- colnames(res_df_freq)[grepl(regex_cond2_pattern, colnames(res_df_freq))]
    
    cat(paste("DEBUG: Final Identified Condition 1 LFQ.intensity columns (for '", condition1_name, "'):", paste(ibaq_cols_cond1_current_file, collapse=", "), "\n"))
    cat(paste("DEBUG: Final Identified Condition 2 LFQ.intensity columns (for '", condition2_name, "'):", paste(ibaq_cols_cond2_current_file, collapse=", "), "\n"))
    
    if (length(ibaq_cols_cond1_current_file) < 3 || length(ibaq_cols_cond2_current_file) < 3) {
      cat(paste("WARNING: Not enough LFQ.intensity columns (need at least 3) for conditions '",
                condition1_name, "' and '", condition2_name,
                "' for Bayesian analysis in file", file_name, ". Skipping Bayesian filter for this method.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_Bayes.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(data.frame(), file = output_file_path)
      cat(paste("Empty Bayesian output file created (due to missing input) at:", output_file_path, "\n"))
      proteins_passing_bayes_filter[[method_name]] <- character(0)
      next
    }
    
    protein_id_col <- get_id_col_name(method_name)
    if (!(protein_id_col %in% colnames(res_df_freq))) {
      cat(paste("WARNING: Protein ID column '", protein_id_col, "' not found in data for", method_name, ". Skipping Bayesian analysis.\n"))
      output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_Bayes.txt")
      output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
      write_tsv(data.frame(), file = output_file_path)
      proteins_passing_bayes_filter[[method_name]] <- character(0)
      next
    }
    
    proteins_to_analyze_bayes <- res_df_freq %>%
      select(all_of(c(protein_id_col, ibaq_cols_cond1_current_file, ibaq_cols_cond2_current_file)))
    
    cat(paste("INFO: Running Bayesian analysis for", nrow(proteins_to_analyze_bayes), "proteins using method", method_name, "...\n"))
    
    # --- HERE IS THE CRITICAL CALL TO future_apply ---
    cat(paste0("   DEBUG: Calling future_apply for ", nrow(proteins_to_analyze_bayes), " proteins. This should now run in parallel.\n"))
    bayesian_results_list <- future_apply(proteins_to_analyze_bayes, 1, function(row) {
      run_bayesian_analysis_single_protein_v7(
        protein_row_data = as.list(row),
        cond1_ibaq_cols = ibaq_cols_cond1_current_file,
        cond2_ibaq_cols = ibaq_cols_cond2_current_file,
        cond1_label = condition1_name, 
        cond2_label = condition2_name, 
        fold_change_threshold_log2 = log2_fc_cutoff
      )
    }, future.seed = TRUE)
    cat(paste0("   DEBUG: future_apply call completed for method ", method_name, ".\n"))
    # --- END OF CRITICAL CALL ---
    
    bayesian_results_df <- bind_rows(bayesian_results_list)
    cat(paste0("   DEBUG: Structure of bayesian_results_df for ", method_name, ":\n"))
    print(str(bayesian_results_df))
    cat(paste0("   DEBUG: Head of bayesian_results_df for ", method_name, ":\n"))
    print(head(bayesian_results_df))
    if (nrow(bayesian_results_df) == 0 || ncol(bayesian_results_df) == 0) {
      cat(paste0("   WARNING: bayesian_results_df for ", method_name, " is empty (0 rows or 0 columns) after bind_rows.\n"))
    }
    
    if (protein_id_col != "Protein.IDs") {
      bayesian_results_df <- bayesian_results_df %>%
        rename("Protein.IDs" = "protein_id")
    } else {
      bayesian_results_df <- bayesian_results_df %>%
        rename("Protein.IDs" = "protein_id")
    }
    
    cat(paste0("   DEBUG: Attempting left_join for method ", method_name, ".\n"))
    cat(paste0("   DEBUG: Columns in res_df_freq: ", paste(colnames(res_df_freq), collapse=", "), "\n"))
    cat(paste0("   DEBUG: Columns in bayesian_results_df: ", paste(colnames(bayesian_results_df), collapse=", "), "\n"))
    
    if ("Protein.IDs" %in% colnames(res_df_freq) && "Protein.IDs" %in% colnames(bayesian_results_df)) {
      combined_results_df <- left_join(res_df_freq, bayesian_results_df, by = "Protein.IDs")
      cat(paste0("   DEBUG: Joined by 'Protein.IDs'.\n"))
    } else if ("Protein.ID" %in% colnames(res_df_freq) && "Protein.IDs" %in% colnames(bayesian_results_df)) {
      combined_results_df <- left_join(res_df_freq, bayesian_results_df, by = c(Protein.ID = "Protein.IDs"))
      cat(paste0("   DEBUG: Joined by 'Protein.ID' (from res_df_freq) and 'Protein.IDs' (from bayesian_results_df).\n"))
    } else {
      cat(paste0("   ERROR: Cannot find a suitable column for joining. res_df_freq has 'Protein.IDs': ", "Protein.IDs" %in% colnames(res_df_freq), ", 'Protein.ID': ", "Protein.ID" %in% colnames(res_df_freq), ". bayesian_results_df has 'Protein.IDs': ", "Protein.IDs" %in% colnames(bayesian_results_df), ".\n"))
      combined_results_df <- data.frame() 
    }
    
    cat(paste0("   DEBUG: Structure of combined_results_df for ", method_name, " BEFORE writing:\n"))
    print(str(combined_results_df))
    cat(paste0("   DEBUG: Head of combined_results_df for ", method_name, " BEFORE writing:\n"))
    print(head(combined_results_df))
    if (nrow(combined_results_df) == 0 || ncol(combined_results_df) == 0) {
      cat(paste0("   WARNING: combined_results_df for ", method_name, " is empty (0 rows or 0 columns) BEFORE writing to file.\n"))
    }
    
    output_file_name <- paste0("W", work_number, "_biolrel_", method_name, "_Bayes.txt")
    output_file_path <- file.path(current_biol_rel_folder_path, output_file_name)
    write_tsv(combined_results_df, file = output_file_path)
    cat(paste("Bayesian results (all proteins) written to:", output_file_path, "\n"))
    
    proteins_bayes_filtered <- combined_results_df %>%
      filter(prob_biol_relevance >= threshold_prob_biol_relevance)
    
    if (nrow(proteins_bayes_filtered) > 0) {
      proteins_passing_bayes_filter[[method_name]] <- unique(proteins_bayes_filtered[[protein_id_col]])
      cat(paste("Method", method_name, ":", length(proteins_passing_bayes_filter[[method_name]]), "proteins passed Bayesian filter threshold for overlap plots.\n"))
    } else {
      proteins_passing_bayes_filter[[method_name]] <- character(0)
      cat(paste("Method", method_name, ": No proteins passed Bayesian filter threshold for overlap plots.\n"))
    }
  }
  
  cat(paste("--- END: Bayesian filter for Work ", work_number, " ---\n"))
  
  # --- 5. Overlap and Visualization ---
  
  all_filtered_proteins <- list(
    FC_tstudent = proteins_passing_fc_filter$tstudent,
    FC_twelch = proteins_passing_fc_filter$twelch,
    FC_limma = proteins_passing_fc_filter$limma,
    FC_deqms = proteins_passing_fc_filter$deqms,
    FC_msstats = proteins_passing_fc_filter$msstats,
    Bayes_tstudent = proteins_passing_bayes_filter$tstudent,
    Bayes_twelch = proteins_passing_bayes_filter$twelch,
    Bayes_limma = proteins_passing_bayes_filter$limma,
    Bayes_deqms = proteins_passing_bayes_filter$deqms,
    Bayes_msstats = proteins_passing_bayes_filter$msstats
  )
  
  all_filtered_proteins <- purrr::discard(all_filtered_proteins, ~ is.null(.) || length(.) == 0)
  
  if (length(all_filtered_proteins) == 0) {
    cat("INFO: No proteins passed any filters across any method. Skipping overlap plots.\n")
    return(invisible(NULL))
  }
  
  if (length(all_filtered_proteins) >= 2 && length(all_filtered_proteins) <= 5) {
    cat(paste0("\nGenerating Venn Diagram for Work ", work_number, "...\n"))
    venn_output_path <- file.path(current_biol_comparison_folder_path, paste0("W", work_number, "_Comparison_Venn.png"))
    
    num_sets <- length(all_filtered_proteins)
    venn_colors <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#FFC700", "#FF6F61", "#6A5ACD", "#20B2AA", "#FFA07A", "#8A2BE2")[1:num_sets]
    
    tryCatch({
      venn.plot <- venn.diagram(
        x = all_filtered_proteins,
        category.names = names(all_filtered_proteins),
        filename = NULL,
        output = TRUE,
        imagetype = "png",
        height = 1200,
        width = 1200,
        resolution = 300,
        compression = "lzw",
        lwd = 2,
        col = venn_colors,
        fill = alpha(venn_colors, 0.3),
        cex = 0.8,
        fontfamily = "sans",
        cat.cex = 0.7,
        cat.fontfamily = "sans",
        cat.default.pos = "outer",
        margin = 0.05
      )
      
      png(filename = venn_output_path, width = 1200, height = 1200, units = "px", res = 300)
      grid.draw(venn.plot)
      dev.off()
      cat(paste("Venn Diagram saved to:", venn_output_path, "\n"))
    }, error = function(e) {
      cat(paste("ERROR generating Venn Diagram:", e$message, "\n"))
    })
  } else {
    cat(paste0("Skipping Venn Diagram for Work ", work_number, " (Number of sets: ", length(all_filtered_proteins), ", Venn supports 2-5 sets).\n"))
  }
  
  if (length(all_filtered_proteins) >= 2) {
    cat(paste0("\nGenerating UpSet Plot for Work ", work_number, "...\n"))
    upset_output_path <- file.path(current_biol_comparison_folder_path, paste0("W", work_number, "_Comparison_UpSet.png"))
    
    tryCatch({
      png(filename = upset_output_path, width = 1200, height = 800, units = "px", res = 300)
      upset(
        fromList(all_filtered_proteins),
        nsets = length(all_filtered_proteins),
        order.by = "freq",
        decreasing = TRUE,
        point.size = 2,
        line.size = 1,
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2),
        mb.ratio = c(0.7, 0.3)
      )
      dev.off()
      cat(paste("UpSet Plot saved to:", upset_output_path, "\n"))
    }, error = function(e) {
      cat(paste("ERROR generating UpSet Plot:", e$message, "\n"))
    })
  } else {
    cat(paste0("Skipping UpSet Plot for Work ", work_number, " (Number of sets: ", length(all_filtered_proteins), ", UpSet requires at least 2 sets).\n"))
  }
  
  cat(paste0("\n--- FINISHED OVERLAP AND FC ANALYSIS FOR WORK_", work_number, " ---\n"))
}

# --- Execute analysis for Work 1 ---
# Call the main function with appropriate work number, condition names, and base path
# IMPORTANT: For Work 1, conditions are "Cold" and "RT" as per LFQ.intensity column names.
process_overlap_analysis(work_number = 1,
                         condition1_name = "Cold",  # Corrected for W1
                         condition2_name = "RT",    # Corrected for W1
                         base_parent_path = base_data_parent_path)

# Example calls for other works based on your description (uncomment and run as needed):
process_overlap_analysis(work_number = 2,
                          condition1_name = "u",
                          condition2_name = "y",
                          base_parent_path = base_data_parent_path)

process_overlap_analysis(work_number = 3,
                          condition1_name = "NaCl",
                          condition2_name = "Control",
                          base_parent_path = base_data_parent_path)

process_overlap_analysis(work_number = 4,
                          condition1_name = "BNF",
                          condition2_name = "C",
                          base_parent_path = base_data_parent_path)

 process_overlap_analysis(work_number = 5,
                          condition1_name = "C",
                          condition2_name = "PM",
                          base_parent_path = base_data_parent_path)


cat("\nScript finished.\n")








# Ensure necessary libraries are installed and loaded
# install.packages("UpSetR")
# install.packages("tidyverse") # For dplyr and other data manipulation functions
library(UpSetR)
library(tidyverse)

# --- Variable and Function Definitions (ADAPTED TO YOUR ENVIRONMENT) ---

# Base path where W1, W2, etc. folders are located
# IMPORTANT: Confirm this path is EXACTLY where 'W1', 'W2', etc. reside
base_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/7_Especial_IJMS/analisis27062025"

# List of work numbers to process (you can adjust this)
work_numbers_to_process <- c(1, 2, 3, 4, 5) # To process all
# Or just a specific number for testing
# work_numbers_to_process <- c(3) # For quick test of W3 only

# List of tests for which you want to generate overlaps
# Ensure these names match your file names (e.g., "tstudent", "bayes", "limma", "deqms", "msstats", "twelch")
tests_for_overlap_analysis <- c("tstudent", "twelch", "msstats", "bayes", "limma", "deqms") 

# Function to get the protein ID column name based on the test
get_id_col_name <- function(test_name) {
  if (test_name == "msstats") {
    return("Protein.ID") # For msstats, it's "Protein.ID"
  } else {
    return("Protein.IDs") # For others, it's "Protein.IDs"
  }
}

# --- Main loop for each work number (W1, W2, etc.) ---
for (work_number in work_numbers_to_process) {
  cat(paste0("\n################################################################################\n"))
  cat(paste0("### Starting processing for Work: W", work_number, " ###\n"))
  cat(paste0("################################################################################\n"))
  
  # Paths for the current work's specific folders
  ht_results_folder_path_current <- file.path(base_path, paste0("W", work_number), paste0("W", work_number, "_H_testing"))
  biol_rel_folder_path_current <- file.path(base_path, paste0("W", work_number), paste0("W", work_number, "_BioRel"))
  biol_comparison_folder_path_current <- file.path(biol_rel_folder_path_current, paste0("W", work_number, "_BioRel_comparison"))
  
  cat(paste0("DEBUG PATHS: HT Folder: ", ht_results_folder_path_current, "\n"))
  cat(paste0("DEBUG PATHS: BioRel Folder: ", biol_rel_folder_path_current, "\n"))
  cat(paste0("DEBUG PATHS: Comparison Folder: ", biol_comparison_folder_path_current, "\n"))
  
  # Create the comparison folder if it doesn't exist
  if (!dir.exists(biol_comparison_folder_path_current)) {
    dir.create(biol_comparison_folder_path_current, recursive = TRUE)
    cat(paste0("Created comparison results folder: ", biol_comparison_folder_path_current, "\n"))
  }
  
  # --- 5. Biologically Relevant Protein Overlap Comparisons and Plots ---
  
  cat(paste("\n### Generating biologically relevant protein overlap comparisons and plots (Work ", work_number, ") ###\n"))
  
  for (test_name in tests_for_overlap_analysis) { # Use the updated list
    
    tryCatch({ 
      cat(paste0("\n--- Processing overlap for test: ", test_name, " ---\n"))
      
      current_test_proteins_list <- list()
      
      # 0. Load the original proteingroups_filtered.txt file (All quantified and filtered proteins)
      file_original_proteingroups <- file.path(ht_results_folder_path_current, paste0("W", work_number, "_proteingroups_filtered.txt"))
      id_col_original_pg <- get_id_col_name(test_name) 
      cat(paste0("DEBUG: Attempting to read 'Total_Proteins' with ID col: '", id_col_original_pg, "' from: ", basename(file_original_proteingroups), " (Full Path: ", file_original_proteingroups, ")\n"))
      
      df_original_pg <- NULL 
      if (file.exists(file_original_proteingroups)) {
        df_original_pg <- tryCatch({
          read.delim(file_original_proteingroups, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
          cat(paste0("WARNING: Could not read ", basename(file_original_proteingroups), " - ", e$message, "\n"))
          return(NULL)
        })
        if (!is.null(df_original_pg) && id_col_original_pg %in% colnames(df_original_pg) && nrow(df_original_pg) > 0) {
          current_test_proteins_list[["Total_Proteins"]] <- as.character(df_original_pg[[id_col_original_pg]])
          cat(paste0("Loaded: Total_Proteins (", test_name, ") (", length(current_test_proteins_list[["Total_Proteins"]]), " proteins)\n"))
        } else {
          cat(paste0("INFO: File ", basename(file_original_proteingroups), " is empty or missing column '", id_col_original_pg, "'. 'Total_Proteins' set is empty.\n"))
        }
      } else {
        cat(paste0("WARNING: File not found: ", basename(file_original_proteingroups), ". 'Total_Proteins' set is empty.\n"))
      }
      
      # 1. Load original hypothesis testing results (significant)
      file_original_ht <- file.path(ht_results_folder_path_current, paste0("W", work_number, "_res_", test_name, ".txt"))
      id_col_original_ht <- get_id_col_name(test_name)
      cat(paste0("DEBUG: Attempting to read 'Original_Significant' with ID col: '", id_col_original_ht, "' from: ", basename(file_original_ht), " (Full Path: ", file_original_ht, ")\n"))
      
      df_original_ht <- NULL 
      if (file.exists(file_original_ht)) {
        df_original_ht <- tryCatch({
          read.delim(file_original_ht, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
          cat(paste0("WARNING: Could not read ", basename(file_original_ht), " - ", e$message, "\n"))
          return(NULL)
        })
        if (!is.null(df_original_ht) && id_col_original_ht %in% colnames(df_original_ht) && nrow(df_original_ht) > 0) {
          current_test_proteins_list[["Original_Significant"]] <- as.character(df_original_ht[[id_col_original_ht]])
          cat(paste0("Loaded: Original_Significant (", test_name, ") (", length(current_test_proteins_list[["Original_Significant"]]), " proteins)\n"))
        } else {
          cat(paste0("INFO: File ", basename(file_original_ht), " is empty or missing column '", id_col_original_ht, "'. 'Original_Significant' set is empty.\n"))
        }
      } else {
        cat(paste0("WARNING: File not found: ", basename(file_original_ht), ". 'Original_Significant' set is empty.\n"))
      }
      
      # 2. Load FC-filtered results
      file_fc_filtered <- file.path(biol_rel_folder_path_current, paste0("W", work_number, "_biolrel_", test_name, "_FC.txt"))
      id_col_fc <- get_id_col_name(test_name)
      cat(paste0("DEBUG: Attempting to read 'Filtered_by_FC' with ID col: '", id_col_fc, "' from: ", basename(file_fc_filtered), " (Full Path: ", file_fc_filtered, ")\n"))
      
      df_fc_filtered <- NULL 
      if (file.exists(file_fc_filtered)) {
        df_fc_filtered <- tryCatch({
          read.delim(file_fc_filtered, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
          cat(paste0("WARNING: Could not read ", basename(file_fc_filtered), " - ", e$message, "\n"))
          return(NULL)
        })
        if (!is.null(df_fc_filtered) && id_col_fc %in% colnames(df_fc_filtered) && nrow(df_fc_filtered) > 0) {
          current_test_proteins_list[["Filtered_by_FC"]] <- as.character(df_fc_filtered[[id_col_fc]])
          cat(paste0("Loaded: Filtered_by_FC (", test_name, ") (", length(current_test_proteins_list[["Filtered_by_FC"]]), " proteins)\n"))
        } else {
          cat(paste0("INFO: File ", basename(file_fc_filtered), " is empty or missing column '", id_col_fc, "'. 'Filtered_by_FC' set is empty.\n"))
        }
      } else {
        cat(paste0("WARNING: File not found: ", basename(file_fc_filtered), ". 'Filtered_by_FC' set is empty.\n"))
      }
      
      # 3. Load Bayesian-filtered results
      file_bayes_filtered <- file.path(biol_rel_folder_path_current, paste0("W", work_number, "_biolrel_", test_name, "_Bayes.txt"))
      id_col_bayes <- get_id_col_name(test_name)
      cat(paste0("DEBUG: Attempting to read 'Filtered_by_Bayes' with ID col: '", id_col_bayes, "' from: ", basename(file_bayes_filtered), " (Full Path: ", file_bayes_filtered, ")\n"))
      
      df_bayes_filtered <- NULL 
      if (file.exists(file_bayes_filtered)) {
        df_bayes_filtered <- tryCatch({
          read.delim(file_bayes_filtered, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
          cat(paste0("WARNING: Could not read ", basename(file_bayes_filtered), " - ", e$message, "\n"))
          return(NULL)
        })
        if (!is.null(df_bayes_filtered) && id_col_bayes %in% colnames(df_bayes_filtered) && nrow(df_bayes_filtered) > 0) {
          current_test_proteins_list[["Filtered_by_Bayes"]] <- as.character(df_bayes_filtered[[id_col_bayes]])
          cat(paste0("Loaded: Filtered_by_Bayes (", test_name, ") (", length(current_test_proteins_list[["Filtered_by_Bayes"]]), " proteins)\n"))
        } else {
          cat(paste0("INFO: File ", basename(file_bayes_filtered), " is empty or missing column '", id_col_bayes, "'. 'Filtered_by_Bayes' set is empty.\n"))
        }
      } else {
        cat(paste0("WARNING: File not found: ", basename(file_bayes_filtered), ". 'Filtered_by_Bayes' set is empty.\n"))
      }
      
      # Clean the list of empty sets for this test's UpSet plot
      current_test_proteins_list_cleaned <- current_test_proteins_list[sapply(current_test_proteins_list, function(x) length(x) > 0)]
      
      cat(paste0("DEBUG: Number of cleaned sets for UpSet plot for ", test_name, ": ", length(current_test_proteins_list_cleaned), "\n"))
      
      if (length(current_test_proteins_list_cleaned) >= 2) { # At least 2 sets are needed for UpSet
        cat(paste0("\n--- Generating UpSet Plot for ", test_name, " ---\n"))
        
        upset_plot_filename <- file.path(biol_comparison_folder_path_current, paste0("UpSet_diagram_", test_name, ".png"))
        cat(paste0("DEBUG: Attempting to create PNG device for: ", upset_plot_filename, "\n"))
        png(upset_plot_filename, width = 2000, height = 1600, res = 300)
        
        set_order <- c("Total_Proteins", "Original_Significant", "Filtered_by_FC", "Filtered_by_Bayes")
        ordered_sets <- current_test_proteins_list_cleaned[intersect(set_order, names(current_test_proteins_list_cleaned))]
        
        tryCatch({ 
          upset_plot <- UpSetR::upset(
            UpSetR::fromList(ordered_sets),
            nsets = length(ordered_sets),
            nintersects = NA,
            order.by = "freq",
            decreasing = TRUE,
            point.size = 3.5,
            line.size = 2,
            mainbar.y.label = paste0("Intersection Size (", test_name, ")"),
            sets.x.label = paste0("Number of Proteins (", test_name, ")"),
            text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
          )
          print(upset_plot)
        }, error = function(e) {
          cat(paste0("ERROR: Failed to generate main UpSet Plot for ", test_name, " (UpSetR call): ", e$message, "\n"))
        })
        dev.off() # Close the PNG device
        cat(paste0("UpSet Plot for ", test_name, " saved to: ", upset_plot_filename, "\n"))
        
        cat(paste0("\n--- Calculating intersections for ", test_name, " (tabular) ---\n"))
        all_intersections_counts_test <- list()
        set_names_for_combinations <- names(ordered_sets)
        for (i in 1:length(set_names_for_combinations)) {
          combinations <- combn(set_names_for_combinations, m = i, simplify = FALSE)
          for (combo in combinations) {
            current_intersection <- ordered_sets[[combo[1]]]
            if (length(combo) > 1) {
              for (j in 2:length(combo)) {
                current_intersection <- intersect(current_intersection, ordered_sets[[combo[j]]])
              }
            }
            combo_name <- paste(sort(combo), collapse = " & ")
            all_intersections_counts_test[[combo_name]] <- length(current_intersection)
          }
        }
        intersections_df_test <- data.frame(
          Intersection = names(all_intersections_counts_test),
          Number_of_Proteins = unlist(all_intersections_counts_test),
          stringsAsFactors = FALSE
        ) %>% arrange(desc(Number_of_Proteins))
        cat(paste0("\n--- Intersection Table for ", test_name, " ---\n"))
        print(intersections_df_test)
        intersection_table_filename <- file.path(biol_comparison_folder_path_current, paste0("Intersections_", test_name, ".txt"))
        cat(paste0("DEBUG: Attempting to write intersection table to: ", intersection_table_filename, "\n"))
        write.table(intersections_df_test, file = intersection_table_filename, sep = "\t", row.names = FALSE, quote = FALSE)
        cat(paste0("Intersection table for ", test_name, " saved to: ", intersection_table_filename, "\n"))
        
      } else {
        cat(paste0("\nCannot generate UpSet Plot for ", test_name, ": At least 2 protein sets are needed for comparison (after cleaning empty sets).\n"))
        cat("Available sets: ", paste(names(current_test_proteins_list_cleaned), collapse = ", "), "\n")
      }
    }, error = function(e) { # Broad tryCatch for the UpSet generation block for current test_name
      cat(paste0("FATAL ERROR: An unhandled error occurred during UpSet plot generation for test '", test_name, "' in work 'W", work_number, "': ", e$message, "\n"))
    }) # End broad tryCatch for current test_name UpSet generation
  } # End of loop for test_name for UpSet plots
  
  
  #--- NEW SECTION: Fold Change Analysis in Intersection Groups (CURRENT WORK) ---
  cat(paste("\n### START: Fold Change Analysis in Intersection Groups (Work ", work_number, ") ###\n"))
  
  # Loop through each test to perform FC analysis in intersection groups
  for (test_name in tests_for_overlap_analysis) { # Use the updated list
    
    tryCatch({ 
      output_comparison_file <- file.path(biol_comparison_folder_path_current, paste0("W", work_number, "_BioRel_", test_name, "_comparison.txt"))
      cat(paste0("DEBUG: Attempting to open sink to: ", output_comparison_file, "\n"))
      sink(output_comparison_file, append = FALSE) 
      cat(paste0("### Fold Change Analysis in Intersection Groups for test: ", toupper(test_name), " ###\n\n"))
      cat("------------------------------------------------------------------------------------------------------\n")
      cat(paste0("\n--- Analyzing FC of intersection groups for test: ", test_name, " ---\n"))
      
      id_col <- get_id_col_name(test_name)
      log2fc_col <- if (test_name %in% c("tstudent", "twelch", "msstats", "bayes")) "log2FC" else "logFC" 
      
      cat(paste0("DEBUG: For test '", test_name, "', ID column is '", id_col, "', Log2FC column is '", log2fc_col, "'\n"))
      
      file_original_ht <- file.path(ht_results_folder_path_current, paste0("W", work_number, "_res_", test_name, ".txt"))
      df_original_ht <- NULL
      if (file.exists(file_original_ht)) {
        df_original_ht <- tryCatch({ read.delim(file_original_ht, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) }, 
                                   error = function(e) { cat(paste0("WARNING: Could not read ", basename(file_original_ht), " for FC analysis - ", e$message, "\n")); return(NULL) })
      } else { cat(paste0("WARNING: File not found: ", basename(file_original_ht), " for FC analysis.\n")) }
      
      file_total_proteins <- file.path(ht_results_folder_path_current, paste0("W", work_number, "_proteingroups_filtered.txt"))
      df_total_proteins <- NULL
      if (file.exists(file_total_proteins)) {
        df_total_proteins <- tryCatch({ read.delim(file_total_proteins, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) }, 
                                      error = function(e) { cat(paste0("WARNING: Could not read ", basename(file_total_proteins), " for total protein IDs - ", e$message, "\n")); return(NULL) })
      } else { cat(paste0("WARNING: File not found: ", basename(file_total_proteins), " for total protein IDs.\n")) }
      
      if (is.null(df_original_ht) || nrow(df_original_ht) == 0 || !(id_col %in% colnames(df_original_ht)) || !(log2fc_col %in% colnames(df_original_ht))) {
        cat("WARNING: Original results file (significant proteins) for this test is empty or missing necessary columns. Skipping FC analysis.\n")
        sink() 
        cat(paste0("Comparison results for ", test_name, " saved to: ", output_comparison_file, "\n"))
        return(NULL)
      }
      
      if (is.null(df_total_proteins) || nrow(df_total_proteins) == 0 || !(id_col %in% colnames(df_total_proteins))) {
        cat("WARNING: Total proteins file is empty or missing the necessary ID column. Skipping FC analysis.\n")
        sink()
        cat(paste0("Comparison results for ", test_name, " saved to: ", output_comparison_file, "\n"))
        return(NULL)
      }
      
      df_original_ht <- as_tibble(df_original_ht)
      df_total_proteins <- as_tibble(df_total_proteins)
      
      ids_total <- as.character(df_total_proteins[[id_col]])
      ids_original_significant <- as.character(df_original_ht[[id_col]])
      
      file_fc_filtered <- file.path(biol_rel_folder_path_current, paste0("W", work_number, "_biolrel_", test_name, "_FC.txt"))
      df_fc_filtered <- NULL
      if (file.exists(file_fc_filtered)) {
        df_fc_filtered <- tryCatch({ read.delim(file_fc_filtered, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) }, 
                                   error = function(e) { cat(paste0("WARNING: Could not read ", basename(file_fc_filtered), " for FC IDs - ", e$message, "\n")); return(NULL) })
      } else { cat(paste0("WARNING: File not found: ", basename(file_fc_filtered), " for FC IDs.\n")) }
      ids_fc <- if (!is.null(df_fc_filtered) && id_col %in% colnames(df_fc_filtered) && nrow(df_fc_filtered) > 0) as.character(df_fc_filtered[[id_col]]) else character(0)
      
      file_bayes_filtered <- file.path(biol_rel_folder_path_current, paste0("W", work_number, "_biolrel_", test_name, "_Bayes.txt"))
      df_bayes_filtered <- NULL
      if (file.exists(file_bayes_filtered)) {
        df_bayes_filtered <- tryCatch({ read.delim(file_bayes_filtered, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) }, 
                                      error = function(e) { cat(paste0("WARNING: Could not read ", basename(file_bayes_filtered), " for Bayes IDs - ", e$message, "\n")); return(NULL) })
      } else { cat(paste0("WARNING: File not found: ", basename(file_bayes_filtered), " for Bayes IDs.\n")) }
      ids_bayes <- if (!is.null(df_bayes_filtered) && id_col %in% colnames(df_bayes_filtered) && nrow(df_bayes_filtered) > 0) as.character(df_bayes_filtered[[id_col]]) else character(0)
      
      cat(paste0("DEBUG: Before df_analysis join. IDs total: ", length(ids_total), ", Original Sig: ", length(ids_original_significant), ", FC: ", length(ids_fc), ", Bayes: ", length(ids_bayes), "\n"))
      
      df_analysis <- df_total_proteins %>%
        select(all_of(id_col)) %>%
        rename(Protein_ID = all_of(id_col)) %>%
        left_join(df_original_ht %>% select(all_of(c(id_col, log2fc_col))) %>%
                    rename(Protein_ID = all_of(id_col), log2FC = all_of(log2fc_col)),
                  by = "Protein_ID") %>%
        mutate(log2FC = as.numeric(log2FC)) %>%
        mutate(abs_log2FC = abs(log2FC))
      
      cat(paste0("DEBUG: After df_analysis join. NROW: ", nrow(df_analysis), ". First few rows:\n"))
      print(head(df_analysis))
      
      df_analysis <- df_analysis %>%
        mutate(
          is_original_significant = Protein_ID %in% ids_original_significant,
          in_FC = Protein_ID %in% ids_fc,
          in_Bayes = Protein_ID %in% ids_bayes,
          Group = case_when(
            is_original_significant & in_FC & in_Bayes ~ "Sig_FC_Bayes",
            is_original_significant & in_FC & !in_Bayes ~ "Sig_FC_only",
            is_original_significant & !in_FC & in_Bayes ~ "Sig_Bayes_only",
            is_original_significant & !in_FC & !in_Bayes ~ "Sig_Discarded_FC_Bayes",
            !is_original_significant & in_FC & in_Bayes ~ "NonSig_FC_Bayes",
            !is_original_significant & in_FC & !in_Bayes ~ "NonSig_FC_only",
            !is_original_significant & !in_FC & in_Bayes ~ "NonSig_Bayes_only",
            !is_original_significant & !in_FC & !in_Bayes ~ "NonSig_Discarded",
            TRUE ~ "Other_Unclassified"
          )
        ) %>%
        filter(Group != "Other_Unclassified")
      
      df_analysis$Group <- factor(df_analysis$Group)
      
      cat("\nDIAGNOSTIC: Structure of df_analysis after classification:\n")
      print(str(df_analysis))
      cat("\n")
      
      if (!inherits(df_analysis, "data.frame") || nrow(df_analysis) == 0 || length(unique(df_analysis$Group)) < 1) {
        cat("WARNING: No valid data or sufficient groups to perform FC comparison analysis for this test.\n")
        cat("Details: df_analysis is of class '", class(df_analysis)[1], "', has ", nrow(df_analysis), " rows, and ", length(unique(df_analysis$Group)), " unique groups.\n", sep="")
        sink()
        cat(paste0("Comparison results for ", test_name, " saved to: ", output_comparison_file, "\n"))
        return(NULL)
      }
      
      cat(paste0("\n--- Generating UpSet Plot for Discarded Significant Proteins (Test: ", test_name, ") ---\n"))
      sig_and_discarded_by_FC <- setdiff(ids_original_significant, ids_fc)
      sig_and_discarded_by_Bayes <- setdiff(ids_original_significant, ids_bayes)
      
      list_discarded_proteins <- list(
        "Sig_Discarded_by_FC" = sig_and_discarded_by_FC,
        "Sig_Discarded_by_Bayes" = sig_and_discarded_by_Bayes
      )
      
      cat(paste0("DEBUG: Counts for discarded sets (", test_name, "):\n"))
      cat(paste0("  Sig_Discarded_by_FC count: ", length(list_discarded_proteins[["Sig_Discarded_by_FC"]]), "\n"))
      cat(paste0("  Sig_Discarded_by_Bayes count: ", length(list_discarded_proteins[["Sig_Discarded_by_Bayes"]]), "\n"))
      
      intersection_discarded <- intersect(list_discarded_proteins[["Sig_Discarded_by_FC"]], list_discarded_proteins[["Sig_Discarded_by_Bayes"]])
      cat(paste0("  Intersection (Sig_Discarded_by_FC & Sig_Discarded_by_Bayes): ", length(intersection_discarded), " proteins\n"))
      only_in_fc_discarded <- setdiff(list_discarded_proteins[["Sig_Discarded_by_FC"]], list_discarded_proteins[["Sig_Discarded_by_Bayes"]])
      cat(paste0("  Only in Sig_Discarded_by_FC: ", length(only_in_fc_discarded), " proteins\n"))
      only_in_bayes_discarded <- setdiff(list_discarded_proteins[["Sig_Discarded_by_Bayes"]], list_discarded_proteins[["Sig_Discarded_by_FC"]])
      cat(paste0("  Only in Sig_Discarded_by_Bayes: ", length(only_in_bayes_discarded), " proteins\n"))
      
      list_discarded_proteins_cleaned <- list_discarded_proteins[sapply(list_discarded_proteins, function(x) length(x) > 0)]
      
      if (length(list_discarded_proteins_cleaned) >= 2) {
        upset_discarded_filename <- file.path(biol_comparison_folder_path_current, paste0("UpSet_discarded_significant_", test_name, ".png"))
        cat(paste0("DEBUG: Attempting to create PNG device for discarded plot: ", upset_discarded_filename, "\n"))
        png(upset_discarded_filename, width = 1500, height = 1200, res = 200)
        set_order_discarded <- c("Sig_Discarded_by_FC", "Sig_Discarded_by_Bayes")
        ordered_sets_discarded <- list_discarded_proteins_cleaned[intersect(set_order_discarded, names(list_discarded_proteins_cleaned))]
        
        tryCatch({ 
          upset_plot_discarded <- UpSetR::upset(
            UpSetR::fromList(ordered_sets_discarded),
            nsets = length(ordered_sets_discarded),
            nintersects = NA,
            order.by = "freq",
            decreasing = TRUE,
            point.size = 3.5,
            line.size = 2,
            mainbar.y.label = paste0("Intersection Size (Significant Discarded, ", test_name, ")"),
            sets.x.label = paste0("Number of Proteins (Significant Discarded, ", test_name, ")"),
            text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
          )
          print(upset_plot_discarded)
        }, error = function(e) {
          cat(paste0("ERROR: Failed to generate discarded UpSet Plot for ", test_name, " (UpSetR call): ", e$message, "\n"))
        })
        dev.off()
        cat(paste0("UpSet Plot of Significant Discarded Proteins for ", test_name, " saved to: ", upset_discarded_filename, "\n"))
      } else {
        cat(paste0("INFO: Cannot generate UpSet Plot of Significant Discarded Proteins for ", test_name, ": Less than 2 non-empty protein sets are available for comparison.\n"))
        cat("  This means there might be no significant proteins discarded by one or both filters, or all discarded proteins are shared.\n")
        cat("  See the 'DEBUG: Counts for discarded sets' above for detailed numbers.\n")
        
        discarded_counts_df <- data.frame(
          Category = names(list_discarded_proteins),
          Count = sapply(list_discarded_proteins, length),
          stringsAsFactors = FALSE
        )
        cat("\nSummary of discarded significant protein counts:\n")
        print(discarded_counts_df)
        discarded_counts_filename <- file.path(biol_comparison_folder_path_current, paste0("Discarded_Significant_Counts_", test_name, ".txt"))
        cat(paste0("DEBUG: Attempting to write discarded counts to: ", discarded_counts_filename, "\n"))
        write.table(discarded_counts_df, 
                    file = discarded_counts_filename, 
                    sep = "\t", row.names = FALSE, quote = FALSE)
        cat(paste0("Summary of discarded significant proteins for ", test_name, " saved to: ", discarded_counts_filename, "\n"))
      }
      
      cat("\n--- Fold Change Statistics by Intersection Group ---\n")
      
      df_fc_analysis <- df_analysis %>% filter(!is.na(log2FC) & !is.infinite(log2FC)) 
      
      cat(paste0("DEBUG: NROW of df_fc_analysis after NA/Inf filter: ", nrow(df_fc_analysis), "\n"))
      
      if (nrow(df_fc_analysis) > 0 && length(unique(df_fc_analysis$Group)) > 1) {
        summary_stats <- df_fc_analysis %>%
          group_by(Group) %>%
          summarise(
            Count = n(),
            Mean_abs_log2FC = mean(abs_log2FC, na.rm = TRUE),
            Median_abs_log2FC = median(abs_log2FC, na.rm = TRUE),
            SD_abs_log2FC = sd(abs_log2FC, na.rm = TRUE),
            Min_abs_log2FC = min(abs_log2FC, na.rm = TRUE),
            Max_abs_log2FC = max(abs_log2FC, na.rm = TRUE)
          ) %>%
          arrange(desc(Count))
        
        cat("\nSummary of Absolute Fold Change by Group:\n")
        print(summary_stats)
        cat("\n")
        
        groups_with_enough_data <- df_fc_analysis %>%
          group_by(Group) %>%
          filter(n() >= 2) %>%
          ungroup()
        
        cat(paste0("DEBUG: Number of groups with enough data for ANOVA: ", length(unique(groups_with_enough_data$Group)), "\n"))
        
        if (length(unique(groups_with_enough_data$Group)) >= 2) {
          tryCatch({
            groups_with_enough_data$Group <- factor(groups_with_enough_data$Group)
            anova_result <- aov(abs_log2FC ~ Group, data = groups_with_enough_data)
            cat("\nANOVA Results (Absolute Fold Change vs. Group):\n")
            print(summary(anova_result))
            
            if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
              cat("\nPost-Hoc Tests (Tukey HSD) if ANOVA is significant:\n")
              tukey_hsd <- TukeyHSD(anova_result)
              print(tukey_hsd)
            }
          }, error = function(e) {
            cat(paste0("WARNING: Could not perform ANOVA or Tukey HSD for ", test_name, " due to: ", e$message, "\n"))
            cat("This might be because there is not enough variation or data in the groups for the analysis.\n")
          })
        } else {
          cat("Not enough groups or data within groups (minimum 2 obs. per group) to perform ANOVA for ", test_name, ".\n")
        }
        
        boxplot_filename <- file.path(biol_comparison_folder_path_current, paste0("Boxplot_absFC_by_Group_", test_name, ".png"))
        cat(paste0("DEBUG: Attempting to create PNG device for boxplot: ", boxplot_filename, "\n"))
        png(boxplot_filename, width = 1200, height = 1000, res = 200)
        tryCatch({ 
          boxplot_plot <- ggplot(df_fc_analysis, aes(x = Group, y = abs_log2FC, fill = Group)) +
            geom_boxplot() +
            labs(title = paste0("Distribution of Absolute Log2FC by Intersection Group (", test_name, ")"),
                 x = "Group",
                 y = "Absolute Log2FC") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          print(boxplot_plot)
        }, error = function(e) {
          cat(paste0("ERROR: Failed to generate boxplot for ", test_name, " (ggplot2 call): ", e$message, "\n"))
        })
        dev.off()
        cat(paste0("Boxplot of Absolute Fold Change by Group for ", test_name, " saved to: ", boxplot_filename, "\n"))
        
      } else {
        cat("No sufficient Fold Change data or groups to perform statistical analysis for this test.\n")
      }
      
      sink() 
      cat(paste0("Comparison results for ", test_name, " saved to: ", output_comparison_file, "\n"))
      
    }, error = function(e) { 
      cat(paste0("FATAL ERROR: An unhandled error occurred during FC analysis for test '", test_name, "' in work 'W", work_number, "': ", e$message, "\n"))
      cat("Attempting to close sink if open and continue to next test/work.\n")
      if (sink.number() > 0) { 
        tryCatch({ sink() }, error = function(e_sink) { cat(paste0("Error closing sink: ", e_sink$message, "\n")) })
      }
    }) 
  } 
  
  cat(paste("\n### END: Fold Change Analysis in Intersection Groups (Work ", work_number, ") ###\n"))
  cat(paste0("################################################################################\n"))
} 
cat(paste0("### Processing completed for all specified works. ###\n"))
cat(paste0("################################################################################\n"))


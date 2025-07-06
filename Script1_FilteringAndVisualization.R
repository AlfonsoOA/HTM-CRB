# --- INSTALLATION AND LOADING OF LIBRARIES (Global) ---
# Install packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car")
}
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan") # For NMDS and ANOSIM
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") # For heatmaps
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot") # For theme_cowplot and no secondary axes
}
if (!requireNamespace("showtext", quietly = TRUE)) {
  install.packages("showtext") # For custom fonts
}


# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(vegan)
library(pheatmap)
library(cowplot) # For clean themes
library(showtext) # For custom fonts

# Add Open Sans font for plotting
# Open Sans is a Google Font, so showtext can download and use it automatically.
font_add_google("Open Sans", "Open Sans")
showtext_auto()

# Define global plot settings
# AUMENTO AGRESIVO DE LOS TAMAÑOS DE FUENTE
base_font_size <- 24 # Aumentado significativamente
title_font_size <- 30 # Aumentado
axis_title_font_size <- 28 # Aumentado
axis_text_font_size <- 24 # Aumentado
dpi_resolution <- 600 # Mantener alta resolución para nitidez al reducir

# --- MAIN FUNCTION TO PROCESS AND ANALYZE A WORK ---
process_and_analyze_proteingroups <- function(work_name, base_path, condition_patterns, peptide_col_prefix = "Peptides.", intensity_col_prefix = "LFQ.intensity.") {
  
  cat(paste0("\n--- PROCESSING WORK: ", work_name, " ---\n"))
  
  # --- INITIAL CONFIGURATION AND FOLDER PREPARATION ---
  main_folder_path <- file.path(base_path, work_name)
  output_folder_name <- paste0("1_ProteinGroups_Analysis_", work_name)
  output_folder_path <- file.path(main_folder_path, output_folder_name)
  
  # Create the folder if it doesn't exist
  if (!dir.exists(output_folder_path)) {
    dir.create(output_folder_path, recursive = TRUE)
    cat(paste("Folder created:", output_folder_path, "\n"))
  } else {
    cat(paste("Folder already exists:", output_folder_path, "\n"))
  }
  
  # Define input and output file paths
  proteinGroups_file <- file.path(main_folder_path, paste0(work_name, "_proteinGroups.txt"))
  output_file_step3 <- file.path(output_folder_path, paste0("1_proteingroups_filtered_", work_name, ".txt"))
  output_file_variability <- file.path(output_folder_path, paste0("2_proteingroups_variability_", work_name, ".txt"))
  boxplot_cv_file <- file.path(output_folder_path, paste0("CV_Distribution_Boxplot_", work_name, ".png"))
  hist_log2fc_file <- file.path(output_folder_path, paste0("Log2FC_Distribution_Histogram_", work_name, ".png"))
  heatmap_file <- file.path(output_folder_path, paste0("Heatmap_LFQ_intensity_", work_name, ".png"))
  nmds_file <- file.path(output_folder_path, paste0("NMDS_LFQ_intensity_", work_name, ".png"))
  anosim_file <- file.path(output_folder_path, paste0("ANOSIM_Results_", work_name, ".txt"))
  
  # --- PROTEINGROUPS CLEANING AND FILTERING ---
  cat(paste0("\n--- STARTING PROTEINGROUPS CLEANING AND FILTERING (", work_name, ") ---\n"))
  
  # 1. Read the proteinGroups.txt file
  if (!file.exists(proteinGroups_file)) {
    # Usar stop() aquí es apropiado porque el archivo base no existe.
    stop(paste0("ERROR: The file '", basename(proteinGroups_file), "' not found for ", work_name, ". Please check the path."))
  } else {
    cat(paste("Reading file:", proteinGroups_file, "\n"))
  }
  protein_groups <- read.delim(proteinGroups_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Check if the file is empty or has issues reading
  if (nrow(protein_groups) == 0) {
    # Usar stop() aquí también es apropiado porque no hay datos que procesar.
    stop(paste0("ERROR: The file '", basename(proteinGroups_file), "' for ", work_name, " was read but is empty. Cannot proceed."))
  }
  
  # 3. Separate proteins by IDs and duplicate information
  protein_groups_separated <- protein_groups %>%
    separate_rows(`Protein.IDs`, sep = ";")
  
  # 4. Remove contaminants and decoys
  protein_groups_filtered <- protein_groups_separated %>%
    filter(!grepl("^CON__", `Protein.IDs`)) %>%
    filter(!grepl("^REV_", `Protein.IDs`))
  
  # --- KEY ADJUSTMENT HERE: Dynamically identify LFQ intensity and Peptides columns ---
  # Condition patterns are defined at the start of the function
  
  # Build regular expressions for peptide and intensity columns
  all_peptide_cols_pattern <- paste0("^", peptide_col_prefix, "(", paste(unlist(condition_patterns), collapse = "|"), ")", ".*")
  all_intensity_cols_pattern <- paste0("^", intensity_col_prefix, "(", paste(unlist(condition_patterns), collapse = "|"), ")", ".*")
  
  replicas_peptides_cols <- grep(all_peptide_cols_pattern, colnames(protein_groups_filtered), value = TRUE)
  replicas_intensity_cols <- grep(all_intensity_cols_pattern, colnames(protein_groups_filtered), value = TRUE)
  
  # --- IDENTIFIED COLUMNS VERIFICATION (DEBUG) ---
  cat("\n--- Peptide Columns found by RegEx ---\n")
  print(replicas_peptides_cols)
  cat("\n--- LFQ Intensity Columns found by RegEx ---\n")
  print(replicas_intensity_cols)
  cat("\n--------------------------------------------------\n")
  
  if (length(replicas_peptides_cols) == 0 || length(replicas_intensity_cols) == 0) {
    # Este es un caso crítico: si no hay columnas para analizar, NO se puede continuar.
    stop(paste0("ERROR: No peptide or LFQ intensity columns found for the specified conditions in work ", work_name, ". Cannot proceed with any analysis."))
  }
  
  # Ensure all relevant columns are numeric
  protein_groups_filtered <- protein_groups_filtered %>%
    mutate(across(all_of(replicas_peptides_cols), as.numeric)) %>%
    mutate(across(all_of(replicas_intensity_cols), as.numeric))
  
  # 5. Set intensity to 0 for proteins identified by a single peptide PER REPLICA
  protein_groups_filtered <- protein_groups_filtered %>%
    rowwise() %>%
    mutate(
      across(all_of(replicas_intensity_cols), ~ {
        peptide_col_name <- gsub(intensity_col_prefix, peptide_col_prefix, cur_column())
        if (peptide_col_name %in% colnames(protein_groups_filtered)) {
          if (!is.na(get(peptide_col_name)) && get(peptide_col_name) == 1) {
            return(0)
          }
        }
        return(.) # Return original value if condition not met or peptide column doesn't exist
      })
    ) %>%
    ungroup()
  
  # 5b. Filter to remove proteins with no identification in any replica (based on LFQ intensity)
  # Esto es un filtro importante: si una proteína tiene 0 o NA en TODAS sus réplicas, se elimina.
  protein_groups_filtered <- protein_groups_filtered %>%
    filter(rowSums(select(., all_of(replicas_intensity_cols)), na.rm = TRUE) > 0)
  
  # 6. Write the filtered file to the created folder
  write.table(protein_groups_filtered, file = output_file_step3, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(paste0("Filtered file '", basename(output_file_step3), "' created in: ", output_folder_path, "\n"))
  cat(paste("Number of proteins after filtering (", work_name, "):", nrow(protein_groups_filtered), "\n"))
  
  # --- VARIABILITY ANALYSIS ON FILTERED PROTEINGROUPS ---
  cat(paste0("\n--- STARTING VARIABILITY ANALYSIS (", work_name, ") ---\n"))
  
  # 1. Read the filtered proteinGroups file
  if (!file.exists(output_file_step3)) {
    # Este caso debería ser raro si el save anterior funcionó, pero lo mantenemos.
    stop(paste0("ERROR: Filtered file '", basename(output_file_step3), "' not found. This should not happen. Skipping variability analysis, heatmap, or NMDS/ANOSIM.\n"))
  } else {
    cat(paste("Reading filtered file:", output_file_step3, "\n"))
  }
  protein_groups <- read.delim(output_file_step3, sep = "\t", stringsAsFactors = FALSE)
  
  # Check if data frame is empty after filtering
  if (nrow(protein_groups) == 0) {
    # Este es un punto crítico: si no hay proteínas después del filtrado, no se puede hacer nada.
    stop(paste0("ERROR: Filtered data for work ", work_name, " is empty after cleaning. Cannot proceed with variability analysis, heatmap, or NMDS/ANOSIM."))
  }
  
  # 2. Identify LFQ intensity columns for each condition
  # Assume the first condition is the "denominator" and the second is the "numerator" for Log2FC
  cond1_pattern <- condition_patterns[[1]] # Denominator
  cond2_pattern <- condition_patterns[[2]] # Numerator
  
  lfq_cols_cond1 <- grep(paste0("^", intensity_col_prefix, cond1_pattern, ".*"), colnames(protein_groups), value = TRUE)
  lfq_cols_cond2 <- grep(paste0("^", intensity_col_prefix, cond2_pattern, ".*"), colnames(protein_groups), value = TRUE)
  
  # --- LFQ INTENSITY COLUMNS VERIFICATION FOR VARIABILITY (DEBUG) ---
  cat(paste0("\n--- LFQ Intensity Columns for ", cond1_pattern, " (variability analysis) ---\n"))
  print(lfq_cols_cond1)
  cat(paste0("\n--- LFQ Intensity Columns for ", cond2_pattern, " (variability analysis) ---\n"))
  print(lfq_cols_cond2)
  cat("\n--------------------------------------------------------------\n")
  
  if (length(lfq_cols_cond1) == 0 || length(lfq_cols_cond2) == 0) {
    # Este es otro punto crítico. Si no se encuentran las columnas de intensidad esperadas, no se puede proceder.
    stop(paste0("ERROR: No LFQ Intensity columns found for one or both conditions (", cond1_pattern, ", ", cond2_pattern, ") in work ", work_name, ". Variability analysis, heatmap, or NMDS/ANOSIM will not be performed."))
  }
  
  # Function to get LFQ intensity values for a condition
  get_lfq_values <- function(row, cols) {
    as.numeric(row[cols])
  }
  
  # Function to calculate Coefficient of Variation (CV) in percentage
  calculate_cv_percent <- function(row, cols) {
    intensities <- as.numeric(row[cols])
    # Consider only values > 0 for CV calculation. NAs are handled by na.rm=TRUE in mean/sd.
    valid_intensities <- intensities[intensities > 0 & !is.na(intensities)]
    if (length(valid_intensities) > 1) { # Need at least 2 values to calculate standard deviation
      mean_intensity <- mean(valid_intensities, na.rm = TRUE)
      sd_intensity <- sd(valid_intensities, na.rm = TRUE)
      if (!is.na(mean_intensity) && mean_intensity > 0) {
        return((sd_intensity / mean_intensity) * 100)
      } else {
        return(NA) # Return NA if mean is 0 or NA
      }
    } else {
      return(NA) # Return NA if no valid values (less than 2)
    }
  }
  
  # Calculate CVs
  cv_cond1_percent <- apply(protein_groups, 1, calculate_cv_percent, cols = lfq_cols_cond1)
  cv_cond2_percent <- apply(protein_groups, 1, calculate_cv_percent, cols = lfq_cols_cond2)
  
  # Calculate means and Log2FC (Cond2/Cond1)
  mean_lfq_cond1 <- rowMeans(protein_groups[, lfq_cols_cond1], na.rm = TRUE)
  mean_lfq_cond2 <- rowMeans(protein_groups[, lfq_cols_cond2], na.rm = TRUE)
  
  # Handle zero values for Log2FC by adding a small epsilon
  epsilon <- 0
  # Calculate epsilon only if there are positive values in the LFQ intensity dataset
  # This epsilon helps avoid -Inf or Inf from log2(0)
  if(any(protein_groups[, c(lfq_cols_cond1, lfq_cols_cond2)] > 0, na.rm = TRUE)) {
    epsilon <- min(protein_groups[, c(lfq_cols_cond1, lfq_cols_cond2)][protein_groups[, c(lfq_cols_cond1, lfq_cols_cond2)] > 0], na.rm = TRUE) / 100
  }
  
  mean_lfq_cond1_adjusted <- ifelse(mean_lfq_cond1 == 0, epsilon, mean_lfq_cond1)
  mean_lfq_cond2_adjusted <- ifelse(mean_lfq_cond2 == 0, epsilon, mean_lfq_cond2)
  
  log2fc <- log2(mean_lfq_cond2_adjusted / mean_lfq_cond1_adjusted)
  
  # Perform statistical tests PER PROTEIN
  shapiro_cond1 <- apply(protein_groups, 1, function(row) {
    values <- get_lfq_values(row, lfq_cols_cond1)
    values_finite <- values[is.finite(values) & values > 0]
    if (length(values_finite) >= 3 && length(values_finite) <= 5000) return(shapiro.test(values_finite)$p.value) else return(NA)
  })
  shapiro_cond2 <- apply(protein_groups, 1, function(row) {
    values <- get_lfq_values(row, lfq_cols_cond2)
    values_finite <- values[is.finite(values) & values > 0]
    if (length(values_finite) >= 3 && length(values_finite) <= 5000) return(shapiro.test(values_finite)$p.value) else return(NA)
  })
  levene_p_value <- apply(protein_groups, 1, function(row) {
    cond1_values <- get_lfq_values(row, lfq_cols_cond1)
    cond2_values <- get_lfq_values(row, lfq_cols_cond2)
    combined_values <- c(cond1_values[is.finite(cond1_values) & cond1_values > 0], cond2_values[is.finite(cond2_values) & cond2_values > 0])
    group <- factor(c(rep(cond1_pattern, sum(is.finite(cond1_values) & cond1_values > 0)), rep(cond2_pattern, sum(is.finite(cond2_values) & cond2_values > 0))))
    
    if (length(combined_values) >= 2 && length(levels(group)) == 2) {
      if (sum(group == cond1_pattern) > 0 && sum(group == cond2_pattern) > 0) {
        # Check for sufficient unique values per group for Levene's test
        if (length(unique(cond1_values[is.finite(cond1_values) & cond1_values > 0])) >= 2 &&
            length(unique(cond2_values[is.finite(cond2_values) & cond2_values > 0])) >= 2) {
          return(leveneTest(combined_values, group)$`Pr(>F)`[1])
        } else {
          return(NA) # Not enough unique values in one or both groups
        }
      } else {
        return(NA) # If one of the groups is empty, test cannot be performed
      }
    } else {
      return(NA)
    }
  })
  
  # Create variability data frame
  variability_data <- data.frame(
    Protein.IDs = protein_groups$Protein.IDs,
    CV_Cond1_Percent = cv_cond1_percent,
    CV_Cond2_Percent = cv_cond2_percent,
    Log2FC = log2fc,
    Shapiro_Wilk_p_Cond1 = shapiro_cond1,
    RShapiro_Wilk_p_Cond2 = shapiro_cond2,
    Levene_p_value = levene_p_value
  )
  colnames(variability_data)[2] <- paste0("CV_", cond1_pattern, "_Percent")
  colnames(variability_data)[3] <- paste0("CV_", cond2_pattern, "_Percent")
  colnames(variability_data)[4] <- paste0("Log2FC_", cond2_pattern, "_vs_", cond1_pattern)
  colnames(variability_data)[5] <- paste0("Shapiro_Wilk_p_", cond1_pattern)
  colnames(variability_data)[6] <- paste0("Shapiro_Wilk_p_", cond2_pattern)
  
  # Calculate global variability
  median_cv_cond1 <- median(variability_data[[paste0("CV_", cond1_pattern, "_Percent")]], na.rm = TRUE)
  median_cv_cond2 <- median(variability_data[[paste0("CV_", cond2_pattern, "_Percent")]], na.rm = TRUE)
  log2fc_finite <- variability_data[[paste0("Log2FC_", cond2_pattern, "_vs_", cond1_pattern)]][is.finite(variability_data[[paste0("Log2FC_", cond2_pattern, "_vs_", cond1_pattern)]])]
  sd_abs_log2fc <- sd(abs(log2fc_finite), na.rm = TRUE)
  
  cat("\n--- Global Intrareplicate Variability (Median CV%) ---\n")
  cat(paste0("Median CV (", cond1_pattern, "): ", round(median_cv_cond1, 2), "%\n"))
  cat(paste0("Median CV (", cond2_pattern, "): ", round(median_cv_cond2, 2), "%\n"))
  
  cat("\n--- Global Interreplicate Variability (SD of abs(Log2FC)) ---\n")
  cat(paste0("Standard Deviation of |Log2FC (", cond2_pattern, " vs ", cond1_pattern, ")|: ", round(sd_abs_log2fc, 2), "\n"))
  
  # Perform global statistical tests per condition
  all_lfq_cond1 <- unlist(protein_groups[, lfq_cols_cond1])
  all_lfq_cond2 <- unlist(protein_groups[, lfq_cols_cond2])
  all_finite_cond1 <- all_lfq_cond1[is.finite(all_lfq_cond1) & all_lfq_cond1 > 0]
  all_finite_cond2 <- all_lfq_cond2[is.finite(all_lfq_cond2) & all_lfq_cond2 > 0]
  
  shapiro_global_cond1_p <- ifelse(length(all_finite_cond1) >= 3 && length(all_finite_cond1) <= 5000, shapiro.test(all_finite_cond1)$p.value, NA)
  shapiro_global_cond2_p <- ifelse(length(all_finite_cond2) >= 3 && length(all_finite_cond2) <= 5000, shapiro.test(all_finite_cond2)$p.value, NA)
  
  global_levene_p <- NA
  if (length(all_finite_cond1) > 0 && length(all_finite_cond2) > 0) {
    combined_all <- c(all_finite_cond1, all_finite_cond2)
    group_all <- factor(c(rep(cond1_pattern, length(all_finite_cond1)), rep(cond2_pattern, length(all_finite_cond2))))
    if (length(levels(group_all)) == 2) {
      # Ensure enough unique values for Levene's test
      if (length(unique(all_finite_cond1)) >= 2 && length(unique(all_finite_cond2)) >= 2) {
        global_levene_p <- leveneTest(combined_all, group_all)$`Pr(>F)`[1]
      }
    }
  }
  
  cat("\n--- Global Statistical Tests per Condition ---\n")
  cat(paste0("Shapiro-Wilk p-value global (", cond1_pattern, "): ", round(shapiro_global_cond1_p, 4), "\n"))
  cat(paste0("Shapiro-Wilk p-value global (", cond2_pattern, "): ", round(shapiro_global_cond2_p, 4), "\n"))
  cat(paste0("Levene p-value global (", cond1_pattern, " vs ", cond2_pattern, "): ", round(global_levene_p, 4), "\n"))
  
  # Visualizations (boxplot and histogram)
  cv_data_boxplot <- data.frame(
    CV_Percent = c(variability_data[[paste0("CV_", cond1_pattern, "_Percent")]], variability_data[[paste0("CV_", cond2_pattern, "_Percent")]]),
    Condition = factor(rep(c(cond1_pattern, cond2_pattern), each = nrow(variability_data)))
  )
  cv_data_boxplot <- na.omit(cv_data_boxplot)
  
  if (nrow(cv_data_boxplot) > 0) {
    boxplot_cv <- ggplot(cv_data_boxplot, aes(x = Condition, y = CV_Percent)) +
      geom_boxplot() +
      labs(title = paste0("Coefficient of Variation Distribution (%) - ", work_name),
           x = "Condition",
           y = "CV (%)") +
      theme_cowplot(font_size = base_font_size, font_family = "Open Sans") +
      theme(plot.title = element_text(size = title_font_size, family = "Open Sans", hjust = 0.5),
            axis.title = element_text(size = axis_title_font_size, family = "Open Sans"),
            axis.text = element_text(size = axis_text_font_size, family = "Open Sans"),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            axis.ticks.length = unit(.15, "cm"),
            axis.text.x = element_text(margin = margin(t = 5)),
            axis.text.y = element_text(margin = margin(r = 5)))
    
    print(boxplot_cv)
    ggsave(file.path(output_folder_path, basename(boxplot_cv_file)), plot = boxplot_cv, width = 6, height = 4, dpi = dpi_resolution)
    cat(paste("CV distribution boxplot saved to:", file.path(output_folder_path, basename(boxplot_cv_file)), "\n"))
  } else {
    cat("No sufficient data to generate CV boxplot after omitting NA.\n")
  }
  
  log2fc_no_inf <- variability_data[[paste0("Log2FC_", cond2_pattern, "_vs_", cond1_pattern)]][is.finite(variability_data[[paste0("Log2FC_", cond2_pattern, "_vs_", cond1_pattern)]])]
  
  if (length(log2fc_no_inf) > 0) {
    hist_log2fc <- ggplot(data.frame(Log2FC = log2fc_no_inf), aes(x = Log2FC)) +
      geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = paste0("Log2 Fold Change Distribution (", cond2_pattern, "/", cond1_pattern, ") - ", work_name),
           x = paste0("Log2 Fold Change (", cond2_pattern, "/", cond1_pattern, ")"),
           y = "Frequency") +
      theme_cowplot(font_size = base_font_size, font_family = "Open Sans") +
      theme(plot.title = element_text(size = title_font_size, family = "Open Sans", hjust = 0.5),
            axis.title = element_text(size = axis_title_font_size, family = "Open Sans"),
            axis.text = element_text(size = axis_text_font_size, family = "Open Sans"),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            axis.ticks.length = unit(.15, "cm"),
            axis.text.x = element_text(margin = margin(t = 5)),
            axis.text.y = element_text(margin = margin(r = 5)))
    
    print(hist_log2fc)
    ggsave(file.path(output_folder_path, basename(hist_log2fc_file)), plot = hist_log2fc, width = 6, height = 4, dpi = dpi_resolution)
    cat(paste("Log2FC distribution histogram saved to:", file.path(output_folder_path, basename(hist_log2fc_file)), "\n"))
  } else {
    cat("No sufficient data to generate Log2FC histogram after omitting Inf/NA.\n")
  }
  
  # Save results to file
  cat("\n--- Statistical Test Results ---\n", file = output_file_variability, append = FALSE)
  cat("  - P-value > 0.05: No significant evidence to reject the null hypothesis that data are normally distributed.\n", file = output_file_variability, append = TRUE)
  cat("  - P-value <= 0.05: Significant evidence to reject the null hypothesis of normality (data may not be normally distributed).\n", file = output_file_variability, append = TRUE)
  cat("\nLevene's Test (Equality of Variances):\n", file = output_file_variability, append = TRUE)
  cat("  - P-value > 0.05: No significant evidence to reject the null hypothesis that variances between groups are equal.\n", file = output_file_variability, append = TRUE)
  cat("  - P-value <= 0.05: Significant evidence to reject the null hypothesis of equal variances (variances between groups may be different).\n", file = output_file_variability, append = TRUE)
  cat(paste0("\n--- Protein-wise Variability Results (", work_name, ") ---\n"), file = output_file_variability, append = TRUE)
  write.table(variability_data, file = output_file_variability, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  
  cat("\n--- Global Variability Measures ---\n", file = output_file_variability, append = TRUE)
  cat(paste0("Median CV (", cond1_pattern, "): ", round(median_cv_cond1, 2), "%\n"), file = output_file_variability, append = TRUE)
  cat(paste0("Median CV (", cond2_pattern, "): ", round(median_cv_cond2, 2), "%\n"), file = output_file_variability, append = TRUE)
  cat(paste0("Standard Deviation of |Log2FC (", cond2_pattern, " vs ", cond1_pattern, ")|: ", round(sd_abs_log2fc, 2), "\n"), file = output_file_variability, append = TRUE)
  
  cat("\n--- Global Statistical Tests per Condition ---\n", file = output_file_variability, append = TRUE)
  cat(paste0("Shapiro-Wilk p-value global (", cond1_pattern, "): ", round(shapiro_global_cond1_p, 4), "\n"), file = output_file_variability, append = TRUE)
  cat(paste0("Shapiro-Wilk p-value global (", cond2_pattern, "): ", round(shapiro_global_cond2_p, 4), "\n"), file = output_file_variability, append = TRUE)
  cat(paste0("Levene p-value global (", cond1_pattern, " vs ", cond2_pattern, "): ", round(global_levene_p, 4), "\n"), file = output_file_variability, append = TRUE)
  
  cat(paste0("Variability measures and statistical tests (protein-wise and global) saved to: ", output_file_variability, "\n"))
  
  # --- ADDITIONAL ANALYSES: HEATMAP AND NMDS/ANOSIM ---
  cat(paste0("\n--- STARTING HEATMAP and NMDS/ANOSIM (", work_name, ") ---\n"))
  
  # Prepare data for heatmap and NMDS
  # MODIFICACIÓN CLAVE AQUÍ: Aseguramos la selección y conversión a matriz de forma robusta
  # -------------------------------------------------------------------------------------
  
  # 1. Seleccionar *solo* las columnas de intensidad LFQ para el análisis de visualización.
  # Usamos `dplyr::select` con `all_of()` para una selección más segura y explícita.
  # El resultado será un tibble.
  lfq_data_selected <- protein_groups %>%
    dplyr::select(all_of(replicas_intensity_cols))
  
  # 2. Convertir el tibble seleccionado a una matriz numérica.
  # Esto es fundamental para asegurar que las operaciones posteriores (log2, min)
  # se realicen sobre un tipo de dato que R espera.
  lfq_data_numeric_matrix <- as.matrix(lfq_data_selected)
  
  # 3. Convertir ceros a NA en esta matriz (para que log2(0) no sea -Inf directamente)
  lfq_data_numeric_matrix[lfq_data_numeric_matrix == 0] <- NA
  
  # 4. Aplicar log2. Los NAs seguirán siendo NAs.
  lfq_data_log2 <- log2(lfq_data_numeric_matrix)
  
  # --- Imputation (temporary for visualization/analysis) ---
  # Find the minimum positive finite value across *all* relevant LFQ data
  # (usando los datos originales antes de la transformación log2 y la conversión a NA).
  #
  # Aquí volvemos a extraer los datos originales de LFQ intensity para calcular el epsilon,
  # pero asegurándonos de que sea un vector numérico plano usando unlist().
  all_original_lfq_values_vector <- unlist(protein_groups %>%
                                             dplyr::select(all_of(replicas_intensity_cols)))
  
  # Filtrar por valores positivos y finitos
  positive_finite_values <- all_original_lfq_values_vector[all_original_lfq_values_vector > 0 & is.finite(all_original_lfq_values_vector)]
  
  # Calcular el mínimo de esos valores
  min_val_global_for_imputation <- min(positive_finite_values, na.rm = TRUE)
  
  # If all original LFQ values are 0 or NA, then min_val_global_for_imputation will be Inf
  if (is.infinite(min_val_global_for_imputation) || is.na(min_val_global_for_imputation)) {
    cat(paste0("WARNING: All LFQ intensity values for ", work_name, " are zero or NA after filtering. Skipping Heatmap, NMDS, and ANOSIM.\n"))
    return(variability_data) # Retornar y continuar con el siguiente trabajo
  }
  
  # Calculate epsilon for imputation (a small value below the minimum)
  epsilon_for_imputation <- log2(min_val_global_for_imputation / 100) # log2 of a value much smaller than the smallest detectable
  
  # Impute NA and -Inf values in the log2-transformed data
  lfq_data_log2_imputed <- lfq_data_log2 # Start with the log2 data
  lfq_data_log2_imputed[is.na(lfq_data_log2_imputed)] <- epsilon_for_imputation
  lfq_data_log2_imputed[is.infinite(lfq_data_log2_imputed)] <- epsilon_for_imputation # This handles -Inf from log2(0)
  
  # Ahora, `lfq_data_log2_imputed` ya es una matriz con todos los valores numéricos finitos.
  # Esta es la matriz que se pasa a pheatmap y a vegan.
  lfq_data_matrix <- lfq_data_log2_imputed # Renombrado para claridad, ya es una matriz
  
  # -------------------------------------------------------------------------------------
  # --- HEATMAP ---
  cat(paste0("\nGenerating Heatmap for ", work_name, "...\n"))
  
  # Check for sufficient data for heatmap and variance (pheatmap needs variance to cluster)
  # Check if there's any variance after imputation. If all values are now identical, var() will be 0.
  if (nrow(lfq_data_matrix) > 1 && ncol(lfq_data_matrix) > 1 && var(as.vector(lfq_data_matrix), na.rm = TRUE) > 0) {
    tryCatch({
      # Create column annotations (conditions)
      sample_conditions <- gsub("(\\d+|\\d+\\.\\d+)", "", colnames(lfq_data_matrix)) # Remove numbers to get condition group
      sample_conditions_df <- data.frame(Condition = factor(sample_conditions))
      rownames(sample_conditions_df) <- colnames(lfq_data_matrix)
      
      pheatmap(lfq_data_matrix,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = FALSE, # Do not show protein names if too many
               main = paste0("LFQ Intensity (Log2) Heatmap - ", work_name),
               annotation_col = sample_conditions_df,
               filename = heatmap_file,
               width = 8, height = 10, # Adjust these as needed for resolution/number of proteins
               dpi = dpi_resolution,
               gaps_col = cumsum(table(sample_conditions_df$Condition)[unique(sample_conditions_df$Condition)])[-length(unique(sample_conditions_df$Condition))], # Add gaps between conditions
               # Ajustar fuentes dentro de pheatmap (se escalan con el base_font_size definido globalmente)
               main_gp = grid::gpar(fontfamily = "Open Sans", fontsize = title_font_size),
               annotation_legend_param = list(
                 labels_gp = grid::gpar(fontfamily = "Open Sans", fontsize = axis_text_font_size),
                 title_gp = grid::gpar(fontfamily = "Open Sans", fontsize = axis_title_font_size)
               ),
               labels_row = NULL # Explicitly set to NULL to ensure no row names
      )
      cat(paste("Heatmap saved to:", heatmap_file, "\n"))
    }, error = function(e) {
      cat(paste0("WARNING: Could not generate Heatmap for ", work_name, ". Error: ", e$message, "\n"))
      # No return(NULL) here; we continue to NMDS if possible
    })
  } else {
    cat(paste0("WARNING: Not enough clean data or variance to generate Heatmap for work ", work_name, " after imputation. Skipping Heatmap.\n"))
  }
  
  # --- NMDS and ANOSIM ---
  cat(paste0("\nGenerating NMDS and ANOSIM for ", work_name, "...\n"))
  
  # The distance matrix for NMDS must be without NA. Use imputed data.
  # Transpose the matrix so samples are rows and proteins are columns
  nmds_data <- t(lfq_data_log2_imputed) # Usamos la matriz ya procesada
  
  # Create the group vector for ANOSIM
  # Make sure group names are clean for the plot legend
  groups <- factor(gsub("(\\d+|\\d+\\.\\d+)", "", rownames(nmds_data)))
  
  # Check if there are enough samples and at least two unique groups with more than one sample per group
  if (nrow(nmds_data) > 1 && ncol(nmds_data) > 1 && length(unique(groups)) > 1 && min(table(groups)) >= 2) { # min 2 samples per group
    tryCatch({
      # Check if all columns (proteins) have zero variance, which would make distance calculation meaningless
      if (all(apply(nmds_data, 2, var, na.rm = TRUE) == 0)) {
        cat(paste0("WARNING: All protein intensities are identical or without variance across samples for NMDS in ", work_name, " after imputation. Skipping NMDS/ANOSIM.\n"))
        # No return(NULL) here; let the function finish.
      } else {
        # Calculate Euclidean distance (or other, Bray-Curtis is common in ecology)
        dist_matrix <- vegdist(nmds_data, method = "euclidean")
        
        # NMDS
        nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100, autotransform = FALSE, expand = FALSE) # Increased trymax
        
        # ANOSIM
        anosim_result <- anosim(dist_matrix, groups)
        
        # Extract ANOSIM R and p-value
        anosim_R_value <- round(anosim_result$statistic, 3) # Redondear R a 3 decimales
        anosim_p_value <- format.pval(anosim_result$signif, digits = 3, eps = 0.001) # Formatear p-value (e.g., < 0.001)
        
        # Crear cadena de texto para los resultados de ANOSIM
        anosim_text <- paste0("ANOSIM R = ", anosim_R_value, ", p = ", anosim_p_value)
        
        # Plot NMDS
        png(filename = nmds_file, width = 7, height = 6, units = "in", res = dpi_resolution, family = "Open Sans")
        plot(nmds_result, type = "n", main = paste0("NMDS of LFQ Intensity (Log2) - ", work_name),
             cex.main = title_font_size / base_font_size * 1.0, # Ajuste fino si es necesario
             cex.axis = axis_text_font_size / base_font_size * 1.0,
             cex.lab = axis_title_font_size / base_font_size * 1.0)
        
        # Plot points and spiders
        points(nmds_result, display = "sites", col = as.numeric(groups), pch = 19, cex = 3.0) # AUMENTADO el tamaño del punto
        ordispider(nmds_result, groups, col = as.numeric(groups), lwd = 3) # AUMENTADO el ancho de la línea
        
        # Add legend
        legend("topright", legend = levels(groups), col = 1:length(levels(groups)), pch = 19, cex = 2.0, bty = "n", # AUMENTADO cex de la leyenda
               text.font = 2)
        
        # Add labels for samples, adjust position to avoid overlap
        text(nmds_result, display = "sites", labels = rownames(nmds_data), cex = 1.8, pos = 3, font = 2) # AUMENTADO cex de las etiquetas de los puntos
        
        # AÑADIR RESULTADOS DE ANOSIM AL PLOT
        # Determinar los límites del plot para posicionamiento (usamos las coordenadas actuales del plot)
        plot_usr <- par("usr") # user coordinates: [x1, x2, y1, y2]
        
        # Posicionar en la esquina inferior izquierda, ligeramente desplazado del borde
        # Usamos los límites del plot directamente y un porcentaje de esos límites para el margen
        text_x_pos <- plot_usr[1] + (plot_usr[2] - plot_usr[1]) * 0.02 # 2% desde la izquierda
        text_y_pos <- plot_usr[3] + (plot_usr[4] - plot_usr[3]) * 0.02 # 2% desde abajo
        
        text(x = text_x_pos, y = text_y_pos,
             labels = anosim_text,
             adj = c(0, 0), # Alinear texto a la esquina inferior izquierda de la posición
             cex = axis_text_font_size / base_font_size * 1.5, # Ligeramente más grande que el texto del eje para destacarlo
             font = 2, # Fuente en negrita
             col = "black") # Color del texto
        
        dev.off()
        cat(paste("NMDS plot saved to:", nmds_file, "\n"))
        
        # Save ANOSIM results to a text file (this part is unchanged)
        cat("\n--- ANOSIM Results ---\n", file = anosim_file, append = FALSE)
        capture.output(summary(anosim_result), file = anosim_file, append = TRUE)
        cat(paste0("\nANOSIM results saved to:", anosim_file, "\n"))
      }
    }, error = function(e) {
      cat(paste0("WARNING: Could not perform NMDS/ANOSIM for ", work_name, ". Error: ", e$message, "\n"))
      # No return(NULL) here; let the function finish.
    })
  } else {
    cat(paste0("WARNING: Not enough samples (minimum 2 per group) or conditions for NMDS/ANOSIM for work ", work_name, ".\n"))
  }
  
  cat(paste0("\n--- ANALYSIS COMPLETED FOR WORK: ", work_name, " ---\n"))
  
  return(variability_data) # Always return variability_data if successfully created
}

# --- EXECUTION FOR EACH WORK ---
base_path_global <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/7_Especial_IJMS/analisis27062025"

# Define conditions for each work
# Each element in the list contains a vector with the patterns of the two conditions
# The order matters: the first pattern is considered the "denominator" and the second the "numerator"
# for Log2FC (Log2FC = Numerator / Denominator).
work_conditions <- list(
  W1 = c("RT", "Cold"), # (Cold/RT) -> Numerator: Cold, Denominator: RT
  W2 = c("u", "y"),     # (y/u)
  W3 = c("Control", "NaCl"), # (NaCl/Control)
  W4 = c("C", "BNF"),   # (BNF/C)
  W5 = c("C", "PM")     # (PM/C)
)

# Iterate over each work and execute the function
for (work in names(work_conditions)) {
  conditions <- work_conditions[[work]]
  
  # Llamada a la función principal y manejo de errores/retornos NULL
  # 'result' capturará lo que la función devuelva (variability_data o NULL)
  result <- tryCatch({
    process_and_analyze_proteingroups(
      work_name = work,
      base_path = base_path_global,
      condition_patterns = conditions
    )
  }, error = function(e) {
    # Si ocurre un error CATASTRÓFICO (ej. archivo no encontrado, datos de entrada vacíos,
    # o columnas críticas no halladas) dentro de la función process_and_analyze_proteingroups,
    # lo capturamos aquí, mostramos un mensaje y devolvemos NULL.
    cat(paste0("\n--- ERROR CATASTROFICO DURANTE EL PROCESAMIENTO DE LA OBRA: ", work, " ---\n"))
    cat(paste0("Mensaje de error: ", e$message, "\n"))
    return(NULL) # Devolvemos NULL para indicar que esta obra falló completamente
  })
  
  # Verificamos si la función devolvió NULL (indicando que no pudo completar el procesamiento)
  if (is.null(result)) {
    cat(paste0("Saltando el procesamiento adicional para la obra: ", work, " debido a problemas críticos anteriores.\n"))
    next # Saltamos a la siguiente iteración del bucle 'for'
  }
  
  # Si la función completó y devolvió 'variability_data', 'result' contendrá esos datos.
  # Aquí podrías añadir código para, por ejemplo, almacenar o combinar los resultados
  # de 'variability_data' de todas las obras que se procesaron correctamente.
  # Por ahora, simplemente el bucle continuará.
}

cat("\n--- PROCESAMIENTO DE TODAS LAS OBRAS COMPLETADO ---\n")
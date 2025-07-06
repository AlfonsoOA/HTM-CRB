# --- GLOBAL CONFIGURATION OF PATHS AND LIBRARIES ---
# Define the base path for all projects (W1, W2, etc.)
# Ensure this path exists and contains the W1, W2, etc. folders
base_analysis_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/7_Especial_IJMS/analisis27062025/"

# --- INSTALL AND LOAD LIBRARIES ---
# Ensure these libraries are installed: install.packages(c("dplyr", "purrr", "stringr"))
library(dplyr)
library(purrr)
library(stringr) # For string manipulation like str_replace

# --- AUXILIARY FUNCTION DEFINITIONS ---

# Helper function to get the protein ID column name based on the method
# This function must be globally accessible for all contexts where it's used.
get_id_col_name <- function(method) {
  if (method == "msstats") {
    return("Protein.ID")
  } else {
    return("Protein.IDs")
  }
}

# --- Function to separate UP/DOWN proteins for ClueGO and save them ---
# file_path_input: Full path to the input file (FC or Bayes results)
# output_cluego_folder: Path to the specific output folder for ClueGO (e.g., W1_ClueGO)
# work_number: The current work number (e.g., 1, 2, 3, etc.)
separate_up_down_and_save_cluego <- function(file_path_input, output_cluego_folder, work_number) {
  
  cat(paste0("\nProcessing file for UP/DOWN separation: ", basename(file_path_input), "\n"))
  
  # Read the results file with error handling
  res_df <- tryCatch({
    read.delim(file_path_input, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    cat(paste("Error reading file:", basename(file_path_input), ". Message:", e$message, ". Skipping UP/DOWN separation.\n"))
    return(NULL) # Return NULL if there's an error
  })
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    cat(paste("INFO: File", basename(file_path_input), "is empty or could not be read. No proteins to separate UP/DOWN.\n"))
    return(NULL)
  }
  
  # --- Identify Protein ID Column Flexibly ---
  protein_id_col_name <- NULL
  if ("Protein.IDs" %in% colnames(res_df)) {
    protein_id_col_name <- "Protein.IDs"
  } else if ("Protein.ID" %in% colnames(res_df)) { # For msstats
    protein_id_col_name <- "Protein.ID"
  }
  
  # --- Identify Fold Change Column Flexibly ---
  # Priority: log2FC_mean_bayes (for our generated Bayes files), then log2FC, then logFC, then log2FC_mean (for original WX_res_bayes.txt)
  fc_col_name <- NULL
  if ("log2FC_mean_bayes" %in% colnames(res_df)) {
    fc_col_name <- "log2FC_mean_bayes"
  } else if ("log2FC" %in% colnames(res_df)) {
    fc_col_name <- "log2FC"
  } else if ("logFC" %in% colnames(res_df)) {
    fc_col_name <- "logFC"
  } else if ("log2FC_mean" %in% colnames(res_df)) { # For WX_res_bayes.txt or WX_resbiolrel_bayes.txt
    fc_col_name <- "log2FC_mean"
  }
  
  if (is.null(fc_col_name) || is.null(protein_id_col_name)) {
    cat(paste("WARNING: File", basename(file_path_input), "does not contain a valid Fold Change column ('log2FC', 'logFC', 'log2FC_mean', 'log2FC_mean_bayes') or a valid Protein ID column ('Protein.IDs', 'Protein.ID'). Skipping.\n"))
    return(NULL)
  }
  
  # --- Determine Base Name for Output Files ---
  # Remove prefixes and extensions to get a clean base name for ClueGO
  base_name <- basename(file_path_input)
  
  # Example: W1_biolrel_tstudent_FC.txt -> W1_tstudent_FC
  base_name <- str_replace(base_name, paste0("^W", work_number, "_biolrel_"), paste0("W", work_number, "_")) 
  # Example: W1_res_bayes.txt -> W1_bayes
  base_name <- str_replace(base_name, paste0("^W", work_number, "_res_"), paste0("W", work_number, "_"))
  # Example: W1_resbiolrel_bayes.txt -> W1_bayes (renaming for consistency with original Bayes)
  base_name <- str_replace(base_name, paste0("^W", work_number, "_resbiolrel_"), paste0("W", work_number, "_")) 
  
  base_name <- str_replace(base_name, "\\.txt$", "") # Remove .txt extension
  
  # Filter and save upregulated (UP) proteins
  # Use !!sym() to dynamically reference the column
  up_data <- res_df %>%
    filter(!!sym(fc_col_name) > 0 | (is.infinite(!!sym(fc_col_name)) & !!sym(fc_col_name) > 0)) %>% 
    select(all_of(c(protein_id_col_name, fc_col_name))) %>% # *** MODIFIED: Select both ID and FC column ***
    distinct() # Ensure unique IDs (based on both columns now, though Protein ID is primary)
  
  # Rename the Protein ID column to "Protein.IDs" for ClueGO consistency
  # And rename the FC column to "log2FC" for consistency in output
  if (protein_id_col_name != "Protein.IDs" || fc_col_name != "log2FC") {
    up_data <- up_data %>% rename(Protein.IDs = !!sym(protein_id_col_name), log2FC = !!sym(fc_col_name))
  }
  
  output_up_name <- file.path(output_cluego_folder, paste0(base_name, "_ClueGO_up.txt"))
  # *** MODIFIED: write.table with col.names = TRUE and quote = FALSE ***
  write.table(up_data, output_up_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) 
  cat(paste("File created:", basename(output_up_name), "with", nrow(up_data), "upregulated proteins (including FC).\n"))
  
  # Filter and save downregulated (DOWN) proteins
  down_data <- res_df %>%
    filter(!!sym(fc_col_name) < 0 | (is.infinite(!!sym(fc_col_name)) & !!sym(fc_col_name) < 0)) %>%
    select(all_of(c(protein_id_col_name, fc_col_name))) %>% # *** MODIFIED: Select both ID and FC column ***
    distinct() # Ensure unique IDs
  
  # Rename the Protein ID column to "Protein.IDs" for consistency
  # And rename the FC column to "log2FC" for consistency in output
  if (protein_id_col_name != "Protein.IDs" || fc_col_name != "log2FC") {
    down_data <- down_data %>% rename(Protein.IDs = !!sym(protein_id_col_name), log2FC = !!sym(fc_col_name))
  }
  
  output_down_name <- file.path(output_cluego_folder, paste0(base_name, "_ClueGO_down.txt"))
  # *** MODIFIED: write.table with col.names = TRUE and quote = FALSE ***
  write.table(down_data, output_down_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) 
  cat(paste("File created:", basename(output_down_name), "with", nrow(down_data), "downregulated proteins (including FC).\n"))
}

# --- Function to generate and save intersection tables (global) ---
# proteins_list: List of protein vectors to compare (each element is a set of IDs)
# regulation_type: "UP" or "DOWN"
# bio_criterion: "FC" or "Bayes"
# output_base_path: Folder where the table will be saved (e.g., W1_ClueGO)
generate_intersection_table <- function(proteins_list, regulation_type, bio_criterion, output_base_path) {
  # Ensure the protein list is not empty and contains at least 2 non-empty sets
  proteins_list_filtered <- purrr::keep(proteins_list, ~ length(.) > 0)
  num_sets <- length(proteins_list_filtered)
  
  if (num_sets > 1) {
    file_name <- file.path(output_base_path,
                           paste0("Intersection_Table_Proteins_", regulation_type, "_", bio_criterion, ".txt"))
    
    cat(paste0("Calculating intersections for ", regulation_type, " proteins (", bio_criterion, ")...\n"))
    
    # Create a matrix to store intersection results
    set_names <- names(proteins_list_filtered)
    intersection_matrix <- matrix(NA, nrow = num_sets, ncol = num_sets,
                                  dimnames = list(set_names, set_names))
    
    # Populate the matrix
    for (i in 1:num_sets) {
      for (j in 1:num_sets) { # Iterate through the entire matrix
        if (i == j) { # Diagonal: set size
          intersection_matrix[i, j] <- length(proteins_list_filtered[[i]])
        } else { # Pairwise intersections
          intersection_matrix[i, j] <- length(intersect(proteins_list_filtered[[i]], proteins_list_filtered[[j]]))
        }
      }
    }
    
    # Convert to dataframe for better visualization and writing
    intersections_df <- as.data.frame(intersection_matrix)
    intersections_df <- cbind(Method = rownames(intersections_df), intersections_df)
    
    # Save the table to a text file
    write.table(intersections_df, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(paste0("Intersection table for ", regulation_type, " proteins (", bio_criterion, ") saved to: ", file_name, "\n"))
    
  } else {
    cat(paste0("Not enough '", regulation_type, "' datasets to generate a meaningful intersection table (", bio_criterion, ", ", num_sets, " found).\n"))
  }
}

# --- MAIN FUNCTION TO PROCESS CLUEGO FILES AND INTERSECTION TABLES ---

process_cluego_analysis <- function(work_number, base_path) {
  
  cat(paste0("\n--- STARTING CLUEGO FILE PREPARATION AND INTERSECTION TABLES FOR WORK_W", work_number, " ---\n"))
  
  # --- 1. Define Paths for the Current Work ---
  current_main_folder_path <- file.path(base_path, paste0("W", work_number))
  current_biol_rel_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_BioRel"))
  current_h_testing_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_H_testing")) # Needed for WX_res_bayes.txt if not copied
  current_cluego_folder_path <- file.path(current_main_folder_path, paste0("W", work_number, "_ClueGO")) # Output folder
  
  # Create the ClueGO folder if it doesn't exist
  if (!dir.exists(current_cluego_folder_path)) {
    dir.create(current_cluego_folder_path, recursive = TRUE)
    cat(paste("ClueGO files folder created at:", current_cluego_folder_path, "\n"))
  } else {
    cat(paste("ClueGO files folder already exists at:", current_cluego_folder_path, "\n"))
  }
  
  # --- 2. Prepare List of Input Files ---
  
  # Files from WX_BioRel folder (tstudent, twelch, limma, deqms, msstats for FC and Bayes)
  biorel_files <- list.files(current_biol_rel_folder_path,
                             pattern = paste0("^W", work_number, "_biolrel_(tstudent|twelch|limma|deqms|msstats)_(FC|Bayes)\\.txt$"),
                             full.names = TRUE)
  
  # The newly copied WX_resbiolrel_bayes.txt file in WX_BioRel
  copied_res_bayes_biorel <- file.path(current_biol_rel_folder_path, paste0("W", work_number, "_resbiolrel_bayes.txt"))
  
  # The original WX_res_bayes.txt file from WX_H_testing (in case it's still needed)
  original_res_bayes_htesting <- file.path(current_h_testing_folder_path, paste0("W", work_number, "_res_bayes.txt"))
  
  # Combine all paths of files to process for UP/DOWN (only existing ones)
  files_for_cluego_separation <- character(0)
  
  if (length(biorel_files) > 0) {
    files_for_cluego_separation <- c(files_for_cluego_separation, biorel_files)
  }
  # Add the new copied Bayes file from BioRel
  if (file.exists(copied_res_bayes_biorel)) {
    files_for_cluego_separation <- c(files_for_cluego_separation, copied_res_bayes_biorel)
  }
  # Add the original Bayes from H_testing (if it exists and is not already covered by the copied version)
  # It's crucial to check if the content of copied_res_bayes_biorel is truly a copy of original_res_bayes_htesting
  # For simplicity, we'll include both if they exist, and the unique() call below will handle actual duplicates.
  if (file.exists(original_res_bayes_htesting)) {
    files_for_cluego_separation <- c(files_for_cluego_separation, original_res_bayes_htesting)
  }
  
  # Ensure no duplicates if, for some reason, a file has been processed already or is redundant
  files_for_cluego_separation <- unique(files_for_cluego_separation)
  
  
  # --- 3. Separate UP and DOWN proteins for ClueGO ---
  cat(paste("\n### Separating UP/DOWN proteins for ClueGO (Work W", work_number, ") ###\n"))
  
  if (length(files_for_cluego_separation) == 0) {
    cat(paste("WARNING: No relevant files found for UP/DOWN separation for ClueGO in Work W", work_number, ".\n"))
  } else {
    for (file_path_cluego in files_for_cluego_separation) {
      if (file.exists(file_path_cluego)) { 
        separate_up_down_and_save_cluego(file_path_cluego, current_cluego_folder_path, work_number)
      } else {
        cat(paste0("WARNING: File expected at ", file_path_cluego, " does not exist. Skipping.\n"))
      }
    }
  }
  cat(paste("\n### END: UP/DOWN protein separation for ClueGO completed (Work W", work_number, ") ###\n"))
  
  # --- 4. Generate Intersection Tables for UP/DOWN proteins ---
  cat(paste("\n### Generating Intersection Tables for UP/DOWN proteins (Work W", work_number, ") ###\n"))
  
  # List all _ClueGO_up.txt and _ClueGO_down.txt files that were just generated
  cluego_up_files <- list.files(current_cluego_folder_path, pattern = "_ClueGO_up\\.txt$", full.names = TRUE)
  cluego_down_files <- list.files(current_cluego_folder_path, pattern = "_ClueGO_down\\.txt$", full.names = TRUE)
  
  # --- Prepare data and generate Intersection Tables for UP proteins ---
  if (length(cluego_up_files) > 0) {
    # Files for the "FC" table: includes those ending in FC_ClueGO_up.txt AND _bayes_ClueGO_up.txt (which now also covers _resbiolrel_bayes.txt output)
    up_fc_cluego_files <- cluego_up_files[grepl("_FC_ClueGO_up\\.txt$|_bayes_ClueGO_up\\.txt$", cluego_up_files)]
    # Files for the "Bayes" table: includes those ending in Bayes_ClueGO_up.txt (filtered) AND _bayes_ClueGO_up.txt (original/copied)
    up_bayes_cluego_files <- cluego_up_files[grepl("_Bayes_ClueGO_up\\.txt$|_bayes_ClueGO_up\\.txt$", cluego_up_files)]
    
    # Process for UP (FC)
    if (length(up_fc_cluego_files) > 0) {
      list_up_fc_proteins <- purrr::map(up_fc_cluego_files, function(file_path) {
        df_temp <- tryCatch({
          read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character") # *** MODIFIED: header = TRUE ***
        }, error = function(e) {
          data.frame(Protein.IDs = character(0), log2FC = numeric(0), stringsAsFactors = FALSE) # *** MODIFIED: include log2FC column ***
        })
        if (nrow(df_temp) == 0) character(0) else df_temp$Protein.IDs %>% unique()
      })
      names(list_up_fc_proteins) <- basename(up_fc_cluego_files) %>% 
        str_replace(paste0("^W", work_number, "_"), "") %>% 
        str_replace("_ClueGO_up\\.txt$", "")
      generate_intersection_table(list_up_fc_proteins, "UP", "FC", current_cluego_folder_path)
    } else {
      cat("No UP files with FC criterion found to generate the intersection table.\n")
    }
    
    # Process for UP (Bayes)
    if (length(up_bayes_cluego_files) > 0) {
      list_up_bayes_proteins <- purrr::map(up_bayes_cluego_files, function(file_path) {
        df_temp <- tryCatch({
          read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character") # *** MODIFIED: header = TRUE ***
        }, error = function(e) {
          data.frame(Protein.IDs = character(0), log2FC = numeric(0), stringsAsFactors = FALSE) # *** MODIFIED: include log2FC column ***
        })
        if (nrow(df_temp) == 0) character(0) else df_temp$Protein.IDs %>% unique()
      })
      names(list_up_bayes_proteins) <- basename(up_bayes_cluego_files) %>% 
        str_replace(paste0("^W", work_number, "_"), "") %>% 
        str_replace("_ClueGO_up\\.txt$", "")
      generate_intersection_table(list_up_bayes_proteins, "UP", "Bayes", current_cluego_folder_path)
    } else {
      cat("No UP files with Bayes criterion found to generate the intersection table.\n")
    }
  } else {
    cat("No UP protein files found to generate any intersection table.\n")
  }
  
  # --- Prepare data and generate Intersection Tables for DOWN proteins ---
  if (length(cluego_down_files) > 0) {
    # Files for the "FC" table: includes those ending in FC_ClueGO_down.txt AND _bayes_ClueGO_down.txt
    down_fc_cluego_files <- cluego_down_files[grepl("_FC_ClueGO_down\\.txt$|_bayes_ClueGO_down\\.txt$", cluego_down_files)]
    # Files for the "Bayes" table: includes those ending in Bayes_ClueGO_down.txt AND _bayes_ClueGO_down.txt
    down_bayes_cluego_files <- cluego_down_files[grepl("_Bayes_ClueGO_down\\.txt$|_bayes_ClueGO_down\\.txt$", cluego_down_files)]
    
    # Process for DOWN (FC)
    if (length(down_fc_cluego_files) > 0) {
      list_down_fc_proteins <- purrr::map(down_fc_cluego_files, function(file_path) {
        df_temp <- tryCatch({
          read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character") # *** MODIFIED: header = TRUE ***
        }, error = function(e) {
          data.frame(Protein.IDs = character(0), log2FC = numeric(0), stringsAsFactors = FALSE) # *** MODIFIED: include log2FC column ***
        })
        if (nrow(df_temp) == 0) character(0) else df_temp$Protein.IDs %>% unique()
      })
      names(list_down_fc_proteins) <- basename(down_fc_cluego_files) %>% 
        str_replace(paste0("^W", work_number, "_"), "") %>% 
        str_replace("_ClueGO_down\\.txt$", "")
      generate_intersection_table(list_down_fc_proteins, "DOWN", "FC", current_cluego_folder_path)
    } else {
      cat("No DOWN files with FC criterion found to generate the intersection table.\n")
    }
    
    # Process for DOWN (Bayes)
    if (length(down_bayes_cluego_files) > 0) {
      list_down_bayes_proteins <- purrr::map(down_bayes_cluego_files, function(file_path) {
        df_temp <- tryCatch({
          read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character") # *** MODIFIED: header = TRUE ***
        }, error = function(e) {
          data.frame(Protein.IDs = character(0), log2FC = numeric(0), stringsAsFactors = FALSE) # *** MODIFIED: include log2FC column ***
        })
        if (nrow(df_temp) == 0) character(0) else df_temp$Protein.IDs %>% unique()
      })
      names(list_down_bayes_proteins) <- basename(down_bayes_cluego_files) %>% 
        str_replace(paste0("^W", work_number, "_"), "") %>% 
        str_replace("_ClueGO_down\\.txt$", "")
      generate_intersection_table(list_down_bayes_proteins, "DOWN", "Bayes", current_cluego_folder_path)
    } else {
      cat("No DOWN files with Bayes criterion found to generate the intersection table.\n")
    }
  } else {
    cat("No DOWN protein files found to generate any intersection table.\n")
  }
  
  cat(paste("\n### END: Intersection Tables generation for UP/DOWN proteins completed (Work W", work_number, ") ###\n"))
  
  cat(paste0("\n--- CLUEGO FILE PREPARATION AND INTERSECTION TABLES PROCESS COMPLETED FOR WORK_W", work_number, " ---\n"))
} # End of process_cluego_analysis function

# --- EXECUTE THE FUNCTION FOR EACH DATASET ---

cat("\n--- STARTING CLUEGO FILE PREPARATION AND INTERSECTION TABLES PROCESS FOR ALL DATASETS ---\n")

# List of work numbers to process (you can adjust this)
work_numbers_to_process <- c(1, 2, 3, 4, 5)

for (work_num in work_numbers_to_process) {
  process_cluego_analysis(
    work_number = work_num,
    base_path = base_analysis_path
  )
}

cat("\n--- CLUEGO FILE PREPARATION AND INTERSECTION TABLES PROCESS COMPLETED FOR ALL DATASETS ---\n")
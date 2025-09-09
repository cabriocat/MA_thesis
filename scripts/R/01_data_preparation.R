# Data Preparation and Model Setup Script
#
# This script handles the initial setup of the R analysis pipeline:
# - Loading required libraries
# - Reading preprocessed EEG data from Python pipeline
# - Setting up experimental parameters
# - Creating cumulative datasets for different repetition ranges
#
# Author: [Your Name]

# Load required libraries
load_required_libraries <- function() {
  cat("Loading required libraries...\n")
  
  # Core libraries
  suppressPackageStartupMessages({
    library(glue)
    library(lme4)
    library(lmerTest)
    library(RColorBrewer)
    library(ggsignif)
    library(emmeans)
    library(broom.mixed)
    library(tidyverse)
  })
  
  # Try to load svglite for SVG export
  if (!require("svglite", quietly = TRUE)) {
    cat("Attempting to install svglite package...\n")
    tryCatch({
      install.packages("svglite", dependencies = TRUE)
      library(svglite)
      cat("✓ svglite installed successfully\n")
    }, error = function(e) {
      cat("⚠️ Warning: svglite package could not be installed.\n")
      cat("SVG export may not work. Consider installing manually: install.packages('svglite')\n")
    })
  } else {
    cat("✓ svglite loaded\n")
  }
  
  cat("✓ All libraries loaded successfully\n\n")
}

# Set up experimental parameters
setup_parameters <- function() {
  cat("Setting up experimental parameters...\n")
  
  # Define electrode ROI
  electrode <- c("FC1", "FCz", "FC2", "FCC1h", "FCC2h", 
                 "C1", "Cz", "C2", "CCP1h", "CCP2h", 
                 "CP1", "CPz", "CP2", "CPP1h", "CPP2h")
  
  # Define conditions
  conditions <- c("animal", "food", "tool", "commun", "emotion", "social")
  
  # Define plot colors (Okabe-Ito palette)
  plot_colors <- c(
    "animal" = "#388E3C",
    "commun" = "#d9a00f", 
    "emotion" = "#616161",
    "food" = "#1976D2",
    "social" = "#ff00ff",
    "tool" = "#D32F2F"
  )
  
  cat("✓ Parameters set up successfully\n")
  cat(sprintf("  - ROI electrodes: %d channels\n", length(electrode)))
  cat(sprintf("  - Conditions: %s\n", paste(conditions, collapse = ", ")))
  cat("\n")
  
  return(list(
    electrode = electrode,
    conditions = conditions,
    plot_colors = plot_colors
  ))
}

# Load preprocessed data from Python pipeline
load_data <- function(data_path = NULL) {
  cat("Loading preprocessed EEG data...\n")
  
  # Default data path if not provided
  if (is.null(data_path)) {
    # Assume we're in the scripts/R directory, go up to project root
    project_root <- file.path("..", "..")
    data_path <- file.path(project_root, "data", "2_preprocessed", "285-345ms", 
                          "combined_task-nouns_285-345ms.csv")
  }
  
  # Check if file exists
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s\n", data_path),
         "Please run the Python preprocessing pipeline first or check the path.")
  }
  
  # Load data
  tryCatch({
    df <- read_csv(data_path, show_col_types = FALSE)
    cat("✓ Data loaded successfully\n")
    cat(sprintf("  - Data shape: %d rows × %d columns\n", nrow(df), ncol(df)))
    cat(sprintf("  - Subjects: %d\n", length(unique(df$subject))))
    cat(sprintf("  - Repetitions: %s\n", paste(sort(unique(df$repetition)), collapse = ", ")))
    cat(sprintf("  - Categories: %s\n", paste(sort(unique(df$category)), collapse = ", ")))
    cat("\n")
    
    return(df)
  }, error = function(e) {
    stop(sprintf("Error loading data: %s\n", e$message))
  })
}

# Create cumulative datasets for different repetition ranges
create_cumulative_datasets <- function(df, electrode, conditions) {
  cat("Creating cumulative repetition datasets...\n")
  
  datasets <- list()
  
  # Dataset 1: Only repetition 1
  cat("  - Creating dataset for repetition 1...\n")
  datasets$df_1 <- df %>%
    filter(category %in% conditions,
           repetition == 1,
           channel %in% electrode) %>%
    group_by(item, subject, category, repetition) %>%
    summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
    mutate(
      category = factor(category, levels = conditions)
    ) %>%
    { contrasts(.$category) <- contr.sum(length(conditions)); . }
  
  # Datasets 2-6: Cumulative repetitions
  for (max_rep in 2:6) {
    dataset_name <- sprintf("df_1_%d", max_rep)
    cat(sprintf("  - Creating dataset for repetitions 1-%d...\n", max_rep))
    
    datasets[[dataset_name]] <- df %>%
      filter(category %in% conditions,
             repetition %in% 1:max_rep,
             channel %in% electrode) %>%
      group_by(item, subject, category, repetition) %>%
      summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
      mutate(
        repetition_c = scale(repetition, scale = FALSE)[, 1],
        category = structure(
          factor(category, levels = conditions),
          contrasts = contr.sum(length(conditions))
        )
      )
  }
  
  cat("✓ All cumulative datasets created successfully\n")
  cat(sprintf("  - Total datasets: %d\n", length(datasets)))
  
  # Print dataset summaries
  for (name in names(datasets)) {
    cat(sprintf("  - %s: %d rows\n", name, nrow(datasets[[name]])))
  }
  cat("\n")
  
  return(datasets)
}

# Validate data structure
validate_data <- function(df, expected_cols = c("subject", "item", "repetition", "category", "channel", "voltage")) {
  cat("Validating data structure...\n")
  
  # Check required columns
  missing_cols <- setdiff(expected_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for missing values
  na_counts <- df %>% summarise(across(everything(), ~ sum(is.na(.))))
  na_cols <- names(na_counts)[na_counts > 0]
  if (length(na_cols) > 0) {
    cat("⚠️ Warning: Missing values found in columns:\n")
    for (col in na_cols) {
      cat(sprintf("  - %s: %d missing values\n", col, na_counts[[col]]))
    }
  }
  
  # Check data ranges
  cat("Data validation summary:\n")
  cat(sprintf("  - Voltage range: %.3f to %.3f μV\n", min(df$voltage, na.rm = TRUE), max(df$voltage, na.rm = TRUE)))
  cat(sprintf("  - Repetitions: %s\n", paste(sort(unique(df$repetition)), collapse = ", ")))
  cat(sprintf("  - Categories: %s\n", paste(sort(unique(df$category)), collapse = ", ")))
  
  cat("✓ Data validation completed\n\n")
}

# Main setup function
setup_data_and_models <- function(data_path = NULL) {
  cat("=========================================\n")
  cat("    R DATA SETUP AND MODEL PREPARATION   \n")
  cat("=========================================\n\n")
  
  # Load libraries
  load_required_libraries()
  
  # Setup parameters
  params <- setup_parameters()
  
  # Load data
  df <- load_data(data_path)
  
  # Validate data
  validate_data(df)
  
  # Create cumulative datasets
  datasets <- create_cumulative_datasets(df, params$electrode, params$conditions)
  
  cat("=========================================\n")
  cat("           SETUP COMPLETE!               \n")
  cat("=========================================\n\n")
  
  cat("Available objects:\n")
  cat("  - df: Main dataframe\n")
  cat("  - datasets: List of cumulative datasets (df_1, df_1_2, ..., df_1_6)\n")
  cat("  - params: Experimental parameters (electrode, conditions, plot_colors)\n\n")
  
  return(list(
    df = df,
    datasets = datasets,
    params = params
  ))
}

# If running this script directly
if (sys.nframe() == 0) {
  cat("Running data setup script...\n\n")
  setup_results <- setup_data_and_models()
  
  # Make objects available in global environment
  df <<- setup_results$df
  datasets <<- setup_results$datasets
  params <<- setup_results$params
  
  cat("Setup complete! Objects are now available in your environment.\n")
}

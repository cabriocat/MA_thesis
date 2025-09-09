#!/usr/bin/env Rscript
#
# Data Preparation and Setup
# ==========================
# 
# This script handles the initial data loading, preparation, and setup
# for statistical analysis of ERP data.
#
# Tasks performed:
# - Load the combined CSV data from Python preprocessing
# - Set up experimental parameters and condition mappings
# - Create cumulative datasets (1, 1-2, 1-3, ..., 1-6 repetitions)
# - Apply proper contrast coding for statistical analysis

# Load required libraries
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

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR DATA
# ============================================================================

# Path to the combined CSV file created by Python script 03_data_preparation.py
# This should be in your data/2_preprocessed/285-345ms/ directory
DATA_CSV_PATH <- "/Users/johannberger/Documents/thesis/data/2_preprocessed/285-345ms/combined_task-nouns_285-345ms.csv"

# ============================================================================

# Install and load svglite for SVG export (with error handling)
if (!require("svglite", quietly = TRUE)) {
  tryCatch({
    install.packages("svglite", dependencies = TRUE)
    library(svglite)
    cat("✓ Successfully installed and loaded svglite\n")
  }, error = function(e) {
    cat("⚠ Warning: svglite package could not be installed.\n")
    cat("SVG export may not work. Consider installing manually:\n")
    cat("install.packages('svglite')\n")
  })
}

#' Setup experimental parameters
#'
#' @return List containing experimental parameters
setup_parameters <- function() {
  
  # Define electrode ROI (frontocentral region)
  electrode <- c(
    "FC1", "FCz", "FC2", "FCC1h", "FCC2h", 
    "C1", "Cz", "C2", "CCP1h", "CCP2h", 
    "CP1", "CPz", "CP2", "CPP1h", "CPP2h"
  )
  
  # Define semantic categories
  conditions <- c("animal", "food", "tool", "commun", "emotion", "social")
  
  # Color palette (matching Python Okabe-Ito colors)
  plot_colors <- c(
    "animal"  = "#388E3C",
    "commun"  = "#d9a00f", 
    "emotion" = "#616161",
    "food"    = "#1976D2",
    "social"  = "#ff00ff",
    "tool"    = "#D32F2F"
  )
  
  return(list(
    electrode = electrode,
    conditions = conditions,
    plot_colors = plot_colors
  ))
}

#' Load and validate the preprocessed data
#'
#' @param data_path Path to the combined CSV file
#' @return Loaded and validated dataframe
load_data <- function(data_path) {
  
  cat("Loading data from:", data_path, "\n")
  
  if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path)
  }
  
  df <- read_csv(data_path, show_col_types = FALSE)
  
  # Validate data structure
  expected_cols <- c("subject", "item", "repetition", "category", "channel", "voltage")
  missing_cols <- setdiff(expected_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Basic data summary
  cat("Data loaded successfully:\n")
  cat("  Shape:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("  Subjects:", length(unique(df$subject)), "\n")
  cat("  Categories:", paste(sort(unique(df$category)), collapse = ", "), "\n")
  cat("  Repetitions:", paste(sort(unique(df$repetition)), collapse = ", "), "\n")
  cat("  Channels:", length(unique(df$channel)), "\n")
  
  return(df)
}

#' Create cumulative datasets for different repetition counts
#'
#' @param df Raw data dataframe
#' @param electrode ROI electrode list
#' @param conditions Condition list
#' @return List of cumulative dataframes
create_cumulative_datasets <- function(df, electrode, conditions) {
  
  cat("Creating cumulative datasets...\n")
  
  # Base filtering function
  create_dataset <- function(rep_range, dataset_name) {
    
    cat("  Creating", dataset_name, "dataset (reps", paste(rep_range, collapse = "-"), ")\n")
    
    df_filtered <- df %>%
      filter(
        category %in% conditions,
        repetition %in% rep_range,
        channel %in% electrode
      ) %>%
      group_by(item, subject, category, repetition) %>%
      summarise(mean_voltage = mean(voltage), .groups = "drop")
    
    # Set up contrasts for statistical analysis
    if (length(rep_range) == 1) {
      # For single repetition, only need category factor
      df_filtered <- df_filtered %>%
        mutate(
          category = factor(category, levels = conditions)
        )
      contrasts(df_filtered$category) <- contr.sum(length(conditions))
      
    } else {
      # For multiple repetitions, add centered repetition variable
      df_filtered <- df_filtered %>%
        mutate(
          repetition_c = scale(repetition, scale = FALSE)[, 1],
          category = structure(
            factor(category, levels = conditions),
            contrasts = contr.sum(length(conditions))
          )
        )
    }
    
    cat("    Final shape:", nrow(df_filtered), "observations\n")
    return(df_filtered)
  }
  
  # Create datasets for each cumulative repetition count
  datasets <- list(
    df_1   = create_dataset(1,     "df_1"),
    df_1_2 = create_dataset(1:2,   "df_1_2"), 
    df_1_3 = create_dataset(1:3,   "df_1_3"),
    df_1_4 = create_dataset(1:4,   "df_1_4"),
    df_1_5 = create_dataset(1:5,   "df_1_5"),
    df_1_6 = create_dataset(1:6,   "df_1_6")
  )
  
  cat("✓ Created", length(datasets), "cumulative datasets\n")
  return(datasets)
}

#' Validate datasets and report statistics
#'
#' @param datasets List of cumulative datasets
#' @param conditions List of conditions
validate_datasets <- function(datasets, conditions) {
  
  cat("\nDataset validation:\n")
  cat(sprintf("%-10s %-15s %-15s %-20s\n", "Dataset", "N_obs", "N_subjects", "Categories"))
  cat(paste(rep("-", 65), collapse = ""), "\n")
  
  for (name in names(datasets)) {
    df <- datasets[[name]]
    n_obs <- nrow(df)
    n_subjects <- length(unique(df$subject))
    categories <- paste(sort(unique(df$category)), collapse = ",")
    
    cat(sprintf("%-10s %-15d %-15d %-20s\n", name, n_obs, n_subjects, categories))
    
    # Check for missing categories
    missing_cats <- setdiff(conditions, unique(df$category))
    if (length(missing_cats) > 0) {
      cat("  ⚠ Warning: Missing categories:", paste(missing_cats, collapse = ", "), "\n")
    }
  }
  
  cat("\n✓ Dataset validation complete\n")
}

#' Main execution function
main <- function() {
  
  cat("==============================================================\n")
  cat("DATA PREPARATION AND SETUP\n")
  cat("==============================================================\n")
  
  # Setup parameters
  params <- setup_parameters()
  cat("Experimental setup:\n")
  cat("  ROI electrodes:", length(params$electrode), "channels\n")
  cat("  Conditions:", paste(params$conditions, collapse = ", "), "\n")
  
  # Load data using configured path
  df <- load_data(DATA_CSV_PATH)
  
  # Create cumulative datasets
  datasets <- create_cumulative_datasets(df, params$electrode, params$conditions)
  
  # Validate datasets
  validate_datasets(datasets, params$conditions)
  
  cat("\n==============================================================\n")
  cat("DATA PREPARATION COMPLETE\n") 
  cat("==============================================================\n")
  cat("Ready for statistical analysis\n")
  cat("Available datasets: ", paste(names(datasets), collapse = ", "), "\n")
  
  # Return everything for use in subsequent scripts
  return(list(
    params = params,
    raw_data = df,
    datasets = datasets
  ))
}

# Execute main function if script is run directly
if (!interactive()) {
  results <- main()
}

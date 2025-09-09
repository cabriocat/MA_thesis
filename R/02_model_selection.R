#!/usr/bin/env Rscript
#  
# Model Selection
# ===============
#
# This script performs systematic model selection for each cumulative dataset
# to determine the best-fitting statistical model structure.
#
# Tasks performed:
# - Fit multiple candidate models for each dataset
# - Compare models using likelihood ratio tests
# - Select best model based on AIC/BIC and significance tests
# - Document model selection rationale

# Load required libraries
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(tidyverse)
})

# Source data preparation script
if (file.exists("01_data_preparation.R")) {
  source("01_data_preparation.R")
} else {
  cat("Warning: 01_data_preparation.R not found. Run data preparation first.\n")
}

#' Fit candidate models for a dataset
#'
#' @param df Dataset to fit models to
#' @param dataset_name Name of the dataset for reporting
#' @return List of fitted models
fit_candidate_models <- function(df, dataset_name) {
  
  cat("Fitting candidate models for", dataset_name, "...\n")
  
  models <- list()
  
  # Model 1: Category only (baseline)
  tryCatch({
    models$mdl1 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item), data = df)
    cat("  ✓ Model 1: Category only\n")
  }, error = function(e) {
    cat("  ✗ Model 1 failed:", e$message, "\n")
    models$mdl1 <<- NULL
  })
  
  # For datasets with repetition_c variable, fit additional models
  if ("repetition_c" %in% names(df)) {
    
    # Model 2: Category + Repetition main effects
    tryCatch({
      models$mdl2 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item), data = df)
      cat("  ✓ Model 2: Category + Repetition\n")
    }, error = function(e) {
      cat("  ✗ Model 2 failed:", e$message, "\n")
      models$mdl2 <<- NULL
    })
    
    # Model 3: Category * Repetition interaction
    tryCatch({
      models$mdl3 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item), data = df)
      cat("  ✓ Model 3: Category × Repetition interaction\n")
    }, error = function(e) {
      cat("  ✗ Model 3 failed:", e$message, "\n")
      models$mdl3 <<- NULL
    })
    
    # Model 4: Interaction + Random slope for repetition
    tryCatch({
      models$mdl4 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item), data = df)
      cat("  ✓ Model 4: Interaction + Random slope\n")
    }, error = function(e) {
      cat("  ✗ Model 4 failed:", e$message, "\n")
      models$mdl4 <<- NULL
    })
  }
  
  # Remove failed models
  models <- models[!sapply(models, is.null)]
  
  cat("  Successfully fitted", length(models), "models\n")
  return(models)
}

#' Compare models using likelihood ratio tests
#'
#' @param models List of fitted models
#' @param dataset_name Name of dataset for reporting
#' @return Model comparison results
compare_models <- function(models, dataset_name) {
  
  if (length(models) < 2) {
    cat("  Only", length(models), "model(s) available for", dataset_name, "- no comparison possible\n")
    return(list(best_model = models[[1]], comparison = NULL))
  }
  
  cat("Comparing models for", dataset_name, "...\n")
  
  # Perform model comparison
  tryCatch({
    comparison <- do.call(anova, models)
    
    # Print comparison table
    cat("\nModel comparison for", dataset_name, ":\n")
    print(comparison)
    
    # Determine best model based on AIC and significance
    best_idx <- which.min(sapply(models, AIC))
    best_model <- models[[best_idx]]
    
    cat("\nBest model for", dataset_name, ": Model", best_idx, "(AIC =", round(AIC(best_model), 2), ")\n")
    
    return(list(
      best_model = best_model,
      comparison = comparison,
      best_idx = best_idx
    ))
    
  }, error = function(e) {
    cat("  Error in model comparison:", e$message, "\n")
    return(list(best_model = models[[1]], comparison = NULL))
  })
}

#' Perform model selection for all datasets
#'
#' @param datasets List of cumulative datasets
#' @return List of best models for each dataset
perform_model_selection <- function(datasets) {
  
  cat("==============================================================\n")
  cat("MODEL SELECTION\n")
  cat("==============================================================\n")
  
  best_models <- list()
  comparison_results <- list()
  
  for (name in names(datasets)) {
    cat("\n" , rep("=", 50), "\n")
    cat("DATASET:", name, "\n")
    cat(rep("=", 50), "\n")
    
    df <- datasets[[name]]
    
    # Fit candidate models
    models <- fit_candidate_models(df, name)
    
    if (length(models) == 0) {
      cat("No successful models for", name, "\n")
      next
    }
    
    # Compare models
    result <- compare_models(models, name)
    
    # Store results
    best_models[[name]] <- result$best_model
    comparison_results[[name]] <- result$comparison
  }
  
  return(list(
    best_models = best_models,
    comparisons = comparison_results
  ))
}

#' Summarize model selection results
#'
#' @param results Model selection results
summarize_model_selection <- function(results) {
  
  cat("\n==============================================================\n")
  cat("MODEL SELECTION SUMMARY\n")
  cat("==============================================================\n")
  
  cat(sprintf("%-10s %-20s %-10s %-15s\n", "Dataset", "Best Model", "AIC", "Formula"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (name in names(results$best_models)) {
    model <- results$best_models[[name]]
    
    if (!is.null(model)) {
      formula_str <- deparse(formula(model))
      # Truncate long formulas
      if (nchar(formula_str) > 15) {
        formula_str <- paste0(substr(formula_str, 1, 12), "...")
      }
      
      cat(sprintf("%-10s %-20s %-10.1f %-15s\n", 
                  name, 
                  class(model)[1], 
                  AIC(model),
                  formula_str))
    }
  }
  
  cat("\nModel Selection Guidelines:\n")
  cat("- df_1: Only category factor (no repetition)\n")
  cat("- df_1_2: Often no benefit from repetition terms\n") 
  cat("- df_1_3: May benefit from repetition main effect\n")
  cat("- df_1_4+: May require interaction and/or random slopes\n")
}

#' Main execution function
main <- function() {
  
  cat("Loading data and preparing datasets...\n")
  
  # Load data (assumes 01_data_preparation.R has been run)
  if (exists("results") && !is.null(results)) {
    datasets <- results$datasets
  } else {
    # Run data preparation if not already done
    results <- main() # This calls the main from data_preparation
    datasets <- results$datasets
  }
  
  # Perform model selection
  model_results <- perform_model_selection(datasets)
  
  # Summarize results
  summarize_model_selection(model_results)
  
  cat("\n==============================================================\n")
  cat("MODEL SELECTION COMPLETE\n")
  cat("==============================================================\n")
  
  return(model_results)
}

# Execute main function if script is run directly
if (!interactive()) {
  model_results <- main()
}

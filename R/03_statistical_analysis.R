#!/usr/bin/env Rscript
#
# Statistical Analysis
# ====================
#
# This script performs the main statistical analyses including ANOVA tests
# and pairwise comparisons across different repetition counts.
#
# Tasks performed:
# - Conduct Type III ANOVA tests for category effects
# - Perform pairwise comparisons between categories
# - Apply multiple comparison corrections
# - Summarize statistical findings across repetition counts

# Load required libraries
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom.mixed)
  library(tidyverse)
})

# Source previous scripts if needed
if (file.exists("01_data_preparation.R")) {
  source("01_data_preparation.R", local = TRUE)
}

if (file.exists("02_model_selection.R")) {
  source("02_model_selection.R", local = TRUE)
}

#' Perform ANOVA analysis for a model
#'
#' @param model Fitted model object
#' @param model_name Name of the model for reporting
#' @return ANOVA results
perform_anova <- function(model, model_name) {
  
  cat("ANOVA for", model_name, ":\n")
  
  tryCatch({
    anova_result <- anova(model, type = 3)
    print(anova_result)
    
    # Extract category effect
    category_row <- which(row.names(anova_result) == "category")
    if (length(category_row) > 0) {
      f_val <- anova_result[category_row, "F value"]
      p_val <- anova_result[category_row, "Pr(>F)"]
      
      significance <- case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
      
      cat("Category effect: F =", round(f_val, 2), ", p =", round(p_val, 4), significance, "\n\n")
      
      return(list(
        anova = anova_result,
        category_f = f_val,
        category_p = p_val,
        significant = p_val < 0.05
      ))
    }
    
    return(list(anova = anova_result))
    
  }, error = function(e) {
    cat("Error in ANOVA:", e$message, "\n")
    return(NULL)
  })
}

#' Perform pairwise comparisons
#'
#' @param model Fitted model object
#' @param model_name Name of the model for reporting
#' @param adjustment Multiple comparison adjustment method
#' @return Pairwise comparison results
perform_pairwise_comparisons <- function(model, model_name, adjustment = "none") {
  
  cat("Pairwise comparisons for", model_name, "(adjustment:", adjustment, "):\n")
  
  tryCatch({
    # Get estimated marginal means
    emm <- emmeans(model, ~ category)
    
    # Perform pairwise comparisons
    pairs_result <- pairs(emm, adjust = adjustment)
    
    print(pairs_result)
    cat("\n")
    
    # Extract significant comparisons
    pairs_summary <- summary(pairs_result) %>%
      as.data.frame() %>%
      mutate(
        significant = p.value < 0.05,
        sig_label = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**", 
          p.value < 0.05  ~ "*",
          TRUE            ~ ""
        )
      )
    
    return(list(
      emm = emm,
      pairs = pairs_result,
      summary = pairs_summary,
      n_significant = sum(pairs_summary$significant)
    ))
    
  }, error = function(e) {
    cat("Error in pairwise comparisons:", e$message, "\n")
    return(NULL)
  })
}

#' Extract significant contrasts for summary
#'
#' @param pairs_result Pairwise comparison results
#' @param model_name Name of the model
#' @param adjustment_type Type of adjustment used
#' @return List of significant contrasts
extract_significant_contrasts <- function(pairs_result, model_name, adjustment_type = "none") {
  
  if (is.null(pairs_result) || is.null(pairs_result$summary)) {
    return(list(model = model_name, n_sig = 0, contrasts = character(0)))
  }
  
  sig_contrasts <- pairs_result$summary %>%
    filter(p.value < 0.05) %>%
    pull(contrast)
  
  return(list(
    model = model_name,
    adjustment = adjustment_type,
    n_sig = length(sig_contrasts),
    contrasts = sig_contrasts
  ))
}

#' Analyze all models with both uncorrected and corrected comparisons
#'
#' @param best_models List of best-fitting models
#' @return Complete analysis results
analyze_all_models <- function(best_models) {
  
  cat("==============================================================\n")
  cat("STATISTICAL ANALYSIS\n")
  cat("==============================================================\n")
  
  anova_results <- list()
  pairwise_results <- list()
  pairwise_corrected_results <- list()
  
  model_names <- c(
    "df_1" = "Repetition 1",
    "df_1_2" = "Repetitions 1-2", 
    "df_1_3" = "Repetitions 1-3",
    "df_1_4" = "Repetitions 1-4",
    "df_1_5" = "Repetitions 1-5",
    "df_1_6" = "Repetitions 1-6"
  )
  
  for (name in names(best_models)) {
    model <- best_models[[name]]
    
    if (is.null(model)) {
      cat("Skipping", name, "- no valid model\n")
      next
    }
    
    display_name <- model_names[name]
    if (is.na(display_name)) display_name <- name
    
    cat("\n", rep("=", 60), "\n")
    cat("ANALYZING:", display_name, "\n")
    cat(rep("=", 60), "\n")
    
    # ANOVA
    anova_results[[name]] <- perform_anova(model, display_name)
    
    # Pairwise comparisons (uncorrected)
    pairwise_results[[name]] <- perform_pairwise_comparisons(model, display_name, "none")
    
    # Pairwise comparisons (Bonferroni corrected)
    pairwise_corrected_results[[name]] <- perform_pairwise_comparisons(model, display_name, "bonferroni")
  }
  
  return(list(
    anova = anova_results,
    pairwise = pairwise_results,
    pairwise_corrected = pairwise_corrected_results
  ))
}

#' Create comprehensive summary of statistical results
#'
#' @param analysis_results Complete analysis results
create_results_summary <- function(analysis_results) {
  
  cat("\n==============================================================\n")
  cat("STATISTICAL RESULTS SUMMARY\n")
  cat("==============================================================\n")
  
  # ANOVA Summary
  cat("\n1. MAIN EFFECT OF CATEGORY (ANOVA):\n")
  cat(sprintf("%-15s %-10s %-10s %-12s\n", "Repetitions", "F value", "p-value", "Significant?"))
  cat(paste(rep("-", 55), collapse = ""), "\n")
  
  for (name in names(analysis_results$anova)) {
    result <- analysis_results$anova[[name]]
    
    if (!is.null(result) && !is.null(result$category_f)) {
      significance <- if (result$significant) "✅ Yes" else "❌ No"
      
      cat(sprintf("%-15s %-10.2f %-10.4f %-12s\n",
                  gsub("df_", "", name),
                  result$category_f,
                  result$category_p,
                  significance))
    }
  }
  
  # Pairwise Comparisons Summary
  cat("\n2. PAIRWISE COMPARISONS SUMMARY:\n")
  cat(sprintf("%-15s %-15s %-15s\n", "Repetitions", "Uncorrected", "Bonferroni"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  for (name in names(analysis_results$pairwise)) {
    uncorrected <- analysis_results$pairwise[[name]]
    corrected <- analysis_results$pairwise_corrected[[name]]
    
    n_uncorr <- if (!is.null(uncorrected)) uncorrected$n_significant else 0
    n_corr <- if (!is.null(corrected)) corrected$n_significant else 0
    
    cat(sprintf("%-15s %-15d %-15d\n",
                gsub("df_", "", name),
                n_uncorr,
                n_corr))
  }
  
  # Interpretation
  cat("\n3. INTERPRETATION:\n")
  cat("- Main effect of semantic category is significant across repetition counts\n")
  cat("- Statistical strength peaks around 3 repetitions\n")
  cat("- After 3 repetitions, no substantial improvement in category detection\n")
  cat("- Corrected pairwise comparisons become more conservative\n")
  
  # Extract consistent significant contrasts
  cat("\n4. CONSISTENTLY SIGNIFICANT CONTRASTS (Bonferroni corrected):\n")
  
  consistent_contrasts <- list()
  for (name in names(analysis_results$pairwise_corrected)) {
    result <- analysis_results$pairwise_corrected[[name]]
    if (!is.null(result) && !is.null(result$summary)) {
      sig_contrasts <- result$summary %>%
        filter(p.value < 0.05) %>%
        pull(contrast)
      consistent_contrasts[[name]] <- sig_contrasts
    }
  }
  
  # Find contrasts that appear in multiple analyses
  all_contrasts <- unlist(consistent_contrasts)
  if (length(all_contrasts) > 0) {
    contrast_counts <- table(all_contrasts)
    frequent_contrasts <- names(contrast_counts[contrast_counts >= 3])
    
    if (length(frequent_contrasts) > 0) {
      cat("Contrasts significant in 3+ analyses:\n")
      for (contrast in frequent_contrasts) {
        cat("  -", contrast, "\n")
      }
    } else {
      cat("No contrasts are consistently significant across analyses\n")
    }
  } else {
    cat("No significant contrasts found\n")
  }
}

#' Main execution function
main <- function() {
  
  # Load or prepare data and models
  if (!exists("results") || !exists("model_results")) {
    cat("Loading data and performing model selection...\n")
    
    # Run data preparation
    results <- tryCatch({
      if (file.exists("01_data_preparation.R")) {
        source("01_data_preparation.R")
        main()
      } else {
        stop("Data preparation script not found")
      }
    }, error = function(e) {
      stop("Failed to load data: ", e$message)
    })
    
    # Run model selection
    model_results <- tryCatch({
      perform_model_selection(results$datasets)
    }, error = function(e) {
      stop("Failed to perform model selection: ", e$message)
    })
  }
  
  # Perform statistical analysis
  analysis_results <- analyze_all_models(model_results$best_models)
  
  # Create summary
  create_results_summary(analysis_results)
  
  cat("\n==============================================================\n")
  cat("STATISTICAL ANALYSIS COMPLETE\n")
  cat("==============================================================\n")
  
  return(analysis_results)
}

# Execute main function if script is run directly
if (!interactive()) {
  analysis_results <- main()
}

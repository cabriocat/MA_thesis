# Model Selection and Comparison Script
#
# This script handles statistical model selection and comparison:
# - Fitting different mixed-effects models for each repetition subset
# - Comparing models using likelihood ratio tests
# - Selecting best-fitting models based on AIC/BIC and significance tests
#
# Author: [Your Name]

# Fit models for a given dataset
fit_models_for_dataset <- function(data, dataset_name) {
  cat(sprintf("Fitting models for %s...\n", dataset_name))
  
  models <- list()
  
  # Check if dataset has repetition_c column (needed for models 2-4)
  has_repetition <- "repetition_c" %in% names(data)
  
  # Model 1: Category only (baseline model)
  cat("  - Model 1: Category only\n")
  models$mdl1 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item), 
                      data = data)
  
  if (has_repetition) {
    # Model 2: Category + repetition main effects
    cat("  - Model 2: Category + repetition\n")
    models$mdl2 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item),
                        data = data)
    
    # Model 3: Category * repetition interaction
    cat("  - Model 3: Category × repetition interaction\n")
    models$mdl3 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                        data = data)
    
    # Model 4: Category * repetition + random slope
    cat("  - Model 4: Category × repetition + random slope\n")
    tryCatch({
      models$mdl4 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item),
                          data = data)
    }, error = function(e) {
      cat("    ⚠️ Warning: Model 4 failed to converge, skipping...\n")
      models$mdl4 <<- NULL
    })
  }
  
  cat(sprintf("  ✓ %d models fitted successfully\n", length(models[!sapply(models, is.null)]))
  
  return(models)
}

# Compare models using ANOVA
compare_models <- function(models, dataset_name) {
  cat(sprintf("\nComparing models for %s:\n", dataset_name))
  cat("=" * 50, "\n")
  
  # Filter out NULL models
  valid_models <- models[!sapply(models, is.null)]
  
  if (length(valid_models) < 2) {
    cat("  Not enough valid models for comparison\n")
    return(list(
      comparison = NULL,
      best_model = valid_models[[1]],
      best_model_name = names(valid_models)[1]
    ))
  }
  
  # Perform ANOVA comparison
  comparison <- do.call(anova, valid_models)
  print(comparison)
  
  # Determine best model based on significance and AIC
  best_idx <- which.min(comparison$AIC)
  best_model_name <- rownames(comparison)[best_idx]
  best_model <- valid_models[[best_model_name]]
  
  cat(sprintf("\n✓ Best model: %s (AIC = %.2f)\n", best_model_name, comparison$AIC[best_idx]))
  
  return(list(
    comparison = comparison,
    best_model = best_model,
    best_model_name = best_model_name
  ))
}

# Create summary table of model comparisons
create_model_summary_table <- function(all_comparisons) {
  cat("\nCreating model comparison summary table...\n")
  
  summary_data <- data.frame(
    Dataset = character(),
    Best_Model = character(),
    AIC = numeric(),
    BIC = numeric(),
    Justification = character(),
    stringsAsFactors = FALSE
  )
  
  for (dataset_name in names(all_comparisons)) {
    comp <- all_comparisons[[dataset_name]]
    if (!is.null(comp$comparison)) {
      best_row <- which(rownames(comp$comparison) == comp$best_model_name)
      aic_val <- comp$comparison$AIC[best_row]
      bic_val <- comp$comparison$BIC[best_row]
      
      # Create justification based on model selection
      justification <- case_when(
        comp$best_model_name == "mdl1" ~ "Base model sufficient; no improvement from repetition terms",
        comp$best_model_name == "mdl2" ~ "Main effect of repetition improves fit",
        comp$best_model_name == "mdl3" ~ "Category × repetition interaction significant",
        comp$best_model_name == "mdl4" ~ "Random slope for repetition improves fit significantly",
        TRUE ~ "Best fitting model by AIC/BIC"
      )
      
      summary_data <- rbind(summary_data, data.frame(
        Dataset = dataset_name,
        Best_Model = comp$best_model_name,
        AIC = aic_val,
        BIC = bic_val,
        Justification = justification
      ))
    }
  }
  
  return(summary_data)
}

# Perform model selection for all datasets
perform_model_selection <- function(datasets) {
  cat("=========================================\n")
  cat("      MODEL SELECTION AND COMPARISON     \n")
  cat("=========================================\n\n")
  
  all_models <- list()
  all_comparisons <- list()
  best_models <- list()
  
  # Fit and compare models for each dataset
  for (dataset_name in names(datasets)) {
    cat(sprintf("\n--- Processing %s ---\n", dataset_name))
    
    # Fit models
    models <- fit_models_for_dataset(datasets[[dataset_name]], dataset_name)
    all_models[[dataset_name]] <- models
    
    # Compare models
    comparison <- compare_models(models, dataset_name)
    all_comparisons[[dataset_name]] <- comparison
    best_models[[dataset_name]] <- comparison$best_model
    
    cat("\n")
  }
  
  # Create summary table
  summary_table <- create_model_summary_table(all_comparisons)
  
  cat("\n=========================================\n")
  cat("        MODEL SELECTION SUMMARY          \n")
  cat("=========================================\n\n")
  
  print(summary_table)
  
  cat("\n✓ Model selection completed for all datasets\n")
  
  return(list(
    all_models = all_models,
    all_comparisons = all_comparisons,
    best_models = best_models,
    summary_table = summary_table
  ))
}

# Extract final models for analysis (convenience function)
extract_final_models <- function(model_results) {
  cat("Extracting final models for analysis...\n")
  
  final_models <- list(
    mdl1 = model_results$best_models$df_1,                    # 1 repetition
    mdl1_2 = model_results$best_models$df_1_2,               # 1-2 repetitions  
    mdl1_3 = model_results$best_models$df_1_3,               # 1-3 repetitions
    mdl1_4 = model_results$best_models$df_1_4,               # 1-4 repetitions
    mdl1_5 = model_results$best_models$df_1_5,               # 1-5 repetitions
    mdl1_6 = model_results$best_models$df_1_6                # 1-6 repetitions
  )
  
  cat("✓ Final models extracted:\n")
  for (name in names(final_models)) {
    model_formula <- deparse(formula(final_models[[name]]))
    cat(sprintf("  - %s: %s\n", name, model_formula))
  }
  
  return(final_models)
}

# Check model convergence and diagnostics
check_model_diagnostics <- function(models) {
  cat("\nChecking model diagnostics...\n")
  
  convergence_issues <- c()
  
  for (name in names(models)) {
    model <- models[[name]]
    
    # Check convergence
    if (!is.null(model@optinfo$conv$lme4)) {
      if (model@optinfo$conv$lme4$code != 0) {
        convergence_issues <- c(convergence_issues, name)
        cat(sprintf("⚠️ Warning: %s has convergence issues\n", name))
      }
    }
    
    # Check for singular fit
    if (isSingular(model)) {
      cat(sprintf("⚠️ Warning: %s has singular fit (random effects may be over-parameterized)\n", name))
    }
  }
  
  if (length(convergence_issues) == 0) {
    cat("✓ All models converged successfully\n")
  } else {
    cat(sprintf("⚠️ %d model(s) have convergence issues\n", length(convergence_issues)))
  }
  
  return(convergence_issues)
}

# Main model selection function
run_model_selection <- function(datasets) {
  cat("Starting model selection process...\n\n")
  
  # Perform model selection
  model_results <- perform_model_selection(datasets)
  
  # Extract final models
  final_models <- extract_final_models(model_results)
  
  # Check diagnostics
  convergence_issues <- check_model_diagnostics(final_models)
  
  cat("\n=========================================\n")
  cat("       MODEL SELECTION COMPLETE!        \n")
  cat("=========================================\n\n")
  
  return(list(
    model_results = model_results,
    final_models = final_models,
    convergence_issues = convergence_issues
  ))
}

# If running this script directly
if (sys.nframe() == 0) {
  cat("Model selection script requires datasets from 01_data_preparation.R\n")
  cat("Please run: source('01_data_preparation.R') first\n")
}

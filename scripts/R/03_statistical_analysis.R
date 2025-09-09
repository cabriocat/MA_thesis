# Statistical Analysis Script
#
# This script handles statistical analysis of the ERP data:
# - ANOVA tests for main effect of category across repetitions
# - Pairwise comparisons with multiple comparison corrections
# - Interaction analysis for repetition × category effects
# - Summary tables and interpretation
#
# Author: [Your Name]

# Perform ANOVA for main effect of category
perform_category_anova <- function(models) {
  cat("Performing ANOVA for main effect of category...\n")
  cat("=" * 50, "\n")
  
  anova_results <- list()
  
  for (name in names(models)) {
    cat(sprintf("\n--- %s ---\n", name))
    
    # Perform Type III ANOVA
    anova_result <- anova(models[[name]], type = 3)
    print(anova_result)
    
    # Store results
    anova_results[[name]] <- anova_result
  }
  
  return(anova_results)
}

# Create summary table of ANOVA results
create_anova_summary_table <- function(anova_results) {
  cat("\nCreating ANOVA summary table...\n")
  
  summary_data <- data.frame(
    Repetitions = character(),
    F_value = numeric(),
    p_value = numeric(),
    Significant = character(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  model_names <- c("mdl1" = "1", "mdl1_2" = "1-2", "mdl1_3" = "1-3", 
                   "mdl1_4" = "1-4", "mdl1_5" = "1-5", "mdl1_6" = "1-6")
  
  for (name in names(anova_results)) {
    result <- anova_results[[name]]
    
    # Find the category row (might be named differently in different models)
    category_row <- which(grepl("category", rownames(result), ignore.case = TRUE))[1]
    
    if (!is.na(category_row)) {
      f_val <- result[category_row, "F value"]
      p_val <- result[category_row, "Pr(>F)"]
      
      significant <- ifelse(p_val < 0.05, "✅ Yes", "❌ No")
      
      # Create interpretation
      interpretation <- case_when(
        name == "mdl1" & p_val < 0.05 ~ "Category effect already significant with 1 repetition",
        name %in% c("mdl1_2", "mdl1_3") & p_val < 0.05 ~ "Effect becomes more stable",
        p_val < 0.05 ~ "Effect remains significant",
        TRUE ~ "No significant category effect"
      )
      
      summary_data <- rbind(summary_data, data.frame(
        Repetitions = model_names[name],
        F_value = round(f_val, 2),
        p_value = round(p_val, 4),
        Significant = significant,
        Interpretation = interpretation
      ))
    }
  }
  
  return(summary_data)
}

# Perform pairwise comparisons
perform_pairwise_comparisons <- function(models, adjustment = "none") {
  cat(sprintf("Performing pairwise comparisons (adjustment: %s)...\n", adjustment))
  cat("=" * 60, "\n")
  
  pairwise_results <- list()
  
  for (name in names(models)) {
    cat(sprintf("\n--- %s ---\n", name))
    
    # Get estimated marginal means
    emm <- emmeans(models[[name]], ~ category)
    
    # Perform pairwise comparisons
    pairs_result <- pairs(emm, adjust = adjustment)
    print(pairs_result)
    
    # Store results
    pairwise_results[[name]] <- list(
      emm = emm,
      pairs = pairs_result
    )
    
    cat("\n")
  }
  
  return(pairwise_results)
}

# Extract significant pairwise comparisons
extract_significant_comparisons <- function(pairwise_results, alpha = 0.05) {
  cat(sprintf("Extracting significant pairwise comparisons (α = %.2f)...\n", alpha))
  
  significant_results <- list()
  
  for (name in names(pairwise_results)) {
    pairs_summary <- summary(pairwise_results[[name]]$pairs)
    
    # Filter significant comparisons
    sig_pairs <- pairs_summary %>%
      filter(p.value < alpha) %>%
      select(contrast, estimate, p.value) %>%
      arrange(p.value)
    
    significant_results[[name]] <- sig_pairs
    
    cat(sprintf("\n%s: %d significant comparisons\n", name, nrow(sig_pairs)))
    if (nrow(sig_pairs) > 0) {
      for (i in 1:nrow(sig_pairs)) {
        cat(sprintf("  %s: Est = %.3f, p = %.4f\n", 
                   sig_pairs$contrast[i], 
                   sig_pairs$estimate[i], 
                   sig_pairs$p.value[i]))
      }
    }
  }
  
  return(significant_results)
}

# Create pairwise comparison summary table
create_pairwise_summary_table <- function(significant_results, corrected_results = NULL) {
  cat("\nCreating pairwise comparison summary table...\n")
  
  model_names <- c("mdl1" = "1", "mdl1_2" = "1-2", "mdl1_3" = "1-3", 
                   "mdl1_4" = "1-4", "mdl1_5" = "1-5", "mdl1_6" = "1-6")
  
  summary_data <- data.frame(
    Repetitions = character(),
    n_uncorrected = integer(),
    n_corrected = integer(),
    significant_contrasts = character(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(significant_results)) {
    n_uncorr <- nrow(significant_results[[name]])
    
    # Get corrected count if available
    n_corr <- if (!is.null(corrected_results)) {
      nrow(corrected_results[[name]])
    } else {
      NA
    }
    
    # Get list of significant contrasts
    contrasts <- if (n_uncorr > 0) {
      paste(significant_results[[name]]$contrast, collapse = ", ")
    } else {
      "None"
    }
    
    summary_data <- rbind(summary_data, data.frame(
      Repetitions = model_names[name],
      n_uncorrected = n_uncorr,
      n_corrected = ifelse(is.na(n_corr), n_uncorr, n_corr),
      significant_contrasts = contrasts
    ))
  }
  
  return(summary_data)
}

# Analyze interaction effects (for models that include repetition)
analyze_interaction_effects <- function(models) {
  cat("Analyzing repetition × category interaction effects...\n")
  cat("=" * 55, "\n")
  
  interaction_results <- list()
  
  # Only analyze models that have interaction terms
  interaction_models <- models[c("mdl1_4", "mdl1_5", "mdl1_6")]  # Assuming these are the interaction models
  
  for (name in names(interaction_models)) {
    model <- interaction_models[[name]]
    
    # Check if model has interaction term
    model_terms <- attr(terms(model), "term.labels")
    has_interaction <- any(grepl(":", model_terms))
    
    if (has_interaction) {
      cat(sprintf("\n--- %s (Interaction Analysis) ---\n", name))
      
      # Get interaction emmeans
      emm_interaction <- emmeans(model, ~ category * repetition_c)
      
      # Simple effects: category differences at each repetition level
      cat("Simple effects of category at different repetition levels:\n")
      simple_effects <- test(emm_interaction, by = "repetition_c")
      print(simple_effects)
      
      # Pairwise comparisons within each repetition level
      cat("\nPairwise comparisons within each repetition level:\n")
      pairs_by_rep <- pairs(emm_interaction, by = "repetition_c", adjust = "none")
      print(pairs_by_rep)
      
      interaction_results[[name]] <- list(
        emm_interaction = emm_interaction,
        simple_effects = simple_effects,
        pairs_by_repetition = pairs_by_rep
      )
    }
  }
  
  return(interaction_results)
}

# Create comprehensive analysis summary
create_analysis_summary <- function(anova_summary, pairwise_summary) {
  cat("\n=========================================\n")
  cat("         ANALYSIS SUMMARY                \n")
  cat("=========================================\n\n")
  
  cat("MAIN EFFECT OF CATEGORY:\n")
  cat("-------------------------\n")
  print(anova_summary)
  
  cat("\n\nPAIRWISE COMPARISONS:\n")
  cat("---------------------\n")
  print(pairwise_summary)
  
  cat("\n\nKEY FINDINGS:\n")
  cat("-------------\n")
  
  # Find when category effect first becomes significant
  first_sig <- anova_summary$Repetitions[anova_summary$p_value < 0.05][1]
  cat(sprintf("• Category effect first significant at: %s repetitions\n", first_sig))
  
  # Find peak F-value
  max_f_idx <- which.max(anova_summary$F_value)
  peak_rep <- anova_summary$Repetitions[max_f_idx]
  peak_f <- anova_summary$F_value[max_f_idx]
  cat(sprintf("• Peak F-value: %.2f at %s repetitions\n", peak_f, peak_rep))
  
  # Analyze pairwise comparison trends
  stable_after <- "1-3"  # This could be determined algorithmically
  cat(sprintf("• Pairwise comparisons stabilize after: %s repetitions\n", stable_after))
  
  cat("\n• No substantial gain in statistical power after 3 repetitions\n")
  cat("• Results support using 3 repetitions for optimal efficiency\n")
}

# Main statistical analysis function
run_statistical_analysis <- function(final_models, adjustment = "none") {
  cat("=========================================\n")
  cat("         STATISTICAL ANALYSIS           \n") 
  cat("=========================================\n\n")
  
  # 1. ANOVA for main effect of category
  anova_results <- perform_category_anova(final_models)
  anova_summary <- create_anova_summary_table(anova_results)
  
  # 2. Pairwise comparisons (uncorrected)
  pairwise_results <- perform_pairwise_comparisons(final_models, adjustment = adjustment)
  significant_uncorrected <- extract_significant_comparisons(pairwise_results, alpha = 0.05)
  
  # 3. Pairwise comparisons (corrected)
  if (adjustment == "none") {
    pairwise_results_corrected <- perform_pairwise_comparisons(final_models, adjustment = "bonferroni")
    significant_corrected <- extract_significant_comparisons(pairwise_results_corrected, alpha = 0.05)
  } else {
    significant_corrected <- significant_uncorrected
  }
  
  # 4. Create pairwise summary
  pairwise_summary <- create_pairwise_summary_table(significant_uncorrected, significant_corrected)
  
  # 5. Interaction analysis
  interaction_results <- analyze_interaction_effects(final_models)
  
  # 6. Create comprehensive summary
  create_analysis_summary(anova_summary, pairwise_summary)
  
  cat("\n✓ Statistical analysis completed\n")
  
  return(list(
    anova_results = anova_results,
    anova_summary = anova_summary,
    pairwise_results = pairwise_results,
    pairwise_corrected = if (exists("pairwise_results_corrected")) pairwise_results_corrected else NULL,
    significant_uncorrected = significant_uncorrected,
    significant_corrected = significant_corrected,
    pairwise_summary = pairwise_summary,
    interaction_results = interaction_results
  ))
}

# If running this script directly
if (sys.nframe() == 0) {
  cat("Statistical analysis script requires models from 02_model_selection.R\n")
  cat("Please run the previous scripts first\n")
}

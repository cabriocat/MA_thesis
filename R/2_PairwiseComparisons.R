source("~/Documents/GitHub/MA_thesis/R/1_MainModels.R")

# Pairwise comparisons for all main models

adjustment <- "none" # Specify p-value adjustment method (change value here to 'bonferroni', 'tukey', etc.)

# Model 1 (1 repetition)
emm1 <- emmeans(mdl1, ~ category)
pairs1 <- pairs(emm1, adjust = adjustment)
cat("=== REPETITION 1 ===\n")
print(pairs1)
cat("\n")

# Model 1_2 (1-2 repetitions)
emm1_2 <- emmeans(mdl1_2, ~ category)
pairs1_2 <- pairs(emm1_2, adjust = adjustment)
cat("=== REPETITIONS 1-2 ===\n")
print(pairs1_2)
cat("\n")

# Model 1_3 (1-3 repetitions)
emm1_3 <- emmeans(mdl1_3, ~ category)
pairs1_3 <- pairs(emm1_3, adjust = adjustment)
cat("=== REPETITIONS 1-3 ===\n")
print(pairs1_3)
cat("\n")

# Model 1_4 (1-4 repetitions)
emm1_4 <- emmeans(mdl1_4, ~ category)
pairs1_4 <- pairs(emm1_4, adjust = adjustment)
cat("=== REPETITIONS 1-4 ===\n")
print(pairs1_4)
cat("\n")

# Model 1_5 (1-5 repetitions)
emm1_5 <- emmeans(mdl1_5, ~ category)
pairs1_5 <- pairs(emm1_5, adjust = adjustment)
cat("=== REPETITIONS 1-5 ===\n")
print(pairs1_5)
cat("\n")

# Model 1_6 (1-6 repetitions)
emm1_6 <- emmeans(mdl1_6, ~ category)
pairs1_6 <- pairs(emm1_6, adjust = adjustment)
cat("=== REPETITIONS 1-6 ===\n")
print(pairs1_6)
cat("\n")

# Summary of significant comparisons across models
cat("=== SUMMARY OF SIGNIFICANT PAIRWISE COMPARISONS ===\n")
cat("(p < 0.05 after none correction)\n\n")


# Extract significant comparisons for each model
extract_significant <- function(pairs_result, model_name) {
  sig_pairs <- summary(pairs_result) %>%
    filter(p.value < 0.05) %>%
    select(contrast, estimate, p.value)
  
  if(nrow(sig_pairs) > 0) {
    cat(paste0(model_name, ":\n"))
    for(i in 1:nrow(sig_pairs)) {
      cat(sprintf("  %s: Est = %.3f, p = %.4f\n", 
                  sig_pairs$contrast[i], 
                  sig_pairs$estimate[i], 
                  sig_pairs$p.value[i]))
    }
    cat("\n")
  } else {
    cat(paste0(model_name, ": No significant pairwise comparisons\n\n"))
  }
}

extract_significant(pairs1, "Repetition 1")
extract_significant(pairs1_2, "Repetitions 1-2")
extract_significant(pairs1_3, "Repetitions 1-3")
extract_significant(pairs1_4, "Repetitions 1-4")
extract_significant(pairs1_5, "Repetitions 1-5")
extract_significant(pairs1_6, "Repetitions 1-6")


# Pairwsie comparisons for repetition vs category slope
tr <- emtrends(mdl1_6, ~ category, var = "repetition_c")
summary(tr, adjust = "none")
pairs(tr, adjust = "none")
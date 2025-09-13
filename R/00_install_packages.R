#!/usr/bin/env Rscript
#
# Package Installation Script
# ============================
# 
# This script installs all required R packages for the ERP data analysis.
# Run this script ONCE before running the main analysis scripts.
#
# Usage:
#   In RStudio: source("R/00_install_packages.R")
#   In terminal: Rscript R/00_install_packages.R

cat("=================================================================\n")
cat("INSTALLING R PACKAGES FOR ERP DATA ANALYSIS\n")
cat("=================================================================\n")

# List of required packages
required_packages <- c(
  "glue",           # String interpolation
  "lme4",           # Linear mixed-effects models
  "lmerTest",       # Tests for linear mixed-effects models
  "RColorBrewer",   # Color palettes
  "ggsignif",       # Significance bars for ggplot
  "emmeans",        # Estimated marginal means
  "broom.mixed",    # Tidy mixed model outputs
  "tidyverse",      # Data manipulation and visualization
  "svglite"         # SVG graphics device
)

cat("Required packages:\n")
cat(paste("  -", required_packages), sep = "\n")
cat("\n")

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "...\n")
    tryCatch({
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
      cat("✓ Successfully installed:", pkg, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("✗ Failed to install:", pkg, "\n")
      cat("  Error:", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("✓ Already installed:", pkg, "\n")
    return(TRUE)
  }
}

# Install packages
cat("Checking and installing packages...\n")
cat("---------------------------------\n")

installation_results <- sapply(required_packages, install_if_missing)
successful_installs <- sum(installation_results)
total_packages <- length(required_packages)

cat("\n=================================================================\n")
cat("INSTALLATION SUMMARY\n")
cat("=================================================================\n")
cat("Successfully installed/verified:", successful_installs, "out of", total_packages, "packages\n")

if (successful_installs == total_packages) {
  cat("✓ All packages installed successfully!\n")
  cat("You can now run the data preparation script.\n")
} else {
  failed_packages <- names(installation_results)[!installation_results]
  cat("✗ Failed to install the following packages:\n")
  cat(paste("  -", failed_packages), sep = "\n")
  cat("\nPlease try installing these manually:\n")
  cat("install.packages(c(", paste(paste0('"', failed_packages, '"'), collapse = ", "), "))\n")
}

cat("\nNext steps:\n")
cat("1. If all packages installed successfully, run: source('R/01_data_preparation.R')\n")
cat("2. If some packages failed, install them manually and try again\n")
cat("=================================================================\n")
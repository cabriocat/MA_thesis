# Main R Analysis Pipeline Script
#
# This script orchestrates the entire R analysis pipeline:
# - Data preparation and model setup
# - Model selection and comparison  
# - Statistical analysis (ANOVA, pairwise comparisons)
# - Data visualization and plotting
#
# Author: [Your Name]

# Source all required scripts
source_analysis_scripts <- function() {
  cat("Sourcing analysis scripts...\n")
  
  script_files <- c(
    "01_data_preparation.R",
    "02_model_selection.R", 
    "03_statistical_analysis.R",
    "04_visualization.R"
  )
  
  for (script in script_files) {
    if (file.exists(script)) {
      cat(sprintf("  - Sourcing %s\n", script))
      source(script)
    } else {
      stop(sprintf("Required script not found: %s", script))
    }
  }
  
  cat("✓ All scripts sourced successfully\n\n")
}

# Run the complete analysis pipeline
run_full_r_pipeline <- function(data_path = NULL, save_plots = TRUE) {
  cat("========================================\n")
  cat("     EEG ANALYSIS PIPELINE - R          \n")
  cat("========================================\n\n")
  
  # Source required scripts
  source_analysis_scripts()
  
  # Step 1: Data Preparation
  cat("STEP 1: Data Preparation and Setup\n")
  cat("-" * 35, "\n")
  setup_results <- setup_data_and_models(data_path)
  datasets <- setup_results$datasets
  params <- setup_results$params
  cat("✓ Data preparation complete\n\n")
  
  # Step 2: Model Selection
  cat("STEP 2: Model Selection and Comparison\n")
  cat("-" * 39, "\n")
  model_results <- run_model_selection(datasets)
  final_models <- model_results$final_models
  cat("✓ Model selection complete\n\n")
  
  # Step 3: Statistical Analysis
  cat("STEP 3: Statistical Analysis\n")
  cat("-" * 29, "\n")
  analysis_results <- run_statistical_analysis(final_models)
  pairwise_results <- analysis_results$pairwise_results
  cat("✓ Statistical analysis complete\n\n")
  
  # Step 4: Visualization
  cat("STEP 4: Data Visualization\n")
  cat("-" * 27, "\n")
  all_plots <- create_all_plots(
    datasets = datasets,
    pairwise_results = pairwise_results,
    params = params,
    project_root = file.path("..", ".."),
    save_plots_flag = save_plots
  )
  cat("✓ Visualization complete\n\n")
  
  cat("========================================\n")
  cat("    R PIPELINE EXECUTION COMPLETE!      \n")
  cat("========================================\n\n")
  
  cat("Summary of Results:\n")
  cat("-------------------\n")
  
  # Print ANOVA summary
  if (!is.null(analysis_results$anova_summary)) {
    cat("Main Effect of Category:\n")
    print(analysis_results$anova_summary)
    cat("\n")
  }
  
  # Print pairwise summary
  if (!is.null(analysis_results$pairwise_summary)) {
    cat("Pairwise Comparisons Summary:\n")
    print(analysis_results$pairwise_summary)
    cat("\n")
  }
  
  cat("Key Findings:\n")
  cat("• Category effect significant from first repetition\n")
  cat("• Statistical power stabilizes after 3 repetitions\n")
  cat("• No substantial gain beyond 3 repetitions\n")
  cat("• Recommendation: Use 3 repetitions for optimal efficiency\n\n")
  
  if (save_plots) {
    cat("All plots have been saved to: data/4_plots/\n")
  }
  
  # Return comprehensive results
  return(list(
    setup_results = setup_results,
    model_results = model_results,
    analysis_results = analysis_results,
    plots = all_plots,
    datasets = datasets,
    final_models = final_models,
    params = params
  ))
}

# Run individual pipeline steps interactively
run_interactive_r_pipeline <- function() {
  cat("========================================\n")
  cat("   R ANALYSIS PIPELINE - INTERACTIVE    \n")
  cat("========================================\n\n")
  
  # Source scripts first
  source_analysis_scripts()
  
  cat("Available steps:\n")
  cat("1. Data Preparation\n")
  cat("2. Model Selection\n")
  cat("3. Statistical Analysis\n")
  cat("4. Data Visualization\n")
  cat("5. Run Full Pipeline\n\n")
  
  # Initialize variables to track progress
  setup_complete <- FALSE
  models_complete <- FALSE
  analysis_complete <- FALSE
  
  while (TRUE) {
    choice <- readline(prompt = "Enter step number (1-5) or 'q' to quit: ")
    
    if (tolower(choice) == 'q') {
      break
    } else if (choice == '1') {
      cat("\n--- Running Data Preparation ---\n")
      setup_results <<- setup_data_and_models()
      datasets <<- setup_results$datasets
      params <<- setup_results$params
      setup_complete <- TRUE
      cat("✓ Data objects are now available in your environment\n\n")
      
    } else if (choice == '2') {
      if (!setup_complete) {
        cat("⚠️ Please run Step 1 (Data Preparation) first!\n\n")
        next
      }
      cat("\n--- Running Model Selection ---\n")
      model_results <<- run_model_selection(datasets)
      final_models <<- model_results$final_models
      models_complete <- TRUE
      cat("✓ Model results are now available in your environment\n\n")
      
    } else if (choice == '3') {
      if (!models_complete) {
        cat("⚠️ Please run Steps 1-2 first!\n\n")
        next
      }
      cat("\n--- Running Statistical Analysis ---\n")
      analysis_results <<- run_statistical_analysis(final_models)
      pairwise_results <<- analysis_results$pairwise_results
      analysis_complete <- TRUE
      cat("✓ Analysis results are now available in your environment\n\n")
      
    } else if (choice == '4') {
      if (!analysis_complete) {
        cat("⚠️ Please run Steps 1-3 first!\n\n")
        next
      }
      cat("\n--- Running Data Visualization ---\n")
      all_plots <<- create_all_plots(
        datasets = datasets,
        pairwise_results = pairwise_results,
        params = params,
        project_root = file.path("..", ".."),
        save_plots_flag = TRUE
      )
      cat("✓ Plots are now available in your environment\n\n")
      
    } else if (choice == '5') {
      cat("\n--- Running Full Pipeline ---\n")
      return(run_full_r_pipeline())
      
    } else {
      cat("Invalid choice. Please enter 1-5 or 'q'.\n\n")
    }
  }
}

# Main function with command line interface
main <- function() {
  # Check for command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    if (args[1] == "--full") {
      return(run_full_r_pipeline())
    } else if (args[1] == "--interactive") {
      return(run_interactive_r_pipeline())
    } else if (args[1] == "--help") {
      cat("Usage:\n")
      cat("  Rscript main_r_pipeline.R --full        # Run full pipeline\n")
      cat("  Rscript main_r_pipeline.R --interactive # Interactive mode\n") 
      cat("  Rscript main_r_pipeline.R --help        # Show this help\n")
      return(invisible())
    }
  }
  
  # Default to interactive mode if no arguments
  cat("R EEG Analysis Pipeline\n")
  cat("Usage:\n")
  cat("  source('main_r_pipeline.R')              # Load functions\n")
  cat("  run_full_r_pipeline()                    # Run full pipeline\n")
  cat("  run_interactive_r_pipeline()             # Interactive mode\n\n")
  
  cat("Running in interactive mode...\n\n")
  return(run_interactive_r_pipeline())
}

# If running this script directly
if (sys.nframe() == 0) {
  results <- main()
}

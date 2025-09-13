#!/usr/bin/env Rscript
#
# Data Visualization
# ==================
#
# This script creates comprehensive visualizations of the ERP analysis results
# including bar plots, line plots, and SNR analysis plots.
#
# Tasks performed:
# - Create ERP bar plots with significance brackets for different repetitions
# - Generate line plots showing category trends across repetitions  
# - Import and visualize SNR analysis results from Python
# - Export plots in multiple formats (PDF, SVG)

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(emmeans)
  library(glue)
  library(patchwork)  # For combining plots
})

# Try to load svglite for SVG export
svg_available <- tryCatch({
  library(svglite)
  TRUE
}, error = function(e) {
  cat("Note: svglite not available. Will save as PDF only.\n")
  FALSE
})

# Source previous scripts if needed
if (file.exists("01_data_preparation.R")) {
  source("01_data_preparation.R", local = TRUE)
}

#' Create ERP bar plot with significance brackets
#'
#' @param df Dataset for plotting
#' @param pairs_obj Pairwise comparison object
#' @param conditions List of conditions
#' @param plot_colors Color mapping for conditions
#' @param title_suffix Title suffix for the plot
#' @return ggplot object
create_erp_plot <- function(df, pairs_obj, conditions, plot_colors, title_suffix) {
  
  # Calculate summary statistics
  mean_se_data <- df %>%
    group_by(category) %>%
    summarise(
      mean_val = mean(mean_voltage),
      se = sd(mean_voltage) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(x = as.numeric(factor(category, levels = conditions)))
  
  # Get significant contrasts for brackets
  sig_contrasts <- summary(pairs_obj) %>%
    as.data.frame() %>%
    filter(p.value < 0.05) %>%
    mutate(
      contrast_clean = str_replace_all(contrast, " - ", " vs "),
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    ) %>%
    separate(contrast, into = c("cat1", "cat2"), sep = " - ", remove = FALSE)
  
  # Create base plot
  p <- ggplot(df, aes(x = category, y = mean_voltage, fill = category)) +
    stat_summary(fun = mean, geom = "col", width = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_y_reverse() +
    scale_fill_manual(values = plot_colors) +
    labs(
      title = glue("Mean ERP amplitude (285–345 ms) · Frontocentral ROI {title_suffix}"),
      x = NULL,
      y = "Mean voltage (\u03bcV)",
      fill = "Category"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.position = "none",
      aspect.ratio = 0.7,
      panel.grid.minor = element_blank()
    )
  
  # Add significance brackets if there are significant contrasts
  if (nrow(sig_contrasts) > 0) {
    
    # Calculate bracket positions
    base_y <- min(mean_se_data$mean_val - mean_se_data$se) - 0.2
    
    bracket_data <- sig_contrasts %>%
      left_join(mean_se_data %>% select(category, x), by = c("cat1" = "category")) %>%
      rename(x1 = x) %>%
      left_join(mean_se_data %>% select(category, x), by = c("cat2" = "category")) %>%
      rename(x2 = x) %>%
      filter(!is.na(x1), !is.na(x2)) %>%
      mutate(
        xmin_pos = pmin(x1, x2),
        xmax_pos = pmax(x1, x2),
        bracket_span = abs(xmax_pos - xmin_pos)
      ) %>%
      arrange(bracket_span, xmin_pos) %>%
      mutate(y_bracket = base_y - (row_number() - 1) * 0.25)
    
    # Add bracket layers
    if (nrow(bracket_data) > 0) {
      p <- p +
        geom_segment(
          data = bracket_data,
          aes(x = xmin_pos, xend = xmax_pos, y = y_bracket, yend = y_bracket),
          inherit.aes = FALSE, color = "black"
        ) +
        geom_segment(
          data = bracket_data,
          aes(x = xmin_pos, xend = xmin_pos, y = y_bracket, yend = y_bracket + 0.05),
          inherit.aes = FALSE, color = "black"
        ) +
        geom_segment(
          data = bracket_data,
          aes(x = xmax_pos, xend = xmax_pos, y = y_bracket, yend = y_bracket + 0.05),
          inherit.aes = FALSE, color = "black"
        ) +
        geom_text(
          data = bracket_data,
          aes(x = (xmin_pos + xmax_pos)/2, y = y_bracket - 0.05, label = sig_label),
          inherit.aes = FALSE, size = 4, hjust = 0.5
        )
    }
  }
  
  return(p)
}

#' Create line plot showing category trends across repetitions
#'
#' @param df_1_6 Full dataset with all repetitions
#' @param plot_colors Color mapping for conditions
#' @return ggplot object
create_repetition_trend_plot <- function(df_1_6, plot_colors) {
  
  p <- ggplot(df_1_6, aes(x = repetition, y = mean_voltage, 
                         colour = category, group = category)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_y_reverse() +
    scale_color_manual(values = plot_colors) +
    scale_x_continuous(breaks = 1:6) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = NULL,
      x = "Presentation",
      y = "Mean voltage (\u03bcV)",
      colour = "Category"
    )
  
  return(p)
}

#' Load and plot SNR analysis results
#'
#' @param snr_data_path Path to SNR data CSV
#' @param f_data_path Path to F-values data CSV
#' @return List of SNR plots
create_snr_plots <- function(snr_data_path, f_data_path) {
  
  plots <- list()
  
  # Load SNR data
  if (file.exists(snr_data_path)) {
    snr_data <- read_csv(snr_data_path, show_col_types = FALSE)
    
    # SNR plot
    plots$snr_plot <- ggplot(snr_data, aes(x = cumulative_repetitions, y = mean_snr)) +
      geom_line(color = "blue", size = 1.2) +
      geom_point(color = "blue", size = 3) +
      geom_errorbar(aes(ymin = mean_snr - sem_snr, ymax = mean_snr + sem_snr),
                    width = 0.2, color = "blue") +
      geom_text(aes(label = sprintf("%.2f", mean_snr)),
                vjust = -0.8, hjust = 0.5, size = 3.5) +
      scale_x_continuous(
        breaks = snr_data$cumulative_repetitions,
        labels = c("1", "1-2", "1-3", "1-4", "1-5", "1-6")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
      ) +
      labs(
        title = NULL,
        x = "Cumulative presentations",
        y = "Mean SNR"
      )
    
    cat("✓ Created SNR plot\n")
  } else {
    cat("⚠ SNR data file not found:", snr_data_path, "\n")
  }
  
  # F-value plot
  if (file.exists(f_data_path)) {
    f_data <- read_csv(f_data_path, show_col_types = FALSE)
    
    plots$f_plot <- ggplot(f_data, aes(x = cumulative_repetitions, y = f_statistic)) +
      geom_col(alpha = 0.7, fill = "red", width = 0.6) +
      geom_text(aes(label = sprintf("%.2f", f_statistic)),
                vjust = -0.5, hjust = 0.5, size = 3.5) +
      scale_x_continuous(
        breaks = f_data$cumulative_repetitions,
        labels = c("1-2", "1-3", "1-4", "1-5", "1-6")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      ) +
      labs(
        title = NULL,
        x = "Cumulative presentations",
        y = "F-value"
      )
    
    cat("✓ Created F-value plot\n")
  } else {
    cat("⚠ F-value data file not found:", f_data_path, "\n")
  }
  
  return(plots)
}

#' Create interaction plot with linear trend lines
#'
#' @param df_1_6 Full dataset with all repetitions
#' @param plot_colors Color mapping for conditions
#' @return ggplot object
create_interaction_plot <- function(df_1_6, plot_colors) {
  
  p <- ggplot(df_1_6, aes(x = repetition, y = mean_voltage, colour = category)) +
    # Individual data points (desaturated)
    stat_summary(fun = mean, geom = "point", size = 1.5, alpha = 0.4) +
    # Linear trend lines for each category
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.5, alpha = 0.8) +
    scale_y_reverse() +
    scale_color_manual(values = plot_colors) +
    scale_x_continuous(breaks = 1:6) +
    theme_minimal(base_size = 17) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Presentation",
      y = "Mean voltage (\u03bcV)",
      colour = "Category"
    )
  
  return(p)
}

#' Save plot in multiple formats
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param output_dir Output directory
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot_multiple_formats <- function(plot, filename, output_dir, width = 10, height = 7) {
  
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Always save as PDF
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  cat("✓ Saved:", pdf_path, "\n")
  
  # Save as SVG if available
  if (svg_available) {
    tryCatch({
      svg_path <- file.path(output_dir, paste0(filename, ".svg"))
      ggsave(svg_path, plot = plot, width = width, height = height, 
             units = "in", device = "svg")
      cat("✓ Saved:", svg_path, "\n")
    }, error = function(e) {
      cat("⚠ SVG save failed for", filename, ":", e$message, "\n")
    })
  }
}

#' Generate all visualizations
#'
#' @param datasets List of cumulative datasets
#' @param pairwise_results Pairwise comparison results
#' @param params Experimental parameters
#' @param output_dir Output directory for plots
generate_all_plots <- function(datasets, pairwise_results, params, output_dir) {
  
  cat("==============================================================\n")
  cat("GENERATING VISUALIZATIONS\n")
  cat("==============================================================\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate ERP bar plots for key repetition counts
  key_datasets <- c("df_1_3", "df_1_6")
  
  for (name in key_datasets) {
    if (name %in% names(datasets) && name %in% names(pairwise_results)) {
      
      cat("Creating ERP plot for", name, "...\n")
      
      title_suffix <- case_when(
        name == "df_1_3" ~ "— Repetitions 1-3",
        name == "df_1_6" ~ "— Repetitions 1-6",
        TRUE ~ paste("—", name)
      )
      
      erp_plot <- create_erp_plot(
        df = datasets[[name]],
        pairs_obj = pairwise_results[[name]]$pairs,
        conditions = params$conditions,
        plot_colors = params$plot_colors,
        title_suffix = title_suffix
      )
      
      # Remove title for cleaner look
      erp_plot <- erp_plot + labs(title = NULL)
      
      filename <- paste0("erp_plot_", gsub("df_", "", name))
      save_plot_multiple_formats(erp_plot, filename, output_dir)
    }
  }
  
  # Create repetition trend plot
  if ("df_1_6" %in% names(datasets)) {
    cat("Creating repetition trend plot...\n")
    
    trend_plot <- create_repetition_trend_plot(datasets$df_1_6, params$plot_colors)
    save_plot_multiple_formats(trend_plot, "repetition_trend_plot", output_dir)
  }
  
  # Create interaction plot
  if ("df_1_6" %in% names(datasets)) {
    cat("Creating interaction plot...\n")
    
    interaction_plot <- create_interaction_plot(datasets$df_1_6, params$plot_colors)
    save_plot_multiple_formats(interaction_plot, "interaction_plot", output_dir, width = 12, height = 8)
  }
  
  # Create SNR plots (if data available)
  cat("Creating SNR analysis plots...\n")
  
  snr_data_path <- "data/3_analysis/cumulative_snr_data.csv"
  f_data_path <- "data/3_analysis/f_values_data.csv"
  
  snr_plots <- create_snr_plots(snr_data_path, f_data_path)
  
  if ("snr_plot" %in% names(snr_plots)) {
    save_plot_multiple_formats(snr_plots$snr_plot, "snr_plot", output_dir)
  }
  
  if ("f_plot" %in% names(snr_plots)) {
    save_plot_multiple_formats(snr_plots$f_plot, "f_plot", output_dir)
  }
  
  cat("\n✓ All visualizations created and saved to:", output_dir, "\n")
}

#' Main execution function
main <- function() {
  
  # Load or prepare required data
  if (!exists("results")) {
    cat("Loading data preparation results...\n")
    if (file.exists("01_data_preparation.R")) {
      source("01_data_preparation.R")
      results <- main()
    } else {
      stop("Data preparation script not found")
    }
  }
  
  if (!exists("analysis_results")) {
    cat("Loading analysis results...\n")
    if (file.exists("03_statistical_analysis.R")) {
      source("03_statistical_analysis.R")
      analysis_results <- main()
    } else {
      stop("Statistical analysis script not found")
    }
  }
  
  # Set output directory
  output_dir <- "data/4_plots"
  
  # Generate all plots
  generate_all_plots(
    datasets = results$datasets,
    pairwise_results = analysis_results$pairwise,
    params = results$params,
    output_dir = output_dir
  )
  
  cat("\n==============================================================\n")
  cat("DATA VISUALIZATION COMPLETE\n")
  cat("==============================================================\n")
  
  return(list(output_dir = output_dir))
}

# Execute main function if script is run directly
if (!interactive()) {
  plot_results <- main()
}

# Data Visualization and Plotting Script
#
# This script handles creation of publication-ready plots:
# - ERP bar plots with significance brackets for different repetition numbers
# - Line plots showing trends across repetitions
# - SNR analysis plots from Python data export
# - Interaction plots for repetition × category effects
# - Comprehensive plot themes and styling
#
# Author: [Your Name]

# Load additional plotting libraries
load_plotting_libraries <- function() {
  suppressPackageStartupMessages({
    if (!require("patchwork", quietly = TRUE)) {
      install.packages("patchwork")
      library(patchwork)
    }
  })
}

# Define plot themes
setup_plot_themes <- function() {
  cat("Setting up plot themes...\n")
  
  # Base theme
  base_theme <- theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # ERP-specific theme
  erp_theme <- base_theme +
    theme(
      aspect.ratio = 0.7,
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none"
    )
  
  # SNR plot theme
  snr_theme <- base_theme +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
    )
  
  # Line plot theme
  line_theme <- theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = "right"
    )
  
  return(list(
    base_theme = base_theme,
    erp_theme = erp_theme,
    snr_theme = snr_theme,
    line_theme = line_theme
  ))
}

# Core ERP plotting function with significance brackets
create_erp_plot <- function(df, pairwise_results, conditions, plot_colors, title_suffix, themes) {
  cat(sprintf("Creating ERP plot for %s...\n", title_suffix))
  
  # Summary stats with numeric x positions
  mean_se_data <- df %>%
    group_by(category) %>%
    summarise(
      mean_val = mean(mean_voltage),
      se = sd(mean_voltage) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(x = as.numeric(factor(category, levels = conditions)))
  
  # Significant contrasts
  pairs_summary <- summary(pairwise_results$pairs)
  sig_contrasts <- pairs_summary %>%
    as.data.frame() %>%
    filter(p.value < 0.05) %>%
    mutate(contrast_clean = str_replace_all(contrast, " - ", " vs ")) %>%
    separate(contrast, into = c("cat1", "cat2"), sep = " - ", remove = FALSE) %>%
    mutate(
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**", 
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  
  # Create bracket data for significance annotations
  bracket_data <- if (nrow(sig_contrasts) > 0) {
    base_y <- min(mean_se_data$mean_val - mean_se_data$se) - 0.2
    
    sig_contrasts %>%
      left_join(mean_se_data %>% select(category, mean_val, se, x),
                by = c("cat1" = "category")) %>%
      rename_with(~ paste0(.x, "_1"), c("mean_val", "se", "x")) %>%
      left_join(mean_se_data %>% select(category, mean_val, se, x),
                by = c("cat2" = "category")) %>%
      rename_with(~ paste0(.x, "_2"), c("mean_val", "se", "x")) %>%
      mutate(
        xmin_pos = x_1,
        xmax_pos = x_2,
        bracket_span = abs(xmax_pos - xmin_pos)
      ) %>%
      arrange(bracket_span, xmin_pos) %>%
      mutate(y_bracket = base_y - (row_number() - 1) * 0.25)
  } else {
    NULL
  }
  
  # Create the plot
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
    themes$erp_theme
  
  # Add significance brackets if any
  if (!is.null(bracket_data) && nrow(bracket_data) > 0) {
    p <- p +
      geom_segment(data = bracket_data,
                   aes(x = xmin_pos, xend = xmax_pos, y = y_bracket, yend = y_bracket),
                   inherit.aes = FALSE, color = "black") +
      geom_segment(data = bracket_data,
                   aes(x = xmin_pos, xend = xmin_pos, y = y_bracket, yend = y_bracket + 0.05),
                   inherit.aes = FALSE, color = "black") +
      geom_segment(data = bracket_data,
                   aes(x = xmax_pos, xend = xmax_pos, y = y_bracket, yend = y_bracket + 0.05),
                   inherit.aes = FALSE, color = "black") +
      geom_text(data = bracket_data,
                aes(x = (xmin_pos + xmax_pos)/2, y = y_bracket - 0.05, label = sig_label),
                inherit.aes = FALSE, size = 4, hjust = 0.5)
  }
  
  return(p)
}

# Create ERP plots for all repetition numbers
create_all_erp_plots <- function(datasets, pairwise_results, params, themes) {
  cat("Creating ERP plots for all repetition numbers...\n")
  
  plots <- list()
  
  # Define which datasets and pairwise results to plot
  plot_specs <- list(
    list(data = "df_1", pairs = "mdl1", title = "— Repetition 1"),
    list(data = "df_1_2", pairs = "mdl1_2", title = "— Repetitions 1-2"),
    list(data = "df_1_3", pairs = "mdl1_3", title = "— Repetitions 1-3"),
    list(data = "df_1_4", pairs = "mdl1_4", title = "— Repetitions 1-4"),
    list(data = "df_1_5", pairs = "mdl1_5", title = "— Repetitions 1-5"),
    list(data = "df_1_6", pairs = "mdl1_6", title = "— Repetitions 1-6")
  )
  
  for (spec in plot_specs) {
    if (spec$data %in% names(datasets) && spec$pairs %in% names(pairwise_results)) {
      plots[[spec$pairs]] <- create_erp_plot(
        df = datasets[[spec$data]],
        pairwise_results = pairwise_results[[spec$pairs]],
        conditions = params$conditions,
        plot_colors = params$plot_colors,
        title_suffix = spec$title,
        themes = themes
      )
    }
  }
  
  cat(sprintf("✓ Created %d ERP plots\n", length(plots)))
  return(plots)
}

# Create repetition trend line plot
create_repetition_trend_plot <- function(datasets, params, themes) {
  cat("Creating repetition trend plot...\n")
  
  # Use the full dataset (1-6 repetitions)
  df_full <- datasets$df_1_6
  
  trend_plot <- ggplot(df_full, aes(x = repetition, y = mean_voltage, 
                                   colour = category, group = category)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_y_reverse() +
    scale_color_manual(values = params$plot_colors) +
    scale_x_continuous(breaks = 1:6) +
    themes$line_theme +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    labs(
      title = NULL,
      x = "Presentation number",
      y = "Mean voltage (\u03bcV)",
      colour = "Category"
    )
  
  cat("✓ Repetition trend plot created\n")
  return(trend_plot)
}

# Create interaction plot with linear trend lines
create_interaction_plot <- function(datasets, params, themes) {
  cat("Creating interaction plot...\n")
  
  df_full <- datasets$df_1_6
  
  interaction_plot <- ggplot(df_full, aes(x = repetition, y = mean_voltage, colour = category)) +
    # Individual data points (desaturated)
    stat_summary(fun = mean, geom = "point", size = 1.5, alpha = 0.4) +
    # Linear trend lines for each category
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.5, alpha = 0.8) +
    scale_y_reverse() +
    scale_color_manual(values = params$plot_colors) +
    scale_x_continuous(breaks = 1:6) +
    themes$line_theme +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    labs(
      x = "Presentation number",
      y = "Mean voltage (\u03bcV)",
      colour = "Category"
    )
  
  cat("✓ Interaction plot created\n")
  return(interaction_plot)
}

# Load and plot SNR data from Python analysis
create_snr_plots <- function(project_root = NULL, themes) {
  cat("Creating SNR analysis plots...\n")
  
  # Default path
  if (is.null(project_root)) {
    project_root <- file.path("..", "..")
  }
  
  # Load SNR data
  snr_data_path <- file.path(project_root, "data", "3_analysis", "cumulative_snr_data.csv")
  f_data_path <- file.path(project_root, "data", "3_analysis", "f_values_data.csv")
  
  plots <- list()
  
  # Check if SNR data exists
  if (file.exists(snr_data_path)) {
    snr_data <- read_csv(snr_data_path, show_col_types = FALSE)
    
    plots$snr_plot <- ggplot(snr_data, aes(x = cumulative_repetitions, y = mean_snr)) +
      geom_line(color = "blue", size = 1.2) +
      geom_point(color = "blue", size = 3) +
      geom_errorbar(aes(ymin = mean_snr - sem_snr, ymax = mean_snr + sem_snr), 
                    width = 0.2, color = "blue") +
      geom_text(aes(label = sprintf("%.2f", mean_snr)), 
                vjust = -0.8, hjust = 0.5, size = 3.5) +
      scale_x_continuous(breaks = snr_data$cumulative_repetitions,
                         labels = c("1", "1-2", "1-3", "1-4", "1-5", "1-6")) +
      themes$snr_theme +
      labs(
        title = "Cumulative SNR Analysis",
        x = "Cumulative presentations",
        y = "Mean SNR"
      )
    
    cat("✓ SNR plot created\n")
  } else {
    cat("⚠️ SNR data not found, skipping SNR plot\n")
  }
  
  # Check if F-value data exists
  if (file.exists(f_data_path)) {
    f_data <- read_csv(f_data_path, show_col_types = FALSE)
    
    plots$f_plot <- ggplot(f_data, aes(x = cumulative_repetitions, y = f_statistic)) +
      geom_col(alpha = 0.7, fill = "red", width = 0.6) +
      geom_text(aes(label = sprintf("%.2f", f_statistic)), 
                vjust = -0.5, hjust = 0.5, size = 3.5) +
      scale_x_continuous(breaks = f_data$cumulative_repetitions,
                         labels = c("1-2", "1-3", "1-4", "1-5", "1-6")) +
      themes$base_theme +
      theme(panel.grid.major.x = element_blank()) +
      labs(
        title = "F-values from Paired t-tests",
        x = "Cumulative presentations", 
        y = "F-value"
      )
    
    cat("✓ F-value plot created\n")
  } else {
    cat("⚠️ F-value data not found, skipping F-value plot\n")
  }
  
  return(plots)
}

# Save plots to files
save_plots <- function(plots, output_dir = NULL, format = "pdf") {
  cat("Saving plots...\n")
  
  if (is.null(output_dir)) {
    output_dir <- file.path("..", "..", "data", "4_plots")
  }
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save each plot
  for (plot_name in names(plots)) {
    filename <- file.path(output_dir, paste0(plot_name, ".", format))
    
    tryCatch({
      if (format == "svg" && require("svglite", quietly = TRUE)) {
        ggsave(filename, plot = plots[[plot_name]], width = 10, height = 7, 
               units = "in", device = "svg")
      } else {
        ggsave(filename, plot = plots[[plot_name]], width = 10, height = 7, 
               units = "in", device = "pdf")
      }
      cat(sprintf("✓ Saved: %s\n", filename))
    }, error = function(e) {
      cat(sprintf("⚠️ Failed to save %s: %s\n", plot_name, e$message))
    })
  }
  
  cat(sprintf("✓ Plot saving completed. Files saved to: %s\n", output_dir))
}

# Create key plots for publication
create_key_plots <- function(datasets, pairwise_results, params, themes) {
  cat("Creating key plots for publication...\n")
  
  key_plots <- list()
  
  # Key ERP plots (1-3 and 1-6 repetitions)
  if ("df_1_3" %in% names(datasets) && "mdl1_3" %in% names(pairwise_results)) {
    key_plots$erp_plot_1_3 <- create_erp_plot(
      df = datasets$df_1_3,
      pairwise_results = pairwise_results$mdl1_3,
      conditions = params$conditions,
      plot_colors = params$plot_colors,
      title_suffix = "— Repetitions 1-3",
      themes = themes
    ) + labs(title = NULL)
  }
  
  if ("df_1_6" %in% names(datasets) && "mdl1_6" %in% names(pairwise_results)) {
    key_plots$erp_plot_1_6 <- create_erp_plot(
      df = datasets$df_1_6,
      pairwise_results = pairwise_results$mdl1_6,
      conditions = params$conditions,
      plot_colors = params$plot_colors,
      title_suffix = "— Repetitions 1-6",
      themes = themes
    ) + labs(title = NULL)
  }
  
  # Repetition trend plot
  key_plots$repetition_trend <- create_repetition_trend_plot(datasets, params, themes)
  
  # Interaction plot
  key_plots$interaction_plot <- create_interaction_plot(datasets, params, themes)
  
  cat(sprintf("✓ Created %d key plots\n", length(key_plots)))
  return(key_plots)
}

# Main plotting function
create_all_plots <- function(datasets, pairwise_results, params, project_root = NULL, save_plots_flag = TRUE) {
  cat("=========================================\n")
  cat("       DATA VISUALIZATION PLOTS         \n")
  cat("=========================================\n\n")
  
  # Load plotting libraries
  load_plotting_libraries()
  
  # Setup themes
  themes <- setup_plot_themes()
  
  # Create all ERP plots
  erp_plots <- create_all_erp_plots(datasets, pairwise_results, params, themes)
  
  # Create key plots
  key_plots <- create_key_plots(datasets, pairwise_results, params, themes)
  
  # Create SNR plots
  snr_plots <- create_snr_plots(project_root, themes)
  
  # Combine all plots
  all_plots <- c(erp_plots, key_plots, snr_plots)
  
  # Display key plots
  cat("\nDisplaying key plots:\n")
  if ("repetition_trend" %in% names(all_plots)) {
    print(all_plots$repetition_trend)
  }
  if ("erp_plot_1_3" %in% names(all_plots)) {
    print(all_plots$erp_plot_1_3)
  }
  if ("erp_plot_1_6" %in% names(all_plots)) {
    print(all_plots$erp_plot_1_6)
  }
  if ("snr_plot" %in% names(all_plots)) {
    print(all_plots$snr_plot)
  }
  if ("interaction_plot" %in% names(all_plots)) {
    print(all_plots$interaction_plot)
  }
  
  # Save plots if requested
  if (save_plots_flag) {
    output_dir <- if (!is.null(project_root)) {
      file.path(project_root, "data", "4_plots")
    } else {
      file.path("..", "..", "data", "4_plots") 
    }
    save_plots(all_plots, output_dir)
  }
  
  cat("\n✓ All plots created successfully\n")
  
  return(all_plots)
}

# If running this script directly
if (sys.nframe() == 0) {
  cat("Plotting script requires data from previous analysis scripts\n")
  cat("Please run the previous scripts first\n")
}

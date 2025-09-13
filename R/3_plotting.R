source("~/Documents/GitHub/MA_thesis/R/2_PairwiseComparisons.R")

category_levels <- c("animal", "food", "tool", "commun", "emotion", "social")
df_1_3$category <- factor(df_1_3$category, levels = category_levels)
df_1_6$category <- factor(df_1_6$category, levels = category_levels)


# Define output directory
out_dir <- "/Users/johannberger/Documents/GitHub/MA_thesis/plots/erp"

# Make sure the folder exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ========== CREATE PLOTS ==========
p_1_3 <- ggplot(df_1_3, aes(x = category, y = mean_voltage, fill = category)) +
  stat_summary(fun = mean, geom = "col", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = plot_colors) +
  scale_y_reverse() +
  labs(
    title = NULL,
    x = NULL,
    y = "Mean voltage (\u03bcV)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

p_1_6 <- ggplot(df_1_6, aes(x = category, y = mean_voltage, fill = category)) +
  stat_summary(fun = mean, geom = "col", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = plot_colors) +
  scale_y_reverse() +
  labs(
    title = NULL,
    x = NULL,
    y = "Mean voltage (\u03bcV)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

p_trend <- ggplot(df_1_6, aes(x = repetition, y = mean_voltage, colour = category, group = category)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_y_reverse() +
  scale_color_manual(values = plot_colors) +
  labs(
    title = NULL,
    x = "Repetition",
    y = "Mean voltage (\u03bcV)",
    colour = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.direction = "horizontal")


interaction_plot <- ggplot(df_1_6, aes(x = repetition, y = mean_voltage, colour = category)) +
  stat_summary(fun = mean, geom = "point", size = 1.5, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.5, alpha = 0.8) +
  scale_y_reverse() +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = 1:6) +
  theme_minimal() +   # <-- set theme first
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(
    x = "Presentation",
    y = expression("Mean voltage ("*mu*"V)"),
    colour = "Category"
  )

# Update labels once
p_1_3 <- p_1_3 +
  labs(y = expression("Mean voltage ("*mu*"V)"),
       title = NULL)

p_1_6 <- p_1_6 +
  labs(y = expression("Mean voltage ("*mu*"V)"),
       title = NULL)

p_trend <- p_trend +
  labs(y = expression("Mean voltage ("*mu*"V)"),
       title = NULL)

p_int <- interaction_plot +
  labs(y = expression("Mean voltage ("*mu*"V)"),
       title = NULL)

ggsave(file.path(out_dir, "barplot_1_3.pdf"), p_1_3,
       device = "pdf", width = 12, height = 8)
ggsave(file.path(out_dir, "barplot_1_6.pdf"), p_1_6,
       device = "pdf", width = 12, height = 8)
ggsave(file.path(out_dir, "trendplot.pdf"), p_trend,
       device = "pdf", width = 12, height = 8)
ggsave(file.path(out_dir, "interaction.pdf"), p_int,
       device = "pdf", width = 12, height = 8)

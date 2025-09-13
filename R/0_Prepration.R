# Loading in the data

# Libraries, source function and plotting parameters
library(glue)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(emmeans)
library(broom.mixed)
library(tidyverse)
library(svglite)

setwd("/Users/johannberger/Documents/GitHub/MA_thesis/data/2_preprocessed/285-345ms")
df <- read_csv("combined_task-nouns_285-345ms.csv")
electrode <- c("FC1", "FCz", "FC2", "FCC1h", "FCC2h", "C1", "Cz", "C2", "CCP1h", "CCP2h", "CP1", "CPz", "CP2", "CPP1h", "CPP2h") # nolint: line_length_linter.
conditions <- c("animal", "food", "tool", "commun", "emotion","social") # nolint

plot_colors = c(
  "animal" = "#388E3C",
  "commun" = "#d9a00f",
  "emotion" = "#616161",
  "food" = "#1976D2",
  "social" = "#ff00ff",
  "tool" = "#D32F2F"
)

# Define data frames
df_1 <- df %>%
  filter(category %in% conditions,
         repetition == 1,
         channel %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    category = factor(category, levels = conditions)
  ) %>%
  { contrasts(.$category) <- contr.sum(length(conditions)); . }


df_1_2 <- df %>%
  filter(category   %in% conditions,
         repetition %in% 1:2,
         channel    %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    repetition_c = scale(repetition, scale = FALSE)[, 1],
    category     = structure(
      factor(category, levels = conditions),
      contrasts = contr.sum(length(conditions))
    )
  )

df_1_3 <- df %>%
  filter(category   %in% conditions,
         repetition %in% 1:3,
         channel    %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    repetition_c = scale(repetition, scale = FALSE)[, 1],
    category     = structure(
      factor(category, levels = conditions),
      contrasts = contr.sum(length(conditions))
    )
  )

df_1_4 <- df %>%
  filter(category   %in% conditions,
         repetition %in% 1:4,
         channel    %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    repetition_c = scale(repetition, scale = FALSE)[, 1],
    category     = structure(
      factor(category, levels = conditions),
      contrasts = contr.sum(length(conditions))
    )
  )

df_1_5 <- df %>%
  filter(category   %in% conditions,
         repetition %in% 1:5,
         channel    %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    repetition_c = scale(repetition, scale = FALSE)[, 1],
    category     = structure(
      factor(category, levels = conditions),
      contrasts = contr.sum(length(conditions))
    )
  )

df_1_6 <- df %>%
  filter(category   %in% conditions,
         repetition %in% 1:6,
         channel    %in% electrode) %>%
  group_by(item, subject, category, repetition) %>%
  summarise(mean_voltage = mean(voltage), .groups = "drop") %>%
  mutate(
    repetition_c = scale(repetition, scale = FALSE)[, 1],
    category     = structure(
      factor(category, levels = conditions),
      contrasts = contr.sum(length(conditions))
    )
  )

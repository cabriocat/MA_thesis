source("~/Documents/GitHub/MA_thesis/R/0_Prepration.R")

mdl1_df1 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item), # This Model is used
                 data = df_1)


# REPETITION 1 + 2
mdl1_df1_2 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item), # This Model is used
                   data = df_1_2)
mdl2_df1_2 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item),
                   data = df_1_2)
mdl3_df1_2 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                   data = df_1_2)
mdl4_df1_2 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item),
                   data = df_1_2)
anova(mdl1_df1_2, mdl2_df1_2, mdl3_df1_2, mdl4_df1_2)


# REPETITION 1 + 2 + 3
mdl1_df1_3 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item),
                   data = df_1_3)
mdl2_df1_3 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item), # This Model is used
                   data = df_1_3)
mdl3_df1_3 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                   data = df_1_3)
mdl4_df1_3 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item),
                   data = df_1_3)
anova(mdl1_df1_3, mdl2_df1_3, mdl3_df1_3, mdl4_df1_3)


# REPETITION 1 + 2 + 3 + 4
mdl1_df1_4 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item),
                   data = df_1_4)
mdl2_df1_4 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item),
                   data = df_1_4)
mdl3_df1_4 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                   data = df_1_4)
mdl4_df1_4 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item), # This Model is used
                   data = df_1_4)
anova(mdl1_df1_4, mdl2_df1_4, mdl3_df1_4, mdl4_df1_4)


# REPETITION 1 + 2 + 3 + 4 + 5
mdl1_df1_5 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item),
                   data = df_1_5)
mdl2_df1_5 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item),
                   data = df_1_5)
mdl3_df1_5 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                   data = df_1_5)
mdl4_df1_5 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item), # This Model is used
                   data = df_1_5)
anova(mdl1_df1_5, mdl2_df1_5, mdl3_df1_5, mdl4_df1_5)


# REPETITION 1 + 2 + 3 + 4 + 5 + 6
mdl1_df1_6 <- lmer(mean_voltage ~ category + (1 | subject) + (1 | item),
                   data = df_1_6)
mdl2_df1_6 <- lmer(mean_voltage ~ category + repetition_c + (1 | subject) + (1 | item),
                   data = df_1_6)
mdl3_df1_6 <- lmer(mean_voltage ~ category * repetition_c + (1 | subject) + (1 | item),
                   data = df_1_6)
mdl4_df1_6 <- lmer(mean_voltage ~ category * repetition_c + (1 + repetition_c || subject) + (1 | item), # This Model is used
                   data = df_1_6)
anova(mdl1_df1_6, mdl2_df1_6, mdl3_df1_6, mdl4_df1_6)


# Model Selection
# Main Models
mdl1 <- mdl1_df1
mdl1_2 <- mdl1_df1_2 # Best model for 1-2 repetitions
mdl1_3 <- mdl2_df1_3 # Best model for 1-3 repetitions
mdl1_4 <- mdl4_df1_4 # Best model for 1-4 repetitions
mdl1_5 <- mdl4_df1_5 # Best model for 1-5 repetitions
mdl1_6 <- mdl4_df1_6 # Best model for 1-6 repetitions

# ANOVA for each model
anova(mdl1, type = 3)
anova(mdl1_2, type = 3)
anova(mdl1_3, type = 3)
anova(mdl1_4, type = 3)
anova(mdl1_5, type = 3)
anova(mdl1_6, type = 3)
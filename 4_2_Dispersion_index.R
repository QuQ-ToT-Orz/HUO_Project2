#### 6 Dispersion Index Analysis ####
# drop + inclusive + covered-segment-bin normalization + event-time mu + 
# bin_size=120 + inclusive lambda + adj_marks uses mu_window #
gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
source("function_dispersion.R")

#### 1 Load Data ####
load(file = paste("./data/runlength/", "event_analysis_old.RData", sep = ""))
load(file = paste("./data/count/", "data_analysis_old.RData", sep = ""))

load(file = paste("./data/runlength/", "event_analysis_new.RData", sep = ""))
load(file = paste("./data/mims/", "data_analysis_new.RData", sep = ""))

circadian_bin_size <- 120
window_sizes <- c(15, 30, 60, 90, 120)

#### 2 Compute for All Participants ####

# Filter active events
active_events <- event_analysis %>%
  filter(categories != "sedentary") %>%
  group_by(SEQN, WEEKDAY) %>%
  mutate(
    # sqrt for hip
    mark_sqrt = (original_value),
    daily_mean_sqrt = mean(mark_sqrt),
    mark_norm = mark_sqrt / daily_mean_sqrt
  ) %>%
  ungroup()

daily_subject_metrics <- active_events %>%
  group_by(SEQN, WEEKDAY) %>%
  summarise(
    Q_day = mean(mark_norm^2) / mean(mark_norm),
    n_events = n(),
    sum_marks = sum(mark_norm),
    obs_period = max(start) - min(start) + 1,
    .groups = "drop"
  ) %>%
  mutate(
    obs_hours = obs_period / 60,
    lambda_count_day = n_events / obs_hours,
    lambda_marks_day = sum_marks / obs_hours
  )

person_metrics <- daily_subject_metrics %>%
  group_by(SEQN) %>%
  summarise(
    Q = mean(Q_day, na.rm = TRUE),
    lambda_count = mean(lambda_count_day, na.rm = TRUE),
    lambda_marks = mean(lambda_marks_day, na.rm = TRUE),
    .groups = "drop"
  )

person_events_list <- split(active_events, active_events$SEQN, drop = TRUE)
cat("n_participants:", length(person_events_list), "\n")

dispersion_results <- lapply(person_events_list, function(person_events) {
  person_disp <- compute_dispersion(split(person_events, person_events$WEEKDAY, drop = TRUE), window_sizes, circadian_bin_size)
  person_disp$SEQN <- person_events$SEQN[[1]]
  person_disp
})

dispersion_df <- bind_rows(dispersion_results) %>%
  left_join(select(person_metrics, SEQN, Q), by = "SEQN") %>%
  mutate(
    n_raw = branching_ratio(D_raw),
    n_adj = branching_ratio(D_adj),
    n_marks = 1 - sqrt(Q / D_marks),
    n_adj_marks = 1 - sqrt(Q / D_adj_marks)
  )

#### 3 Window Size Sensitivity ####
window_sensitivity <- dispersion_df %>%
  group_by(window_size) %>%
  summarise(
    mean_D_raw = mean(D_raw, na.rm = TRUE),
    sd_D_raw = sd(D_raw, na.rm = TRUE),
    mean_D_marks = mean(D_marks, na.rm = TRUE),
    sd_D_marks = sd(D_marks, na.rm = TRUE),
    mean_D_adj = mean(D_adj, na.rm = TRUE),
    sd_D_adj = sd(D_adj, na.rm = TRUE),
    mean_D_adj_marks = mean(D_adj_marks, na.rm = TRUE),
    sd_D_adj_marks = sd(D_adj_marks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff_count = mean_D_raw - mean_D_adj,
    diff_marks = mean_D_marks - mean_D_adj_marks
  )

print(window_sensitivity)

p_windows <- ggplot(window_sensitivity, aes(x = window_size)) +
  # D_raw with error ribbon
  geom_ribbon(aes(ymin = mean_D_raw - sd_D_raw, ymax = mean_D_raw + sd_D_raw),
              fill = "red", alpha = 0.2) +
  geom_line(aes(y = mean_D_raw, color = "D_raw"), linewidth = 1) +
  geom_point(aes(y = mean_D_raw, color = "D_raw"), size = 3) +
  # D_marks with error ribbon
  geom_ribbon(aes(ymin = mean_D_marks - sd_D_marks, ymax = mean_D_marks + sd_D_marks),
              fill = "orange", alpha = 0.2) +
  geom_line(aes(y = mean_D_marks, color = "D_marks"), linewidth = 1) +
  geom_point(aes(y = mean_D_marks, color = "D_marks"), size = 3) +
  # D_adj with error ribbon
  geom_ribbon(aes(ymin = mean_D_adj - sd_D_adj, ymax = mean_D_adj + sd_D_adj),
              fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean_D_adj, color = "D_adj"), linewidth = 1) +
  geom_point(aes(y = mean_D_adj, color = "D_adj"), size = 3) +
  # D_adj_marks with error ribbon
  geom_ribbon(aes(ymin = mean_D_adj_marks - sd_D_adj_marks, ymax = mean_D_adj_marks + sd_D_adj_marks),
              fill = "green", alpha = 0.2) +
  geom_line(aes(y = mean_D_adj_marks, color = "D_adj_marks"), linewidth = 1) +
  geom_point(aes(y = mean_D_adj_marks, color = "D_adj_marks"), size = 3) +
  scale_color_manual(values = c("D_raw" = "red", "D_marks" = "orange", "D_adj" = "blue", "D_adj_marks" = "darkgreen")) +
  labs(x = "Window Size (min)", y = "Dispersion Index", color = "") +
  theme_minimal()
print(p_windows)
ggsave("Output/dispersion/hip_windows.pdf", p_windows, width = 10, height = 6)
ggsave("Output/dispersion/wrist_windows.pdf", p_windows, width = 10, height = 6)

#### 4 Create Summary and Save ####
primary_window <- 30

dispersion_summary <- dispersion_df %>%
  filter(window_size == primary_window) %>%
  left_join(select(person_metrics, SEQN, lambda_count, lambda_marks), by = "SEQN") %>%
  mutate(
    # Immigration rates using Hawkes relationship: μ* = λ(1-n)
    # 4 versions to match 4 versions of n
    mu_star_raw = lambda_count * (1 - n_raw),           # count, no circadian adj
    mu_star_adj = lambda_count * (1 - n_adj),           # count, circadian adj
    mu_star_marks = lambda_marks * (1 - n_marks),       # marks, no circadian adj
    mu_star_adj_marks = lambda_marks * (1 - n_adj_marks) # marks, circadian adj
  ) %>%
  select(SEQN, Q, D_raw, D_marks, D_adj, D_adj_marks,
         n_raw, n_marks, n_adj, n_adj_marks,
         lambda_count, lambda_marks,
         mu_star_raw, mu_star_adj, mu_star_marks, mu_star_adj_marks)

# Save results
dir_path <- "./data/dispersion/"
save(dispersion_df, file = paste0(dir_path, "dispersion_df_old.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_old.RData"))

save(dispersion_df, file = paste0(dir_path, "dispersion_df_new.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_new.RData"))

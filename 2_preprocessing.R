#### 0 Packages ####
gc()
rm(list = ls())
R.version

library(dplyr)
library(tidyr)
library(ggplot2)
# remotes::install_github('coolbutuseless/purler')
library(purler)
library(tibble)
library(stringr)
library(purrr)
source("./function_runlength.R")
source("./function_plots.R")

process_daily_data <- function(data, start_time = 0, end_time = 1440) {
  # Step 1: Apply gap processing - keep longest continuous period
  data_after_gaps <- data %>%
    group_by(SEQN, WEEKDAY) %>%
    do({
      filtered_data <- .
      times <- filtered_data$start

      # Find gap positions (> 60 minutes)
      gap_positions <- which(diff(times) > 60)

      if (length(gap_positions) > 0) {
        # Define segment boundaries (start and end indices)
        segment_starts <- c(1, gap_positions + 1)
        segment_ends <- c(gap_positions, nrow(filtered_data))

        # Calculate length of each segment
        segment_lengths <- segment_ends - segment_starts + 1

        # Find the longest segment
        longest_segment_idx <- which.max(segment_lengths)
        start_idx <- segment_starts[longest_segment_idx]
        end_idx <- segment_ends[longest_segment_idx]

        # Keep only the longest continuous period
        filtered_data <- filtered_data[start_idx:end_idx, ]
      }
      # If no gaps, keep all data

      filtered_data
    }) %>%
    ungroup()

  # Step 2: Filter data within time window and ensure category balance per day
  daily_balanced <- data_after_gaps %>%
    filter(
      start >= start_time,
      start <= end_time
    ) %>%
    group_by(SEQN, WEEKDAY) %>%
    summarise(
      sedentary_count = sum(categories == "sedentary"),
      other_count = sum(categories != "sedentary"),
      obs_hours = (max(start) - min(start)) / 60,
      sedentary_rate = sedentary_count / obs_hours,
      active_rate = other_count / obs_hours,
      .groups = "drop"
    ) %>%
    # Filter each day by event rates per hour
    # min 1 event/hour, max 50 events/hour for active
    filter(
      sedentary_rate >= 1,
      active_rate >= 1,
      active_rate <= 50
    )

  # Now calculate weekly summary from balanced days only
  weekly_summary <- daily_balanced %>%
    group_by(SEQN) %>%
    summarise(
      total_days = n_distinct(WEEKDAY),
      sedentary_count = sum(sedentary_count),
      other_count = sum(other_count),
      .groups = "drop"
    ) %>%
    # Filter to ensure subject has all 7 days after day-level filtering
    filter(total_days == 7)


  # Step 3: Filter the original data to include only valid subject-days
  valid_subjects <- weekly_summary$SEQN

  daily_data <- data_after_gaps %>%
    filter(SEQN %in% valid_subjects) %>%
    filter(
      start >= start_time,
      start <= end_time
    ) %>%
    # normalize time
    group_by(SEQN, WEEKDAY) %>%
    mutate(
      start_normalized = start - min(start), # First event of each day becomes t=0
      time_since_start = start - start_time
    ) %>%
    ungroup()

  # Step 4: Calculate activity timing tags
  activity_timing_summary <- daily_data %>%
    group_by(SEQN, WEEKDAY) %>%
    summarise(
      total_events = n(),
      starts_before_7am = sum(start < 480), # Before 7am = 420 minutes
      starts_after_11pm = sum(start >= 1320), # After 11pm = 1380 minutes
      .groups = "drop"
    ) %>%
    mutate(
      early_activer = starts_before_7am > (total_events * 0.1),
      late_activer = starts_after_11pm > (total_events * 0.1)
    )

  # Step 5: Calculate sleep regularity during core sleep time (1am-5am)
  sleep_metrics <- daily_data %>%
    group_by(SEQN, WEEKDAY) %>%
    summarise(
      sleep_events = sum(start >= 60 & start < 300),
      h1am_2am = sum(start >= 60 & start < 120),
      h2am_3am = sum(start >= 120 & start < 180),
      h3am_4am = sum(start >= 180 & start < 240),
      h4am_5am = sum(start >= 240 & start < 300),
      .groups = "drop"
    ) %>%
    mutate(
      hourly_window_mean = (h1am_2am + h2am_3am + h3am_4am + h4am_5am) / 4,
      hourly_window_sd = sqrt(((h1am_2am - hourly_window_mean)^2 +
        (h2am_3am - hourly_window_mean)^2 +
        (h3am_4am - hourly_window_mean)^2 +
        (h4am_5am - hourly_window_mean)^2) / 4),
      hourly_window_cv = ifelse(hourly_window_mean == 0, 0, hourly_window_sd / hourly_window_mean),
      night_shifter_status = ifelse(sleep_events >= quantile(sleep_events, 0.9, na.rm = TRUE), "Night_Shifter", "Normal_Circadian"),
      sleep_regularity = case_when(
        hourly_window_mean == 0 ~ "Consolidated",
        hourly_window_cv <= 0.3 ~ "Fragmented",
        hourly_window_cv <= 0.8 ~ "Somewhat Consolidated",
        TRUE ~ "Consolidated"
      )
    )

  # Step 6: Create subject-level tags based on subject-weekday calculations
  subject_tags <- activity_timing_summary %>%
    group_by(SEQN) %>%
    summarise(
      days_early_activer = sum(early_activer),
      days_late_activer = sum(late_activer),
      .groups = "drop"
    ) %>%
    left_join(
      sleep_metrics %>%
        group_by(SEQN) %>%
        summarise(days_consolidated_sleep = sum(sleep_regularity == "Consolidated"), .groups = "drop"),
      by = "SEQN"
    )

  # Merge tags with daily_data
  daily_data <- daily_data %>%
    left_join(subject_tags, by = "SEQN")

  return(daily_data)
}
#### 1 Load the data and functions ####
load("../2025/data/count/Act_Analysis_old.RData")
load("../2025/data/count/Flags_Analysis_old.RData")

load("../2025/data/mims/Act_Analysis_new.RData")
load("../2025/data/mims/Flags_Analysis_new.RData")

set.seed(123)
sample_id <- Act_Analysis %>%
  group_by(Gender, Race, MobilityProblem) %>%
  summarise(ID = list(unique(SEQN)), .groups = "keep") %>%
  mutate(sampled = lapply(ID, function(x) sample(x, size = round(length(x) * 0.15))))
sample_id <- sample_id %>%
  unnest(cols = c(sampled))
weekdays <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
id <- c(head(sample_id, n = 1)$sampled, tail(sample_id, n = 1)$sampled)
id

#### 2 Run-length encoding ####
# "hip"
# 1 - Simple minute-level events - activity intepretation
events <- simple_events(Act_Analysis, activity_cutoffs = "hip")
# filtered_events <- filter_short_bouts(events, min_bout_length = 2)
# event_analysis <- process_daily_data(filtered_events, start_time = 0, end_time = 1440)
filtered_events <- process_daily_data(events, start_time = 0, end_time = 1440)
n_after_daily_process_hip <- length(unique(filtered_events$SEQN))  # Capture for flowchart
event_analysis <- clean_events_spikes(filtered_events, activity_cutoffs = "hip")
# 2 - Generate RLE data with manual binning
# 2.1 - runlength_single_manual - transition intepretation
rle_data <- runlength_single_manual(Act_Analysis, activity_cutoffs = "hip")
result_rle <- semi_join(rle_data, events, by = c("SEQN", "WEEKDAY", "start"))
# categorize_activity_bouts(rle_data, activity_cutoffs = "hip")
# # 2.2 - runlength_single → all thresholds → filter to max → 1 event/minute - activity intepretation
# single_rle <- runlength_single_manual(Act_Analysis[1:100,], activity_cutoffs = "hip")
# combined_df(single_rle)
# 3 - Generate RLE data with quantile binning
runlength_single(Act_Analysis[1:100, ],
  n_bins = 10,
  bin_method = "quantile",
  global_bins = FALSE
)
runlength(Act_Analysis[1:100, ],
  n_bins = 10,
  bin_method = "quantile",
  global_bins = FALSE
)
# 4 - runlength_single (TRUE/FALSE) → all thresholds → filter to max (or original values) - multivariate_marked modelling

# save
dir_path <- "../2025/data/runlength/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

save(filtered_events, file = paste(dir_path, "filtered_events_old.RData", sep = ""))
save(event_analysis, file = paste(dir_path, "event_analysis_old.RData", sep = ""))

# "wrist"
# 1 - Simple minute-level events - activity intepretation
events <- simple_events(Act_Analysis, activity_cutoffs = "wrist")
# filtered_events <- filter_short_bouts(events, min_bout_length = 2)
# event_analysis <- process_daily_data(filtered_events, start_time = 0, end_time = 1440)
filtered_events <- process_daily_data(events, start_time = 0, end_time = 1440)
n_after_daily_process_wrist <- length(unique(filtered_events$SEQN))  # Capture for flowchart
event_analysis <- clean_events_spikes(filtered_events, activity_cutoffs = "wrist")
# 2 - Generate RLE data with manual binning
# 2.1 - runlength_single_manual - transition intepretation
rle_data <- runlength_single_manual(Act_Analysis, activity_cutoffs = "wrist")
result_rle <- semi_join(rle_data, filtered_events, by = c("SEQN", "WEEKDAY", "start"))
# categorize_activity_bouts(rle_data, activity_cutoffs = "wrist")
# # 2.2 - runlength_single → all thresholds → filter to max → 1 event/minute - activity intepretation
# single_rle <- runlength_single_manual(Act_Analysis, activity_cutoffs = "wrist")
# combined_df(single_rle)
# 3 - Generate RLE data with quantile binning
runlength_single(Act_Analysis[1:100, ],
  n_bins = 10,
  bin_method = "quantile",
  global_bins = FALSE
)
runlength(Act_Analysis[1:100, ],
  n_bins = 10,
  bin_method = "quantile",
  global_bins = FALSE
)
# 4 - runlength_single (TRUE/FALSE) → all thresholds → filter to max (or original values) - multivariate_marked modelling

# save
dir_path <- "../2025/data/runlength/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
save(filtered_events, file = paste(dir_path, "filtered_events_new.RData", sep = ""))
save(event_analysis, file = paste(dir_path, "event_analysis_new.RData", sep = ""))

#### 3 run length encoding (rle) ####
# Reconstruct Act_Analysis format from event_analysis
reconstruct_act_analysis <- function(event_data, original_act_data) {
  # Get unique subject-day combinations from event_data
  subject_days <- event_data %>%
    select(SEQN, WEEKDAY) %>%
    distinct()

  # Start with original Act_Analysis structure for these subject-days
  result <- original_act_data %>%
    inner_join(subject_days, by = c("SEQN", "WEEKDAY"))

  # Zero out all MIN columns first
  min_cols <- paste0("MIN", 1:1440)
  result[min_cols] <- 0

  # Convert event_data to wide format matching Act_Analysis structure
  event_wide <- event_data %>%
    mutate(min_col = paste0("MIN", start)) %>%
    select(SEQN, WEEKDAY, min_col, original_value) %>%
    pivot_wider(names_from = min_col, values_from = original_value, values_fill = 0)

  # Merge the activity values back into result
  for (col in intersect(names(event_wide), min_cols)) {
    result[[col]] <- result[[col]] + event_wide[[col]][match(
      paste(result$SEQN, result$WEEKDAY),
      paste(event_wide$SEQN, event_wide$WEEKDAY)
    )]
  }

  result
}


filter_rle_data <- function(data, min_events_per_category = 5) {
  # Step 1: Filter days to ensure both categories have sufficient events per day
  valid_days <- data %>%
    group_by(SEQN, WEEKDAY) %>%
    summarise(
      sedentary_count = sum(categories == "sedentary"),
      other_count = sum(categories != "sedentary"),
      .groups = "drop"
    ) %>%
    filter(
      sedentary_count >= min_events_per_category,
      other_count >= min_events_per_category
    )

  # Step 2: Keep only subjects who have all 7 days of valid data
  subjects_with_7_days <- valid_days %>%
    group_by(SEQN) %>%
    summarise(n_days = n(), .groups = "drop") %>%
    filter(n_days == 7) %>%
    pull(SEQN)

  # Step 3: Return the original data filtered for these subjects only
  data %>%
    filter(SEQN %in% subjects_with_7_days)
}

# Reconstruct filtered Act_Analysis from quality-controlled events
filtered_act_analysis <- reconstruct_act_analysis(event_analysis, Act_Analysis)
# Apply runlength_single_manual to filtered data to get pure alternating structure
result_rle <- runlength_single_manual(filtered_act_analysis, activity_cutoffs = "hip")
# Apply additional filtering to ensure sufficient events per category per day
rle_analysis <- filter_rle_data(result_rle, min_events_per_category = 5)
# Capture counts BEFORE intersect (for hip)
n_after_spike_clean_hip <- length(unique(event_analysis$SEQN))
n_after_rle_filter_hip <- length(unique(rle_analysis$SEQN))
# Align event_analysis and rle_analysis to have same subjects
common_subjects <- intersect(unique(event_analysis$SEQN), unique(rle_analysis$SEQN))
event_analysis <- event_analysis %>% filter(SEQN %in% common_subjects)
rle_analysis <- rle_analysis %>% filter(SEQN %in% common_subjects)
n_final_hip <- length(common_subjects)

result_rle <- runlength_single_manual(filtered_act_analysis, activity_cutoffs = "wrist")
# Apply additional filtering to ensure sufficient events per category per day
rle_analysis <- filter_rle_data(result_rle, min_events_per_category = 5)
# Capture counts BEFORE intersect (for wrist)
n_after_spike_clean_wrist <- length(unique(event_analysis$SEQN))
n_after_rle_filter_wrist <- length(unique(rle_analysis$SEQN))
# Align event_analysis and rle_analysis to have same subjects
common_subjects <- intersect(unique(event_analysis$SEQN), unique(rle_analysis$SEQN))
event_analysis <- event_analysis %>% filter(SEQN %in% common_subjects)
rle_analysis <- rle_analysis %>% filter(SEQN %in% common_subjects)
n_final_wrist <- length(common_subjects)

# Check subject counts
cat("Original subjects:", length(unique(Act_Analysis$SEQN)), "\n")
cat("Valid subjects:", length(unique(event_analysis$SEQN)), "\n")
cat("Valid subjects:", length(unique(rle_analysis$SEQN)), "\n")

flowchart_counts_step2 <- list(
  # Step 2 pipeline counts - common
  n_step1_subjects = length(unique(Act_Analysis$SEQN)),

  # Hip counts (captured before intersect)
  n_after_daily_process_hip = n_after_daily_process_hip,
  n_after_spike_clean_hip = n_after_spike_clean_hip,
  n_after_rle_filter_hip = n_after_rle_filter_hip,
  n_final_hip = n_final_hip
)

flowchart_counts_step2 <- list(
  # Step 2 pipeline counts - common
  n_step1_subjects = length(unique(Act_Analysis$SEQN)),

  # Wrist counts (captured before intersect)
  n_after_daily_process_wrist = n_after_daily_process_wrist,
  n_after_spike_clean_wrist = n_after_spike_clean_wrist,
  n_after_rle_filter_wrist = n_after_rle_filter_wrist,
  n_final_wrist = n_final_wrist
)

# Calculate exclusions at each step for HIP
flowchart_counts_step2$excluded_daily_hip <- flowchart_counts_step2$n_step1_subjects - flowchart_counts_step2$n_after_daily_process_hip
flowchart_counts_step2$excluded_spikes_hip <- flowchart_counts_step2$n_after_daily_process_hip - flowchart_counts_step2$n_after_spike_clean_hip
flowchart_counts_step2$excluded_rle_hip <- flowchart_counts_step2$n_after_spike_clean_hip - flowchart_counts_step2$n_after_rle_filter_hip
flowchart_counts_step2$excluded_align_hip <- flowchart_counts_step2$n_after_rle_filter_hip - flowchart_counts_step2$n_final_hip

# Calculate exclusions at each step for WRIST
flowchart_counts_step2$excluded_daily_wrist <- flowchart_counts_step2$n_step1_subjects - flowchart_counts_step2$n_after_daily_process_wrist
flowchart_counts_step2$excluded_spikes_wrist <- flowchart_counts_step2$n_after_daily_process_wrist - flowchart_counts_step2$n_after_spike_clean_wrist
flowchart_counts_step2$excluded_rle_wrist <- flowchart_counts_step2$n_after_spike_clean_wrist - flowchart_counts_step2$n_after_rle_filter_wrist
flowchart_counts_step2$excluded_align_wrist <- flowchart_counts_step2$n_after_rle_filter_wrist - flowchart_counts_step2$n_final_wrist

print(flowchart_counts_step2) 

# save
save(rle_analysis, file = paste(dir_path, "rle_analysis_old.RData", sep = ""))
save(event_analysis, file = paste(dir_path, "event_analysis_old.RData", sep = ""))

save(rle_analysis, file = paste(dir_path, "rle_analysis_new.RData", sep = ""))
save(event_analysis, file = paste(dir_path, "event_analysis_new.RData", sep = ""))

#### 4 data viz with examples ####
load("../2025/data/count/Act_Analysis_old.RData")
load("../2025/data/count/Flags_Analysis_old.RData")
load("../2025/data/runlength/rle_analysis_old.RData")
load("../2025/data/runlength/event_analysis_old.RData")

load("../2025/data/mims/Act_Analysis_new.RData")
load("../2025/data/mims/Flags_Analysis_new.RData")
load("../2025/data/runlength/rle_analysis_new.RData")
load("../2025/data/runlength/event_analysis_new.RData")
result_rle <- rle_analysis
seqn <- levels(Act_Analysis$SEQN)[1]
seqn = 40583 # hip
seqn = 64161 # wrist
weekday <- levels(Act_Analysis$WEEKDAY)[2]

#### 4_1 data viz - stationary ####
vec <- as.numeric(Act_Analysis[Act_Analysis$SEQN == seqn & Act_Analysis$WEEKDAY == weekday, paste0("MIN", 1:1440)]) # 10am to 20pm

par(mfrow = c(3, 1))
plot(ts(vec))
# obtain acf and pacf below
acf_result <- acf(vec)
pacf_result <- pacf(vec)
par(mfrow = c(1, 1))

par(mfrow = c(3, 1))
detrended_diff <- diff(vec)
plot(ts(detrended_diff))
# obtain acf and pacf below
acf_result <- acf(detrended_diff)
pacf_result <- pacf(detrended_diff)
par(mfrow = c(1, 1))
# ACF cutting off sharply
# PACF tailing off slowly

#### 4_2 data viz - sleep/wear/non-wear ####
plot_activity_and_runlength(Act_Analysis, seqn, weekday, result_rle)
plot_events(Act_Analysis, seqn, weekday, result_rle)
# Function to plot 24-hour data with waking periods highlighted
plot_subject_flag(seqn, weekday, "hip")
plot_subject_flag(seqn, weekday, "wrist")

plot_daily_or_weekly(Act_Analysis,
  id = seqn, data_type = "Raw", weekday = weekday,
  dataset = "daily_long", period = "daily"
) +
  labs(y = "acceleration")
plot_daily_or_weekly(Act_Analysis,
  id = seqn, data_type = "Raw",
  dataset = "daily_long", period = "daily"
) +
  labs(y = "acceleration")

# ?inverse.rle()

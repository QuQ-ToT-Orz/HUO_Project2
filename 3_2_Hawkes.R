#### doParallel ####
#### 0 Packages ####
gc()
rm(list = ls())
R.version

library(purrr)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(stringr)
library(HUO)
source("./function_process.R")

fit_list_hawkes_marked <- function(runs_df, n_cores, penalty_coef, single_day = FALSE, data_type) {
  cat(sprintf("Starting marked Hawkes process fitting\n"))

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)

  parallel::clusterExport(cl, c(
    "process_run_marked", "check_conv", "transform_hawkes_data",
    "fit_hawkes", "fit_hawkes_multi_series"
  ))
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(HUO)
  })

  if (single_day) {
    combinations <- runs_df %>%
      ungroup() %>%
      dplyr::select(WEEKDAY, SEQN) %>%
      unique() %>%
      mutate(
        WEEKDAY = as.character(WEEKDAY),
        SEQN = as.character(SEQN)
      )
  } else {
    # use multi-series fitting
    combinations <- runs_df %>%
      ungroup() %>%
      dplyr::select(SEQN) %>%
      unique() %>%
      mutate(
        SEQN = as.character(SEQN)
      )
  }

  total_combinations <- nrow(combinations)
  fits_list <- list()

  for (i in 1:total_combinations) {
    if (single_day) {
      weekday <- combinations$WEEKDAY[i]
      seqn <- combinations$SEQN[i]
      runs_df_single <- runs_df %>%
        filter(SEQN == seqn, WEEKDAY == weekday)
      fit_name <- paste0("fit_", seqn, "_weekday_", weekday)

      fit <- process_run_marked(
        seqn = seqn,
        weekday = weekday,
        runs_df = runs_df_single,
        penalty_coef = penalty_coef,
        single_day = TRUE,
        data_type = data_type
      )

      fits_list[[fit_name]] <- fit
    } else {
      # Multi-series fitting for all 7 days together
      seqn <- combinations$SEQN[i]
      runs_df_single <- runs_df %>%
        filter(SEQN == seqn)

      fit <- process_run_marked(
        seqn = seqn,
        runs_df = runs_df_single,
        penalty_coef = penalty_coef,
        single_day = FALSE,
        data_type = data_type
      )

      fits_list[[paste0("fit_", seqn, "_all_days")]] <- fit
    }

    if (i %% max(ceiling(total_combinations / 10), 1) == 0) {
      cat(sprintf("Progress: %d/%d combinations\n", i, total_combinations))
    }
  }

  cat("Processing completed.\n")
  return(fits_list)
}

#### 1 hip ####
load("../2025/data/count/Act_Analysis_old.RData")
load("../2025/data/count/Flags_Analysis_old.RData")
load("../2025/data/count/data_analysis_old.RData")
load("../2025/data/runlength/event_analysis_old.RData")

### 1-1 Unmarked ####
daily_data <- event_analysis %>% mutate(original_value = 1)

fits_list_all <- fit_list_hawkes_marked(daily_data,
  single_day = FALSE,
  penalty_coef = 0, n_cores = 40, data_type = "hip"
)
convergence_results_all <- check_convergence(fits_list_all)
fits_hawkes_unmarked_all <- process_fits_data(fits_list_all)
View(fits_hawkes_unmarked_all)

dir_path <- "../2025/data/hawkes/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
save(fits_list_all, file = paste(dir_path, "fits_unmarked_all_old.RData", sep = ""))

### 1-2 Marked ####
daily_data <- event_analysis

fits_list_all <- fit_list_hawkes_marked(daily_data,
  single_day = FALSE,
  penalty_coef = 0, n_cores = 40, data_type = "hip"
)
convergence_results_all <- check_convergence(fits_list_all)
fits_hawkes_marked_all <- process_fits_data(fits_list_all)
View(fits_hawkes_marked_all)

dir_path <- "../2025/data/hawkes/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
save(fits_list_all, file = paste(dir_path, "fits_marked_all_old.RData", sep = ""))

#### 2 wrist ####
load("../2025/data/mims/Act_Analysis_new.RData")
load("../2025/data/mims/Flags_Analysis_new.RData")
load("../2025/data/mims/data_analysis_new.RData")
load("../2025/data/runlength/event_analysis_new.RData")

### 2-1 Unmarked ####
daily_data <- event_analysis %>% mutate(original_value = 1)

fits_list_all <- fit_list_hawkes_marked(daily_data,
  single_day = FALSE,
  penalty_coef = 0, n_cores = 40, data_type = "wrist"
)
convergence_results_all <- check_convergence(fits_list_all)
fits_hawkes_unmarked_all <- process_fits_data(fits_list_all)

dir_path <- "../2025/data/hawkes/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
save(fits_list_all, file = paste(dir_path, "fits_unmarked_all_new.RData", sep = ""))

### 2-2 Marked ####
daily_data <- event_analysis

fits_list_all <- fit_list_hawkes_marked(daily_data,
  single_day = FALSE,
  penalty_coef = 0, n_cores = 40, data_type = "wrist"
)
convergence_results_all <- check_convergence(fits_list_all)
fits_hawkes_marked_all <- process_fits_data(fits_list_all)

dir_path <- "../2025/data/hawkes/"
for (path in c(dir_path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
save(fits_list_all, file = paste(dir_path, "fits_marked_all_new.RData", sep = ""))

#### 3 Compare unmarked vs marked ####
# Load data (choose old or new)
load("../2025/data/hawkes/fits_unmarked_all_old.RData")
fits_unmarked <- fits_list_all
load("../2025/data/hawkes/fits_marked_all_old.RData")
fits_marked <- fits_list_all

load("../2025/data/hawkes/fits_unmarked_all_new.RData")
fits_unmarked <- fits_list_all
load("../2025/data/hawkes/fits_marked_all_new.RData")
fits_marked <- fits_list_all

# Process fits and add identifiers
df_unmarked <- process_fits_data(fits_unmarked)
df_marked <- process_fits_data(fits_marked)
df_unmarked$model <- "unmarked"
df_marked$model <- "marked"
df_unmarked$n <- df_unmarked$alpha / df_unmarked$beta
df_marked$n <- df_marked$alpha / df_marked$beta

df_unmarked$n <- plogis(df_unmarked$logit_abratio)                
df_marked$n <- plogis(df_marked$logit_abratio) 

# Create wide format for paired comparison
df_wide <- merge(
  df_unmarked[, c("id", "mu", "alpha", "beta", "n", "ks_D")],
  df_marked[, c("id", "mu", "alpha", "beta", "n", "ks_D")],
  by = "id", suffixes = c("_unmarked", "_marked")
)

cor(df_wide$mu_unmarked, df_wide$n_unmarked)
cor(df_wide$mu_marked, df_wide$n_marked)
t.test(df_wide$alpha_marked, df_wide$alpha_unmarked)
t.test(df_wide$beta_marked, df_wide$beta_unmarked)
t.test(df_wide$n_marked, df_wide$n_unmarked)
library(ggplot2)
library(effsize)

# Calculate paired Cohen's d for each parameter
params <- c("mu", "alpha", "beta", "n", "ks_D")
effect_results <- lapply(params, function(p) {
  d <- cohen.d(df_wide[[paste0(p, "_unmarked")]], df_wide[[paste0(p, "_marked")]],
               paired = TRUE, na.rm = TRUE)$estimate
  list(param = p, cohens_d = d)
})

# Create labels with effect size
sig_labels <- sapply(effect_results, function(x) {
  d <- x$cohens_d
  effect <- ifelse(abs(d) < 0.2, "negligible",
    ifelse(abs(d) < 0.5, "small",
      ifelse(abs(d) < 0.8, "medium", "large")
    )
  )
  paste0(x$param, "\n(d=", round(d, 3), ", ", effect, ")")
})
names(sig_labels) <- params

# Reshape data for plotting
df_plot <- bind_rows(
  df_unmarked %>% select(id, model, mu, alpha, beta, n, ks_D),
  df_marked %>% select(id, model, mu, alpha, beta, n, ks_D)
) %>%
  pivot_longer(
    cols = c(mu, alpha, beta, n, ks_D),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(parameter = factor(parameter, levels = params))

p_hist <- ggplot(df_plot, aes(x = value, fill = model)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  facet_wrap(~parameter, scales = "free", ncol = 2, labeller = labeller(parameter = sig_labels)) +
  scale_fill_manual(values = c("unmarked" = "steelblue", "marked" = "coral")) +
  labs(
    title = "Distribution of Hawkes Parameters: Unmarked vs Marked",
    subtitle = "Paired Cohen's d: |d|<0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, >0.8 large",
    x = "Value",
    y = "Count",
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )

print(p_hist)

p_violin <- ggplot(df_plot, aes(x = model, y = value, fill = model)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8, outlier.size = 0.5) +
  facet_wrap(~parameter, scales = "free", ncol = 3, labeller = labeller(parameter = sig_labels)) +
  scale_fill_manual(values = c("unmarked" = "steelblue", "marked" = "coral")) +
  labs(
    title = "Violin Plot of Hawkes Parameters: Unmarked vs Marked",
    subtitle = "Paired Cohen's d: |d|<0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, >0.8 large",
    x = "Model",
    y = "Value",
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.subtitle = element_text(size = 8)
  )

print(p_violin)

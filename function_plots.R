# Common theme for all plots
common_theme <- theme_bw() + theme(
  legend.title = element_blank(),
  legend.key.height = unit(0.5, "line"),
  text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  legend.text = element_text(size = 7)
)

plot_daily <- function(df, weekday, data_type) {
  ggplot(df, aes(x = time, y = value, group = SEQN, color = SEQN)) +
    geom_line(linewidth = 0.5) +
    scale_x_continuous(
      breaks = c(1, 481, 721, 1021, 1201),
      labels = c("12 AM", "8 AM", "12 PM", "5 PM", "8 PM")
    ) +
    scale_color_manual(values = c("red", "deepskyblue")) +
    common_theme +
    theme(legend.position = c(0.16, 0.95)) +
    labs(x = "", y = "", title = paste(data_type, "data: Day", weekday))
}

plot_weekly <- function(df, data_type) {
  ggplot(df, aes(
    x = interaction(WEEKDAY, time, lex.order = TRUE),
    y = value, group = SEQN, colour = SEQN
  )) +
    geom_line() +
    geom_vline(xintercept = 1440 * (0:6), linetype = "dotted", linewidth = 1) +
    annotate("text",
      x = 100 + 1440 * (0:6), y = -0.10,
      label = levels(df$WEEKDAY), size = 3
    ) +
    common_theme +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(x = "", y = "", title = paste(data_type, "data: Weekly"))
}

plot_daily_or_weekly <- function(df, id, weekday = NULL, data_type, dataset, period) {
  df_formatted <- reformat(df, dataset, period, id, weekday)
  if (is.null(weekday)) {
    plot_weekly(df_formatted, data_type)
  } else {
    plot_daily(df_formatted, weekday, data_type)
  }
}

reformat <- function(df, dataset = c("daily_long", "weekly_wide"),
                     period = c("daily", "weekly"), id = NULL, weekday = NULL) {
  dataset <- match.arg(dataset)
  period <- match.arg(period)

  if (!is.null(id)) df <- df %>% filter(SEQN %in% id)

  if (dataset == "daily_long") {
    if (period == "weekly") {
      df %>%
        select(SEQN, WEEKDAY, starts_with("MIN")) %>%
        pivot_longer(cols = starts_with("MIN"), names_to = "MIN", values_to = "value") %>%
        mutate(Accelerometry = paste(WEEKDAY, MIN, sep = "_")) %>%
        arrange(SEQN, WEEKDAY) %>%
        pivot_wider(names_from = Accelerometry, values_from = value)
    } else {
      df %>%
        filter(if (!is.null(weekday)) WEEKDAY == weekday else TRUE) %>%
        select(SEQN, WEEKDAY, starts_with("MIN")) %>%
        pivot_longer(cols = starts_with("MIN"), names_to = "MIN", values_to = "value") %>%
        mutate(time = as.numeric(str_extract(MIN, "\\d+")))
    }
  } else {
    if (period == "daily") {
      df %>%
        pivot_longer(cols = -SEQN, names_to = "Accelerometry", values_to = "value") %>%
        separate(Accelerometry, into = c("WEEKDAY", "MIN"), sep = "_") %>%
        mutate(
          WEEKDAY = factor(WEEKDAY, levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")),
          time = as.numeric(str_extract(MIN, "\\d+"))
        ) %>%
        filter(if (!is.null(weekday)) WEEKDAY == weekday else TRUE)
    } else {
      df
    }
  }
}

#' plot functions 
#' After doing run-length transformation
#' 
plot_activity_and_runlength <- function(Act_Analysis, seqn, weekday, rle) {
  # Get activity data for specific subject and day
  act_day <- Act_Analysis[Act_Analysis$SEQN == seqn & Act_Analysis$WEEKDAY == weekday, ]

  if (nrow(act_day) == 0) {
    stop("No activity data found for subject ", seqn, " on weekday ", weekday)
  }

  # Extract minute-by-minute activity data
  activity <- as.numeric(act_day[1, paste0("MIN", 1:1440)])
  time_minutes <- 1:1440
  time_hours <- time_minutes / 60

  # Set up two-panel plot
  par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))

  # Top panel: Original activity data
  plot(time_hours, activity,
    type = "l",
    main = paste("Subject", seqn, "- Weekday", weekday, "- Activity Data"),
    xlab = "",
    ylab = "Acceleration",
    col = "gray60",
    lwd = 1
  )

  # Bottom panel: Run length data
  # Filter run data for this subject and day
  subject_run_data <- rle[rle$SEQN == seqn & rle$WEEKDAY == weekday, ]

  # Create colors based on activity categories
  color_map <- c("sedentary" = "green", "active" = "blue", "light_moderate" = "#FFB3B3", "high" = "yellow")
  run_colors <- color_map[subject_run_data$categories]
  run_colors[is.na(run_colors)] <- "gray" # Default for missing categories

  plot(subject_run_data$start / 60, subject_run_data$grid_value,
    type = "h",
    xlim = c(0, 24),
    main = paste("Subject", seqn, "- Weekday", weekday, "- Run Length Data"),
    xlab = "Hour of Day",
    ylab = "Acceleration",
    col = run_colors,
    lwd = 2
  )

  # Add legend for activity categories
  unique_categories <- unique(subject_run_data$categories)
  legend_colors <- color_map[unique_categories[!is.na(unique_categories)]]
  legend_labels <- unique_categories[!is.na(unique_categories)]
  legend("topright",
    legend = legend_labels,
    col = legend_colors,
    lty = 1, lwd = 2,
    cex = 0.8
  )

  # Reset plot layout
  par(mfrow = c(1, 1))
}

# Function to plot 24-hour data with waking periods highlighted
plot_subject_flag <- function(seqn, weekday, data_type = c("wrist", "hip")) {
  data_type <- match.arg(data_type)

  # Get data for specific subject and day
  act_day <- Act_Analysis[Act_Analysis$SEQN == seqn & Act_Analysis$WEEKDAY == weekday, ]
  flags_day <- Flags_Analysis[Flags_Analysis$SEQN == seqn & Flags_Analysis$WEEKDAY == weekday, ]

  if (nrow(act_day) == 0 || nrow(flags_day) == 0) {
    stop("No data found for subject ", seqn, " on weekday ", weekday)
  }

  # Extract minute-by-minute data
  activity <- as.numeric(act_day[1, paste0("MIN", 1:1440)])
  flags <- as.numeric(flags_day[1, paste0("MIN", 1:1440)])

  # Create time axis (hours)
  time_hours <- (1:1440) / 60

  # Create the plot
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

  # Plot 1: Activity data
  plot(time_hours, activity,
    type = "l",
    main = paste("Subject", seqn, "- Weekday", weekday),
    ylab = "Acceleration",
    col = "gray60", lwd = 1
  )

  # Plot 2: Flags with different colors based on data type
  if (data_type == "hip") {
    # Hip data: 0=non-wear, 1=wear
    flag_colors <- c("red", "green")
    point_colors <- flag_colors[flags + 1] # 0->1, 1->2
    ylab_text <- "Flag (0=Non-wear, 1=Wear)"

    unique_flags <- sort(unique(flags))
    legend_colors <- flag_colors[unique_flags + 1]
    legend_labels <- c("Non-wear", "Wear")[unique_flags + 1]
  } else {
    # Wrist data: 1=wake, 2=sleep, 3=non-wear, and potentially others
    unique_flags <- sort(unique(flags))
    max_flag <- max(unique_flags)

    # Create color mapping
    point_colors <- rep("black", length(flags)) # default color
    legend_colors <- c()
    legend_labels <- c()

    for (flag_val in unique_flags) {
      if (flag_val == 1) {
        color <- "green"
        label <- "Wake"
      } else if (flag_val == 2) {
        color <- "blue"
        label <- "Sleep"
      } else if (flag_val == 3) {
        color <- "yellow"
        label <- "Non-wear"
      } else {
        color <- "red"
        label <- "Unknown"
      }
      point_colors[flags == flag_val] <- color
      legend_colors <- c(legend_colors, color)
      legend_labels <- c(legend_labels, label)
    }

    ylab_text <- "Flag (1=Wake, 2=Sleep, 3=Non-wear, 4=Unknown)"
  }

  plot(time_hours, flags,
    type = "p",
    main = paste("Subject", seqn, "- Weekday", weekday, paste0("(", data_type, " data)")),
    ylab = ylab_text,
    col = point_colors,
    pch = 16
  )

  # Add legend
  legend("right",
    legend = legend_labels,
    col = legend_colors,
    pch = 16,
    cex = 0.8
  )
}

plot_events <- function(Act_Analysis, seqn, weekday, event_data = NULL) {
  # Get activity data for specific subject and day
  act_day <- Act_Analysis[Act_Analysis$SEQN == seqn & Act_Analysis$WEEKDAY == weekday, ]

  if (nrow(act_day) == 0) {
    stop("No activity data found for subject ", seqn, " on weekday ", weekday)
  }

  # Extract minute-by-minute activity data
  activity <- as.numeric(act_day[1, paste0("MIN", 1:1440)])
  time_minutes <- 1:1440
  time_hours <- time_minutes / 60

  # Set up three-panel plot: raw data + event types + marked values
  n_panels <- ifelse(is.null(event_data), 1, 3)
  par(mfrow = c(n_panels, 1), mar = c(3, 4, 3, 2))

  # Top panel: Original activity data
  plot(time_hours, activity,
    type = "l",
    main = paste("Subject", seqn, "- Weekday", weekday, "- Raw Activity Data"),
    xlab = "",
    ylab = "Acceleration",
    col = "gray60",
    lwd = 1
  )

  # Second panel: event types/categories (multivariate aspect)
  if (!is.null(event_data)) {
    subject_event_data <- event_data[event_data$SEQN == seqn & event_data$WEEKDAY == weekday, ]

    if (nrow(subject_event_data) > 0) {
      # Create colors based on categories
      category_colors <- c("sedentary" = "blue", "active" = "green", "light_moderate" = "orange", "high" = "red")
      event_colors <- category_colors[subject_event_data$categories]
      event_colors[is.na(event_colors)] <- "gray"

      # Plot event types colored by categories
      plot(subject_event_data$start / 60, subject_event_data$run_id,
        type = "h",
        xlim = c(0, 24),
        main = paste("Subject", seqn, "- Weekday", weekday, "- Event Types (Categories)"),
        xlab = "",
        ylab = "Run ID",
        col = event_colors,
        lwd = 2
      )

      # Add legend for categories
      unique_categories <- unique(subject_event_data$categories)
      legend_colors <- category_colors[unique_categories]
      legend("topright",
        legend = unique_categories,
        col = legend_colors,
        lty = 1, lwd = 2,
        cex = 0.8
      )
    } else {
      plot(0, 0,
        type = "n", xlim = c(0, 24), ylim = c(0, 1),
        main = paste("Subject", seqn, "- Weekday", weekday, "- No Event Data"),
        xlab = "", ylab = "No Data"
      )
    }

    # Third panel: marked values with categories
    if (nrow(subject_event_data) > 0) {
      # Use same colors as categories for consistency
      mark_colors <- category_colors[subject_event_data$categories]
      mark_colors[is.na(mark_colors)] <- "gray"

      # Plot original values colored by categories
      plot(subject_event_data$start / 60, subject_event_data$original_value,
        type = "h",
        xlim = c(0, 24),
        main = paste("Subject", seqn, "- Weekday", weekday, "- Marked Hawkes Events"),
        xlab = "Hour of Day",
        ylab = "Original Value",
        col = mark_colors,
        lwd = 2
      )

      # Add legend for categories
      unique_categories_marked <- unique(subject_event_data$categories)
      legend_colors_marked <- category_colors[unique_categories_marked]
      legend("topright",
        legend = unique_categories_marked,
        col = legend_colors_marked,
        lty = 1, lwd = 2,
        cex = 0.8
      )
    } else {
      plot(0, 0,
        type = "n", xlim = c(0, 24), ylim = c(0, 1),
        main = paste("Subject", seqn, "- Weekday", weekday, "- No Marked Data"),
        xlab = "Hour of Day", ylab = "No Data"
      )
    }
  }

  # Reset plot layout
  par(mfrow = c(1, 1))
}
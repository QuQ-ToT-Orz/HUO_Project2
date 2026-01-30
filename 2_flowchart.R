#### Full Preprocessing Flowchart with Count Extraction ####
# Includes steps from 1_data_*.R and 2_Preprocessing.R

library(DiagrammeR)
library(dplyr)

# ==============================================================================
# PART 1: Extract counts from 1_data_*.R preprocessing
# ==============================================================================

#' Extract exclusion counts from Step 1 preprocessing (1_data_*.R)
#' Run this WITHIN 1_data_wrist.R or 1_data_hip.R at the appropriate lines
#'
#' @param AllAct The AllAct dataframe before any exclusions
#' @param table_dat The table_dat dataframe after age exclusion
#' @param nms_rm Vector of SEQNs to remove (bad accelerometer data)
#' @return List of counts for each exclusion criterion
extract_step1_counts <- function(AllAct, table_dat, nms_rm) {

  # Count before any exclusions (unique participants)
  n_initial <- length(unique(AllAct$SEQN))

  # Count after age exclusion
  n_after_age <- nrow(table_dat)

  # Individual exclusion criteria counts
  counts <- list(
    n_initial = n_initial,
    n_after_age = n_after_age,

    # Missing BMI
    n_missing_bmi = sum(is.na(table_dat$BMI_cat), na.rm = TRUE),

    # Missing education
    n_missing_edu = sum(is.na(table_dat$EducationAdult), na.rm = TRUE),

    # Bad accelerometer data (<7 days with ≥10h wear)
    n_bad_accel = sum(table_dat$SEQN %in% nms_rm, na.rm = TRUE),

    # Missing/invalid mortality or accidental death
    n_missing_mort = sum(
      (!table_dat$eligstat %in% 1) |
        is.na(table_dat$mortstat) |
        is.na(table_dat$permth_exm) |
        (table_dat$ucod_leading %in% "004"),
      na.rm = TRUE
    ),

    # Alive with <5 years follow-up
    n_short_followup = sum(
      table_dat$mortstat == 0 & (table_dat$permth_exm / 12 < 5),
      na.rm = TRUE
    ),

    # Missing lab measures (SYS, total cholesterol, HDL)
    n_missing_lab = sum(
      is.na(table_dat$SYS) | is.na(table_dat$LBXTC) | is.na(table_dat$LBDHDD),
      na.rm = TRUE
    )
  )

  return(counts)
}

# ==============================================================================
# PART 2: Extract counts from 2_Preprocessing.R
# ==============================================================================

# Full Step 2 Pipeline:
# 1. simple_events()           → events (minute-level activity events)
# 2. process_daily_data()      → filtered_events (gap removal, event rate filtering, 7-day requirement)
# 3. clean_events_spikes()     → event_analysis (spike cleaning)
# 4. reconstruct_act_analysis() → filtered_act_analysis (reconstruct Act_Analysis format)
# 5. runlength_single_manual() → result_rle (run-length encoding)
# 6. filter_rle_data()         → rle_analysis (≥5 events/category/day, 7-day requirement)
# 7. intersect()               → common_subjects (align event_analysis & rle_analysis)

#' Extract counts from Step 2 preprocessing (2_Preprocessing.R)
#'
#' @param Act_Analysis The Act_Analysis dataframe from Step 1
#' @param events The events dataframe after simple_events()
#' @param filtered_events The filtered_events dataframe after process_daily_data()
#' @param event_analysis The event_analysis dataframe after clean_events_spikes()
#' @param rle_analysis The rle_analysis dataframe after filter_rle_data()
#' @param common_subjects The common_subjects vector after intersect()
#' @return List of counts for Step 2 filtering
extract_step2_counts <- function(Act_Analysis, events, filtered_events,
                                  event_analysis, rle_analysis, common_subjects) {

  counts <- list(
    # After Step 1 (input to Step 2)
    n_step1_subjects = length(unique(Act_Analysis$SEQN)),
    n_step1_days = nrow(Act_Analysis),

    # After simple_events (minute-level events)
    n_after_simple_events = length(unique(events$SEQN)),

    # After process_daily_data (gap removal + event rate filtering + 7-day req)
    n_after_daily_processing = length(unique(filtered_events$SEQN)),

    # After clean_events_spikes
    n_after_spike_cleaning = length(unique(event_analysis$SEQN)),

    # After filter_rle_data (≥5 events/category/day + 7-day req)
    n_after_rle_filter = length(unique(rle_analysis$SEQN)),

    # Final aligned subjects (intersect of event_analysis and rle_analysis)
    n_final = length(common_subjects)
  )

  # Calculate exclusions at each step
  counts$excluded_simple_events <- counts$n_step1_subjects - counts$n_after_simple_events
  counts$excluded_daily_processing <- counts$n_after_simple_events - counts$n_after_daily_processing
  counts$excluded_spike_cleaning <- counts$n_after_daily_processing - counts$n_after_spike_cleaning
  counts$excluded_rle_filter <- counts$n_after_spike_cleaning - counts$n_after_rle_filter
  counts$excluded_alignment <- counts$n_after_rle_filter - counts$n_final

  return(counts)
}

# ==============================================================================
# PART 3: Create the full PRISMA flowchart
# ==============================================================================

create_full_flowchart <- function(
    dataset_name = "Dataset",
    # Step 1 counts (from 1_data_*.R)
    n_initial = NA,
    n_after_age = NA,
    n_missing_bmi = NA,
    n_missing_edu = NA,
    n_bad_accel = NA,
    n_missing_mort = NA,
    n_short_followup = NA,
    n_missing_lab = NA,
    n_step1_final = NA,
    # Step 2 counts (from 2_Preprocessing.R)
    n_after_daily_process = NA,  # After process_daily_data()
    n_after_spike_clean = NA,    # After clean_events_spikes()
    n_after_rle_filter = NA,     # After filter_rle_data()
    n_final = NA) {              # After intersect() alignment

  # Calculate totals
  n_total_step1_excluded <- n_after_age - n_step1_final
  n_total_step2_excluded <- n_step1_final - n_final

  # Step 2 sub-exclusions
  n_excluded_daily <- n_step1_final - n_after_daily_process
  n_excluded_spikes <- n_after_daily_process - n_after_spike_clean
  n_excluded_rle <- n_after_spike_clean - n_after_rle_filter
  n_excluded_align <- n_after_rle_filter - n_final

  flowchart <- grViz(sprintf('
    digraph full_flowchart {
      graph [layout = dot, rankdir = TB, splines = ortho, nodesep = 0.5, ranksep = 0.6]
      node [shape = rectangle, style = filled, fontname = "Helvetica", fontsize = 10]
      edge [color = "#555555"]

      # ==================== IDENTIFICATION ====================
      subgraph cluster_id {
        label = "Identification"
        style = rounded
        color = "#1565C0"
        fontcolor = "#1565C0"
        fontsize = 12
        bgcolor = "#E3F2FD"

        id_box [
          label = "%s\\n\\nParticipants with\\naccelerometer data\\nn = %s",
          fillcolor = "#BBDEFB"
        ]
      }

      # ==================== STEP 1: DATA CLEANING ====================
      subgraph cluster_step1 {
        label = "Screening"
        style = rounded
        color = "#7B1FA2"
        fontcolor = "#7B1FA2"
        fontsize = 12
        bgcolor = "#F3E5F5"

        age_box [
          label = "After excluding\\nmissing age / age ≥85\\nn = %s",
          fillcolor = "#E1BEE7"
        ]

        step1_final [
          label = "Meeting all\\ninclusion criteria\\nn = %s",
          fillcolor = "#CE93D8"
        ]
      }

      # Step 1 Exclusion Box
      excl_step1 [
        label = <<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1">
          <TR><TD><B>Excluded* (n = %s)</B></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Missing BMI: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Missing education: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• &lt;7 days with ≥10h wear time: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Ineligible for mortality linkage,</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">  missing mortality, or accidental death: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Alive with &lt;5 years follow-up: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Missing lab measures (SBP/TC/HDL): %s</FONT></TD></TR>
        </TABLE>>,
        shape = rectangle,
        fillcolor = "#FFF3E0",
        color = "#E65100"
      ]

      # ==================== STEP 2: EVENT QUALITY CONTROL ====================
      subgraph cluster_step2 {
        label = "Eligibility"
        style = rounded
        color = "#00695C"
        fontcolor = "#00695C"
        fontsize = 12
        bgcolor = "#E0F2F1"

        step2_input [
          label = "Participants for\\nevent-level analysis\\nn = %s",
          fillcolor = "#B2DFDB"
        ]

        step2_final [
          label = "After event\\nquality control\\nn = %s",
          fillcolor = "#80CBC4"
        ]
      }

      # Step 2 Exclusion Box
      excl_step2 [
        label = <<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1">
          <TR><TD><B>Excluded (n = %s)</B></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Gaps &gt;60min in daily recording,</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">  sedentary or active rate &lt;1/hour,</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">  active rate &gt;50/hour,</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">  or &lt;7 valid days: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Activity spike artifacts: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• &lt;5 sedentary or &lt;5 active</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">  events per day: %s</FONT></TD></TR>
          <TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9">• Misaligned event/bout data: %s</FONT></TD></TR>
        </TABLE>>,
        shape = rectangle,
        fillcolor = "#FFF3E0",
        color = "#E65100"
      ]

      # ==================== FINAL SAMPLE ====================
      subgraph cluster_final {
        label = "Included in Analysis"
        style = rounded
        color = "#2E7D32"
        fontcolor = "#2E7D32"
        fontsize = 12
        bgcolor = "#E8F5E9"

        final_box [
          label = "Final Analytical Sample\\n\\nn = %s participants\\n%s person-days",
          fillcolor = "#A5D6A7",
          fontsize = 11
        ]
      }

      # ==================== EDGES ====================
      id_box -> age_box
      age_box -> excl_step1 [style = dashed, color = "#E65100"]
      age_box -> step1_final
      step1_final -> step2_input
      step2_input -> excl_step2 [style = dashed, color = "#E65100"]
      step2_input -> step2_final
      step2_final -> final_box
    }
  ',
  # Dataset name
  dataset_name,
  # Identification
  format(n_initial, big.mark = ","),
  # After age
  format(n_after_age, big.mark = ","),
  # Step 1 final
  format(n_step1_final, big.mark = ","),
  # Step 1 exclusions
  format(n_total_step1_excluded, big.mark = ","),
  format(n_missing_bmi, big.mark = ","),
  format(n_missing_edu, big.mark = ","),
  format(n_bad_accel, big.mark = ","),
  format(n_missing_mort, big.mark = ","),
  format(n_short_followup, big.mark = ","),
  format(n_missing_lab, big.mark = ","),
  # Step 2 input
  format(n_step1_final, big.mark = ","),
  # Step 2 final
  format(n_final, big.mark = ","),
  # Step 2 exclusions
  format(n_total_step2_excluded, big.mark = ","),
  format(n_excluded_daily, big.mark = ","),
  format(n_excluded_spikes, big.mark = ","),
  format(n_excluded_rle, big.mark = ","),
  format(n_excluded_align, big.mark = ","),
  # Final
  format(n_final, big.mark = ","),
  format(n_final * 7, big.mark = ",")
  ))

  return(flowchart)
}

# ==============================================================================
# PART 4: Example usage - Replace with your actual counts
# ==============================================================================

# ----- HIP DATA (NHANES 2003-2006) -----
# Run these lines inside 1_data_hip.R after the relevant sections:
#
# After line 631 (nparticipants <- dim(table_dat)[1]):
#   n_initial_hip <- nparticipants
#
# After line 636 (table_dat <- subset(table_dat, !is.na(Age))):
#   n_after_age_hip <- nrow(table_dat)
#
# After line 693 (data_analysis <- subset(table_dat, Exclude == 0, select = -Exclude)):
#   n_step1_final_hip <- nrow(data_analysis)
#   step1_counts_hip <- extract_step1_counts(AllAct, table_dat, nms_rm)
#
# In 2_Preprocessing.R after line 342:
#   step2_counts_hip <- extract_step2_counts(Act_Analysis, filtered_events, event_analysis, rle_analysis)

# ----- WRIST DATA (NHANES 2011-2014) -----
# Same approach for 1_data_wrist.R

# ==============================================================================
# Example flowchart with placeholder values
# ==============================================================================

# Hip/Count data flowchart
flowchart_hip <- create_full_flowchart(
  dataset_name = "NHANES 2003-2006\\n(Hip-worn, Activity Counts)",
  # Step 1 (1_data_hip.R)
  n_initial = 14631,
  n_after_age = 14331,
  n_missing_bmi = 96,
  n_missing_edu = 6001,
  n_bad_accel = 10437,
  n_missing_mort = 5083,
  n_short_followup = 0,
  n_missing_lab = 1976,
  n_step1_final = 2607,
  # Step 2 (2_Preprocessing.R)
  n_after_daily_process = 2571,  # After process_daily_data()
  n_after_spike_clean = 2571,    # After clean_events_spikes()
  n_after_rle_filter = 2564,     # After filter_rle_data()
  n_final = 2564                 # After intersect() alignment
)

# Wrist/MIMS data flowchart
flowchart_wrist <- create_full_flowchart(
  dataset_name = "NHANES 2011-2014\\n(Wrist-worn, MIMS)",
  # Step 1 (1_data_wrist.R)
  n_initial = 6905,
  n_after_age = 6905,
  n_missing_bmi = 82,
  n_missing_edu = 1564,
  n_bad_accel = 2554,
  n_missing_mort = 1318,
  n_short_followup = 18,
  n_missing_lab = 932,
  n_step1_final = 3208,
  # Step 2 (2_Preprocessing.R)
  n_after_daily_process = 3120,  # After process_daily_data()
  n_after_spike_clean = 3120,    # After clean_events_spikes()
  n_after_rle_filter = 3120,     # After filter_rle_data()
  n_final = 3120                 # After intersect() alignment
)

# Display
print(flowchart_hip)
print(flowchart_wrist)

# ==============================================================================
# PART 5: Code already added to preprocessing scripts
# ==============================================================================

# The flowchart count extraction code has been added to the preprocessing files:
#
# 1_data_hip.R (line 696-720):
#    - flowchart_counts_hip with: n_initial, n_excluded_age, n_after_age,
#      n_missing_bmi, n_missing_edu, n_bad_accel, n_missing_mort,
#      n_short_followup, n_missing_lab, n_step1_final
#
# 1_data_wrist.R (line 585-609):
#    - flowchart_counts_wrist with same fields as hip
#
# 2_Preprocessing.R:
#    - Line 189: n_after_daily_process_hip
#    - Line 229: n_after_daily_process_wrist
#    - Line 327-328: n_after_spike_clean_hip, n_after_rle_filter_hip
#    - Line 333: n_final_hip
#    - Line 339-340: n_after_spike_clean_wrist, n_after_rle_filter_wrist
#    - Line 345: n_final_wrist
#    - Line 354-383: flowchart_counts_step2 with all hip and wrist counts
#
# To generate the flowchart with actual values, run the preprocessing scripts
# first, then load the saved counts and call create_full_flowchart()

# ==============================================================================
# PART 6: Save flowchart
# ==============================================================================

save_flowchart <- function(flowchart, filename) {
  if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
    install.packages("DiagrammeRsvg")
  }
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    install.packages("rsvg")
  }

  svg_code <- DiagrammeRsvg::export_svg(flowchart)

  if (grepl("\\.pdf$", filename)) {
    rsvg::rsvg_pdf(charToRaw(svg_code), filename, width = 900, height = 1100)
  } else if (grepl("\\.png$", filename)) {
    rsvg::rsvg_png(charToRaw(svg_code), filename, width = 900, height = 1100)
  } else if (grepl("\\.svg$", filename)) {
    writeLines(svg_code, filename)
  }
  message("Saved: ", filename)
}

# Uncomment to save:
# save_flowchart(flowchart_hip, "flowchart_hip.pdf")
# save_flowchart(flowchart_wrist, "flowchart_wrist.pdf")

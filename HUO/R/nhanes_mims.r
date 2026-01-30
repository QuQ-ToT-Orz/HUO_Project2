# Declare global variables for R CMD check (dplyr/data.table NSE)
utils::globalVariables(c(
  "PAXDAYM", "PAXQFM", "SEQN", "PAXDAYWM", "PAXSSNMP", "PAXTRANM",
  "PAXPREDM", "MIN", "PAXMTSM", "..selected_vars"
))

#' Process NHANES 2011-2014 MIMS accelerometer data
#'
#' @description
#' This function processes NHANES 2011-2014 accelerometer data in MIMS format,
#' converting from minute-level data to wide format with proper quality control.
#'
#' @param PAXMIN Data frame with MIMS accelerometer data containing variables:
#'   SEQN, PAXDAYM, PAXDAYWM, PAXMTSM, PAXPREDM, PAXTRANM, PAXQFM, PAXSSNMP
#'
#' @details
#' This function:
#' 1. Filters out first/last days and bad quality data
#' 2. Handles transition minutes by forward-filling predictions
#' 3. Ensures complete 7-day, 1440-minute datasets
#' 4. Creates separate wide-format datasets for activity counts and flags
#'
#' The PAXPREDM variable codes:
#' - 1: wake/wear
#' - 2: sleep/wear
#' - 3: non-wear
#' - 4: unknown
#'
#' @return
#' List with two elements:
#' - accel: Wide format data frame with MIMS values (MIN1-MIN1440)
#' - flags: Wide format data frame with wear/sleep/non-wear flags (MIN1-MIN1440)
#'
#' @importFrom dplyr filter arrange select group_by mutate ungroup n_distinct n row_number
#' @importFrom tidyr fill
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#'
#' @export
processNHANESAccelerometer <- function(PAXMIN) {
  # Process the data first
  PAXMIN <- PAXMIN %>%
    # Step 1: Filter out first/last days and bad quality data
    filter(PAXDAYM != 1, PAXDAYM != 9) %>%
    filter(PAXQFM == 0) %>%
    arrange(SEQN, PAXDAYWM, PAXSSNMP) %>%
    select(-PAXSSNMP, -PAXDAYM) %>%
    # Step 2: Handle transition minutes
    group_by(SEQN, PAXDAYWM) %>%
    mutate(
      PAXPREDM = ifelse(PAXTRANM == 1, NA, PAXPREDM)
    ) %>%
    fill(PAXPREDM, .direction = "down") %>%
    mutate(PAXPREDM = ifelse(is.na(PAXPREDM), 4, PAXPREDM)) %>%
    select(-PAXTRANM) %>%
    # Step 3: Check completeness
    filter(n() == 1440) %>%
    mutate(MIN = row_number()) %>%
    group_by(SEQN) %>%
    filter(n_distinct(PAXDAYWM) == 7) %>%
    ungroup()

  # Step 4: Create two separate wide format datasets with modified column names
  accel <- PAXMIN %>%
    pivot_wider(
      id_cols = c(SEQN, PAXDAYWM, PAXQFM),
      names_from = MIN,
      values_from = PAXMTSM,
      names_prefix = "MIN"
      # names_prefix = "_",
      # names_glue = "{.value}{names_prefix}{sprintf('%d', MIN)}"
    ) %>%
    as.data.frame()

  flags <- PAXMIN %>%
    pivot_wider(
      id_cols = c(SEQN, PAXDAYWM, PAXQFM),
      names_from = MIN,
      values_from = PAXPREDM,
      names_prefix = "MIN"
    ) %>%
    as.data.frame()

  # Return both datasets as a list
  return(list(
    accel = accel,
    flags = flags
  ))
}


#' Extract NHANES laboratory data for specific variables
#'
#' @description
#' This function extracts specific laboratory variables from NHANES data files
#' for a given year, with options to exclude certain tables and variables.
#'
#' @param year Numeric year for NHANES cycle (e.g., 2011 for 2011-2012 cycle)
#' @param exclude_tables Character vector of table names to exclude
#' @param extractAll Logical indicating whether to extract all variables
#' @param varnames Character vector of specific variable names to extract
#'
#' @details
#' This function uses the nhanesA package to download and process laboratory
#' data from NHANES. It automatically handles table merging and variable
#' naming conflicts by prefixing with table names.
#'
#' @return
#' List containing:
#' - data: Merged data frame with requested variables
#' - summary: Summary information about the extraction
#' - variable_info: Information about variables and their source tables
#'
#' @importFrom nhanesA nhanesTables nhanes
#'
#' @export
extract_nhanes_data <- function(year,
                                exclude_tables = character(0),
                                extractAll = FALSE,
                                varnames = c("LBXTC", "LBDHDD")) {
  # Get all LAB tables for the year
  tables <- nhanesTables("LAB", year)
  tables <- tables[!startsWith(tables[, 1], "S"), ]

  # Remove excluded tables from the tables dataframe
  if (length(exclude_tables) > 0) {
    tables <- tables[!tables[, 1] %in% exclude_tables, ]
    if (nrow(tables) == 0) {
      stop("All tables were excluded")
    }
  }

  # Create a named vector of descriptions
  table_descriptions <- setNames(tables[, 2], tables[, 1])

  # Initialize list to store extracted data
  extracted_data_list <- list()

  # Initialize list to store variable information
  var_info_list <- list()

  # Initialize list to track skipped tables
  skipped_tables <- character(0)

  # Initialize list to track used tables
  used_tables <- character(0)

  # Loop through each table
  for (table_name in tables[, 1]) {
    # Get data from current table
    message(sprintf(
      "Extracting data from %s (%s)...",
      table_name, table_descriptions[table_name]
    ))

    # Extract data
    current_data <- nhanes(table_name)
    if (is.null(current_data) || nrow(current_data) == 0) {
      stop(sprintf("Dataset %s not found or empty", table_name))
    }

    # Check if SEQN exists in the table
    if (!"SEQN" %in% names(current_data)) {
      message(sprintf("Skipping table %s - No SEQN found", table_name))
      skipped_tables <- c(skipped_tables, table_name)
      next # Skip to next table
    }

    # If extractAll is FALSE and varnames is not empty, filter columns
    if (!extractAll && length(varnames) > 0) {
      # Always keep SEQN
      cols_to_keep <- c("SEQN", intersect(names(current_data), varnames))
      if (length(cols_to_keep) <= 1) { # Only SEQN
        message(sprintf("Skipping table %s - No requested variables found", table_name))
        skipped_tables <- c(skipped_tables, table_name)
        next # Skip to next table
      }
      current_data <- current_data[, cols_to_keep, drop = FALSE]
      used_tables <- c(used_tables, table_name) # Track that we used this table
    } else {
      used_tables <- c(used_tables, table_name) # Track that we used this table
    }

    # Store variable names and their source table
    var_info <- data.frame(
      Variable = names(current_data),
      SourceTable = rep(table_name, length(names(current_data))),
      TableDescription = rep(table_descriptions[table_name], length(names(current_data))),
      stringsAsFactors = FALSE
    )
    var_info_list[[table_name]] <- var_info

    # Add combined table identifier to prevent column name conflicts
    # but keep SEQN as is
    non_seqn_cols <- names(current_data) != "SEQN"
    names(current_data)[non_seqn_cols] <- paste0(
      table_name, ".",
      names(current_data)[non_seqn_cols]
    )

    # Store in list
    extracted_data_list[[table_name]] <- current_data
  }

  # Check if we have any valid tables
  if (length(extracted_data_list) == 0) {
    stop("No valid tables with SEQN found")
  }

  # Combine all variable information
  all_variables <- do.call(rbind, var_info_list)

  # Merge all datasets using SEQN as key
  message("Merging all datasets...")
  merged_data <- Reduce(
    function(x, y) base::merge(x, y, by = "SEQN", all = TRUE),
    extracted_data_list
  )

  # Create summary of the merged dataset
  summary_info <- list(
    cycle = paste0(year, "-", year + 1),
    n_tables_total = nrow(tables),
    n_tables_included = length(extracted_data_list),
    n_tables_skipped = length(skipped_tables),
    skipped_tables = skipped_tables,
    used_tables = used_tables,
    excluded_tables = exclude_tables,
    n_variables = ncol(merged_data),
    n_observations = nrow(merged_data),
    variables_by_table = table(all_variables$SourceTable),
    missing_by_variable = sapply(merged_data, function(x) sum(is.na(x)))
  )

  # Print summary based on extractAll parameter
  if (!extractAll) {
    message("\nUsed tables for variable extraction:")
    for (table in used_tables) {
      message(sprintf("- %s (%s)", table, table_descriptions[table]))
    }
  } else {
    if (length(skipped_tables) > 0) {
      message("\nSkipped tables (no SEQN):")
      for (table in skipped_tables) {
        message(sprintf("- %s (%s)", table, table_descriptions[table]))
      }
    }
  }

  # Print summary of excluded tables
  if (length(exclude_tables) > 0) {
    message("\nExcluded tables:")
    for (table in exclude_tables) {
      message(sprintf("- %s", table))
    }
  }

  return(list(
    data = merged_data,
    summary = summary_info,
    variable_info = all_variables
  ))
}


#' Read selected variables from NHANES XPT files
#'
#' @description
#' This function reads specific variables from NHANES XPT files, designed for
#' processing large accelerometer data files efficiently.
#'
#' @param file_path Path to the XPT file
#'
#' @details
#' This function is optimized for reading large NHANES accelerometer files by
#' selecting only essential variables and using data.table for efficiency.
#'
#' @return
#' Data table with selected variables: SEQN, PAXDAYM, PAXDAYWM, PAXMTSM,
#' PAXPREDM, PAXTRANM, PAXQFM, PAXSSNMP
#'
#' @importFrom foreign read.xport
#' @importFrom data.table setDT
#'
#' @export
read_selected_xpt <- function(file_path) {
  # Define variables to extract
  selected_vars <- c(
    "SEQN", # Participant ID
    "PAXDAYM", # Day (e.g., 1,2,3)
    "PAXDAYWM", # Day (e.g., Mon, Tue, Wed)
    "PAXMTSM", # MIMS
    "PAXPREDM", # 1-wake/2-sleep/3-non-wear/4-unknown
    "PAXTRANM", # transition for wake/sleep/non-wear
    "PAXQFM", # QC, Values >0 indicate that this minute is invalid based on the QC review
    "PAXSSNMP"
  )

  message("Reading XPT file...")

  # Read XPT file
  data <- foreign::read.xport(file_path)

  # Convert to data.table and select only needed columns
  data.table::setDT(data)
  data <- data[, ..selected_vars]

  return(data)
}
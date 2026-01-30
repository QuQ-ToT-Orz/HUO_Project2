#' Process NHANES 2011-2012 and 2013-2014 mortality data for waves G and H
#'
#' @description
#' This function creates a clean mortality dataset which can be combined with data from the
#' NHANES 2011-2012/2013-2014 waves (G and H).
#'
#' @param waves Character vector indicating the waves. Defaults to a vector with "G" and "H",
#' corresponding to the 2011-2012 and 2013-2014 waves.
#'
#' @param mort_release_yr Numeric value indicating the year associated with the raw mortality
#' data to be processed. The default, 2019, corresponds to the most recent raw mortality data
#' for these waves.
#'
#' @param localpath Character scalar describing the location where the raw data are stored.
#' If NULL, the function will look in package data directory for the requested raw mortality data.
#' Defaults to NULL.
#'
#' @details
#' This function has been adapted for NHANES 2011-2014 waves (G and H) and tested on the 2019
#' release. The raw data comes in the form of fixed-width format files with specific column positions.
#'
#' @return
#' This function returns a list with elements corresponding to each wave specified. Each element
#' is a data frame with mortality information including SEQN, eligibility status, mortality status,
#' follow-up times, and cause of death information.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom readr read_fwf fwf_cols
#'
#' @export
process_mort <- function(waves = c("G", "H"), mort_release_yr = 2019, localpath = NULL) {
  # Input validation
  if (!all(waves %in% LETTERS)) {
    stop("One or more waves invalid")
  }
  if (!is.numeric(mort_release_yr)) {
    stop("mort_release_yr must be numeric")
  }

  # Sort waves and convert to years
  waves <- sort(waves)
  waves_num <- vapply(waves, function(x) {
    seq(1999, 2049, by = 2)[which(LETTERS == x)]
  }, numeric(1))

  # Create file names
  waves_mort <- paste0(
    "NHANES_", waves_num, "_", waves_num + 1,
    "_MORT_", mort_release_yr, "_PUBLIC.dat"
  )

  if (is.null(localpath)) {
    localpath <- system.file(file.path("extdat", "mort"), package = "HUO")
  }
  filepath <- file.path(localpath, waves_mort)

  # Check if files exist
  dne <- !file.exists(filepath)
  if (any(dne)) {
    stop(cat(
      "Specified files: \n ",
      paste0(filepath[which(dne)], collapse = ", \n "),
      "\n not found."
    ))
  }

  # Initialize results list
  ret <- list()
  pb <- txtProgressBar(min = 0, max = length(waves_mort), style = 3)

  # Process each wave
  for (i in seq_along(waves_mort)) {
    out_name <- paste0("Mortality_", mort_release_yr, "_", waves[i])

    # Read data using the new format for 2019 release
    out <- read_fwf(
      file = filepath[i],
      col_types = "iiiiiiii",
      fwf_cols(
        seqn = c(1, 6),
        eligstat = c(15, 15),
        mortstat = c(16, 16),
        ucod_leading = c(17, 19),
        diabetes = c(20, 20),
        hyperten = c(21, 21),
        permth_int = c(43, 45),
        permth_exm = c(46, 48)
      ),
      na = c("", ".")
    )

    ret[[out_name]] <- data.frame(out)
    setTxtProgressBar(pb, i)
  }

  close(pb)
  return(ret)
}


#' Merge non-accelerometry data for NHANES waves G and H
#'
#' @description
#' This function retrieves and merges covariate data from one or more NHANES data files
#' across waves G and H (2011-2014). Variables are merged using the NHANES unique subject
#' identifier (SEQN).
#'
#' @param waves character vector with entries of (capitalized) letter of the alphabet
#' corresponding to the NHANES wave of interest. Defaults to a vector containing "G" and "H"
#' corresponding to the NHANES 2011-2012 and 2013-2014 waves.
#'
#' @param varnames character vector indicating which column names are to be searched for.
#' Will check all .XPT files in located in the directory specified by localpath.
#' If extractAll = TRUE, then this argument is effectively ignored.
#'
#' @param localpath file path where covariate data are saved. Covariate data must be in .XPT format.
#'
#' @param extractAll logical argument indicating whether all columns of all .XPT files in the
#' search path should be returned. Defaults to FALSE.
#'
#' @param exclude_files character vector of file patterns to exclude from processing.
#'
#' @details
#' This function has been adapted for NHANES 2011-2014 waves (G and H). It searches all .XPT
#' files which match the NHANES naming convention and merges variables using SEQN.
#'
#' @return
#' This function returns a list with elements equal to the number of waves specified.
#' Each element is named Covariate_* where * corresponds to the wave letter.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom haven read_xpt
#'
#' @export
process_covar <- function(waves = c("G", "H"),
                          varnames = c(
                            "SDDSRVYR", "WTMEC2YR", "WTINT2YR",
                            "SDMVPSU", "SDMVSTRA",
                            "RIDAGEMN", "RIDAGEEX", "RIDAGEYR", "RIDRETH1", "RIAGENDR",
                            "BMXWT", "BMXHT", "BMXBMI", "DMDEDUC2",
                            "ALQ101", "ALQ110", "ALQ120Q", "ALQ120U", "ALQ130",
                            "SMQ020", "SMD030", "SMQ040",
                            "MCQ220", "MCQ160F", "MCQ160B", "MCQ160C",
                            "PFQ049", "PFQ054", "PFQ057", "PFQ059", "PFQ061B", "PFQ061C",
                            "DIQ010",
                            "LBXTC", "LBXHDD", "LBDHDD",
                            "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"
                          ),
                          localpath = NULL,
                          extractAll = FALSE,
                          exclude_files = c()) {
  if (!all(waves %in% LETTERS)) stop("One or more waves invalid.")
  stopifnot(length(waves) >= 1)
  stopifnot(is.vector(waves))

  waves <- sort(waves)
  cnt <- ifelse("SEQN" %in% varnames, 1, 0)
  varnames <- unique(c("SEQN", varnames))

  if (is.null(localpath)) {
    localpath <- system.file(file.path("extdat", "covar"), package = "HUO")
  }

  files_full <- list.files(localpath)

  # Convert exclude_files to uppercase for case-insensitive matching
  exclude_files <- toupper(exclude_files)

  if (any(!sapply(paste0("DEMO_", waves, ".XPT"), function(pattern) {
    any(grepl(pattern, files_full, ignore.case = TRUE))
  }))) {
    warning(cat(paste0(
      "One or more demographic files were not found in the data directory (",
      paste0("DEMO_", waves, ".XPT", collapse = ", "),
      "). There is no guarantee all participants for a particular wave will be included in the returned object."
    )))
  }

  ret <- rep(list(NULL), length(waves))
  names(ret) <- paste0("Covariate_", waves)
  pb <- txtProgressBar(min = 0, max = length(waves), style = 3)

  for (i in seq_along(waves)) {
    cohort <- waves[i]

    pathExt <- paste0("_", cohort, ".", sep = "")
    # find only those files associated with a particular NHANES wave - case insensitive for XPT extension
    files <- files_full[grepl(paste0(pathExt, "xpt$", collapse = ""), files_full, ignore.case = TRUE)]

    # Exclude files that contain any of the patterns
    files <- files[!sapply(toupper(files), function(x) any(sapply(exclude_files, function(pattern) grepl(pattern, x))))]

    if (length(files) == 0) {
      message(cat(paste0("\n No data associated with wave ", waves[i], " was found. \n")))
      setTxtProgressBar(pb, i)
      next
    }

    # find which files contain the variables requested by the user
    covarMats <- lapply(files, function(x) {
      mat <- try(read_xpt(file.path(localpath, x)))
      if (inherits(mat, "try-error")) {
        return(NULL)
      }
      if (!"SEQN" %in% colnames(mat)[1]) {
        return(NULL)
      }
      if (!extractAll) {
        if (length(setdiff(intersect(colnames(mat), varnames), "SEQN")) == 0) {
          return(NULL)
        }
        mat <- mat[, colnames(mat) %in% varnames, drop = FALSE]
      }
      if (!any(is.null(dim(mat)), all(colnames(mat) %in% "SEQN"))) {
        return(mat)
      } else {
        return(NULL)
      }
    })

    covarMats <- covarMats[!vapply(covarMats, is.null, logical(1))]
    matchedNames <- lapply(covarMats, colnames)
    numMatched <- length(unlist(matchedNames))

    if (numMatched == 0) {
      message(cat(paste0("\n No variables specified by the varnames argument was found for wave ", waves[i], "\n")))
      setTxtProgressBar(pb, i)
      next
    }

    if (numMatched > 0 & !extractAll) {
      message(
        cat(paste(
          "\n For", cohort, "cohort,",
          (numMatched - length(matchedNames) + cnt),
          "Covariates Found of", (length(varnames) - 1 + cnt), "specified.",
          "Missing covariates:",
          paste(setdiff(varnames, unlist(matchedNames)), collapse = ", ")
        ))
      )
    }

    ids <- lapply(covarMats, function(x) x[["SEQN"]])
    uids <- as.integer(sort(unique(unlist(ids))))
    rep_SEQN <- vapply(ids, function(x) any(duplicated(x)), logical(1))
    rep_inx <- which(rep_SEQN)
    notrep_inx <- which(!rep_SEQN)

    totalCols <- sum(vapply(covarMats, ncol, numeric(1))) - length(covarMats) + 1
    CovarMat <- matrix(NA, ncol = totalCols, nrow = length(uids))
    colnames(CovarMat) <- c("SEQN", unlist(sapply(covarMats, function(x) colnames(x)[-1])))
    CovarMat[, "SEQN"] <- uids
    CovarMat <- data.frame(CovarMat)
    invisible(lapply(covarMats[notrep_inx], function(x) CovarMat[, colnames(x)[-1]] <<- x[match(CovarMat$SEQN, x$SEQN), -1, drop = FALSE]))

    if (length(rep_inx) != 0) {
      invisible(lapply(covarMats[rep_inx], function(x) {
        n_vars <- length(colnames(x))
        for (j in 2:n_vars) {
          max_reps <- max(table(x$SEQN), na.rm = TRUE)
          mat_tmp <- matrix(NA, ncol = max_reps, nrow = nrow(CovarMat))

          for (i in 1:nrow(CovarMat)) {
            if (!CovarMat$SEQN[i] %in% x$SEQN) next
            ji_inx <- which(x$SEQN == CovarMat$SEQN[i])
            mat_tmp[i, 1:length(ji_inx)] <- x[ji_inx, j, drop = TRUE]

            rm(list = c("ji_inx"))
          }
          CovarMat[colnames(x)[j]] <<- I(mat_tmp)
          rm(list = c("max_reps", "mat_tmp"))
        }
      }))
      message(
        cat(paste(
          "\n Variables with repeated observations per subject found for the following variables:",
          paste(sapply(covarMats[rep_inx], function(x) colnames(x)[-1]), collapse = ","),
          "Note that these variables will be stored with class AsIs() objects in resulting data frames. ",
          "See ?I for details on AsIs class. \n"
        ))
      )
    }

    out_name <- paste0("Covariate_", cohort)
    ret[[out_name]] <- CovarMat
    rm(list = c("CovarMat", "out_name"))
    setTxtProgressBar(pb, i)
  }

  ret
}


#' Check valid wake days for NHANES 2011-2014 accelerometer data
#'
#' @description
#' This function checks accelerometer data for valid days based on wake wear time criteria
#' and data quality flags. Adapted for NHANES 2011-2014 MIMS data format.
#'
#' @param act Activity data frame with accelerometer data including SEQN, PAXDAYWM, PAXQFM
#' and MIN1-MIN1440 columns.
#'
#' @param flags Flag data frame with wear/non-wear flags. Should have same structure as act
#' but with wear/sleep/non-wear status flags (1=wake wear, 2=sleep wear, 3=non-wear, 4=unknown).
#'
#' @param threshold_lower Minimum number of wake wear minutes required for a valid day.
#' Defaults to 600 (10 hours).
#'
#' @param rm_PAXQFM Logical indicating whether to exclude days with quality control flags.
#' Defaults to TRUE.
#'
#' @details
#' This function is specifically designed for NHANES 2011-2014 MIMS data which uses different
#' flag coding than earlier NHANES waves. Wake wear time is calculated as minutes where
#' flag == 1 (wake wear).
#'
#' @return
#' Numeric vector of indices for valid days that meet the wear time and quality criteria.
#'
#' @export
check_valid_wake_days <- function(act, flags, threshold_lower = 600, rm_PAXQFM = TRUE) {
  # Input validation
  stopifnot(all(is.data.frame(act), is.data.frame(flags)))
  stopifnot(all(colnames(act) == colnames(flags)))
  stopifnot(all(c("SEQN", "PAXDAYWM", "PAXQFM", paste0("MIN", 1:1440)) %in% colnames(act)))

  # Get activity and flag columns
  act_cols <- which(colnames(act) %in% paste0("MIN", 1:1440))

  # Verify that act and flags dataframes are identical except for activity/flag columns
  if (!identical(act[, -act_cols], flags[, -act_cols])) {
    stop("One or more columns of the act and flags do not match. Please double check that these two dataframes are identical except for the activity count and wear/non-wear columns")
  }

  # Create binary matrix for wake wear time (1 = wake wear, 0 = everything else)
  # Convert 2/3/4 to 0 and keep 1 as 1
  wake_wear_binary <- as.matrix(flags[, act_cols]) == 1

  # Calculate wake wear time (sum of periods where flag == 1)
  wake_minutes <- rowSums(wake_wear_binary, na.rm = TRUE)

  # Check PAXQFM validity
  stopifnot(all(is.finite(act$PAXQFM)))

  # Combine criteria for valid days
  flag_nonwear <- wake_minutes < threshold_lower

  # Build condition string similar to reference function
  cond <- c("flag_nonwear")
  if (rm_PAXQFM) {
    cond <- c(cond, "!(act$PAXQFM %in% 0)") # PAXQFM should be 0 for valid
  }

  cond <- paste(cond, collapse = "|")

  # Get indices of days with sufficient wake wear time
  # valid_indices <- which(wake_minutes >= threshold_lower)
  valid_indices <- eval(parse(text = paste("which(!(", cond, "))")))

  # Optional: Print summary of wear time distribution
  cat("Summary of wake wear minutes per day:\n")
  print(summary(wake_minutes))
  cat("\nNumber of valid days (>=10 hours wake wear):", length(valid_indices), "\n")

  return(valid_indices)
}

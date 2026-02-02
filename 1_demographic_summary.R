#### 0. Setup ####
library(dplyr)
library(gtsummary)
#### 1. Load Data ####
dir_path <- "./data/"
# Load the wrist-worn accelerometer data (2011-2014)
load(file.path(dir_path, "mims/data_analysis_new.RData"))
data_wrist <- data_analysis
rm(data_analysis) # Clean up to prevent confusion

# Load the hip-worn accelerometer data (2003-2006)
load(file.path(dir_path, "count/data_analysis_old.RData"))
data_hip <- data_analysis
rm(data_analysis) # Clean up

# Load activity metrics (daily level)
load(file.path(dir_path, "mims/Act_Analysis_new.RData"))
Act_Analysis_wrist <- Act_Analysis
rm(Act_Analysis)

load(file.path(dir_path, "count/Act_Analysis_old.RData"))
Act_Analysis_hip <- Act_Analysis
rm(Act_Analysis)

# Load filtered subjects from 2_preprocessing.R
# These contain only subjects who passed all quality control steps
load(file.path(dir_path, "runlength/event_analysis_new.RData"))
filtered_seqn_wrist <- unique(event_analysis$SEQN)
rm(event_analysis)

load(file.path(dir_path, "runlength/event_analysis_old.RData"))
filtered_seqn_hip <- unique(event_analysis$SEQN)
rm(event_analysis)

# Filter demographic data to include only subjects from 2_preprocessing.R
data_wrist <- data_wrist %>% filter(SEQN %in% filtered_seqn_wrist)
data_hip <- data_hip %>% filter(SEQN %in% filtered_seqn_hip)

# Filter activity data to match
Act_Analysis_wrist <- Act_Analysis_wrist %>% filter(SEQN %in% filtered_seqn_wrist)
Act_Analysis_hip <- Act_Analysis_hip %>% filter(SEQN %in% filtered_seqn_hip)

cat("Wrist cohort subjects after preprocessing:", nrow(data_wrist), "\n")
cat("Hip cohort subjects after preprocessing:", nrow(data_hip), "\n")

# Aggregate activity measures from daily to subject level
activity_summary_wrist <- Act_Analysis_wrist %>%
    group_by(SEQN) %>%
    summarise(
        TAC = mean(TAC, na.rm = TRUE),
        MVPA = mean(MVPA, na.rm = TRUE),
        Peak30 = mean(Peak30, na.rm = TRUE),
        ASTP = mean(ASTP, na.rm = TRUE),
        SBout = mean(SBout, na.rm = TRUE),
        ABout = mean(ABout, na.rm = TRUE),
        ST = mean(ST, na.rm = TRUE),
        WT = mean(WT, na.rm = TRUE),
        .groups = "drop"
    )

activity_summary_hip <- Act_Analysis_hip %>%
    group_by(SEQN) %>%
    summarise(
        TAC = mean(TAC, na.rm = TRUE),
        MVPA = mean(MVPA, na.rm = TRUE),
        Peak30 = mean(Peak30, na.rm = TRUE),
        ASTP = mean(ASTP, na.rm = TRUE),
        SBout = mean(SBout, na.rm = TRUE),
        ABout = mean(ABout, na.rm = TRUE),
        ST = mean(ST, na.rm = TRUE),
        WT = mean(WT, na.rm = TRUE),
        .groups = "drop"
    )

# Merge activity metrics with demographic data
data_wrist <- data_wrist %>%
    left_join(activity_summary_wrist, by = "SEQN")

data_hip <- data_hip %>%
    left_join(activity_summary_hip, by = "SEQN")

#### 2. Prepare Data for Summary ####
# Select the relevant columns for the demographic table
hip_summary_data <- data_hip %>%
    select(
        Age,
        Gender,
        Race,
        EducationAdult,
        BMI,
        SmokeCigs,
        DrinkStatus,
        SYS,
        CHF,
        CHD,
        Cancer,
        Stroke,
        Diabetes,
        MobilityProblem,
        # Activity metrics
        WT,
        TAC,
        MVPA,
        Peak30,
        ASTP,
        ST,
        SBout,
        ABout
    )

wrist_summary_data <- data_wrist %>%
    select(
        Age,
        Gender,
        Race,
        EducationAdult,
        BMI,
        SmokeCigs,
        DrinkStatus,
        SYS,
        CHF,
        CHD,
        Cancer,
        Stroke,
        Diabetes,
        MobilityProblem,
        # Activity metrics
        WT,
        TAC,
        MVPA,
        Peak30,
        ASTP,
        ST,
        SBout,
        ABout
    )


#### 3. Generate Demographic Summary Tables ####
# --- Table for Hip-Worn Cohort (2003-2006) ---
hip_table <- hip_summary_data %>%
    gtsummary::tbl_summary(
        label = list(
            Age ~ "Age (years)",
            Gender ~ "Gender",
            Race ~ "Race/Ethnicity",
            EducationAdult ~ "Education Level",
            BMI ~ "Body Mass Index (kg/m²)",
            SmokeCigs ~ "Smoking Status",
            DrinkStatus ~ "Drinking Status",
            SYS ~ "Systolic Blood Pressure (mmHg)",
            CHF ~ "Congestive Heart Failure",
            CHD ~ "Coronary Heart Disease",
            Cancer ~ "History of Cancer",
            Stroke ~ "History of Stroke",
            Diabetes ~ "History of Diabetes",
            MobilityProblem ~ "Mobility Problem",
            # Activity metrics
            WT ~ "Wear Time (min/day)",
            TAC ~ "Total Activity Count",
            MVPA ~ "MVPA (min/day)",
            Peak30 ~ "Peak 30-min Activity",
            ASTP ~ "Active-to-Sedentary Transition Prob.",
            ST ~ "Sedentary Time (min/day)",
            SBout ~ "Sedentary Bout Duration (min)",
            ABout ~ "Active Bout Duration (min)"
        ),
        statistic = list(
            all_continuous() ~ "{mean} ({sd})",
            all_categorical() ~ "{n} ({p}%)"
        ),
        missing = "no"
    )

hip_table
# --- Table for Wrist-Worn Cohort (2011-2014) ---
wrist_table <- wrist_summary_data %>%
    tbl_summary(
        label = list(
            Age ~ "Age (years)",
            Gender ~ "Gender",
            Race ~ "Race/Ethnicity",
            EducationAdult ~ "Education Level",
            BMI ~ "Body Mass Index (kg/m²)",
            SmokeCigs ~ "Smoking Status",
            DrinkStatus ~ "Alcohol Consumption",
            SYS ~ "Systolic Blood Pressure (mmHg)",
            CHF ~ "Congestive Heart Failure",
            CHD ~ "Coronary Heart Disease",
            Cancer ~ "History of Cancer",
            Stroke ~ "History of Stroke",
            Diabetes ~ "History of Diabetes",
            MobilityProblem ~ "Mobility Problem",
            # Activity metrics
            WT ~ "Wear Time (min/day)",
            TAC ~ "Total Activity Count",
            MVPA ~ "MVPA (min/day)",
            Peak30 ~ "Peak 30-min Activity",
            ASTP ~ "Active-to-Sedentary Transition Prob.",
            ST ~ "Sedentary Time (min/day)",
            SBout ~ "Sedentary Bout Duration (min)",
            ABout ~ "Active Bout Duration (min)"
        ),
        statistic = list(
            all_continuous() ~ "{mean} ({sd})",
            all_categorical() ~ "{n} ({p}%)"
        ),
        missing = "no"
    ) %>%
    modify_caption("Table 1b. Demographic Summary for Wrist-Worn Accelerometer Cohort (2011-2014)")
wrist_table

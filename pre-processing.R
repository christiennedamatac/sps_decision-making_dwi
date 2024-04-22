# PRE-PROCESSING SPS, ECONOMIC PREFERENCE, AND DEMOGRAPHIC DATA

setwd("/Volumes/project/3022060.01")

########################################################################
# Copy specific subject files to be analyzed
########################################################################
# At this stage, N = 413
library(fs)

# Define directory paths
working_dir <- 'data_download/data_2024-04-09'
out_dir <- 'analysis/data'

# Define required files in each subject's dir in working_dir
req_files <- c('demographics_age', 'demographics_sex', 'handedness', 'demographics_education_1', 'big-5_1', 'demographics_income_3', 'sensory-processing_3', 'online_task_ambiguity_processed_3', 'online_task_risk_processed_3', 'online_task_time_processed_3', 'online_task_trust_processed_3')

# Copy and process the MRI data later
#req_files <- c('mri_t1_1', 'mri_diffusion_1')

# List all subject directories 
subject_directories <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Ensure that output dir exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Iterate over each subject dir and look for required files
for (subject_dir in subject_directories) {
  # Find existing files in the subject dir
  existing_files <- list.files(path = subject_dir, full.names = TRUE)
  
  # Check if subject dir has any required files
  if (any(basename(req_files) %in% basename(existing_files))) {
    # Subject directory has some required files
    complete_subjects <- c(complete_subjects, basename(subject_dir))
    
    # Create subject dir in output dir if it doesn't exist
    subject_out_dir <- file.path(out_dir, basename(subject_dir))
    if (!dir.exists(subject_out_dir)) {
      dir.create(subject_out_dir, recursive = TRUE)
    }
    
    # Copy required files from subject dir to output dir if they exist
    for (file_name in req_files) {
      if (basename(file_name) %in% basename(existing_files)) {
        file.copy(from = file.path(subject_dir, file_name), 
                  to = file.path(subject_out_dir, basename(file_name)), 
                  overwrite = FALSE)
      }
    }
  } else {
    # Output missing information
    message(paste(basename(subject_dir), "missing:", paste(req_files, collapse = ", ")))
  }
}


########################################################################
# Copy whatever files are missing from HBS assessment 3
########################################################################
working_dir <- 'data_download/data_2024-04-09'
out_dir <- 'analysis/data'
files <- c('demographics_income_', 'online_task_ambiguity_processed_', 'online_task_risk_processed_', 'online_task_time_processed_')

# Function to check if a file exists in a directory
file_exists <- function(directory, file_prefix) {
  file_path <- file.path(directory, paste0(file_prefix, '.csv'))
  file.exists(file_path)
}

# List all directories within working_dir
subject_directories <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)

# Iterate over each subject directory and look for required files
for (subject_dir in subject_directories) {
  for (file in files) {
    file_3_missing <- !file_exists(subject_dir, paste0(file, '3'))
    file_2_missing <- !file_exists(subject_dir, paste0(file, '2'))
    file_1_missing <- !file_exists(subject_dir, paste0(file, '1'))
    
    if (file_3_missing) {
      if (file_2_missing && file_exists(working_dir, paste0(file, '2'))) {
        file.copy(file.path(working_dir, subject_dir, paste0(file, '2', '.csv')), file.path(out_dir, subject_dir))
      } else if (file_exists(working_dir, paste0(file, '1'))) {
        file.copy(file.path(working_dir, subject_dir, paste0(file, '1', '.csv')), file.path(out_dir, subject_dir))
      }
    }
  }
  
  # Different condition for 'big-5' files
  big_5_1_missing <- !file_exists(subject_dir, 'big-5_1')
  big_5_2_missing <- !file_exists(subject_dir, 'big-5_2')
  big_5_3_missing <- !file_exists(subject_dir, 'big-5_3')
  
  if (big_5_1_missing) {
    if (big_5_2_missing && file_exists(working_dir, 'big-5_2')) {
      file.copy(file.path(working_dir, subject_dir, 'big-5_2.csv'), file.path(out_dir, subject_dir))
    } else if (file_exists(working_dir, 'big-5_3')) {
      file.copy(file.path(working_dir, subject_dir, 'big-5_3.csv'), file.path(out_dir, subject_dir))
    }
  }
}


########################################################################
# Check that each subject dir in analysis dir has all 11 required files 
########################################################################
req_files <- c('demographics_age', 'demographics_sex', 'handedness', 'demographics_education_', 'big-5_', 'demographics_income_', 'sensory-processing_3', 'online_task_ambiguity_processed_', 'online_task_risk_processed_', 'online_task_time_processed_', 'online_task_trust_processed_')

working_dir <- 'analysis/data'

# List all directories within working_dir
subject_directories <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Function to check if all required files exist in a directory
check_files_exist <- function(subject_dir, req_files) {
  # Extract directory name
  dir_name <- basename(subject_dir)
  
  # Initialize a vector to store missing files
  missing_files <- character(0)
  
  # Iterate over required files
  for (req_file in req_files) {
    # Find files in the directory that partially match the required file
    matching_files <- list.files(subject_dir, pattern = req_file)
    
    # If no matching files found, add to missing files vector
    if (length(matching_files) == 0) {
      missing_files <- c(missing_files, req_file)
    }
  }
  
  # Print the status for the directory if it has missing files
  if (length(missing_files) > 0) {
    cat("Missing files in", dir_name, ":", paste(missing_files, collapse = ", "), "\n")
  }
}

# Iterate over each subject directory and look for required files
for (subject_dir in subject_directories) {
  check_files_exist(subject_dir, req_files)
}


########################################################################
# Copy and paste files missing from HBS assessment 3
########################################################################
working_dir <- 'data_download/data_2024-04-09'
out_dir <- 'analysis/data'

# 4 outcomes missing
ambiguity_risk_time_trust <- c('HBU1C7F00B8E3B78579C', 'HBU70C31D20413EDD686', 'HBU80C014278320BC5AF', 'HBUACC40BBEF7AF5931C', 'HBUACD35D4A769D07EFE', 'HBUF037A68F610E1D3EB')

# 3 outcomes missing 
ambiguity_risk_time <- c('HBU74100672298DBC45C', 'HBU9E913DD1A608951D7')

# income missing
demographics_income <- c('HBU16D020EB3E6DC2E1E', 'HBU6264BF7DC3F9EDB81', 'HBU70C31D20413EDD686', 'HBU90466FD5B5CAE14B4', 'HBUB40FAE2450C1438A5', 'HBUD2E4AC96D3E67CABE', 'HBUEAA43518E1553E346', 'HBUF8762898381D3DA02', 'HBU80C014278320BC5AF')

# missing files
missing_files_4 <- c('online_task_ambiguity_processed_2', 'online_task_risk_processed_2', 'online_task_time_processed_2', 'online_task_trust_processed_2')
missing_files_3 <- c('online_task_ambiguity_processed_2', 'online_task_risk_processed_2', 'online_task_time_processed_2')

# copy files from working_dir into out_dir
for (subject in ambiguity_risk_time_trust) {
  for (file in missing_files_4) {
    file_path <- paste0(working_dir, '/', subject, '/', file)
    out_path <- paste0(out_dir, '/', subject, '/', file)
    file.copy(from = file_path, to = out_path, overwrite = FALSE)
  }
}

for (subject in ambiguity_risk_time) {
  for (file in missing_files_3) {
    file_path <- paste0(working_dir, '/', subject, '/', file)
    out_path <- paste0(out_dir, '/', subject, '/', file)
    file.copy(from = file_path, to = out_path, overwrite = FALSE)
  }
}

for (subject in demographics_income) {
  file <- 'demographics_income_2'
  file_path <- paste0(working_dir, '/', subject, '/', file)
  out_path <- paste0(out_dir, '/', subject, '/', file)
  file.copy(from = file_path, to = out_path, overwrite = FALSE)
}


########################################################################
# Exclude subject who are missing data from main variables
########################################################################
# 3 subjects with no online economic preference task data
exclude_list <- c('HBU7A944716853CD1C65', 'HBUA496755A17B382863', 'HBUAC444923D61B631E5')
working_dir <- 'analysis/data'

for (subject in exclude_list) {
  dir_path <- file.path(working_dir, subject)
  if (file.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)  # Remove directory and subdirectories
    cat(paste("Directory", dir_path, "has been removed.\n"))
  } else {
    cat(paste("Directory", dir_path, "does not exist.\n"))
  }
}

# After this exclusion, N = 410 subjects remain


########################################################################
# Rename subject directories with an ordered prefix for ease (sub-0000)
########################################################################
library(fs)

working_dir <- 'analysis/data'

# List all directories within the main directory
subject_directories <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)
#subject_directories <- ('analysis/data/HBU0A6F67FF578F01743')

# Function to generate subject prefix
generate_subject_prefix <- function(i) {
  sprintf("sub-%04d_", i)
}

# Iterate over each directory and rename them with the prefix
for (i in seq_along(subject_directories)) {
  old_dir <- subject_directories[i]
  new_dir <- file.path(dirname(old_dir), paste0(generate_subject_prefix(i), basename(old_dir)))
  file.rename(old_dir, new_dir)
}


########################################################################
# Add .zip extension, unzip, and delete zip files in each subject's directory
########################################################################

working_dir <- 'analysis/data'

# List all directories within the main directory
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

for (subject in subjects) {
  # Add .zip extension to each file in the subject dir
  files <- list.files(subject, full.names = TRUE)
  zip_files <- paste0(files, ".zip")
  file.rename(files, zip_files)
  
  # Unzip each file
  for (zip_file in zip_files) {
    unzip(zip_file, exdir = subject)
  }
  
  # Delete the .zip file
  file.remove(zip_files)
}


########################################################################
# Concatenate files across variables per subject
########################################################################
library(tidyverse)

# Define the working directory and output directory
working_dir <- 'analysis/data'

subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

for (subject in subjects) {
  # Initialize empty lists to store data
  tsk_ambiguity_data <- list()
  tsk_trust_A_data <- list()
  qst_post_3_SPS_data <- list()
  qst_pre_age_data <- list()
  tsk_risk_data <- list()
  tsk_trust_B_data <- list()
  LIS_education_data <- list()
  qst_pre_handedness_data <- list()
  tsk_time_data <- list()
  BIG_data <- list()
  LIS_income_data <- list()
  qst_pre_sex_data <- list()
  
  # Read and concatenate files for each type of data
  tryCatch({
    tsk_ambiguity_data <- map(tsk_ambiguity, ~ read_csv(.x) %>% select(m50, m10, m90))
  }, error = function(e) {})
  
  tryCatch({
    tsk_trust_A_data <- map(tsk_trust_A, ~ read_csv(.x) %>% select(tg_sent))
  }, error = function(e) {})
  
  tryCatch({
    qst_post_3_SPS_data <- map(qst_post_3_SPS, ~ read_csv(.x) %>% select(SPS_sensibility_internal, SPS_emotional_physiological_reactivity, SPS_sensory_discomfort, SPS_sensory_comfort, SPS_social_affective_sensitivity, SPS_esthetic_sensitivity))
  }, error = function(e) {})
  
  tryCatch({
    qst_pre_age_data <- map(qst_pre_age, ~ read_csv(.x) %>% select(Age))
  }, error = function(e) {})
  
  tryCatch({
    tsk_risk_data <- map(tsk_risk, ~ read_csv(.x) %>% select(riskaversion, prudence))
  }, error = function(e) {})
  
  tryCatch({
    tsk_trust_B_data <- map(tsk_trust_B, ~ read_csv(.x) %>% select(starts_with('tg_return')))
  }, error = function(e) {})
  
  tryCatch({
    LIS_education_data <- map(LIS_education, ~ read_csv(.x) %>% select(LIS07_education_completed))
  }, error = function(e) {})
  
  tryCatch({
    qst_pre_handedness_data <- map(qst_pre_handedness, ~ read_csv(.x) %>% select(VAR101_left_right_handed))
  }, error = function(e) {})
  
  tryCatch({
    tsk_time_data <- map(tsk_time, ~ read_csv(.x) %>% select(patience_outcome))
  }, error = function(e) {})
  
  tryCatch({
    BIG_data <- map(BIG, ~ read_csv(.x) %>% select(BIG_openness_sum, BIG_neuroticism_sum))
  }, error = function(e) {})
  
  tryCatch({
    LIS_income_data <- map(LIS_income, ~ read_csv(.x) %>% select(LIS23_net_household_income))
  }, error = function(e) {})
  
  tryCatch({
    qst_pre_sex_data <- map(qst_pre_sex, ~ read_csv(.x) %>% select(Sex))
  }, error = function(e) {})
  
  # Combine all data
  combined_data <- bind_cols(subject = tsk_ambiguity_data, tsk_trust_A_data, qst_post_3_SPS_data, qst_pre_age_data, tsk_risk_data, tsk_trust_B_data, LIS_education_data, qst_pre_handedness_data, tsk_time_data, BIG_data, LIS_income_data, qst_pre_sex_data)
  
  # Output combined data to CSV
  write_csv(combined_data, file.path(subject, 'concatenated_data.csv'))
}


########################################################################
# Find any subjects whose data files did not concatenate
########################################################################

working_dir <- 'analysis/data'

# List all sub-directories (subjects) in the working directory
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Function to check if a file has less than 1 row of data (does not include header row)
has_less_than_two_rows <- function(file_path) {
  data <- read.csv(file_path)
  return(nrow(data) < 1)
}

# Loop through each subject directory
for (subject in subjects) {
  # Get the file path for concatenated_data.csv
  file_path <- file.path(subject, "concatenated_data.csv")
  
  # Check if the file has more than 2 rows
  if (file.exists(file_path) && has_less_than_two_rows(file_path)) {
    cat(basename(subject),", ")
  }
}


########################################################################
# Insert subject number column 
########################################################################
working_dir <- 'analysis/data'

# List all sub-directories (subjects) in the working directory
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Specify the starting subject number
starting_subject <- 'sub-0054'

# Loop through each subject directory
for (subject in subjects) {
  # Check if the subject directory meets the starting condition
  if (basename(subject) >= starting_subject) {
    # Extract subject number from directory name: everything before the first '_'
    subject_num <- sub("_.*", "", basename(subject))
    
    # Read the concatenated data file for the subject
    data <- read.csv(file.path(subject, 'concatenated_data.csv'))
    
    # Insert a first column called 'subject' with the subject number
    data <- cbind(subject = subject_num, data)
    
    # Write the modified data back to the same file
    write.csv(data, file.path(subject, 'concatenated_data.csv'), row.names = FALSE)
  }
}


########################################################################
# Concatenate files across subjects
########################################################################
library(dplyr)
library(readr)

working_dir <- 'analysis/data'
out_dir <- 'analysis'
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Define the column names for merging
columns <- c("subject",
             "m50", "m10", "m90",
             "tg_sent", 
             "SPS_sensibility_internal", "SPS_emotional_physiological_reactivity", "SPS_sensory_discomfort", "SPS_sensory_comfort", "SPS_social_affective_sensitivity", "SPS_esthetic_sensitivity", 
             "Age", 
             "riskaversion", "prudence", 
             "tg_return00", "tg_return05","tg_return10", "tg_return15", "tg_return20", "tg_return25", "tg_return30", "tg_return35","tg_return40","tg_return45", "tg_return50", 
             "LIS07_education_completed", 
             "VAR101_left_right_handed", 
             "patience_outcome", 
             "BIG_openness_sum","BIG_neuroticism_sum", 
             "LIS23_net_household_income",
             "Sex")

concatenated_all <- NULL

for (i in 1:length(subjects)) {
  subject <- subjects[i]
  print(paste("Processing:", subject))
  data <- read_csv(file.path(subject, 'concatenated_data.csv'))
  
  if (is.null(concatenated_all)) {
    concatenated_all <- data
  } else {
    # Identify common columns for merging
    common_cols <- intersect(names(concatenated_all), names(data))
    
    # Merge data frames using common columns
    concatenated_all <- merge(concatenated_all, data, by = common_cols, all = TRUE)
  }
}

write_csv(concatenated_all, file.path(out_dir, 'concatenated_all.csv'))


########################################################################
# Add another column that I missed earlier: tg_sent
########################################################################
library(dplyr)

data <- concatenated_all

# Create a new column called 'tg_sent'
data$tg_sent <- NA

# Define working directory and list subjects
working_dir <- 'analysis/data'
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Loop through each subject
for (subject in subjects) {
  subject_num <- sub("_.*", "", basename(subject))
  
  # Find files matching the pattern
  subject_files <- list.files(path = subject, pattern = "*_tsk_trust_A_processed.csv", full.names = TRUE)
  
  # Ensure only one file is found, assuming one file per subject with this pattern
  if (length(subject_files) != 1) {
    warning(paste("Expected one file matching pattern for subject", subject_num, "but found", length(subject_files), "files. Skipping this subject."))
    next
  }
  
  # Read the data for the current subject
  subject_data <- read.csv(subject_files[1])
  
  # Get the value under the column called 'tg_sent'
  tg_sent_value <- subject_data$tg_sent
  
  # Put the value into 'data' under the column 'tg_sent' matching the row 'subject_num' in the column 'subject'
  data[data$subject == subject_num, "tg_sent"] <- tg_sent_value
}

write_csv(data, file.path(out_dir, 'concatenated_all.csv'))


########################################################################
# Calculate ambiguity aversion
########################################################################
library(AlgDesign)
library(readr)

# Load data
concatenated_all <- read_csv("analysis/concatenated_all.csv")
data <- concatenated_all

# Initialize variables
initialize_variables <- function(data) {
  # Initialize matrices to store indexes
  AA_0.1 <- numeric(length(data$subject))
  AA_0.5 <- numeric(length(data$subject))
  AA_0.9 <- numeric(length(data$subject))
  ambiguity_aversion <- numeric(length(data$subject))
  a_insensitivity <- numeric(length(data$subject))
  
  return(list(AA_0.1 = AA_0.1, AA_0.5 = AA_0.5, AA_0.9 = AA_0.9, ambiguity_aversion = ambiguity_aversion, a_insensitivity = a_insensitivity))
}

# Test initialize_variables function
initial_values <- initialize_variables(data)
print(head(initial_values$AA_0.1))
print(head(initial_values$AA_0.5))
print(head(initial_values$AA_0.9))
print(head(initial_values$ambiguity_aversion))
print(head(initial_values$a_insensitivity))

# Calculate scaled matching probabilities
calculate_scaled_probabilities <- function(data) {
  scaled_m10 <- vector("list", length(data$subject))
  scaled_m50 <- vector("list", length(data$subject))
  scaled_m90 <- vector("list", length(data$subject))
  
  for (i in 1:nrow(data)) {
    mp <- c(data$m10[i], data$m50[i], data$m90[i]) / 100  # Scaling to [0, 1]
    mp <- na.omit(mp)  # Remove NA values
    
    scaled_m10[[i]] <- ifelse(length(mp) >= 1, mp[1], NA)
    scaled_m50[[i]] <- ifelse(length(mp) >= 2, mp[2], NA)
    scaled_m90[[i]] <- ifelse(length(mp) >= 3, mp[3], NA)
  }
  
  data$scaled_m10 <- unlist(scaled_m10)
  data$scaled_m50 <- unlist(scaled_m50)
  data$scaled_m90 <- unlist(scaled_m90)
  
  return(data)
}

# Test calculate_scaled_probabilities function
data_with_scaled_probabilities <- calculate_scaled_probabilities(data)
print(head(data_with_scaled_probabilities$scaled_m10))
print(head(data_with_scaled_probabilities$scaled_m50))
print(head(data_with_scaled_probabilities$scaled_m90))

# Calculate ambiguity aversion indexes
calculate_AA_indexes <- function(data) {
  # Initialize AA indexes columns
  data$AA_0.1 <- numeric(nrow(data))
  data$AA_0.5 <- numeric(nrow(data))
  data$AA_0.9 <- numeric(nrow(data))
  
  for (i in 1:nrow(data)) {
    mp <- c(data$m10[i], data$m50[i], data$m90[i]) / 100  # Scaling to [0, 1]
    mp <- na.omit(mp)  # Remove NA values
    
    # Calculate AA indexes
    data$AA_0.1[i] <- 0.1 - mp[1]
    data$AA_0.5[i] <- 0.5 - mp[2]
    data$AA_0.9[i] <- 0.9 - mp[3]
  }
  
  return(data)
}

# Test calculate_AA_indexes function
data_with_AA_indexes <- calculate_AA_indexes(data)
print(head(data_with_AA_indexes$AA_0.1))
print(head(data_with_AA_indexes$AA_0.5))
print(head(data_with_AA_indexes$AA_0.9))

# Calculate ambiguity aversion
calculate_ambiguity_aversion <- function(data) {
  for (i in 1:nrow(data)) {
    mp <- c(data$m10[i], data$m50[i], data$m90[i]) / 100  # Scaling to [0, 1]
    mp <- na.omit(mp)  # Remove NA values
    
    lin_fit <- lm(mp ~ I(c(0.1, 0.5, 0.9)), data = data.frame(p = c(0.1, 0.5, 0.9), mp = mp), na.action = na.exclude)
    c_coef <- coef(lin_fit)[1]  # Intercept
    s_coef <- coef(lin_fit)[2]  # Slope
    data$ambiguity_aversion[i] <- 1 - s_coef - 2 * c_coef
  }
  
  return(data)
}

# Test calculate_ambiguity_aversion function
data_with_ambiguity_aversion <- calculate_ambiguity_aversion(data)
print(head(data_with_ambiguity_aversion$ambiguity_aversion))

# Calculate a-insensitivity
calculate_a_insensitivity <- function(data) {
  for (i in 1:nrow(data)) {
    mp <- c(data$m10[i], data$m50[i], data$m90[i]) / 100  # Scaling to [0, 1]
    mp <- na.omit(mp)  # Remove NA values
    
    quad_fit <- lm(mp ~ poly(c(0.1, 0.5, 0.9), 2), data = data.frame(p = c(0.1, 0.5, 0.9), mp = mp), na.action = na.exclude)
    s_quad <- coef(quad_fit)[2]  # Slope
    data$a_insensitivity[i] <- 1 - s_quad
  }
  
  return(data)
}

# Test calculate_a_insensitivity function
data_with_a_insensitivity <- calculate_a_insensitivity(data)
print(head(data_with_a_insensitivity$a_insensitivity))

# Apply the functions to the data
data <- calculate_scaled_probabilities(data)
data <- calculate_AA_indexes(data)
data <- calculate_ambiguity_aversion(data)
data <- calculate_a_insensitivity(data)

# Save 
write_csv(data, file.path('analysis/concatenated_all.csv'))

# Columns to summarize
ambiguity_columns <- c("scaled_m10", "scaled_m50", "scaled_m90", "aa_0.1", "aa_0.5", "aa_0.9", "ambiguity_aversion", "a_insensitivity")

# Initialize an empty dataframe to store summary statistics
summary_df <- data.frame(Variable = character(length(ambiguity_columns)), Mean = numeric(length(ambiguity_columns)), Median = numeric(length(ambiguity_columns)), SD = numeric(length(ambiguity_columns)), Min = numeric(length(ambiguity_columns)), Max = numeric(length(ambiguity_columns)))

# Loop through each column
for (i in seq_along(ambiguity_columns)) {
  column <- ambiguity_columns[i]
  summary_df[i, "Variable"] <- column
  summary_df[i, c("Mean", "Median", "SD", "Min", "Max")] <- c(
    round(mean(data[[column]], na.rm = TRUE), 2),
    round(median(data[[column]], na.rm = TRUE), 2),
    round(sd(data[[column]], na.rm = TRUE), 2),
    round(min(data[[column]], na.rm = TRUE), 2),
    round(max(data[[column]], na.rm = TRUE), 2)
  )
}

# Save
write_csv(summary_df, file.path('analysis/ambiguity_attitude_indexes_stats.csv'))


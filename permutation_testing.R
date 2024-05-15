# PERMUTATION TESTING FOR PLS-PM

library(readr)
library(dplyr)

# Read the CSV file
setwd('/Volumes/project/3022060.01')
print("Reading CSV file...")
concatenated_all_z_reg_neuroticism_openness <- read_csv("analysis/concatenated_all_z_reg_neuroticism_openness.csv")

# Subset to measured variables for the predictor latent variable
data_predictor <- select(concatenated_all_z_reg_neuroticism_openness, subject, SPS_positive, SPS_negative)

# Subset to measured variables for the outcome latent variable
data_outcome <- select(concatenated_all_z_reg_neuroticism_openness, subject, prudence, riskaversion, ambiguity_aversion, tg_sent, tg_return_sum, patience_outcome)

# Drop NA values
data <- na.omit(concatenated_all_z_reg_neuroticism_openness)
data_predictor <- na.omit(data_predictor)
data_outcome <- na.omit(data_outcome)

####################################
# Generate permutations
####################################

# Specify the number of permutations
num_permutations <- 10000

# Initialize list to store permuted dataframes
perm_data <- list()

# Preserve the order of subjects; permutation occurs independently within each column while maintaining subject order
print("Permuting SPSQ dimension score values...")
for (i in 1:num_permutations) {
  print(paste("Permuting data matrix:", i))
  perm_data[[i]] <- data_predictor %>%
    mutate(across(-1, ~ sample(.)))  # Sample values within each column
}

# Save permuted dataframes as CSV files
#print("Saving permuted dataframes...")
#dir.create("analysis/generate_permutations", showWarnings = FALSE)  # Create directory if it doesn't exist
#for (i in 1:num_permutations) {
#  write.csv(perm_data[[i]], file = paste0("analysis/generate_permutations/permuted_data_", i, ".csv"), row.names = FALSE)}
#print("Permutations saved successfully.")

####################################
# Function to run PLS-PM
####################################

library(plspm)

run_pls_pm <- function(data) {
  # INNER MODEL
  sps <- c(0, 0)
  decision_making <- c(1, 0)
  path <- rbind(sps, decision_making)
  colnames(path) <- rownames(path)
  
  # OUTER MODEL
  blocks <- list(c('SPS_negative', 'SPS_positive'),
                 c('prudence', 'riskaversion',
                   'ambiguity_aversion', "tg_sent",
                   "tg_return_sum", 'patience_outcome'))
  
  # SCALING
  scaling <- list(c('num', 'num'),
                  c('ord', 'ord', 'num', 'ord', 'ord', 'ord'))
  
  # MODES
  modes <- c('A', 'A')
  
  # RUN PLS-PM MODEL
  pls_result <- plspm(data,
                      path_matrix = path,
                      blocks = blocks,
                      scaling = scaling,
                      modes = modes,
                      scheme = 'factor') 
  
  return(pls_result)
}

####################################
# MAIN PLS_PM 
####################################

pls_main_output <- run_pls_pm(data)

#############################################################################
# Run PLS-PM on each permutation to generate a null distribution of estimates
#############################################################################

# Initialize a list to store t-values from permutations
t_values_permutations <- list()

# Loop through each permutation
for (i in 635:num_permutations) {
  print(paste("Running PLS-PM on permuted dataframe:", i))
  
  # Initialize t_value_perm as NA in case of an error
  t_value_perm <- NA
  
  # Attempt to run PLS-PM on permuted data
  tryCatch({
    # Merge permuted dataframe with 'data_outcome' based on 'subject'
    merged_data <- merge(perm_data[[i]], data_outcome, by = "subject", all = TRUE)
    
    # Perform PLS-PM on permuted data
    pls_pm_result <- run_pls_pm(merged_data)
    
    # Get t-value from permutation PLS-PM output
    t_value_perm <- pls_pm_result$inner_model$decision_making[2, "t value"]
  }, error = function(e) {
    # Print error message if PLS-PM fails to converge
    cat("Error in permutation", i, ": ", conditionMessage(e), "\n")
  })
  
  # Append t-value to the list if successful
  if (!is.na(t_value_perm)) {
    t_values_permutations[[i]] <- t_value_perm
  }
}

# Save list of t-values to a text file
output_file <- "analysis/plspm/permutation_t_values.txt"
writeLines(as.character(unlist(t_values_permutations)), con = output_file)

# Calculate corrected p-value based on main PLS-PM output
original_p_value <- pls_main_output$inner_model$decision_making[2, "Pr(>|t|)"]
t_value_main <- pls_main_output$inner_model$decision_making[2, "t value"]  
num_significant <- sum(abs(t_value_main) >= abs(t_value_perm))
p_value_corrected <- (num_significant + 1) / (num_permutations + 1)

# Print the corrected p-value
print("Corrected p-value:")
print(p_value_corrected)

# Save the output to a text file
output_file <- "analysis/plspm/permutation_testing.txt"
cat("t:", t_value_main, "\n", file = output_file)
cat("n_permutations:", num_permutations, "\n", file = output_file)
cat("p_uncorrected:", original_p_value, "\n", file = output_file)
cat("p_corrected:", p_value_corrected, "\n", file = output_file)

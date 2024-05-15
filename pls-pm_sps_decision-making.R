# Partial least squares path modeling for SPS and decision-making

library(readr)

setwd('/Volumes/project/3022060.01')
concatenated_all_z_reg_neuroticism_openness <- read_csv("analysis/concatenated_all_z_reg_neuroticism_openness.csv")
data <- na.omit(concatenated_all_z_reg_neuroticism_openness)

####################################################
# PLS-PM
####################################################
library(plspm)
library(readr)

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
pls = plspm(data,
            path_matrix = path,
            blocks = blocks,
            scaling = scaling,
            modes = modes,
            scheme = 'factor',
            boot.val = TRUE, br = 1000) # SEs were stable at 1000 bootstrap samples; did not differ from 2000 samples; no need to sample more

# VIEW RESULTS
summary_pls <- summary(pls)
print(summary_pls)

# SAVE
file_path <- "analysis/plspm/aim1_1000.txt"
file_conn <- file(file_path, "w")
writeLines(capture.output(summary_pls), file_conn)
close(file_conn)

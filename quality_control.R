# QUALITY CONTROL

########################################################################
# Questionnaire Response Keys
########################################################################

# Highest level of education completed (LIS07_education_completed)
# 1 = Basisonderwijs = primary education
# 2 = VMBO = pre-vocational secondary education
# 3 = HAVO/VWO = senior general seconday education/pre-university secondary education 
# 4 = MBO = secondary vocational education
# 5 = HBO = higher professional education
# 6 = WO = academic education
# 7 = Other
# 9 = Did not receive any education

# Handedness (VAR101_left_right_handed)
# 1 = left-handed
# 2 = right-handed

# Net monthly income for the complete household, including income from work, allowances, alimony, and benefits (LIS23_net_household_income)
# 0 = no income
# 1 = 500 euros or less
# 2 = 501 – 1000 euros
# 3 = 1001 – 1500 euros
# 4 = 1501 – 2000 euros
# 5 = 2001 – 2500 euros
# 6 = 2501 – 3000 euros
# 7 = 3001 – 3500 euros
# 8 = 3501 – 4000 euros
# 9 = 4001 – 4500 euros
# 10 = 4501 – 5000 euros
# 11 = 5001 – 7500 euros
# 12 = more than 7500 euros
# 13 = I don’t know
# 14 = I prefer not to answer

# Sex (sex)
# 1 = male
# 2 = female

library(readr)

setwd('/Volumes/project/3022060.01')
concatenated_all <- read_csv("analysis/concatenated_all.csv")
data <- concatenated_all
names(data)

########################################################################
# Histograms for categorical variables
########################################################################
library(ggplot2)
library(cowplot)

columns_categorical <- c("LIS07_education_completed", 
                         "VAR101_left_right_handed", 
                         "LIS23_net_household_income",
                         "sex")

# Create an empty list to store plots
plot_list <- list()

# Loop through each column
for (column in columns_categorical) {
  # Calculate the number of bins
  min_val <- min(data[[column]], na.rm = TRUE)
  max_val <- max(data[[column]], na.rm = TRUE)
  num_bins <- max_val - min_val + 1  # Number of integer values between min and max, including min and max
  
  # Create histogram using ggplot2
  p <- ggplot(data, aes(x = !!rlang::sym(column))) + 
    geom_histogram(binwidth = 1, fill = "#c9c9c9", color = "#000000") +  # Set binwidth to 1 to ensure each integer value has its own bin
    scale_x_continuous(breaks = seq(min_val, max_val, by = 1)) +  # Set breaks to show every integer value
    ggtitle(column) +
    theme_classic() +
    theme(plot.margin = margin(10, 10, 10, 10, "pt"),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))  # Adjust plot margins
  
  # Add the plot to the list
  plot_list[[length(plot_list) + 1]] <- p
}

# Arrange plots in a grid
#combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
#print(combined_plot)

# Save the plot to a PDF file
#ggsave("plots/histograms_categorical_variables.pdf", combined_plot, width = 8.27, height = 11.69)  # A4 size in inches

########################################################################
# Violin, box, and scatterplots for continuous variables
########################################################################
columns_continuous <- c("SPS_positive", "SPS_negative", "SPS_sum",
                        "ambiguity_aversion",
                        "tg_sent", # trust
                        "tg_return_sum", # trustworthiness
                        "riskaversion", "prudence", 
                        "patience_outcome", 
                        "BIG_openness_sum","BIG_neuroticism_sum", 
                        "age")

# Create an empty list to store plots
#plot_list <- list()

# Loop through each column
for (column in columns_continuous) {
  # Calculate mean
  mean_val <- mean(data[[column]], na.rm = TRUE)
  
  # Create violin plot
  p <- ggplot(data, aes(x = 1, y = !!rlang::sym(column))) + 
    geom_jitter(width = 0.2, alpha = 0.2, size = 1, color='#2389CA') +  
    geom_violin(fill = "transparent", color='#000000') +
    geom_boxplot(width = 0.1, fill = "transparent", color = "#000000") +  
    geom_hline(yintercept = mean_val, color = "darkred", linetype = "dashed", size = 0.75) +  # Add mean line
    ggtitle(column) +
    theme_classic() +
    theme(plot.margin = margin(10, 10, 10, 10, "pt"), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)) 
  
  # Add the plot to the list
  plot_list[[length(plot_list) + 1]] <- p
}

# Arrange plots in a grid
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
print(combined_plot)

# Save the plot to a PDF file
ggsave("plots/qc_plots.pdf", combined_plot, width = 11, height = 11)  


########################################################################
# Print summary stats
########################################################################
# Create an empty data frame to store summary statistics
summary_df <- data.frame(Column = character(),
                         Min = numeric(),
                         Q1 = numeric(),
                         Median = numeric(),
                         Mean = numeric(),
                         Q3 = numeric(),
                         Max = numeric(),
                         SD = numeric(),
                         stringsAsFactors = FALSE)

# Iterate over continuous columns to calculate summary statistics
for (column in columns_continuous) {
  # Calculate summary statistics
  column_summary <- summary(data[[column]], na.rm = TRUE)
  column_summary <- round(column_summary, 2)  # Round to 2 decimal points
  column_sd <- round(sd(data[[column]], na.rm = TRUE), 2)  # Round to 2 decimal points
  
  # Append summary statistics to the data frame
  summary_df <- rbind(summary_df, c(column, column_summary, column_sd))
}

# Assign column names to the data frame
colnames(summary_df) <- c("Variable", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD")

# Print the summary data frame
print(summary_df)

# Export the summary data frame to a CSV file
write.csv(summary_df, file = "analysis/summary_statistics.csv", row.names = FALSE)


########################################################################
# No removed outliers
########################################################################
# Variables are all within a certain range or even ordinal such that excluding outliers would be arbitrary

########################################################################
# No transformations
########################################################################
# Using non-parametric methods or robust SEM techniques such as bootstrapping or permutation tests do not rely on assumptions about the distribution of the data and can be more robust to violations of normality. Robust SEM techniques that are designed to handle violations of distributional assumptions: diagonally weighted least squares (DWLS) estimator or the robust maximum likelihood estimator (MLR)


########################################################################
# Standardization: z-score across all variables
########################################################################
library(dplyr)

data <- concatenated_all  

# Function to calculate z-score
calculate_zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Z-score all columns except 'subject'
zscored_data <- data %>%
  mutate(across(-subject, calculate_zscore))
head(zscored_data)

# Save
write.csv(zscored_data, "analysis/concatenated_all_z.csv", row.names = FALSE)


########################################################################
# Exploratory factor analysis: check underlying structure
########################################################################
library(readr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggpubr)
library(grid)
library(rstatix)
library(data.table)
library(ggcorrplot)
library(igraph)
library(ggraph)
library(tidygraph)
library(viridis)
library(Hmisc)

# Function to grab just the legend from a plot 
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

concatenated_all_z <- read_csv("analysis/concatenated_all_z.csv")
data <- as.data.table(concatenated_all_z)
names(data)

# Rename columns
setnames(data, old = c('SPS_sum', 'SPS_positive', 'SPS_negative', 'ambiguity_aversion', 'tg_sent', 'tg_return_sum', 'riskaversion', 'prudence', 'patience_outcome', 'sex', 'age', 'BIG_neuroticism_sum', 'BIG_openness_sum'), 
         new = c('SPS sum', 'SPS positive dimension', 'SPS negative dimension', 'ambiguity aversion index', 'trust', 'trustworthiness', 'risk aversion', 'prudence', 'patience', 'sex (M-, F+)', 'age', 'neuroticism', 'openness')) 
names(data)

# Define columns to be selected
selected_columns <- c('SPS sum', 'SPS positive dimension', 'SPS negative dimension', 'ambiguity aversion index', 'trust', 'trustworthiness', 'risk aversion', 'prudence', 'patience', 'sex (M-, F+)', 'age', 'neuroticism', 'openness')

# Subset data to include only the selected columns using .SD
data_subset <- data[, .SD, .SDcols = selected_columns]
names(data_subset)

# Convert data subset to numeric matrix
data_subset_matrix <- as.matrix(data_subset)

# Compute correlation matrix between selected variables
cor.mat <- cor(data_subset_matrix, method = "spearman", use = "pairwise.complete.obs")

# Compute p-values for correlation matrix
p.mat <- rcorr(as.matrix(data_subset_matrix))$P

# Create correlation plot
pi <- ggcorrplot(cor.mat, 
                 method = "square",
                 type = "lower",
                 outline.color = "#ffffff",
                 lab = TRUE,
                 digits = 2,
                 colors = c("#0048A0","#ffffff","#C10000"), # order: low, middle, high
                 tl.col = "#000000",
                 legend.title = "Spearman\ncorrelation",
                 hc.order = FALSE,
                 theme(legend.key.width = unit(2.5, 'cm'),
                       legend.position = "bottom"))

p1 <- ggcorrplot(cor.mat,
                 p.mat = p.mat,
                 method = "square",
                 type = "lower",
                 outline.color = "#ffffff",
                 lab = TRUE,
                 digits = 2,
                 colors = c("#0048A0","#ffffff","#C10000"), # order: low, middle, high
                 tl.col = "#000000",
                 legend.title = "Spearman\ncorrelation",
                 hc.order = FALSE,
                 lab_size = 3) 

# Pick which legend from which plot 
legend <- get_legend(pi)

# Remove legend from main plot 
p1 <- p1 + theme(legend.position="none") 

########################################################################
# Network analysis 
########################################################################
cor_mat <- as.data.table(cor.mat)
p_mat <- as.data.table(p.mat)

# Vertices (or nodes): name of each variable
va <- names(cor.mat)

# Edges: correlation value (r) and significance (p)
## melt data
id.vars <- names(cor.mat[,c(1)])
measure.vars <- names(cor.mat[,c(1:13)])

cor_melt <- melt(cor.mat, id.vars=c(id.vars), measure.vars=c(measure.vars))
head(cor_melt)

## melt p-values
p_melt <- melt(p.mat, id.vars=c(id.vars), measure.vars=c(measure.vars))
head(p_melt)

# merge vertices and edges together 
ed <- merge(cor_melt, p_melt[,c("Var1","Var2","value")], by=c("Var1","Var2"), all=TRUE) # merge correlations and p-values
head(ed)
colnames(ed) <- c('rowname','variable','r','p') # rename

## delete nonsignificant correlations
ed_sig <- ed[!(ed$p > 0.05),]

# generate graph object
ig <- igraph::graph_from_data_frame(d=ed_sig, vertices=va, directed = FALSE)

# add labels to nodes
tg <- tidygraph::as_tbl_graph(ig) %>% # graph object
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(label=c('SPS\nsum','SPS\nsum', 'SPS\npositive\ndimension', 'SPS\nnegative\ndimension', 'ambiguity\naversion\nindex', 'trust', 'trustworthiness', 'risk\naversion', 'prudence', 'patience', 'sex\n(M-, F+)', 'age', 'neuroticism', 'openness'))

# add questionnaire membership for each vertex 
#V(tg)$questionnaire <- subscale_totals$questionnaire

tg #tidygraph object with node and edge data

# Calculate degree centrality of this undirected network
degree.cent <- centr_degree(tg, mode = "all")
degree.cent$res

# Calculate closeness centrality
closeness.cent <- closeness(tg, mode="all")
closeness.cent

# Plot in a network
library(ggraph)
library(ggplot2)

p2 <- ggraph(tg, layout = 'stress', circular = FALSE) +
  geom_edge_arc(lineend = 'butt', linejoin = 'round', 
                linemitre = 2, 
                strength = 0,
                edge_width = 0.5,
                aes(colour = r)) +
  geom_node_point(size = 7,
                  alpha = 0.6) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(2, "lines"), 
                 size=4, 
                 colour="#000000") +
  theme_graph(background = "white") +
  theme(legend.position = "right") +
  guides(edge_width = 'none',
         edge_alpha = 'none') +
  scale_edge_colour_gradient2(
    'Spearman r',
    low = "#0048A0",
    mid = "#ffffff",
    high = "#C10000",
    midpoint = 0,
    space = "Lab",
    na.value = "#000000",
    guide = "edge_colourbar",
    aesthetics = "edge_colour",
    limits = c(-1, 1)
  )

# Arrange plot and legend vertically
combined_plot <- grid.arrange(p1, legend, p2,
                              ncol = 1, nrow=3,
                              heights = c(1, 0.1, 1))

# Save
ggsave("plots/corr_network_plot.pdf", plot = combined_plot, width = 8.27, height = 11.69) # A4 size


########################################################################
# Internal consistency
########################################################################
# Merge SPSQ item-wise scores across subjects 
library(dplyr)
library(readr)

working_dir <- 'analysis/data'
out_dir <- 'analysis'
subjects <- list.dirs(working_dir, recursive = FALSE, full.names = TRUE)

# Define the column names for merging
columns <- c(subject, names(qst_post_3_SPS))

spsq_all <- NULL

for (i in 1:length(subjects)) {
  subject <- subjects[i]
  subject_num <- sub("_.*", "", basename(subject))
  print(paste("Processing:", subject_num))
  data <- read_csv(file.path(subject, 'qst-post-3_SPS.csv'), show_col_types = FALSE)
  
  # Add subject column to data
  data$subject <- subject_num
  
  if (is.null(spsq_all)) {
    spsq_all <- data
  } else {
    # Identify common columns for merging
    common_cols <- intersect(names(spsq_all), names(data))
    
    # Merge data frames using common columns
    spsq_all <- merge(spsq_all, data, by = common_cols, all = TRUE)
  }
}

write_csv(spsq_all, file.path(out_dir, 'spsq_all.csv'))
names(spsq_all)

# Calculate Cronbach's alpha for the positive and negative dimensions of SPSQ-24
library(psych)

# Define function to calculate Cronbach's alpha for a given set of items
calculate_alpha <- function(items, dataframe) {
  # Subset dataframe with only the columns corresponding to the given items
  subset_data <- dataframe[, items]
  
  # Calculate Cronbach's alpha
  alpha_result <- psych::alpha(subset_data)
  
  return(alpha_result$total$raw_alpha)
}

# Define the lists of items
spsq_sum <- c('SPS13_emotionally_touched_music_art', 'SPS29_notice_subtle_touching_tones_music', 'SPS31_very_movedy_nice_work_of_art', 'SPS02_nervous_to_many_things_at_once', 'SPS24_rushed_too_much_little_time', 'SPS28_upset_when_people_ask_many_things_at_once', 'SPS03_see_sad_eyes_behind_smile', 'SPS04_strikes_tone_voice_not_matching_words', 'SPS09_looking_eyes_telling_truth', 'SPS12_strikes_when_acting_not_afraid', 'SPS22_tell_smile_masking_feelings',
              'SPS05_reversed', #'SPS05_hard_enjoy_little_things_reversed'
              'SPS10_feel_good_with_people_I_love', 'SPS17_enjoy_humour_situations', 'SPS19_enjoy_relaxing_activity', 'SPS30_watching_nice_movie_feels_good', 'SPS06_flashing_lights_bother', 'SPS16_easily_disturbed_light_odors', 'SPS20_loud_noises_irritating', 'SPS25_suffer_bright_light', 'SPS14_immediately_feel_mouth_throat_drier', 'SPS15_hardly_visible_details_attract_attention', 'SPS21_quickly_aware_changes_body', 'SPS23_notice_faints_smells')

positive_dimension <- c('SPS13_emotionally_touched_music_art', 'SPS29_notice_subtle_touching_tones_music', 'SPS31_very_movedy_nice_work_of_art', 'SPS03_see_sad_eyes_behind_smile', 'SPS04_strikes_tone_voice_not_matching_words', 'SPS09_looking_eyes_telling_truth', 'SPS12_strikes_when_acting_not_afraid', 'SPS22_tell_smile_masking_feelings',
                        'SPS05_reversed', #'SPS05_hard_enjoy_little_things_reversed'
                        'SPS10_feel_good_with_people_I_love', 'SPS17_enjoy_humour_situations', 'SPS19_enjoy_relaxing_activity', 'SPS30_watching_nice_movie_feels_good', 'SPS14_immediately_feel_mouth_throat_drier', 'SPS15_hardly_visible_details_attract_attention', 'SPS21_quickly_aware_changes_body', 'SPS23_notice_faints_smells')

negative_dimension <- c('SPS02_nervous_to_many_things_at_once', 'SPS24_rushed_too_much_little_time', 'SPS28_upset_when_people_ask_many_things_at_once', 'SPS06_flashing_lights_bother', 'SPS16_easily_disturbed_light_odors', 'SPS20_loud_noises_irritating', 'SPS25_suffer_bright_light')

# Calculate Cronbach's alpha for each list
alpha_spsq_sum <- calculate_alpha(spsq_sum, spsq_all)
alpha_positive_dimension <- calculate_alpha(positive_dimension, spsq_all)
alpha_negative_dimension <- calculate_alpha(negative_dimension, spsq_all)

# Output the results
cat("Cronbach's alpha for spsq_sum:", alpha_spsq_sum, "\n")
cat("Cronbach's alpha for positive_dimension:", alpha_positive_dimension, "\n")
cat("Cronbach's alpha for negative_dimension:", alpha_negative_dimension, "\n")

# Create a data frame for the results
alpha_results <- data.frame(
  items = c("spsq_sum", "positive_dimension", "negative_dimension"),
  alpha = c(alpha_spsq_sum, alpha_positive_dimension, alpha_negative_dimension)
)

# Write the results to a CSV file
write.csv(alpha_results, 'analysis/cronbachs_alpha.csv', row.names = FALSE)

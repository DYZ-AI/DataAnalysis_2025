#--------------------start-------------------------------
# Get current working directory
current_dir <- getwd()
# Define the path to dataset relative to project root
data_path <- "DataSet_No_Details.csv"

#----------------Read the dataset-------------------------
dataset <- read.csv(data_path)
# Display the structure and types of variables
str(dataset)

# Generate a beautiful summary with histograms for numeric variables
library(skimr)
skim(dataset)

#----------------Data Preparation-------------------------
library(dplyr)
# Removing specific columns
columns_to_remove <- c("h_index_34", "h_index_56", "hormone10_1", "hormone10_2", "an_index_23", "outcome", 
                       "factor_eth", "factor_h", "factor_pcos", "factor_prl")
cleaned_data <- dataset %>% select(-any_of(columns_to_remove))

# Create a new data frame for factors of interest
factor_data <- dataset %>% select(record_id, outcome, factor_eth, factor_h, factor_pcos, factor_prl)

# Display the structure of cleaned data
str(cleaned_data)

# Summarize factor data
summary(factor_data)

#--------------Identifying Missing Values-----------------
# Calculate missing values in the dataset
missing_values <- sum(is.na(cleaned_data))
column_missing_values <- colSums(is.na(cleaned_data))

# Display summary of missing values
skim(cleaned_data)

# Calculate percentage of missing values per column
missing_percentage <- colMeans(is.na(cleaned_data)) * 100
missing_percentage_filtered <- missing_percentage[missing_percentage <= 35]
missing_percentage_filtered

missing_percentage_high <- missing_percentage[missing_percentage > 35]
missing_percentage_high

#----------------Visualizing Missing Data Patterns---------
library(visdat)
vis_miss(cleaned_data)  # Visualizes NA patterns

library(naniar)
gg_miss_var(cleaned_data)  # Barplot of missingness per variable

#----------------Handle Missing Data---------------------
# Remove more columns with excess missing values
columns_to_remove_additional <- c("hormone9", "hormone11", "hormone12", "hormone13", "hormone14")
modified_data <- cleaned_data %>% select(-any_of(columns_to_remove_additional))

# View structure of the modified data
str(modified_data)

#----------------Perform Little's MCAR Test----------------
library(mice)
library(dplyr)

# Cleaning the data to avoid character errors
cleaned_for_mcar <- modified_data %>% 
  select(where(~!all(is.na(.)))) %>%
  mutate(across(where(is.character), as.factor))

# Perform the MCAR test
mcar_test_result <- mcar_test(cleaned_for_mcar)

# Output result and interpretation function
interpret_mcar_test <- function(result) {
  p_val <- result$p.value
  if (p_val > 0.05) {
    message("✅ p-value > 0.05 → Data is likely MCAR. Safe to delete or impute.")
  } else {
    message("❗ p-value <= 0.05 → Data is NOT MCAR. Consider using MAR or MNAR.")
    message("➡️ Suggest using multiple imputation (e.g., pmm / rf).")
  }
}

interpret_mcar_test(mcar_test_result)

#----------------Multiple Imputation using MICE---------------
# Perform multiple imputation using Random Forest method
imputed_data_rf <- mice(modified_data, m = 5, method = 'rf', print = FALSE)
imputed_data_rf_final <- complete(imputed_data_rf)

# Density plot comparison: Original vs Imputed (Random Forest)
library(ggplot2)
ggplot(modified_data, aes(x = hormone10_generated, fill = "Original")) +
  geom_density(alpha = 0.5) +
  geom_density(data = imputed_data_rf_final, aes(x = hormone10_generated, fill = "Imputed"), alpha = 0.5) +
  labs(title = "Density Plot of hormone10_generated: Original vs. Imputed") +
  scale_x_continuous(limits = c(0, 2))

# Perform imputation using Predictive Mean Matching (PMM)
imputed_data_pmm <- mice(modified_data[, !names(modified_data) %in% "New"], method = "pmm")
imputed_data_pmm_final <- complete(imputed_data_pmm)

# Density plot for PMM imputed data
ggplot(modified_data, aes(x = hormone10_generated, fill = "Original")) +
  geom_density(alpha = 0.5) +
  geom_density(data = imputed_data_pmm_final, aes(x = hormone10_generated, fill = "Imputed"), alpha = 0.5) +
  labs(title = "Density Plot of hormone10_generated: Original vs. Imputed") +
  scale_x_continuous(limits = c(0, 2))

#----------------Outlier Detection--------------------------
# Outlier detection using Boxplot
outlier_data <- imputed_data_rf_final %>% select(lipids1, lipids2, lipids3, lipids4, lipids5) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# Create boxplot for outliers
ggplot(outlier_data, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Outlier Detection", x = "Variables", y = "Value") +
  theme_minimal()

# Outlier detection across entire dataset
imputed_data_rf_final %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Boxplots for Outlier Detection")

# LOF Outlier Detection
library(dbscan)
lof_scores <- lof(imputed_data_rf_final, k = 20)
lof_df <- imputed_data_rf_final %>% mutate(lof_score = lof_scores)

# Visualizing LOF scores
ggplot(lof_df, aes(x = lof_score)) +
  geom_histogram(binwidth = 0.05, fill = "#FF7F00", color = "black", alpha = 0.7) +
  labs(title = "Histogram of LOF Scores", x = "LOF Score", y = "Frequency") +
  theme_minimal()

# LOF Scatter Plot based on two variables
ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = lof_score), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "LOF-based Outlier Detection", x = "lipids1", y = "lipids2", color = "LOF Score") +
  theme_minimal()

# Identify Top 5% Outliers based on LOF
threshold <- quantile(lof_scores, 0.95)
lof_df <- lof_df %>% mutate(is_outlier = lof_score > threshold)

# Highlight outliers in scatter plot
ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = is_outlier), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(title = "LOF Outliers (Top 5%)", color = "Outlier") +
  theme_minimal()

# PCA Analysis for Outlier Detection
library(FactoMineR)
library(factoextra)
pca_results <- PCA(imputed_data_rf_final, scale.unit = TRUE, graph = FALSE)
pca_data <- data.frame(pca_results$ind$coord) %>% mutate(lof_score = lof_scores)

# Plot PCA with LOF scores
ggplot(pca_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = lof_score), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "LOF-based Outlier Detection in PCA Space", x = "PC1", y = "PC2") +
  theme_minimal()

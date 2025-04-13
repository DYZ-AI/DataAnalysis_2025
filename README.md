# Missing Data Imputation and Outlier Detection

This repository presents a complete workflow for handling missing values and identifying potential error records in a clinical dataset using R. The techniques include MCAR testing, multiple imputation with MICE, and outlier detection via Local Outlier Factor (LOF) and z-score diagnostics.

The working language is R, based on the scripts provided in the course materials.

Dataset Source: Course materials: Module 1

Source Code: ./DataCleaning/

Output: ./Outputs/graphs/

# Working flows
1. **Data Import & Cleaning**
Imported the dataset using read.csv()
Removed irrelevant variables with more than 35% missing values
Filtered the dataset to retain only variables with <=35% missing data
2. **Missing Data Analysis**
Visualization:
Used the visdat and naniar packages to visualize missing data patterns across variables
Generated bar plots and heatmaps to inspect the distribution and patterns of missing values
MCAR Testing:
Applied mcar_test() from the mice package to evaluate if the missingness mechanism is MCAR
Concluded that the data is not MCAR, suggesting the missing data is likely missing at random (MAR) or missing not at random (MNAR)
3. **Multiple Imputation**
Performed Imputation with the mice package:
Random Forest (method = "rf") for complex, non-linear relationships
Predictive Mean Matching (method = "pmm") for preserving the original data distribution
Results:
Generated imputed datasets for both methods and compared the imputation results
Visualized the distributions of original vs. imputed data
Concluded that PMM is the preferred method for imputation due to its ability to preserve the original data distribution, with effects comparable to the Random Forest method
4.  **Visualization**
Plotted distributions comparing original vs. imputed data:
For both PMM and RF methods, the distributions of imputed data closely matched the original data, validating the imputation process.
Comparison:
Compared the distributions of each imputed feature against the original to ensure the imputation did not alter the data's inherent properties.
Outcome: The PMM method was chosen as the recommended interpolation method because it better preserves the original data distribution.
5.  **Outlier Detection**
Local Outlier Factor (LOF):
Used the dbscan::lof() function to detect multivariate anomalies in the dataset.
Visualized top LOF scores via histograms and scatter plots to identify potential outliers.
Z-Score Method:
Calculated z-scores for each numeric variable to flag extreme values where z-score > 3.
Flagged these extreme values as potential outliers or data entry errors.
6.  **Suspected Error Reporting**
Extracted High-LOF Records for further inspection.
Flagged Extreme Values where the z-score was greater than 3.
Exported Results for manual review:
Suspected outliers were exported to a CSV file for human review to determine if they are actual outliers or data entry mistakes.
Future work involves screening methods to identify suspicious records that require further investigation.


# Insights
PMM preserves the original variable distribution better than RF.
LOF combined with z-score provides a robust method for error detection.
Manual review is essential to differentiate between true outliers and data entry errors.





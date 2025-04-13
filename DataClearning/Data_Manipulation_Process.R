#--------------------start-------------------------------
# Get current working directory
getwd()
# create path relative to project root
#data_path <- "/cloud/project/data/DataSet_No_Details.csv"
data_path <- "DataSet_No_Details.csv"
#----------------read dataset--------------------------
df <- read.csv(data_path)
# Display structure with variable types
str(df)
# Beautiful summary with histograms for numeric variables
# install.packages("skimr")
library(skimr)
skim(df) 
#---------------data set preparation------------------
library(dplyr)
# Delete a few columns 
cols_to_remove <- c("h_index_34", "h_index_56", "hormone10_1", "hormone10_2","an_index_23","outcome","factor_eth","factor_h","factor_pcos","factor_prl")
MD_df <- df %>% select(-any_of(cols_to_remove))
factor_df <- df %>% select (record_id, outcome, factor_eth, factor_h,factor_pcos,factor_prl)
str(MD_df)
summary(factor_df)

#--------------Identify Missing Values-----------------
sum(is.na(MD_df))               # Total NAs in entire dataset
colSums(is.na(MD_df))           # NA counts per column
skim(MD_df)
na_stats <- colMeans(is.na(MD_df)) * 100 # % missing data
na_stats
na_stats_filtered <- na_stats[na_stats <= 35] #  missing data <=35 %
# result as a table
data.frame(
  Column = names(na_stats_filtered),
  NA_Percent = na_stats_filtered,
  row.names = NULL
)

na_stats_filtered_1 <- na_stats[na_stats > 35] # missing data >35 %
# result as a table
data.frame(
  Column = names(na_stats_filtered_1),
  NA_Percent = na_stats_filtered_1,
  row.names = NULL
)

#-------------------Visualizing Missing Data Patterns------------------
# install.packages(visdat)
library(visdat)
vis_miss(MD_df)  # Visualizes NA patterns

# install.packages(naniar)
library(naniar)
gg_miss_var(MD_df)  # Barplot of missingness per variable
#------------------ Analyzing the Impact of Missing Data--------------
# Delete a few columns 
library(dplyr)
cols_to_remove1 <- c("hormone9", "hormone11", "hormone12", "hormone13","hormone14")
handle_MD_df <- MD_df %>% select(-any_of(cols_to_remove1))
str(handle_MD_df)

#------------------Performing Little's MCAR Test----------------------
#Homework!!!
#Hypotheses:
#  H₀ (Null Hypothesis): Data is MCAR.

#H₁ (Alternative Hypothesis): Data is not MCAR (either MAR or MNAR).

#If p-value > 0.05, we fail to reject H₀ (data is likely MCAR).
#If p-value ≤ 0.05, we reject H₀ (data is likely not MCAR).

library(mice)
library(dplyr)

# 防止字符变量报错（mcar_test不支持 character 类型）
handle_MD_df_clean <- handle_MD_df %>%
  select(where(~!all(is.na(.)))) %>%
  mutate(across(where(is.character), as.factor))

# 执行 MCAR 检验
mcar_result <- mcar_test(handle_MD_df_clean)

# 输出结果
print(mcar_result)

interpret_mcar <- function(mcar_result) {
  p <- mcar_result$p.value
  if (p > 0.05) {
    message("✅ p-value > 0.05 → Data is likely MCAR. Safe to delete or impute.")
  } else {
    message("❗ p-value <= 0.05 → Data is NOT MCAR. Assume MAR or MNAR.")
    message("➡️ 建议使用多重插补（mice）方法，如 pmm / rf 等。")
  }
}

interpret_mcar(mcar_result)




#------------------Imputation with MICE-------------------------------
# Install packages if they are not already installed
# install.packages(c("mice", "ggplot2", "naniar"))
# Load the packages
library(mice)
library(ggplot2)
library(naniar)

# Perform Multiple Imputation 
#Random Forest method 
## Perform Multiple Imputation
imputed_handle_MD_df <- mice(handle_MD_df, m=5, method='rf', print=FALSE)
# ------For complex nonlinear relationships between variables------------

imputed_handle_MD_df <- mice(handle_MD_df[, !names(handle_MD_df) %in% "New"], method="rf")  
  imputed_handle_MD_df_final <- complete(imputed_handle_MD_df)  # generate full data
# Density plots 
ggplot(handle_MD_df, aes(x=hormone10_generated, fill="Original")) +
  geom_density(alpha=0.5) +
  geom_density(data=imputed_handle_MD_df_final, aes(x=hormone10_generated, fill="Imputed"), alpha=0.5) +
  labs(title="Density Plot of hormone10_generated: Original vs. Imputed")+
  scale_x_continuous(limits = c(0, 2))

#Predictive Mean Matching 
#------default for numeric data------------
#Homework!!!
imputed_handle_MD_df1 <- mice(handle_MD_df[, !names(handle_MD_df) %in% "New"], method="pmm")  
imputed_handle_MD_df_final1 <- complete(imputed_handle_MD_df)  # generate full data
# Density plots 
ggplot(handle_MD_df, aes(x=hormone10_generated, fill="Original")) +
  geom_density(alpha=0.5) +
  geom_density(data=imputed_handle_MD_df_final1, aes(x=hormone10_generated, fill="Imputed"), alpha=0.5) +
  labs(title="Density Plot of hormone10_generated: Original vs. Imputed")+
  scale_x_continuous(limits = c(0, 2))
                     
# Camparsion
library(ggplot2)
library(dplyr)
library(visdat)

# 自动提取需要插值的变量名（有缺失的数值变量）
vars_to_plot <- handle_MD_df %>%
  select(where(is.numeric)) %>%
  select(where(~ any(is.na(.)))) %>%
  names()

# 定义统一的绘图函数
plot_imputation_comparison <- function(var_name) {
  # 提取三组数据
  original <- handle_MD_df %>%
    select(all_of(var_name)) %>%
    mutate(source = "Original")
  
  rf_imputed <- imputed_handle_MD_df_final %>%
    select(all_of(var_name)) %>%
    mutate(source = "Random Forest")
  
  pmm_imputed <- imputed_handle_MD_df_final1 %>%
    select(all_of(var_name)) %>%
    mutate(source = "PMM")
  
  # 合并数据
  combined <- bind_rows(original, rf_imputed, pmm_imputed)
  
  if (var_name == "hormone10_generated") {
    x_lim <- c(0, 5)
  } else {
    # 自动获取范围
    x_min <- min(combined[[var_name]], na.rm = TRUE)
    x_max <- max(combined[[var_name]], na.rm = TRUE)
    
    # 加入buffer，防止贴边
    x_range <- x_max - x_min
    buffer <- x_range * 0.05
    x_lim <- c(x_min - buffer, x_max - buffer)
  }
  
  # 绘制密度图
  ggplot(combined, aes_string(x = var_name, fill = "source", color = "source")) +
    geom_density(alpha = 0.4, size = 1, na.rm = TRUE) +
    labs(title = paste("Density Comparison of:", var_name),
         x = var_name,
         y = "Density") +
    scale_x_continuous(limits = x_lim) +
    scale_fill_manual(values = c(
      "Original" = "black",
      "PMM" = "#E69F00",
      "Random Forest" = "#56B4E9"
    )) +
    scale_color_manual(values = c(
      "Original" = "black",
      "PMM" = "#E69F00",
      "Random Forest" = "#56B4E9"
    )) +
    theme_minimal() +
    theme(legend.position = "top")
}

for (v in vars_to_plot) {
  print(plot_imputation_comparison(v))     # 显示图形
  # 可选：保存图像（取消注释保存）
  ggsave(paste0("imputation_density_", v, ".png"), 
         plot = plot_imputation_comparison(v),
         width = 6, height = 4)
}



#----------------Outlier Detection Methods------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
outliers_data <- imputed_handle_MD_df_final %>%
  select(lipids1, lipids2, lipids3, lipids4, lipids5) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# build a graph 
ggplot(outliers_data, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Outlier Detection",
       x = "variables",
       y = "value") +
  theme_minimal()
#build a graph for all dataset
imputed_handle_MD_df_final %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Boxplots for Outlier Detection")

library(ggplot2)
library(dbscan)
lof_data <- imputed_handle_MD_df_final %>%
  select(where(is.numeric)) %>%
  na.omit()

lof_score <- lof(lof_data, k = 20)
lof_df <- lof_data %>%
  mutate(lof_score = lof_score)

ggplot(lof_df, aes(x = lof_score)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "#FF7F00",
    color = "black",
    alpha = 0.7
  ) +
  labs(
    title = "Histogram of  LOF Scores",
    x = "LOF Score",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = lof_score), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "LOF-based Outlier Detection",
       x = "lipids1",
       y = "lipids2",
       color = "LOF Score") +
  theme_minimal()
  
  
# 加入是否为异常标签（Top 5%）
threshold <- quantile(lof_score, 0.95)
lof_df <- lof_df %>%
  mutate(is_outlier = lof_score > threshold)

# 高亮异常点的散点图
ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(aes(color = is_outlier), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(title = "LOF Outliers (Top 5%)",
       color = "Outlier") +
  theme_minimal()
  

lof_suspects <- lof_df %>%
  filter(lof_score > 1.5) %>%
  arrange(desc(lof_score))

top_lof_errors <- lof_suspects %>%
  slice(1:10)

print(top_lof_errors)

# 用 z-score 方式标注每列是否极端
zscore_outlier_check <- function(df_row, reference_df, threshold = 3) {
  z_check <- sapply(names(df_row), function(col) {
    mu <- mean(reference_df[[col]], na.rm = TRUE)
    sigma <- sd(reference_df[[col]], na.rm = TRUE)
    abs((df_row[[col]] - mu) / sigma) > threshold
  })
  return(z_check)
}

# 应用于 top_lof_errors 的每一行
error_flags <- apply(top_lof_errors[, -ncol(top_lof_errors)], 1, zscore_outlier_check, reference_df = lof_data)

# 转置 & 转成可视化表格
error_matrix <- as.data.frame(t(error_flags))
rownames(error_matrix) <- paste0("row_", 1:nrow(error_matrix))
error_matrix

ggplot(lof_df, aes(x = lipids1, y = lipids2)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = top_lof_errors, aes(x = lipids1, y = lipids2), color = "red", size = 3) +
  labs(title = "LOF Top 10 Suspected Error Points")

library(FactoMineR)
library(factoextra)

pca_res <- PCA(lof_data, scale.unit = TRUE, graph = FALSE)
pca_df <- data.frame(pca_res$ind$coord) %>%
  mutate(lof_score = lof_scores)

ggplot(pca_df, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = lof_score), size = 2, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "LOF-based Outlier Detection in PCA Space",
       x = "PC1", y = "PC2") +
  theme_minimal()

  

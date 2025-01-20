library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(MASS)
library(caret)
library(pROC)
library(ggpubr)
library(glmnet)
library(missForestPredict)


data_path = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/"
ukbb <- read.csv(paste0(data_path,'UKBB_prostate_Prevalent_HC_1.csv'), row.names = 1)

ukbb$response = as.factor(ukbb$response)


#---------------divide and impute the data ------------------

# Set a seed for reproducibility
set.seed(42)

# Calculate the threshold for non-NA values
threshold <- 0.8 * nrow(ukbb)

# Filter columns based on the threshold
filtered_data <- ukbb[, colSums(!is.na(ukbb)) >= threshold]

n_proteins = dim(filtered_data)[2]-2

# Create a partition index
trainIndex <- createDataPartition(filtered_data$response, p = 0.7, list = FALSE)

# Split the data into training and testing sets
train_data <- filtered_data[trainIndex, ]
test_data <- filtered_data[-trainIndex, ]

# Imputation

imputation_object <- missForestPredict::missForest(train_data[,c(1:n_proteins)], 
                                                   save_models = TRUE, 
                                                   num.threads = 2)

train_data_imputed <- imputation_object$ximp
train_data_imputed$status = train_data$response


test_data_imputed <- missForestPredict::missForestPredict(imputation_object, 
                                                          newdata = test_data[,c(1:n_proteins)])
test_data_imputed$status = test_data$response

#-------------------- GLM PSA ---------------

# Initialize a data frame to store the results
results <- data.frame(Gene=character(), AUC=numeric(), stringsAsFactors=FALSE)

# Fit the GLM model for PSA
glm_model <- glm(status ~ KLK3, data=train_data_imputed, family=binomial)

# Make predictions on the test data
predicted_prob <- predict(glm_model, newdata=test_data_imputed, type="response")

#--------------- save current data ----------------
write.csv(predicted_prob, file = paste0(data_path,"UKBB_KLK3_predicted_prob.csv"))

# Calculate the AUC
roc_obj <- roc(test_data_imputed$status, predicted_prob)
auc_value_PSA <- roc_obj$auc
roc_plot_PSA <- ggroc(roc_obj)

################# glmnet status #############################3


x_train <- as.matrix(train_data_imputed[, !names(train_data_imputed) %in% c("status")])
y_train <- ifelse(train_data_imputed$status == 'Prevalent', 1, 0)

x_test <- as.matrix(test_data_imputed[, !names(test_data_imputed) %in% c("status")])
y_test  <- ifelse(test_data_imputed$status == 'Prevalent', 1, 0)


cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)

# Get the best lambda value
best_lambda <- cv_fit$lambda.min

# Fit the final model with the best lambda
final_model <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = best_lambda)

# Display the coefficients of the final model
print(coef(final_model))

# Make predictions on the test data
predicted_prob <- predict(final_model, x_test, s="lambda.min")

#--------------- save current data ----------------
data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'
write.csv(y_test, file = paste0(data_path,"UKBB_y_test.csv"))
write.csv(predicted_prob, file = paste0(data_path,"UKBB_predicted_prob.csv"))

# Calculate the AUC
roc_curve <- roc(y_test, predicted_prob)
roc_plot <- ggroc(roc_curve)

auc_value <- roc_curve$auc

roc_curve_combined <- ggplot() +
  geom_line(data = roc_plot_PSA$data, aes(x = 1 - specificity, y = sensitivity, color = "PSA blood")) +
  geom_line(data = roc_plot$data, aes(x = 1 - specificity, y = sensitivity, color = "Our biomarkers")) +
  labs(title = "UKBB proteomics", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("PSA blood" = "#D81B60", "Our biomarkers" = "#FFC107")) +
  theme_minimal() +
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC PSA blood:", round(auc_value_PSA, 2)), color = "#D81B60") +
  annotate("text", x = 0.6, y = 0.1, label = paste("AUC Our biomarkers:", round(auc_value, 2)), color = "#FFC107")

roc_curve_combined

# 
# # Convert the sparse matrix to a dense matrix
# dense_matrix <- as.matrix(coef(final_model))
# 
# # Convert the dense matrix to a data frame
# df <- as.data.frame(dense_matrix)
# 
# df_non_zero  <- df[rowSums(df != 0) > 0, , drop = FALSE]
# 
# genes_data <- data.frame(
#   Gene = rownames(df_non_zero)[2:nrow(df_non_zero)],
#   Coefficient = df_non_zero[2:nrow(df_non_zero),1]
# )
# 
# 
# # Plot the coefficients of the selected genes
# genes_plot <- ggplot(genes_data, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
#   geom_bar(stat = "identity", fill = "skyblue") +
#   coord_flip() +
#   labs(title = "Selected Genes and Their Coefficients", x = "Gene", y = "Coefficient") +
#   theme_minimal()
# 
# # Combine the plots
# combined_plot <- ggarrange(roc_curve_combined, genes_plot, ncol = 2, nrow = 1)
# combined_plot <- annotate_figure(combined_plot, top = text_grob("UKBB proteomics", size = 14))
# 
# # Print the combined plot
# print(combined_plot)
# print(paste("AUC:", auc_value))
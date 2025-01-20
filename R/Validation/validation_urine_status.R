library("readxl")
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


# Data from this paper:
#https://www.nature.com/articles/s41467-024-49424-5
# Can be downloaded as supplement 1 and supplement 2


#------------------ Read data -----------------------------
data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'

patient_info = read_excel(paste0(data_path,'Table1bwithpatientcodes.xlsx'), sheet = '1b', skip = 1)
urine_train = read_excel(paste0(data_path,'41467_2024_49424_MOESM5_ESM.xlsx'), sheet = '2a', skip = 1)
urine_validation = read_excel(paste0(data_path,'41467_2024_49424_MOESM5_ESM.xlsx'), sheet = '2b', skip = 2)

features <- c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
              'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
              'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA-E', 'HLA-DRA',
              'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
              'PTGDS', 'HLA-A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
              'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB')
Sample_type = 'uEV-P20'


#------------------ preprocess data -----------------------------

rnames = urine_train$`Gene name`
cnames = colnames(urine_train)

# transpose
urine_train <- transpose(urine_train)

# get row and colnames in order
colnames(urine_train) <- make.names(rnames)
rownames(urine_train) <- make.names(cnames)

rnames = urine_validation$`gene`
cnames = colnames(urine_validation)

# transpose
urine_validation <- transpose(urine_validation)

# get row and colnames in order
colnames(urine_validation) <- make.names(rnames)
rownames(urine_validation) <- make.names(cnames)

features = features[features %in% colnames(urine_train) & features %in% colnames(urine_validation)]

patient_info_Discovery = patient_info[patient_info$`Sample type`==Sample_type & patient_info$DRE == 'Post' & patient_info$Cohort == 'Discovery', ] 

urine_train = urine_train[patient_info_Discovery$`Sample ID`, features]

urine_train$PSA = patient_info_Discovery$`serum PSA (ng/mL)`
urine_train[] <- lapply(urine_train, as.numeric)
urine_train$grade = patient_info_Discovery$`cISUP Grade Group`
urine_train$status = patient_info_Discovery$`Prostate cancer status`

patient_info_Validation = patient_info[patient_info$`Sample type`==Sample_type & patient_info$DRE == 'Post' & patient_info$Cohort == 'Validation 1', ] 

urine_validation = urine_validation[patient_info_Validation$`Sample ID`, features]

urine_validation$PSA = patient_info_Validation$`serum PSA (ng/mL)`
urine_validation[] <- lapply(urine_validation, as.numeric)



urine_validation$grade = patient_info_Validation$`cISUP Grade Group`
urine_validation$status = patient_info_Validation$`Prostate cancer status`
urine_train$grade = as.numeric(urine_train$grade)
urine_train$grade[is.na(urine_train$grade)] <- 0
urine_train$status = as.factor(urine_train$status)



urine_validation$grade = as.numeric(urine_validation$grade)
urine_validation$grade[is.na(urine_validation$grade)] <- 0
urine_validation$status = as.factor(urine_validation$status)


urine_all = rbind(urine_train, urine_validation)



#---------------divide and impute the data ------------------

# Set a seed for reproducibility
set.seed(42)

# Calculate the threshold for non-NA values
threshold <- 0.8 * nrow(urine_all)

# Filter columns based on the threshold
filtered_data <- urine_all[, colSums(!is.na(urine_all)) >= threshold]

n_proteins = dim(filtered_data)[2]-2

# Create a partition index
trainIndex <- createDataPartition(filtered_data$grade, p = 0.7, list = FALSE)

# Split the data into training and testing sets
train_data <- filtered_data[trainIndex, ]
test_data <- filtered_data[-trainIndex, ]

# Imputation

imputation_object <- missForestPredict::missForest(train_data[,c(1:n_proteins)], 
                                                       save_models = TRUE, 
                                                       num.threads = 2)

train_data_imputed <- imputation_object$ximp
train_data_imputed$grade = train_data$grade
train_data_imputed$status = train_data$status


test_data_imputed <- missForestPredict::missForestPredict(imputation_object, 
                                                      newdata = test_data[,c(1:n_proteins)])
test_data_imputed$grade = test_data$grade
test_data_imputed$status = test_data$status



#-------------------- GLM PSA ---------------

# Initialize a data frame to store the results
results <- data.frame(Gene=character(), AUC=numeric(), stringsAsFactors=FALSE)

# Fit the GLM model for PSA
glm_model <- glm(status ~ PSA, data=train_data_imputed, family=binomial)

# Make predictions on the test data
predicted_prob <- predict(glm_model, newdata=test_data_imputed, type="response")

#--------------- save current data ----------------
data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'
write.csv(predicted_prob, file = paste0(data_path,"PSA_predicted_prob.csv"))
write.csv(test_data_imputed$status, file = paste0(data_path,"y_status.csv"))

# Calculate the AUC
roc_obj <- roc(test_data_imputed$status, predicted_prob)
auc_value_PSA <- roc_obj$auc
roc_plot_PSA <- ggroc(roc_obj)

#-------------------- GLM KLK3 ---------------

# Initialize a data frame to store the results
results <- data.frame(Gene=character(), AUC=numeric(), stringsAsFactors=FALSE)

# Fit the GLM model for KLK3
glm_model <- glm(status ~ KLK3, data=train_data_imputed, family=binomial)

# Make predictions on the test data
predicted_prob <- predict(glm_model, newdata=test_data_imputed, type="response")

#--------------- save current data ----------------
write.csv(predicted_prob, file = paste0(data_path,"KLK3_predicted_prob.csv"))

# Calculate the AUC
roc_obj <- roc(test_data_imputed$status, predicted_prob)
auc_value_KLK3 <- roc_obj$auc
roc_plot_KLK3 <- ggroc(roc_obj)

################# glmnet status #############################3


x_train <- as.matrix(train_data_imputed[, !names(train_data_imputed) %in% c("PSA", "grade", "status")])
y_train <- ifelse(train_data_imputed$status == 'Non-Prostate cancer (BPH)', 0, 1)
x_test <- as.matrix(test_data_imputed[, !names(test_data_imputed) %in% c("PSA", "grade", "status")])
y_test  <- ifelse(test_data_imputed$status == 'Non-Prostate cancer (BPH)', 0, 1)


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
write.csv(y_test, file = paste0(data_path,"y_test.csv"))
write.csv(predicted_prob, file = paste0(data_path,"predicted_prob.csv"))
write.csv(test_data_imputed$grade, file = paste0(data_path,"y_grade.csv"))

# Calculate the AUC
roc_curve <- roc(y_test, predicted_prob)
roc_plot <- ggroc(roc_curve)

auc_value <- roc_curve$auc

roc_curve_combined <- ggplot() +
  geom_line(data = roc_plot_PSA$data, aes(x = 1 - specificity, y = sensitivity, color = "PSA blood")) +
  geom_line(data = roc_plot$data, aes(x = 1 - specificity, y = sensitivity, color = "Urine biomarkers")) +
  labs(title = "ROC Curves", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("PSA blood" = "#D81B60", "Urine biomarkers" = "#FFC107")) +
  theme_minimal() +
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC PSA blood:", round(auc_value_PSA, 2)), color = "#D81B60") +
  annotate("text", x = 0.6, y = 0.1, label = paste("AUC Urine biomarkers:", round(auc_value, 2)), color = "#FFC107")




# Convert the sparse matrix to a dense matrix
dense_matrix <- as.matrix(coef(final_model))

# Convert the dense matrix to a data frame
df <- as.data.frame(dense_matrix)

df_non_zero  <- df[rowSums(df != 0) > 0, , drop = FALSE]

genes_data <- data.frame(
  Gene = rownames(df_non_zero)[2:nrow(df_non_zero)],
  Coefficient = df_non_zero[2:nrow(df_non_zero),1]
)


# Plot the coefficients of the selected genes
genes_plot <- ggplot(genes_data, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Selected Genes and Their Coefficients", x = "Gene", y = "Coefficient") +
  theme_minimal()

# Combine the plots
combined_plot <- ggarrange(roc_curve_combined, genes_plot, ncol = 2, nrow = 1)
combined_plot <- annotate_figure(combined_plot, top = text_grob(paste("Urine proteomics", Sample_type), size = 14))

# Print the combined plot
print(combined_plot)
print(paste("AUC:", auc_value))




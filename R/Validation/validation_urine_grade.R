library("readxl")
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(MASS)
library(caret)
library(pROC)
library(ggpubr)
library(Metrics)
library(glmnet)

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

urine_all = urine_all[urine_all$status != 0,]

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




################# glmnet grade #############################

x_train <- as.matrix(train_data_imputed[, !names(train_data_imputed) %in% c("grade", "status", 'PSA')])
y_train <- train_data_imputed$grade
x_test <- as.matrix(test_data_imputed[, !names(test_data_imputed) %in% c("grade", "status", 'PSA')])
y_test  <- test_data_imputed$grade


cv_fit <- cv.glmnet(x_train, y_train, family = "gaussian", alpha = 1)

# Get the best lambda value
best_lambda <- cv_fit$lambda.min

# Fit the final model with the best lambda
final_model <- glmnet(x_train, y_train, family = "gaussian", alpha = 1, lambda = best_lambda)


# Make predictions on the test data
predicted_prob <- predict(final_model, x_test, type = 'response', s="lambda.min")


# Metrics
sse <- sum((predicted_prob - y_test)^2)
sst <- sum((y_test - mean(y_test))^2)
r_squared <- 1 - (sse / sst) 
print(paste("R-squared:", r_squared))

pred_corr <- cor(predicted_prob, y_test)
print(paste("correlation:", pred_corr))


model_rmse <- rmse(y_test, predicted_prob)


# Boxplot 

data_ggplot <- data.frame(TrueGrade = factor(y_test, levels = 0:5), PredictedProb = array(predicted_prob))

data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'
write.csv(data_ggplot, file = paste0(data_path,"Our_model_grade.csv"))

df_extra_data = data.frame("r_squared"=r_squared, "model_rmse" = model_rmse, "pred_corr" = pred_corr)
model_boxplot <- ggplot(data_ggplot, aes(x = TrueGrade, y = PredictedProb)) +
  geom_boxplot(trim = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_smooth(method="lm", col="black") + 
  labs(title = "Candidate biomarkers",
       x = "True Cancer Grade",
       y = "Predicted Severity") +
  ylim(0.5,7.5) + 
  annotate("text", x=1, y=6, label=paste0("R-squared = ", round(r_squared, 3)), hjust=0) +
  annotate("text", x=1, y=5, label=paste0("RMSE = ", round(model_rmse, 3)), hjust=0) +  
  annotate("text", x=1, y=5.5, label=paste0("Correlation = ", round(pred_corr, 3)), hjust=0) +  
  theme_minimal()

model_boxplot

# Convert the sparse matrix to a dense matrix
dense_matrix <- as.matrix(coef(final_model))

# Convert the dense matrix to a data frame
df <- as.data.frame(dense_matrix)

df_non_zero  <- df[rowSums(df != 0) > 0, , drop = FALSE]

genes_data <- data.frame(
  Gene = rownames(df_non_zero)[2:nrow(df_non_zero)],
  Coefficient = df_non_zero[2:nrow(df_non_zero),1]
)

genes_data
############# PSA grade prediction ########################


train_data <- train_data_imputed
test_data  <- test_data_imputed


# Initialize a data frame to store the results
results <- data.frame(Gene=character(), AUC=numeric(), stringsAsFactors=FALSE)

# Fit the GLM model

glm_model <- glm(grade ~ PSA, data=train_data, family='gaussian')


# Make predictions on the test data
predicted_prob <- predict(glm_model, test_data, type = 'response')


# Metrics
sse <- sum((predicted_prob - y_test)^2)
sst <- sum((y_test - mean(y_test))^2)
r_squared <- 1 - (sse / sst) 
print(paste("R-squared:", r_squared))

model_rmse <- rmse(y_test, predicted_prob)


pred_corr <- cor(predicted_prob, y_test)
print(paste("correlation:", pred_corr))


# Boxplot
data_ggplot <- data.frame(TrueGrade = factor(y_test, levels = 0:5), PredictedProb = predicted_prob)
write.csv(data_ggplot, file = paste0(data_path,"PSA_grade.csv"))

df_extra_data = rbind(df_extra_data,list(r_squared,model_rmse,pred_corr))
rownames(df_extra_data) =c("OurModel","PSAModel")
write.csv(df_extra_data, file = paste0(data_path,"EvaluationMeasures_grade.csv"))
write.csv(data_ggplot, file = paste0(data_path,"PSA_grade.csv"))

PSA_boxplot <- ggplot(data_ggplot, aes(x = TrueGrade, y = PredictedProb)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_smooth(method="lm", col="black") + 
  labs(title = "PSA",
       x = "True Cancer Grade",
       y = "Predicted Severity") +
  ylim(0.5,7.5) + 
  annotate("text", x=1, y=6, label=paste0("R-squared = ", round(r_squared, 3)), hjust=0) +
  annotate("text", x=1, y=5, label=paste0("RMSE = ", round(model_rmse, 3)), hjust=0) +  
  annotate("text", x=1, y=5.5, label=paste0("Correlation = ", round(pred_corr, 3)), hjust=0) +  
  
  theme_minimal()

model_boxplot
PSA_boxplot

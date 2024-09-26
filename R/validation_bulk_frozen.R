library(AnnotationDbi)
library(org.Hs.eg.db)
library("readxl")
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(MASS)
library(caret)
library(pROC)
library(ggpubr)
library(DHARMa)


pcatlas = read.table('data/pcatlas_dataset_vst_normalized.txt')
pcatlas_annotation = read_excel('data/pcatlas_dataset_annotations.xlsx')

pcatlas_annotation = pcatlas_annotation[pcatlas_annotation$`Tumor Tissue site` == 'PROSTATE',]




# List of ENSG gene annotations
ensg_ids <- rownames(pcatlas)


# Get gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensg_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")


features <- c('TMEFF2', 'SPON2', 'CRISP3', 'TRGC1', 'MYO6', 'STEAP2', 'RPS8', 'RPS4X',
              'EEF1G', 'KLK3', 'NPY', 'KLK2', 'RPL7A', 'IGKC', 'MT2A', 'NEFH',
              'LGALS1', 'HLA-E', 'C1R', 'HLA-A', 'HLA-DRA', 'IGFBP7', 'DES', 'COL1A2',
              'LGALS3BP', 'HSPB1', 'GSTP1', 'AZGP1', 'IFITM3', 'TAGLN', 'ACTG2',
              'B2M', 'FLNA', 'MYL9', 'CD74', 'TPM2', 'ACTA2', 'S100A6', 'MSMB',
              'TIMP1')

gene_symbols = gene_symbols[gene_symbols %in% features]

pcatlas = pcatlas[names(gene_symbols),]
pcatlas['genes'] = gene_symbols

# Group by 'genes' and calculate the mean for each group
pcatlas <- pcatlas %>%
  group_by(genes) %>%
  summarise(across(starts_with("entry"), mean))

pcatlas = column_to_rownames(pcatlas, var = "genes")
pcatlas_annotation = column_to_rownames(pcatlas_annotation, var = "Entry")

pcatlas = pcatlas[,rownames(pcatlas_annotation)]
pcatlas




#------------------- GLM gene specific ---------------

# transpose
t_pcatlas <- transpose(pcatlas)

# get row and colnames in order
colnames(t_pcatlas) <- make.names(rownames(pcatlas))
rownames(t_pcatlas) <- colnames(pcatlas)

t_pcatlas$SampleType = pcatlas_annotation$`Sample Type`
t_pcatlas$BiopsyAge = as.numeric(pcatlas_annotation$`Biopsy Age`)

t_pcatlas = t_pcatlas[t_pcatlas$SampleType %in% c('NORMAL', 'PRIMARY') ,]
t_pcatlas = t_pcatlas[!is.na(t_pcatlas$BiopsyAge),]

t_pcatlas$SampleType = as.factor(t_pcatlas$SampleType)

p_value_list <- list()
estimate_list <- list()

# Run the logistic regression
for(gene in colnames(t_pcatlas)[1:40]){
  print(gene)
  model <- glm(as.formula(paste("SampleType ~", gene, "+ BiopsyAge")),
               data = t_pcatlas, family = binomial)
  print(summary(model))
  # Get the summary of the model
  summary_model <- summary(model)
  
  # Check if the p-value for the gene is less than 0.05
  p_value <- summary_model$coefficients[gene, "Pr(>|z|)"]
  estimate <- summary_model$coefficients[gene, "Estimate"]
  
  p_value_list[[gene]] <- p_value
  estimate_list[[gene]] <- estimate
  
  
}

glm_GeneSpecific = data.frame('estimate' = unlist(estimate_list),
                            'p_val' = unlist(p_value_list),
                            row.names = colnames(t_pcatlas)[1:40])

glm_GeneSpecific['p_val_adj'] = p.adjust(glm_GeneSpecific$p_val, method = 'hochberg', n = dim(glm_GeneSpecific)[1])
glm_GeneSpecific = glm_GeneSpecific[order(glm_GeneSpecific$estimate),]

glm_GeneSpecific = glm_GeneSpecific[glm_GeneSpecific$p_val_adj < 0.05,]
glm_GeneSpecific

glm_GeneSpecific$gene = rownames(glm_GeneSpecific)

# Melt the data frame for ggplot2
data_melted <- melt(glm_GeneSpecific, id.vars = 'gene')

# Filter to keep only the Estimate values
data_estimate <- subset(data_melted, variable == 'estimate')

data_estimate$gene <- factor(data_estimate$gene, levels = features <- c('TMEFF2', 'SPON2', 'CRISP3', 'TRGC1', 'MYO6', 'STEAP2', 'RPS8', 'RPS4X',
                                                                        'EEF1G', 'KLK3', 'NPY', 'KLK2', 'RPL7A', 'IGKC', 'MT2A', 'NEFH',
                                                                        'LGALS1', 'HLA-E', 'C1R', 'HLA-A', 'HLA-DRA', 'IGFBP7', 'DES', 'COL1A2',
                                                                        'LGALS3BP', 'HSPB1', 'GSTP1', 'AZGP1', 'IFITM3', 'TAGLN', 'ACTG2',
                                                                        'B2M', 'FLNA', 'MYL9', 'CD74', 'TPM2', 'ACTA2', 'S100A6', 'MSMB',
                                                                        'TIMP1'))



# Create the heatmap
ggplot(data_estimate, aes(x = variable, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, 
                       limit = c(min(data_estimate$value), max(data_estimate$value)), 
                       name = "coef") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())




#------------------- GLM all genes ---------------

# transpose
t_pcatlas <- transpose(pcatlas)

# get row and colnames in order
colnames(t_pcatlas) <- make.names(rownames(pcatlas))
rownames(t_pcatlas) <- colnames(pcatlas)

t_pcatlas$SampleType = pcatlas_annotation$`Sample Type`
t_pcatlas$BiopsyAge = as.numeric(pcatlas_annotation$`Biopsy Age`)

t_pcatlas = t_pcatlas[t_pcatlas$SampleType %in% c('NORMAL', 'PRIMARY') ,]
t_pcatlas = t_pcatlas[!is.na(t_pcatlas$BiopsyAge),]

t_pcatlas$SampleType = as.factor(t_pcatlas$SampleType)

p_value_list <- list()
estimate_list <- list()

set.seed(42)  # For reproducibility
trainIndex <- createDataPartition(t_pcatlas$SampleType, p = .7, 
                                  list = FALSE, 
                                  times = 1)
trainData <- t_pcatlas[ trainIndex,]
testData  <- t_pcatlas[-trainIndex,]

initial_model <- glm(SampleType ~ 1, data = trainData, family = binomial(link = 'logit'))

# Perform forward stepwise selection
model.AIC <- stepAIC(initial_model, direction = 'forward', 
                     scope = list(lower = initial_model, upper = glm(SampleType ~ ., data = trainData, family = binomial(link = 'logit'))), 
                     k = 10, steps = 1, trace = 0)

# Extract p-values and continue forward selection until all p-values are below the threshold
pvals <- summary(model.AIC)$coefficients[,4]
aic = model.AIC$aic
aic_old = 0
while(aic_old != aic){
  aic_old = model.AIC$aic
  model.AIC <- stepAIC(model.AIC, direction = 'forward', 
                       scope = list(lower = initial_model, upper = glm(SampleType ~ ., data = trainData, family = binomial(link = 'logit'))), 
                       k = 10, steps = 1, trace = 0)
  pvals <- summary(model.AIC)$coefficients[2:dim(summary(model.AIC)$coefficients)[1],4]
  print(summary(model.AIC))
  aic = model.AIC$aic
}




# Extract selected genes and their coefficients
selected_genes <- summary(model.AIC)$coefficients[2:dim(summary(model.AIC)$coefficient)[1],]
selected_genes <- selected_genes[order(selected_genes[,4]),]  # Order by p-value
selected_genes <- selected_genes[selected_genes[,4] <= 0.05,]  # Keep significant genes


#create ROC plot
predictions <- predict(model.AIC, newdata = testData, type = "response")
roc_curve <- roc(testData$SampleType, predictions)
roc_plot <- ggroc(roc_curve)



# Create a data frame for the selected genes
genes_data <- data.frame(
  Gene = rownames(selected_genes),
  Coefficient = selected_genes[,1]
)

# Plot the coefficients of the selected genes
genes_plot <- ggplot(genes_data, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Selected Genes and Their Coefficients", x = "Gene", y = "Coefficient") +
  theme_minimal()

# Combine the plots
combined_plot <- ggarrange(roc_plot, genes_plot, ncol = 2, nrow = 1)

# Print the combined plot
print(combined_plot)
summary(model.AIC)



# Evaluation of the model
#The residuals follow the glm assumtions
# The plot on the left should be diagonal, the one on the right should be uniform
simulationOutput <- simulateResiduals(fittedModel = model.AIC, plot = T) #The residuals follow the glm assumtions



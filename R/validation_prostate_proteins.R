library("readxl")
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)


# Data from this paper:
#https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(24)00400-2
#For each sample, “_1” represents the lower Gleason score tumor sample;
#“_2” represents the higher Gleason score tumor samples; and “
#_3 represents the Benign part from tumor sample

#------------------ preprocessing -----------------------------
data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'

proteomics = read_excel(paste0(data_path,'mmc3.xlsx'))
proteomics_annotation = read_excel(paste0(data_path,'mmc2.xlsx'))

features <- c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
              'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
              'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA-E', 'HLA-DRA',
              'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
              'PTGDS', 'HLA-A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
              'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB')

proteomics = proteomics[proteomics$GeneName %in% features,]
rownames(proteomics) <- proteomics$GeneName

# transpose
t_proteomics <- transpose(proteomics)

# get row and colnames in order
colnames(t_proteomics) <- make.names(rownames(proteomics))
rownames(t_proteomics) <- colnames(proteomics)


t_proteomics = t_proteomics[3:dim(t_proteomics)[1],]

sample_characteristics = transpose(data.frame(str_split(rownames(t_proteomics), '_')))

t_proteomics['Patient_ID'] = sample_characteristics[2]
t_proteomics['sample'] = sample_characteristics[3]

proteomics_annotation = proteomics_annotation[c('Patient_ID', 'Primary_GleasonScore(G1)', 'Second_GleasonScore(G2)')]
t_proteomics = merge(x = t_proteomics, y = proteomics_annotation, by = "Patient_ID")
t_proteomics['grade'] = rep(0, dim(t_proteomics)[1])

subset_sample = t_proteomics$sample == 1
t_proteomics[subset_sample, 'grade'] = pmax(t_proteomics$`Primary_GleasonScore(G1)`[subset_sample],
                                            t_proteomics$`Second_GleasonScore(G2)`[subset_sample])

subset_sample = t_proteomics$sample == 2
t_proteomics[subset_sample, 'grade'] = pmin(t_proteomics$`Primary_GleasonScore(G1)`[subset_sample],
                                            t_proteomics$`Second_GleasonScore(G2)`[subset_sample])

features <- c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
                               'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
                               'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA-E', 'HLA-DRA',
                               'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
                               'PTGDS', 'HLA-A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
                               'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB', 'grade')

t_proteomics_PID = t_proteomics$Patient_ID
t_proteomics = t_proteomics[,colnames(t_proteomics) %in% features]

t_proteomics <- mutate_all(t_proteomics, function(x) as.numeric(as.character(x)))
t_proteomics2 <- t_proteomics
t_proteomics2$Patient_ID <- t_proteomics_PID
############################
gene = "AMACR"
features_gene <- c("Patient_ID",'grade',gene)
t_proteomics3 <- t_proteomics2[,colnames(t_proteomics2) %in% features_gene]
filtered_data <- t_proteomics3 %>%
  pivot_wider(
    names_from = grade,
    values_from = AMACR,
    names_prefix = "Grade_"
  )
t_proteomics2 <- t_proteomics2[order(t_proteomics2$Patient_ID, t_proteomics2$grade), ]
t_proteomics3 <- t_proteomics2
colnames(t_proteomics3) = paste0(colnames(t_proteomics2),"_high")

df_3_0 = cbind(t_proteomics2[0,],t_proteomics3[0,])
df_4_0 = cbind(t_proteomics2[0,],t_proteomics3[0,])
df_5_0 = cbind(t_proteomics2[0,],t_proteomics3[0,])
df_4_3 = cbind(t_proteomics2[0,],t_proteomics3[0,])
df_5_3 = cbind(t_proteomics2[0,],t_proteomics3[0,])
df_5_4 = cbind(t_proteomics2[0,],t_proteomics3[0,])

for(i in 1:(dim(t_proteomics2)[1]-1)){
  pid = t_proteomics2[i,"Patient_ID"]
  sample_grade = t_proteomics2[i,"grade"]
  for(j in (i+1):dim(t_proteomics3)[1]){
    if(pid != t_proteomics3[j,"Patient_ID_high"]){
      break
    }
    sample_grade2 = t_proteomics3[j,"grade_high"]
    if(sample_grade2 > sample_grade){
      new_row = cbind(t_proteomics2[i,],t_proteomics3[j,])
      if(sample_grade2 == 3){
        df_3_0 = rbind(df_3_0,new_row)
      }else if(sample_grade == 0 & sample_grade2 == 4){
        df_4_0 = rbind(df_4_0,new_row)
      }else if(sample_grade == 0 & sample_grade2 == 5){
        df_5_0 = rbind(df_5_0,new_row)
      }else if(sample_grade == 3 & sample_grade2 == 4){
        df_4_3 = rbind(df_4_3,new_row)
      }else if(sample_grade == 3 & sample_grade2 == 5){
        df_5_3 = rbind(df_5_3,new_row)
      }else if(sample_grade == 4){
        df_5_4 = rbind(df_5_4,new_row)
      }
    }
  }
}

############## grade 3 - 0
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_3_0[,paste0(gene,"_high")]
  group2 <- df_3_0[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics3_0.csv"), row.names = FALSE)

############## grade 4 - 0
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_4_0[,paste0(gene,"_high")]
  group2 <- df_4_0[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics4_0.csv"), row.names = FALSE)

############## grade 5 - 0
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_5_0[,paste0(gene,"_high")]
  group2 <- df_5_0[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics5_0.csv"), row.names = FALSE)

############## grade 4 - 3
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_4_3[,paste0(gene,"_high")]
  group2 <- df_4_3[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics4_3.csv"), row.names = FALSE)

############## grade 5 - 3
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_5_3[,paste0(gene,"_high")]
  group2 <- df_5_3[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics5_3.csv"), row.names = FALSE)

############## grade 5 - 4
results2 <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- df_5_4[,paste0(gene,"_high")]
  group2 <- df_5_4[,gene]
  t_test <- t.test(group1, group2,paired=TRUE)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results2 <- rbind(results2, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results2$adj.p.value <- p.adjust(results2$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes2 <- results2 %>% filter(adj.p.value < 0.05)

significant_genes2$Gene
significant_genes2$Gene <- factor(significant_genes2$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes2, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes2, paste0(data_path,"significant_genes_ProstateProteomics5_4.csv"), row.names = FALSE)


############################
# Perform t-test for each gene
results <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(t_proteomics)[1:(ncol(t_proteomics)-1)]) { # Exclude grade, and status columns
  group1 <- t_proteomics %>% filter(grade > 0) %>% pull(gene)
  group2 <- t_proteomics %>% filter(grade == 0) %>% pull(gene)
  t_test <- t.test(group1, group2)
  logFC <- mean(log2(group1), na.rm=TRUE) - mean(log2(group2), na.rm=TRUE)
  results <- rbind(results, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results$adj.p.value <- p.adjust(results$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes <- results %>% filter(adj.p.value < 0.05)


significant_genes$Gene
significant_genes$Gene <- factor(significant_genes$Gene, levels = features)

# Create a bar plot for significant correlations
ggplot(significant_genes, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate proteins',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(significant_genes, paste0(data_path,"significant_genes_ProstateProteomics.csv"), row.names = FALSE)


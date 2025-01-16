library(AnnotationDbi)
library(org.Hs.eg.db)
library("readxl")
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggpubr)


# Data from this paper:
#https://www.nature.com/articles/s41467-021-26840-5
# Can be downloaded from here:
#https://zenodo.org/records/5546618


#------------------ preprocessing -----------------------------

data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'

pcatlas = read.table(paste0(data_path,'pcatlas_dataset_vst_normalized.txt'))
pcatlas_annotation = read_excel(paste0(data_path,'pcatlas_dataset_annotations.xlsx'))

pcatlas_annotation = pcatlas_annotation[pcatlas_annotation$`Tumor Tissue site` == 'PROSTATE',]


# List of ENSG gene annotations
ensg_ids <- rownames(pcatlas)


# Get gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensg_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")



features <- c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
              'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
              'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA-E', 'HLA-DRA',
              'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
              'PTGDS', 'HLA-A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
              'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB')


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


#-------------------logFC----------------------

set.seed(42)

# transpose
t_pcatlas <- transpose(pcatlas)

# get row and colnames in order
colnames(t_pcatlas) <- make.names(rownames(pcatlas))
rownames(t_pcatlas) <- colnames(pcatlas)

t_pcatlas$SampleType = pcatlas_annotation$`Sample Type`

t_pcatlas = t_pcatlas[t_pcatlas$SampleType %in% c('NORMAL', 'PRIMARY') ,]


# Encode status as a binary variable
t_pcatlas$status  <- ifelse(t_pcatlas$SampleType == 'NORMAL', 0, 1)

data = t_pcatlas

# Encode status as a binary variable

# Perform t-test for each gene
results <- data.frame(Gene=character(), logFC=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

for (gene in colnames(data)[1:(ncol(data)-2)]) { # Exclude grade, and status columns
  group1 <- data %>% filter(status == 1) %>% pull(gene)
  group2 <- data %>% filter(status == 0) %>% pull(gene)
  t_test <- t.test(group1, group2)
  logFC <- mean(group1) - mean(group2)
  results <- rbind(results, data.frame(Gene=gene, logFC=logFC, p.value=t_test$p.value))
}

# Apply FDR correction
results$adj.p.value <- p.adjust(results$p.value, method = "fdr")

# Filter significant genes (adjusted p-value < 0.05)
significant_genes <- results %>% filter(adj.p.value < 0.05)


significant_genes$Gene
significant_genes$Gene <- factor(significant_genes$Gene, levels = c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
                                                                    'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
                                                                    'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA.E', 'HLA.DRA',
                                                                    'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
                                                                    'PTGDS', 'HLA.A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
                                                                    'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB'))
write.csv(significant_genes, paste0(data_path,"significant_genes_ProstateRNA.csv"), row.names = FALSE)

# Create a bar plot for significant correlations
ggplot(significant_genes, aes(x=Gene, y=logFC)) +
  geom_bar(stat="identity") +
  labs(title='Prostate RNA',
       x="Gene",
       y="lofGC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



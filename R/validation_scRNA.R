library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(cowplot)
library(patchwork)

###########################################

# Article for data:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/

###########################################
data_path = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/"
file_name = "dge_E.rds"

dge <- readRDS(paste0(data_path,file_name))

all_features <- rownames(dge)

features <- c('TMEFF2', 'SPON2', 'TRGC1', 'RPS8', 'RPS4X', 'CRISP3', 'AMACR', 'MYO6',
              'STEAP2', 'RPL7A', 'RPS12', 'RPL14', 'RPS27A', 'EEF1G', 'KLK3', 'NPY',
              'KLK2', 'IGKC', 'MT2A', 'LGALS1', 'NEFH', 'C1R', 'HLA-E', 'HLA-DRA',
              'DES', 'IGFBP7', 'COL1A2', 'LGALS3BP', 'HSPB1', 'AZGP1', 'GSTP1',
              'PTGDS', 'HLA-A', 'TAGLN', 'ACTG2', 'B2M', 'FLNA', 'MYL9', 'IFITM3',
              'CD74', 'TPM2', 'ACTA2', 'TIMP1', 'S100A6', 'MSMB')

features <- features[features %in% all_features]

patients <- c("PR5249", "PR5251", "PR5254", "PR5261")
samples <- c("PR5249_T", "PR5251_T", "PR5254_T", "PR5261_T","PR5249_N", "PR5251_N", "PR5254_N", "PR5261_N")

dge@meta.data$orig.ident2 <- dge@meta.data$orig.ident

dge@meta.data$orig.ident2[dge$sample %in% names(dge$sample[dge$orig.ident=="AUG_PB1"])[1:167]] <- "AUG_PB1A"
dge@meta.data$orig.ident2[dge@meta.data$orig.ident2 =="AUG_PB1"] <- "AUG_PB1B"

dge@meta.data$orig.ident2[dge$sample %in% names(dge$sample[dge$orig.ident=="MAY_PB1"])[1:355]] <- "MAY_PB1A"
dge@meta.data$orig.ident2[dge@meta.data$orig.ident2 =="MAY_PB1"] <- "MAY_PB1B"

dge@meta.data$orig.ident2[dge$sample %in% names(dge$sample[dge$orig.ident=="MAY_PB2"])[1:466]] <- "MAY_PB2A"
dge@meta.data$orig.ident2[dge@meta.data$orig.ident2 =="MAY_PB2"] <- "MAY_PB2B"

samples <- sort(unique(dge@meta.data$orig.ident2))

markers_list <- list()

for(sample_var in samples){
  dge@meta.data$group <- ifelse(dge@meta.data$malignancy == "Non-Malignant", "N",
                                ifelse(dge@meta.data$malignancy == "Malignant", "T", NA))
  dge@meta.data$group <- ifelse(dge@meta.data$orig.ident2 %in% c(sample_var), dge@meta.data$group, NA)
  markers <- FindMarkers(dge, ident.1 = "T", ident.2 = "N", group.by = "group", features=features)
  markers_list[[sample_var]] <- markers
}

dge@meta.data$group <- ifelse(dge@meta.data$malignancy == "Non-Malignant", "N",
                              ifelse(dge@meta.data$malignancy == "Malignant", "T", NA))

n_patients <- c("PR5249_N", "PR5251_N", "PR5254_N", "PR5261_N")
t_patients <- c("PR5249_T", "PR5251_T", "PR5254_T", "PR5261_T")

dge@meta.data$group <- ifelse(dge@meta.data$orig.ident2 %in% samples, dge@meta.data$group, NA)
markers <- FindMarkers(dge, ident.1 = "T", ident.2 = "N", group.by = "group", features=features)
markers_list[["All"]] <- markers

p_values_list <- list()
for (id in names(markers_list)) {
  markers <- markers_list[[id]]
  
  # Extract the p-values and create a data frame
  p_values_df <- data.frame(gene = rownames(markers), p_val = ifelse(markers$p_val_adj <= 0.05,markers$avg_log2FC,NaN))
  
  # Set the column name to the ID_coarse
  colnames(p_values_df)[2] <- id
  
  # Add to the list
  p_values_list[[id]] <- p_values_df
}

# Merge the list of data frames by 'gene'
p_values_matrix <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), p_values_list)

p_values_long <- melt(p_values_matrix)

p_values_long$gene <- factor(p_values_long$gene, levels = rev(unique(features)))

p_values_long = p_values_long[p_values_long$variable != "All",] #removing "All" column
# p_values_long$variable != "All" # to check that "All" column is removed

write.csv(p_values_long, file = "InterestingFoldChanges2.csv")

ggplot(p_values_long, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = c(rep("darkred",1), rep("red",1), "white", rep("blue",1),rep("darkblue",1)),  # c("lightblue", "blue", "darkblue")# Color gradient
                       name = "avg_log2FC", limits = c(-10,10)) +  # Title for the legend
  theme_minimal() +  # Clean theme
  labs(x = "Samples", y = "Genes") +  # Axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        axis.text.y = element_text(size = 10))

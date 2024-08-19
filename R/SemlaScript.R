# Loading libraries
library(semla)
library(Seurat)

# path and file name to the RData object
data_path <- "/Spatial_CCC/DATA/"
file_name <- "seurat_st_NCBI569"

# Loading data
load(paste0(data_path,file_name,".RData"))

# Updating object to use semla
seurat_st2 <- UpdateSeuratForSemla(seurat_st)
seurat_st2 <- LoadImages(seurat_st2)

# Use semla and save result
se <- FeatureViewer(seurat_st2)
save(se, file=paste0(data_path,file_name,"_2.RData"))


# Same code but used when we have previously annotated and saved a part of the object
load(paste0(data_path,file_name,"_2.RData"))
se<- UpdateSeuratForSemla(se)
se <- LoadImages(se)
se <- FeatureViewer(se)
save(se, file=paste0(data_path,file_name,"_3.RData"))

# Code used when previous annotations are considered wrong
se@meta.data$Label[se@meta.data$Label == "wd"] = "Old_wd"
se@meta.data$Label[se@meta.data$Label == "md"] = "Old_md"
se@meta.data$Label[se@meta.data$Label == "Benign"] = "Old_Benign"
se@meta.data$Label[se@meta.data$Label == "Normal"] = "Old_Normal"

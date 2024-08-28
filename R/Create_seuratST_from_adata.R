rm(list=ls())
library(SeuratDisk)
library(Seurat)
library(reticulate)
library(zellkonverter)
library(dplyr)
library(semla)
library(reticulate)

use_python("C:\\Users\\dandia\\AppData\\Local\\anaconda3\\envs\\NewEnv", required = TRUE)
scanpy <- import("scanpy")
anndata <- import("anndata")

# example seurat_st_example####
path0 = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/"
file_name = "seurat.rdata"
load(paste0(path0,file_name))
# load 10x h5 data from path2
# seurat_st_example = Load10X_Spatial(
#   data.dir = path0,
#   filename = "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5",
#   assay = "Spatial",
#   slice = "slice1")
seurat_st_example <- subset(seurat_st_example, nCount_Spatial>0)
seurat_st_example
#View(seurat_st_example)

#SpatialPlot(seurat_st_example, features = "nFeature_Spatial",crop = T,alpha = 0.1)
#SpatialPlot(seurat_st_example2, features = "nFeature_Spatial",crop = T,alpha = 0.1)


#read sample of interest '.h5ad file from path1' and create seurat object####
path = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/prostate/NotAnnotated/"
path2 = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/prostate/NotAnnotated/NewFiles/"

all_hest_samples = list.files(path2,pattern = '.h5ad')
all_hest_samples = gsub('.h5ad','',all_hest_samples)

for (sample in all_hest_samples) {
  print(sample)
  path1 = paste0(path2,sample,'.h5ad')
  adata <- anndata$read_h5ad(path1)
  mtx = as.matrix(adata$X)
  mtx = t(mtx)
  mtx = as(mtx, "dgCMatrix")
  metadata = adata$obs
  rownames(mtx) = rownames(adata$var)
  colnames(mtx) = rownames(metadata)
  seurat_st <- CreateSeuratObject(counts = mtx, meta.data = metadata)

  seurat_st@assays$Spatial = seurat_st@assays$RNA
  seurat_st@assays$RNA = NULL
  seurat_st@assays$Spatial@key = 'spatial_'
  seurat_st@meta.data$nCount_Spatial = seurat_st@meta.data$nCount_RNA
  seurat_st@meta.data$nFeature_Spatial = seurat_st@meta.data$nFeature_RNA
  seurat_st@meta.data <- seurat_st@meta.data[,!colnames(seurat_st@meta.data) %in% c('nCount_RNA','nFeature_RNA')]
  seurat_st@active.assay = 'Spatial'
  
  seurat_st@images = seurat_st_example@images
  seurat_st@images$slice1@image = adata$uns[["spatial"]][["ST"]][["images"]][["downscaled_fullres"]]/255
  seurat_st@images[["slice1"]]@scale.factors[["spot"]] = adata$uns[["spatial"]][["ST"]][["scalefactors"]][["spot_diameter_fullres"]]
  seurat_st@images[["slice1"]]@scale.factors[["lowres"]] = adata$uns[["spatial"]][["ST"]][["scalefactors"]][["tissue_downscaled_fullres_scalef"]]
  seurat_st@images[["slice1"]]@scale.factors[["fiducial"]] = NULL
  seurat_st@images[["slice1"]]@scale.factors[["hires"]] = NULL
  seurat_st@images[["slice1"]]@coordinates = data.frame(matrix(ncol = 0, nrow = dim(seurat_st)[2]))
  rownames(seurat_st@images[["slice1"]]@coordinates) = rownames(seurat_st@meta.data)
  if(!is.null(adata$obs[["in_tissue"]])){
    seurat_st@images[["slice1"]]@coordinates[["tissue"]] = as.numeric(adata$obs[["in_tissue"]]) 
  }
  seurat_st@images[["slice1"]]@coordinates[["row"]] = adata$obs[["array_row"]]
  seurat_st@images[["slice1"]]@coordinates[["col"]] = adata$obs[["array_col"]]
  seurat_st@images[["slice1"]]@coordinates[["imagerow"]] = adata$obs[["pxl_row_in_fullres"]]
  seurat_st@images[["slice1"]]@coordinates[["imagecol"]] = adata$obs[["pxl_col_in_fullres"]]
  seurat_st@images[["slice1"]]@spot.radius
  
  pdf(paste0(path,sample,'.pdf'))
  SpatialPlot(seurat_st, features = "nFeature_Spatial",crop = T,alpha = 0.5)
  dev.off()
  
  seurat_semla <- UpdateSeuratForSemla(seurat_st,image_type = c("tissue_lowres"))
  se <- LoadImages(seurat_semla)
  
  st <- GetStaffli(se)
  new_h <- st@image_info$full_width
  st@image_info$full_width <- st@image_info$full_height
  st@image_info$full_height <- new_h
  se@tools$Staffli = st
  save(se, file = paste0(path,sample,'.rdata'))
}

# CODE FOR ONE SINGLE FILE

#load(paste0(wd,"/pancreatic_cancer/Input_data_files/hest_data/seurat/seurat_st_",sample,'.rdata'))
path = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/prostate/NotAnnotated/"
sample = "INT28"
path1 = paste0(path,sample,'.h5ad')
adata <- anndata$read_h5ad(path1)
adata2 <- readH5AD(path1)
View(adata)
#adata.X.toarray()
mtx = as.matrix(adata$X)
mtx = t(mtx)
mtx = as(mtx, "dgCMatrix")
metadata = adata$obs
rownames(mtx) = rownames(adata$var)
colnames(mtx) = rownames(metadata)
seurat_st <- CreateSeuratObject(counts = mtx, meta.data = metadata)

seurat_st@assays$Spatial = seurat_st@assays$RNA
seurat_st@assays$RNA = NULL
seurat_st@assays$Spatial@key = 'spatial_'
seurat_st@meta.data$nCount_Spatial = seurat_st@meta.data$nCount_RNA
seurat_st@meta.data$nFeature_Spatial = seurat_st@meta.data$nFeature_RNA
seurat_st@meta.data <- seurat_st@meta.data[,!colnames(seurat_st@meta.data) %in% c('nCount_RNA','nFeature_RNA')]
seurat_st@active.assay = 'Spatial'

seurat_st@images = seurat_st_example@images
seurat_st@images$slice1@image = adata$uns[["spatial"]][["ST"]][["images"]][["downscaled_fullres"]]/255
seurat_st@images[["slice1"]]@scale.factors[["spot"]] = adata$uns[["spatial"]][["ST"]][["scalefactors"]][["spot_diameter_fullres"]]
seurat_st@images[["slice1"]]@scale.factors[["lowres"]] = adata$uns[["spatial"]][["ST"]][["scalefactors"]][["tissue_downscaled_fullres_scalef"]]
seurat_st@images[["slice1"]]@scale.factors[["fiducial"]] = NULL
seurat_st@images[["slice1"]]@scale.factors[["hires"]] = NULL
seurat_st@images[["slice1"]]@coordinates = data.frame(matrix(ncol = 0, nrow = dim(seurat_st)[2]))
rownames(seurat_st@images[["slice1"]]@coordinates) = rownames(seurat_st@meta.data)
seurat_st@images[["slice1"]]@coordinates[["tissue"]] = as.numeric(adata$obs[["in_tissue"]])
seurat_st@images[["slice1"]]@coordinates[["row"]] = adata$obs[["array_row"]]
seurat_st@images[["slice1"]]@coordinates[["col"]] = adata$obs[["array_col"]]
seurat_st@images[["slice1"]]@coordinates[["imagerow"]] = adata$obs[["pxl_row_in_fullres"]]
seurat_st@images[["slice1"]]@coordinates[["imagecol"]] = adata$obs[["pxl_col_in_fullres"]]
seurat_st@images[["slice1"]]@spot.radius

pdf(paste0(path,sample,'.pdf'))

SpatialPlot(seurat_st, features = "nFeature_Spatial",crop = T,alpha = 0.5)
dev.off()
# Run semla####

seurat_semla <- UpdateSeuratForSemla(seurat_st,image_type = c("tissue_lowres"))
se <- LoadImages(seurat_semla)
se <- FeatureViewer(se, launch.browser = .rs.invokeShinyWindowViewer)
save(se, file = paste0(path,sample,'.rdata'))


st <- GetStaffli(se)
new_h <- st@image_info$full_width
st@image_info$full_width <- st@image_info$full_height
st@image_info$full_height <- new_h
se@tools$Staffli = st
se <- FeatureViewer(se, launch.browser = .rs.invokeShinyWindowViewer)
save(se, file = paste0(path,sample,'.rdata'))

#se@images$slice1@coordinates$imagerow = se@images$slice1@coordinates$imagerow +100
#se@images$slice1@coordinates$imagecol = se@images$slice1@coordinates$imagecol +100

#se@images$slice1@scale.factors$lowres = 0.060 
se <- UpdateSeuratForSemla(se,image_type = c("tissue_lowres"))
se <- LoadImages(se)
se <- FeatureViewer(se)

save(se, file = paste0(path,sample,'.rdata'))

# CODE TO VERIFY THE Rdata IS GOOD TO USE IN SEMLA
adjust_image <- function(se){
  st <- GetStaffli(se)
  new_h <- st@image_info$full_width
  st@image_info$full_width <- st@image_info$full_height
  st@image_info$full_height <- new_h
  se@tools$Staffli = st
  return(se)
}
path2 = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/prostate/NotAnnotated/NewFiles/"
path = "C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/prostate/NotAnnotated/"

all_hest_samples = list.files(path2,pattern = '.h5ad')
all_hest_samples = gsub('.h5ad','',all_hest_samples)

for (sample in all_hest_samples) {
  print(sample)
  load(paste0(path,sample,".RData"))
  se <- UpdateSeuratForSemla(se)
  se <- LoadImages(se)
  se <- adjust_image(se)
  se <- FeatureViewer(se, launch.browser = .rs.invokeShinyWindowViewer)
}


  
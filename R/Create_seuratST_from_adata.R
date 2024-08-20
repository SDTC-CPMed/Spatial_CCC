rm(list=ls())
library(SeuratDisk)
library(Seurat)
library(reticulate)
library(zellkonverter)
library(dplyr)
library(semla)
library(reticulate)
# virtualenv_create("r-reticulate")
# virtualenv_install("r-reticulate", "anndata")
use_virtualenv("r-reticulate")
anndata <- import("anndata")

wd = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects'

# example seurat_st_example####
path0 = paste0(wd,'/Spatial/SpatialInferCNV/10x breastcancer')
# load 10x h5 data from path2
seurat_st_example = Load10X_Spatial(
  data.dir = path0,
  filename = "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1")
seurat_st_example <- subset(seurat_st_example, nCount_Spatial>0)
seurat_st_example
View(seurat_st_example)

#SpatialPlot(seurat_st_example, features = "nFeature_Spatial",crop = T,alpha = 0.1)
#SpatialPlot(seurat_st_example2, features = "nFeature_Spatial",crop = T,alpha = 0.1)

#read sample of interest '.h5ad file from path1' and create seurat object####
all_hest_samples = list.files(paste0(wd,'/Spatial/Spatial_CCC/data/hest_data_prostate/st/'),pattern = '.h5ad')
all_hest_samples = gsub('.h5ad','',all_hest_samples)

for (sample in all_hest_samples) {
  #sample = 'MEND156'
  path1 = paste0(wd,'/Spatial/Spatial_CCC/data/hest_data_prostate/st/',sample,'.h5ad')
  adata <- anndata$read_h5ad(path1)
  #View(adata)
  adata.X.toarray()
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
  
  pdf(paste0(wd,"/Spatial/Spatial_CCC/data/hest_data_prostate/st/seurat_st_spatialplot_",sample,'.pdf'))
  SpatialPlot(seurat_st, features = "nFeature_Spatial",crop = T,alpha = 0.5)
  dev.off()

  # save seurat_st as rdata
  save(seurat_st, file = paste0(wd,"/Spatial/Spatial_CCC/data/hest_data_prostate/seurat/seurat_st_",sample,'.rdata'))
}
  #load(paste0(wd,"/pancreatic_cancer/Input_data_files/hest_data/seurat/seurat_st_",sample,'.rdata'))

# Run semla####
#remotes::install_github("ludvigla/semla")
library(semla)
seurat_semla <- UpdateSeuratForSemla(seurat_st,image_type = c("tissue_lowres"))
se <- FeatureViewer(seurat_semla)


  
  
  
  
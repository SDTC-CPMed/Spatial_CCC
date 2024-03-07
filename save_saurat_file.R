library(Seurat)
library(SeuratDisk)



############ Load data and save for use in python ###########

load('data/seurat_ST1.rda')

write.table(as.matrix(GetAssayData(object = seurat, slot = "counts")), 
            'seurat_ST1_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

write.table(seurat@images$slice1@coordinates[,1:3], 'data/seurat_ST1_obs.csv')
write.table(seurat@images$slice1@coordinates[,4:5], 'data/seurat_ST1_obsm.csv')
write.table(seurat@images$slice1@image[,,1], 'data/seurat_ST1_imageR.csv')
write.table(seurat@images$slice1@image[,,2], 'data/seurat_ST1_imageG.csv')
write.table(seurat@images$slice1@image[,,3], 'data/seurat_ST1_imageB.csv')
write.table(seurat@images$slice1@scale.factors$spot, 'data/seurat_ST1_spot_diameter.csv')
write.table(seurat@images$slice1@scale.factors$fiducial, 'data/seurat_ST1_fiducial_diameter.csv')
write.table(seurat@images$slice1@scale.factors$lowres, 'data/seurat_ST1_hires_scalef.csv')
write.table(seurat@images$slice1@scale.factors$lowres, 'data/seurat_ST1_lowres_scalef.csv')
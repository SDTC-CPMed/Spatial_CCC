library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(scriabin)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)

#https://www.nature.com/articles/s41587-023-01782-z
#https://github.com/BlishLab/scriabin
#https://drive.google.com/drive/folders/1dkGF4kKMbHSm04DKHg6P_y_iTD1pXCH8

############ Load data and save for use in python ###########

load('data/seurat_ST1.rda')

SaveH5Seurat(seurat, filename = "seurat_ST1.h5Seurat")
Convert("seurat_ST1.h5Seurat", assay = 'Spatial', dest = "h5ad")
write.table(seurat@images$slice1@coordinates[,1:3], 'data/seurat_ST1_obs.csv')
write.table(seurat@images$slice1@coordinates[,4:5], 'data/seurat_ST1_obsm.csv')
write.table(seurat@images$slice1@image[,,1], 'data/seurat_ST1_imageR.csv')
write.table(seurat@images$slice1@image[,,2], 'data/seurat_ST1_imageG.csv')
write.table(seurat@images$slice1@image[,,3], 'data/seurat_ST1_imageB.csv')
write.table(seurat@images$slice1@scale.factors$spot, 'data/seurat_ST1_spot_diameter.csv')
write.table(seurat@images$slice1@scale.factors$fiducial, 'data/seurat_ST1_fiducial_diameter.csv')
write.table(seurat@images$slice1@scale.factors$lowres, 'data/seurat_ST1_hires_scalef.csv')
write.table(seurat@images$slice1@scale.factors$lowres, 'data/seurat_ST1_lowres_scalef.csv')


############ Basic processing of Spatial data ###########
plot1 <- VlnPlot(seurat, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

seurat <- SCTransform(seurat, assay = "Spatial", verbose = FALSE)

seurat <- RunPCA(seurat, assay = "SCT", verbose = FALSE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)
p1 <- DimPlot(seurat, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat, label = TRUE, label.size = 3)
p1
p2





########### Connections between clusters 1 and 2 ############# 

# this is the method that I had in mind to use. Instead of using
# cluster 1 and cluster 2 you can use spot in the middle as sender
# and spots around as receivers. 
# Then make a loop going through each spot as sender and hopefully
# the senders and receivers would be somewhat consistent within
# different cell types

seurat_ALRA <- SeuratWrappers::RunALRA(seurat)
scriabin::load_nichenet_database()




variant_genes <- IDVariantGenes(seurat_ALRA, assay = "alra", group.by = "seurat_clusters")
seurat_ALRA@assays$data <- seurat_ALRA@assays$alra #Fixing the new version of packages
gene_signature <- GenerateCellSignature(seurat_ALRA, variant_genes = variant_genes)
active_ligands <- RankActiveLigands(seurat_ALRA, assay = "data", signature_matrix = gene_signature)


i = 1
cell = seurat@images$slice1@coordinates[i,]
interaction = array('non-relevant', dim = dim(seurat@images$slice1@coordinates)[1])

receiver_idx = abs(cell$row - seurat@images$slice1@coordinates$row) <= 1 &
           abs(cell$col - seurat@images$slice1@coordinates$col) <= 2

interaction[receiver_idx] = 'receiver'
interaction[i] = 'sender'


seurat_ALRA[['interaction']] = as.factor(interaction)

SpatialFeaturePlot(seurat_ALRA, features  = 'interaction')

soi <- 1 #sender of interest
roi <- 2 #receiver of interest
TopLigandsByIdent(seurat_ALRA, assay = 'data', active_ligands = active_ligands, 
                  sender = soi, receiver = roi, group.by = "seurat_clusters")


receiver_cells <- colnames(seurat_ALRA)[seurat_ALRA$seurat_clusters==roi]

# calculates the predicted target genes within a set of receiver cells 
PlotLigandTargetAlluvium(seurat_ALRA, signature_matrix = gene_signature,
                         active_ligands = active_ligands, receiver_cells = receiver_cells,
                         ligands_of_interest = c("WNT7A", "WNT7B"))




########## Another tool of Scriabin, not sure how could be used ##############


#find interaction programs
Programs <- FindAllInteractionPrograms(seurat_ALRA, iterate.threshold = 300, group.by = "seurat_clusters",
                                      assay = "alra", r2_cutoff = 0.6, sim_threshold = 0.4)
#test for interaction program significance
Programs_sig <- InteractionProgramSignificance(Programs, n.replicate = 500, min.members = 1)
#keep IP that are significant in at least one cell type
#in this example, all programs are significant in at least one cell type
ip_pvals <- Programs_sig %>% as_tibble() %>%
  dplyr::select(name,ends_with("pval")) %>% unique() %>%
  pivot_longer(!name, names_to = "seurat_clusters", values_to = "pval") %>%
  group_by(name) %>% dplyr::mutate(min_p = min(pval)) %>%
  dplyr::select(name,min_p) %>% unique() %>% 
  dplyr::filter(min_p < 0.05) %>% pull(name)


#score cells by expression of interaction program
Programs <- ScoreInteractionPrograms(seurat_ALRA, Programs_sig)


#a = ScoreInteractionPrograms(panc_id, panc_ip_sig, return.assay = F)

seurat_ALRA_ip_lig <- as.matrix(Programs[["IPligands"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = Programs$seurat_clusters) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))

# A plot of different programs in columns and how much they are expressed as ligands.
Heatmap(seurat_ALRA_ip_lig, show_column_names = F, name = "Ligands", column_names_gp = gpar(fontsize =10))

seurat_ALRA_ip_rec <- as.matrix(Programs[["IPreceptors"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = Programs$seurat_clusters) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))

# A plot of different programs in columns and how much they are expressed as receptors.
Heatmap(seurat_ALRA_ip_rec, show_column_names = F, name = "Receptors")




poi = 'IP-1'

IPFeaturePlot(Programs, ip = poi)

poi_ligand_score = data.frame(Programs@assays$IPligands@data)
poi_ligand_score = poi_ligand_score %>% filter(rownames(poi_ligand_score) %in% c(poi))

Programs[['poi_ligand_score']] = as.numeric(unlist(poi_ligand_score))
p1 = SpatialFeaturePlot(Programs, features  = 'poi_ligand_score')

poi_receptor_score = data.frame(Programs@assays$IPreceptors@data)
poi_receptor_score = poi_receptor_score %>% filter(rownames(poi_receptor_score) %in% c(poi))

Programs[['poi_receptor_score']] = as.numeric(unlist(poi_receptor_score))
p2 = SpatialFeaturePlot(Programs, features  = 'poi_receptor_score')
p1+p2
p2
moi <- reshape2::melt(Programs_sig %>% dplyr::filter(name==poi) %>%
                        select("lr_pair",contains("connectivity"))) %>% arrange(-value)
moi$lr_pair <- factor(moi$lr_pair, levels = unique(moi$lr_pair))
ggplot(moi, aes(x = lr_pair, y = value, color = variable)) + 
  geom_point() + theme_cowplot() + ggpubr::rotate_x_text() + labs(x = NULL, y = "Intramodular\nconnectivity")

ip_by_celltype <- IPCellTypeSummary(Programs, group.by = "seurat_clusters")
ip_by_celltype %>% group_by(sender) %>% top_n(n = 1, wt = additive.score)



